#include "decimate.h"
#include "collapse_edge_seam.h"
#include <igl/edge_flaps.h>
#include <igl/remove_unreferenced.h>
#include <igl/slice_mask.h>
#include <igl/connect_boundary_to_infinity.h>
#include <igl/seam_edges.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <chrono>
#include <igl/hausdorff.h>
#include <igl/seam_edges.h>
#include "cost_and_placement.h"

bool debug = false;

void clean_mesh(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & TC,
	const Eigen::MatrixXi & FT,
	const int nF,
	Eigen::MatrixXd & V_out,
	Eigen::MatrixXi & F_out,
	Eigen::MatrixXd & TC_out,
	Eigen::MatrixXi & FT_out) 
{
	using namespace Eigen;
	using namespace igl;
	MatrixXi F2(nF,3);
	MatrixXi FT2(nF,3);
	int m = 0;
	for(int f = 0;f<nF;f++)
	{
		if( F(f,0) != DUV_COLLAPSE_EDGE_NULL || 
			F(f,1) != DUV_COLLAPSE_EDGE_NULL || 
			F(f,2) != DUV_COLLAPSE_EDGE_NULL)
		{
			assert( FT(f,0) != DUV_COLLAPSE_EDGE_NULL );
			assert( FT(f,1) != DUV_COLLAPSE_EDGE_NULL );
			assert( FT(f,2) != DUV_COLLAPSE_EDGE_NULL );

			F2.row(m) = F.row(f);
			FT2.row(m) = FT.row(f);

			m++;
		}
	}
	
	F2.conservativeResize(m,F2.cols());
	FT2.conservativeResize(m,FT2.cols());
	VectorXi _1;
	remove_unreferenced(V,F2,V_out,F_out,_1);
	remove_unreferenced(TC,FT2,TC_out,FT_out,_1);
}

void prepare_decimate_halfedge_5d(
	const Eigen::MatrixXd & OV,
    const Eigen::MatrixXi & OF,
    const Eigen::MatrixXd & OTC,
    const Eigen::MatrixXi & OFT,
    EdgeMap & seam_edges,
    MapV5d & Vmetrics,
    int & target_num_vertices,
    const int seam_aware_degree,
    // output
    Eigen::MatrixXd & V,
	Eigen::MatrixXi & F,
	Eigen::MatrixXd & TC,
	Eigen::MatrixXi & FT,
    Eigen::VectorXi & EMAP,
	Eigen::MatrixXi & E,
	Eigen::MatrixXi & EF,
	Eigen::MatrixXi & EI,
    PriorityQueue & Q,
	std::vector<PriorityQueue::iterator > & Qit,
	std::vector< placement_info_5d > & C
	)
{
	using namespace Eigen;
	using namespace std;
	using namespace igl;
	
	// Working copies
	V = OV;
	F = OF;
	igl::connect_boundary_to_infinity(OV,OF,V,F);
	const bool has_infinity_vertex = V.row( V.rows()-1 ).minCoeff() == std::numeric_limits<double>::infinity();
	if( has_infinity_vertex )	target_num_vertices++;
	TC = OTC;
	FT = OFT;
	
	// priorityQueue
	edge_flaps(F,E,EMAP,EF,EI);
	
	if( has_infinity_vertex ){
        // Let's add infinity faces to FT and an infinity vertex to TC.
        // We are doing this with intimate knowledge of how connect_boundary_to_infinity()
        // works.
        
        /// 1 Add a texture coordinate at infinity to TC with a zero quadric.
        /// 2 Iterate over newly added faces in F (the bottom ones that aren't in OF).
        /// 3 Find the edge between the two non-infinite vertices.
        /// 4 Find the original face opposite this edge.
        /// 5 Add a corresponding new face to FT referencing the same texture coordinates as
        ///   the same vertices in the opposite face, connected to the texture coordinate
        ///   at infinity.
        
        /// 1
        TC.conservativeResize( OTC.rows() + 1, Eigen::NoChange );
        auto inf = std::numeric_limits<double>::infinity();
        TC.row( OTC.rows() ).setConstant( inf );
        // Add a zero quadric.
        Vmetrics[OV.rows()][OTC.rows()].setZero( 6,6 );
        
        
        // Allocate space for the new faces added by step 5.
        FT.conservativeResize( F.rows(), Eigen::NoChange );
        
        /// 2
        for( int fi = OF.rows(); fi < F.rows(); ++fi ) {
            // In connect_boundary_to_infinity(), the new face always has the infinity
            // vertex last.
            const int fi_vinf = 2;
            assert( F( fi, fi_vinf ) == OV.rows() );
            /// 3
            // The edge is in EMAP across from the infinity vertex in the face.
            const int e = EMAP( fi_vinf*F.rows() + fi );
            // Make sure this is the right edge; one side should point to `fi`.
            assert( EF(e,0) == fi || EF(e,1) == fi );
            /// 4
            // Get the opposite face index.
            const int fi_opp =
                EF(e,0) == fi
                ? EF(e,1)
                : EF(e,0)
                ;
            // Find the index of fi's vertex fi_vinf + 1 in fi_opp.
            int fi_opp_v1 = 0;
            for( ; fi_opp_v1 < 3; ++fi_opp_v1 ) {
                if( F( fi_opp, fi_opp_v1 ) == F( fi, ( fi_vinf + 1 ) % 3 ) ) break;
            }
            assert( fi_opp_v1 < 3 );
            // Since fi and fi_opp have the opposite orientation,
            // the other vertex we want, fi_vinf + 2, is at (fi_opp_v1-1).
            const int fi_opp_v2 = ( fi_opp_v1 - 1 + 3 ) % 3;
            
            /// 5
            FT( fi, fi_vinf ) = OTC.rows();
            FT( fi, ( fi_vinf+1 ) % 3 ) = FT( fi_opp, fi_opp_v1 );
            FT( fi, ( fi_vinf+2 ) % 3 ) = FT( fi_opp, fi_opp_v2 );
        }
    }
    
	Qit.resize(E.rows());
	// If an edge were collapsed, we'd collapse it to these points:
	C.resize( E.rows() );
	cout << "# edges: " << C.size() << endl;
	
	for ( int e=0; e<E.rows(); e++ )
	{
		double cost = -31337;
		placement_info_5d new_placement;
		Bundle b = get_half_edge_bundle( e, E, EF, EI, F, FT );
		cost_and_placement_qslim5d_halfedge(b,V,F,TC,FT,seam_edges,Vmetrics,seam_aware_degree,cost,new_placement);
		C.at(e) = new_placement;
		Qit[e] = Q.insert(std::pair<double,int>(cost,e)).first;
	}
	assert( Q.size() == E.rows() );
	std::cout << "building PriorityQueue succeeds." << std::endl;		
}


bool collapse_one_edge(
	Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXd & TC,
    Eigen::MatrixXi & FT,
    Eigen::VectorXi & EMAP,
	Eigen::MatrixXi & E,
	Eigen::MatrixXi & EF,
	Eigen::MatrixXi & EI,
    EdgeMap & seam_edges,
    MapV5d & Vmetrics,
    const int seam_aware_degree,
	PriorityQueue & Q, 
	std::vector<PriorityQueue::iterator > & Qit, 
	std::vector< placement_info_5d > & C, 
	int & prev_e)
{
	using namespace std;
	using namespace Eigen;
	using namespace igl;
	
	bool success = false;
	int e;
	while(true) {
		if(Q.empty())
		{
			cout << "empty queue" << endl;
			break;
		}
		if(Q.begin()->first == std::numeric_limits<double>::infinity())
		{
			// min cost edge is infinite cost
			cout << "min cost edge is infinite cost" << endl;
			break;
		}

		if(collapse_edge_with_uv(V,F,E,EMAP,EF,EI,TC,FT,seam_edges,Vmetrics,seam_aware_degree,Q,Qit,C,e,debug))
		{
			success = true;
			break;
		} 
		else if(prev_e == e) 
		{
			assert(false && "Edge collapse no progress... bad stopping condition?");
			break;
		}
	}
	prev_e = e;
			
	return success;
}

bool decimate_halfedge_5d(
    const Eigen::MatrixXd & OV,
    const Eigen::MatrixXi & OF,
    const Eigen::MatrixXd & OTC,
    const Eigen::MatrixXi & OFT,
    EdgeMap & seam_edges,
    MapV5d & Vmetrics,
    int target_num_vertices,
    const int seam_aware_degree,
    Eigen::MatrixXd & V_out,
    Eigen::MatrixXi & F_out,
    Eigen::MatrixXd & TC_out,
    Eigen::MatrixXi & FT_out
    )
{
	using namespace Eigen;
	using namespace std;
	using namespace igl;	
	
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd TC;
	Eigen::MatrixXi FT;
	Eigen::VectorXi EMAP;
	Eigen::MatrixXi E;
	Eigen::MatrixXi EF;
	Eigen::MatrixXi EI;
	PriorityQueue Q;
	std::vector<PriorityQueue::iterator > Qit;
	std::vector< placement_info_5d > C;
	prepare_decimate_halfedge_5d(OV,OF,OTC,OFT,seam_edges,Vmetrics,target_num_vertices,seam_aware_degree,
			V,F,TC,FT,EMAP,E,EF,EI,Q,Qit,C);
	
	int prev_e = -1;
	bool clean_finish = false;
	int remain_vertices=V.rows();
	
	int next_output_target = V.rows();
	int suffix = 0;
	while(true)
	{		
		if(Q.empty())
		{
			cout << "empty queue" << endl;
			break;
		}
		if(Q.begin()->first == std::numeric_limits<double>::infinity())
		{
			// min cost edge is infinite cost
			cout << "min cost edge is infinite cost, left nums of vertices: " << remain_vertices << endl;
			break;
		}
		
		bool collapse_success = collapse_one_edge(V,F,TC,FT,EMAP,E,EF,EI,seam_edges,Vmetrics,seam_aware_degree,Q,Qit,C,prev_e);
		if(!collapse_success) {
			clean_finish = false;
			break;
		}
		else if(remain_vertices <= target_num_vertices) {
			clean_finish = true;
			break;
		}
		remain_vertices--;
	}

	// remove all DUV_COLLAPSE_EDGE_NULL faces
	clean_mesh(V,F,TC,FT,OF.rows(),V_out,F_out,TC_out,FT_out);
	return clean_finish;
	
}
    
