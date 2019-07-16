#include <igl/circulation.h>
#include <igl/edge_collapse_is_valid.h>
#include "collapse_edge_seam.h"
#include "decimate.h"
#include "neighbor_faces_and_boundary.h"
#include "detect_foldover.h"
#include "cost_and_placement.h"

bool try_collapse_5d_Edge(
	const int e,
    const placement_info_5d & new_placement, // vertex position, texture coordinate, metric
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
  	Eigen::VectorXi & EMAP,
  	Eigen::MatrixXi & EF,
  	Eigen::MatrixXi & EI,
    Eigen::MatrixXd & TC, // TODO: Texture coordinates
    Eigen::MatrixXi & FT, // TODO: Texture coordinates per face.
    EdgeMap & seam_edges, // TODO: A set of edges in V for vertices which lie on edges which should be preserved.
    MapV5d & Vmetrics, // TODO: The per-vertex data.
    int & a_e1,
    int & a_e2)
{
// Assign this to 0 rather than, say, -1 so that deleted elements will get
	// draw as degenerate elements at vertex 0 (which should always exist and
	// never get collapsed to anything else since it is the smallest index)
	using namespace Eigen;
	using namespace std;

	// Helper function to replace edge and associate information with NULL
	const auto & kill_edge = [&E,&EI,&EF](const int e)
	{
		E(e,0) = DUV_COLLAPSE_EDGE_NULL;
		E(e,1) = DUV_COLLAPSE_EDGE_NULL;
		EF(e,0) = DUV_COLLAPSE_EDGE_NULL;
		EF(e,1) = DUV_COLLAPSE_EDGE_NULL;
		EI(e,0) = DUV_COLLAPSE_EDGE_NULL;
		EI(e,1) = DUV_COLLAPSE_EDGE_NULL;
	};

	const int eflip = E(e,0)>E(e,1);
	// source and destination
	const int s = eflip?E(e,1):E(e,0);
	const int d = eflip?E(e,0):E(e,1);
	const bool collapse_on_seam = contains_edge( seam_edges, s, d );

	// If both endpoints are on seams, but there is no seam between them, reject it.
	if( seam_edges.count( s ) && seam_edges.count( d ) && !collapse_on_seam ) {
	    return false;
	}

	// Link condition
	if(!igl::edge_collapse_is_valid(e,F,E,EMAP,EF,EI) )
	{
		return false;
	}
	// Important to grab neighbors of d before monkeying with edges
	const std::vector<int> nV2Fd = igl::circulation(e,!eflip,EMAP,EF,EI);
	// We need the neighbors of s in case d is a seam "corner".
	const std::vector<int> nV2Fs = igl::circulation(e, eflip,EMAP,EF,EI);

	Bundle bundle = get_half_edge_bundle( e, E, EF, EI, F, FT );
    assert( bundle.size() == 2 );

	// 	sGet s_tc and d_tc by looking at F via EF/EI.
	int s_tc = bundle[0].p[0].tci;
	int d_tc = bundle[0].p[1].tci;
	// Match what eflip does.
	if( bundle[0].p[0].vi == d ) {
	    assert( bundle[0].p[1].vi == s );
		std::swap( s_tc, d_tc );
	} else {
	    assert( bundle[0].p[0].vi == s );
	    assert( bundle[0].p[1].vi == d );
	}

	// test fold-over
	const bool disable_boundary_foldover = true;
	if( disable_boundary_foldover && !collapse_on_seam ) {
		for(int i=1; i<nV2Fd.size()-1; i++) {
			const int f = nV2Fd[i];
			for(int v=0; v<3; v++) {
				if( F(f,v) == d ) {
					RowVectorXd uv = TC.row(FT(f,v));
					RowVectorXd uv1 = TC.row(FT(f,(v+1)%3));
					RowVectorXd uv2 = TC.row(FT(f,(v+2)%3));
					if( !two_points_on_same_side( uv1, uv2, uv, new_placement.tcs[0] ) /*&&
						 contains_edge( seam_edges, F(f,(v+1)%3), F(f,(v+2)%3) )*/) {
						 return false;
					}
				}
			}
		}
		for(int i=1; i<nV2Fs.size()-1; i++) {
			const int f = nV2Fs[i];
			for(int v=0; v<3; v++) {
				if( F(f,v) == s ) {
					RowVectorXd uv = TC.row(FT(f,v));
					RowVectorXd uv1 = TC.row(FT(f,(v+1)%3));
					RowVectorXd uv2 = TC.row(FT(f,(v+2)%3));
					if( !two_points_on_same_side( uv1, uv2, uv, new_placement.tcs[0] ) /*&&
						 contains_edge( seam_edges, F(f,(v+1)%3), F(f,(v+2)%3) )*/) {
						 return false;
					}
				}
			}
		}
	}

	// The following implementation strongly relies on s<d
	assert(s<d && "s should be less than d");
	// Preserve seams.
	// We can handle the case when (s, d) is a seam edge.
	if( collapse_on_seam )
	{
		assert( new_placement.tcs.size() == 2 );
		assert( new_placement.metrics.size() == 2 );
		// move source and destination to midpoint
		V.row(s) = new_placement.p;
		V.row(d) = new_placement.p;
		// Update UV coordinates of the edge endpoints here. If both edge endpoints are on the seam, return false (see above). If one of the edge endpoints is on the seam, we should be able to handle it preserving that endpoint and collapsing the other one.
		int he0_ts = bundle[0].p[0].tci;
		int he0_td = bundle[0].p[1].tci;
        if( bundle[0].p[0].vi == d ) 	std::swap( he0_ts, he0_td );
		TC.row(he0_ts) = new_placement.tcs[0];
		TC.row(he0_td) = new_placement.tcs[0];
		int he1_ts = bundle[1].p[0].tci;
		int he1_td = bundle[1].p[1].tci;
        if( bundle[1].p[0].vi == d ) 	std::swap( he1_ts, he1_td );
		TC.row(he1_ts) = new_placement.tcs[1];
		TC.row(he1_td) = new_placement.tcs[1];
		// Update the per-vertex metric.
        // Move the other d metrics to s.
		Vmetrics[d].erase(he0_td);
		Vmetrics[d].erase(he1_td);
		Vmetrics[s].insert(Vmetrics[d].begin(), Vmetrics[d].end());
		Vmetrics.erase(d);
		Vmetrics[s][he0_ts] = new_placement.metrics[0];
		Vmetrics[s][he1_ts] = new_placement.metrics[1];
	}
	else {
//      if(seam_edges.count( s ) || seam_edges.count( d ))    cout << "try to decimate edge attached to seam." << endl;
		assert( new_placement.tcs.size() == 1 );
		assert( new_placement.metrics.size() == 1 );
		// move source and destination to midpoint
		V.row(s) = new_placement.p;
		V.row(d) = new_placement.p;
		// Update UV coordinates of the edge endpoints here. If both edge endpoints are on the seam, return false (see above). If one of the edge endpoints is on the seam, we should be able to handle it preserving that endpoint and collapsing the other one.
		TC.row(s_tc) = new_placement.tcs[0];
		TC.row(d_tc) = new_placement.tcs[0];
		// Update the per-vertex metric.
		// Move the other d metrics to s.
		assert(bundle[0].p[0] == bundle[1].p[0] || bundle[0].p[0] == bundle[1].p[1]);
		assert(bundle[0].p[1] == bundle[1].p[0] || bundle[0].p[1] == bundle[1].p[1]);
		Vmetrics[d].erase(d_tc);
		Vmetrics[s].insert(Vmetrics[d].begin(), Vmetrics[d].end());
		Vmetrics.erase(d);
		Vmetrics[s][s_tc] = new_placement.metrics[0];
	}

	// finally, reindex faces and edges incident on d. Do this last so asserts
	// make sense.
	//
	// Could actually skip first and last, since those are always the two
	// collpased faces.
	const int m = F.rows();
	// vi, ti at s,d
	VertexBundle s_pair[2];
	VertexBundle d_pair[2];
	// update edge info
	// for each flap
	for(int side = 0;side<2;side++)
	{
		const int f = EF(e,side);
		const int v = EI(e,side);
		const int sign = (eflip==0?1:-1)*(1-2*side);
		// next edge emanating from d
		const int e1 = EMAP(f+m*((v+sign*1+3)%3));
		// prev edge pointing to s
		const int e2 = EMAP(f+m*((v+sign*2+3)%3));
		assert(E(e1,0) == d || E(e1,1) == d);
		assert(E(e2,0) == s || E(e2,1) == s);
		// We can't handle e1, which will be deleted, being a seam edge.
		// We also can't handle both endpoints of e1 being seam vertices,
		// even if there is no edge between them.
		// Add a test at the beginning which rejects collapsing (s,d) if
		//       all three vertices on either face flap are seam vertices.
		// If only e1 is a seam edge, meaning d is a seam vertex, but s is not,
		// we need to use our rename_vertex(d, s) function. We do it at the end.

		// face adjacent to f on e1, also incident on d
		const bool flip1 = EF(e1,1)==f;
		const int f1 = flip1 ? EF(e1,0) : EF(e1,1);
		assert(f1!=f);
		assert(F(f1,0)==d || F(f1,1)==d || F(f1,2) == d);
		// across from which vertex of f1 does e1 appear?
		const int v1 = flip1 ? EI(e1,0) : EI(e1,1);
		assert(F(f1,v1)!=s);
		// find vi and ti at s
		const int e_vi = bundle[side].p[0].vi == s ? 0 : 1;
		s_pair[side] = bundle[side].p[e_vi];
		d_pair[side] = bundle[side].p[1-e_vi];
		// Kill e1
		kill_edge(e1);
		// Kill f
		F(f,0) = DUV_COLLAPSE_EDGE_NULL;
		F(f,1) = DUV_COLLAPSE_EDGE_NULL;
		F(f,2) = DUV_COLLAPSE_EDGE_NULL;
		// Kill the texture coordinates FT face as well
		FT(f,0) = DUV_COLLAPSE_EDGE_NULL;
		FT(f,1) = DUV_COLLAPSE_EDGE_NULL;
		FT(f,2) = DUV_COLLAPSE_EDGE_NULL;
		// map f1's edge on e1 to e2
		assert(EMAP(f1+m*v1) == e1);
		EMAP(f1+m*v1) = e2;
		// side opposite f2, the face adjacent to f on e2, also incident on s
		const int opp2 = (EF(e2,0)==f?0:1);
		assert(EF(e2,opp2) == f);
		EF(e2,opp2) = f1;
		EI(e2,opp2) = v1;
		// remap e2 from d to s
		E(e2,0) = E(e2,0)==d ? s : E(e2,0);
		E(e2,1) = E(e2,1)==d ? s : E(e2,1);
		if(side==0)	a_e1 = e1;
		else		a_e2 = e1;
	}

	// Loop over face neighborhood of d.
	// We always collapse d to s, so these faces should be
	// updated to have their d vertex index point to s.
	for(int i=1; i<nV2Fd.size()-1; i++)
	{
		const int f = nV2Fd[i];
		for(int v=0; v<3; v++)
		{
			if(F(f,v) == d)
			{
				const int flip1 = (EF(EMAP(f+m*((v+1)%3)),0)==f)?1:0;
				const int flip2 = (EF(EMAP(f+m*((v+2)%3)),0)==f)?0:1;
				assert(E(EMAP(f+m*((v+1)%3)),flip1) == d
					|| E(EMAP(f+m*((v+1)%3)),flip1) == s);
				E(EMAP(f+m*((v+1)%3)),flip1) = s;
				assert(E(EMAP(f+m*((v+2)%3)),flip2) == d
					|| E(EMAP(f+m*((v+2)%3)),flip2) == s);
				E(EMAP(f+m*((v+2)%3)),flip2) = s;
				F(f,v) = s;
				// Update FT to point to the other endpoint.
				// If d is on the seam, some F neighbors won't be FT neighbors. Only update FT entries
				// that actually pointed to d_tc. (If "s" were the seam vertex, we wouldn't have to do this check.)
				// If d is a seam vertex, some face(s) incident at d
				// will not have texture coordinate d_tc.
				// Only update the ones that do.

                // Normally F and FTC have the same number of faces.
                // For meshes with boundaries, though, they won't.
                // The outside-of-the-mesh side of the boundary
                // edge will be connected with new faces to a vertex at infinity.
                // These new faces exist in F, but not FT.
                // Skip these new faces.
                // UPDATE: They now do have the same number of faces;
                //         FTC is also augmented with a vertex to infinity.
                // if( f < FT.rows() )
                {
                    assert( seam_edges.count(d) or FT(f,v) == d_tc );
                    if( !collapse_on_seam ) {
                        if( FT(f,v) == d_tc ) FT(f,v) = s_tc;
                    }
                    else {
                        if     (FT(f,v) == d_pair[0].tci )  FT(f,v) = s_pair[0].tci;
                        else if(FT(f,v) == d_pair[1].tci )	FT(f,v) = s_pair[1].tci;
                    }
                }
			}
		}
	}

	// Check for seam corners.
	const bool seam_corner = d_pair[0].tci == d_pair[1].tci;
	if( collapse_on_seam && seam_corner ) {
		for(int i=1; i<nV2Fs.size()-1; i++)
		{
			const int f = nV2Fs[i];
			for(int v=0; v<3; v++)
			{
				if(FT(f,v) == s_pair[1].tci) FT(f,v) = s_pair[0].tci;
			}
		}
	}

	// If 's' was a seam vertex and 'd' wasn't, we don't need to do anything. d was removed.
	// If 'd' was a seam vertex and 's' wasn't, we need to rename 'd' to 's'.
	if( seam_edges.count(d) && !seam_edges.count(s) ) rename_vertex( seam_edges, d, s );
	// If 's' and 'd' were both seam vertices but (s,d) is not a seam edge, we have a problem.
	assert( !( seam_edges.count(d) && seam_edges.count(s) && !contains_edge( seam_edges, d, s ) ) );
	// If (s,d) is a seam_edge, collapse it 'd' to 's'.
	if( contains_edge( seam_edges, d, s ) ) 		::collapse_edge( seam_edges, d, s );

	// Finally, "remove" this edge and its information
	kill_edge(e);

	return true;
}

bool collapse_edge_with_uv(
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & E,
    Eigen::VectorXi & EMAP,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    Eigen::MatrixXd & TC, //  Texture coordinates
    Eigen::MatrixXi & FT, //  Texture coordinates per face.
    EdgeMap & seam_edges, //  A set of indices into V or TC for vertices which lie on edges which should be preserved.
    MapV5d & Vmetrics, //  The per-vertex data.
    int seam_aware_degree,
    std::set<std::pair<double,int> > & Q,
    std::vector<std::set<std::pair<double,int> >::iterator > & Qit,
    std::vector< placement_info_5d > & C,
    int & e,
    bool test)
{
  	using namespace std;
  	using namespace Eigen;
  	using namespace igl;
	if(Q.empty())
	{
		// no edges to collapse
		return false;
	}
  	std::pair<double,int> p = *(Q.begin());
  	if(p.first == std::numeric_limits<double>::infinity())
  	{
    	// min cost edge is infinite cost
    	return false;
	}
	Q.erase(Q.begin());
	e = p.second;
	Qit[e] = Q.end();

	// Get the one-ring of faces as N.
	std::unordered_set<int> N;
	std::vector<int> Ne = circulation(e, true,EMAP,EF,EI);
	N.insert( Ne.begin(), Ne.end() );
	Ne = circulation(e,false,EMAP,EF,EI);
	N.insert( Ne.begin(), Ne.end() );

	int e1,e2;
	const bool collapsed =
    	try_collapse_5d_Edge(e,C.at(e),V,F,E,EMAP,EF,EI,TC,FT,seam_edges,Vmetrics,e1,e2);
	if(collapsed)
	{
		if(test == true ) {
    		cout << "try collapse succeed." << endl;
    	}
		// Erase the two, other collapsed edges
		Q.erase(Qit[e1]);
		Qit[e1] = Q.end();
		Q.erase(Qit[e2]);
		Qit[e2] = Q.end();
		// update local neighbors
		// loop over original face neighbors
		std::unordered_set< int > affected_edges;
		for(auto n : N)
		{
			if(F(n,0) != DUV_COLLAPSE_EDGE_NULL &&
			   F(n,1) != DUV_COLLAPSE_EDGE_NULL &&
			   F(n,2) != DUV_COLLAPSE_EDGE_NULL)
			{
				for(int v = 0;v<3;v++)
				{
					// get edge id
					const int ei = EMAP(v*F.rows()+n);
					// Because faces share edges, recomputing here would recompute
					// every affected edge twice. Instead, build a set< int > of edge
					// indices, and iterate over it once.
					affected_edges.insert( ei );
				}
			}
		}
		if(test == true ) {
    		cout << "first loop succeed." << endl;
    	}
		for( auto ei : affected_edges )
		{
			if( E(ei,0) != DUV_COLLAPSE_EDGE_NULL && E(ei,1) != DUV_COLLAPSE_EDGE_NULL )
			{
				// erase old entry
				Q.erase(Qit[ei]);
				Qit[ei] = Q.end();
				// compute cost and potential placement
				double cost;
				placement_info_5d place;
				Bundle b = get_half_edge_bundle( ei, E, EF, EI, F, FT );
				cost_and_placement_qslim5d_halfedge(b,V,F,TC,FT,seam_edges,Vmetrics,seam_aware_degree,cost,place);
				// Replace in queue
				Qit[ei] = Q.insert(std::pair<double,int>(cost,ei)).first;
				C.at(ei) = place;
			}
		}
		if(test == true ) {
    		cout << "second loop succeed." << endl;
    	}
	} else
	{
		if(test == true ) {
    		cout << "try collapse failed." << endl;
    	}
		// reinsert with infinite weight (the provided cost function must **not**
		// have given this un-collapsable edge inf cost already)
		p.first = std::numeric_limits<double>::infinity();
		Qit[e] = Q.insert(p).first;
	}
	return collapsed;
}
