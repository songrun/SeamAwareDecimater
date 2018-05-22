#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"
// #include "igl/decimate.h"

#include <Eigen/Core>

#include "pythonlike.h"

#include <cstdlib> // exit()
#include <iostream>
#include <cassert>
#include <cstdio> // printf()

#include <igl/seam_edges.h>
#include "decimate.h"
#include "quadric_error_metric.h"
#include <igl/writeDMAT.h>

// An anonymous namespace. This hides these symbols from other modules.
namespace {

void usage( const char* argv0 )
{
    std::cerr << "Usage: " << argv0 << " <path/to/input.obj> num-vertices     <target_number_of_vertices>  [--strict] [<strictness>]" << std::endl;
    std::cerr << "Usage: " << argv0 << " <path/to/input.obj> percent-vertices <target_percent_of_vertices> [--strict] [<strictness>]" << std::endl;
    exit(-1);
}

int count_seam_edge_num(const EdgeMap& seam_vertex_edges)
{
	int count = 0;
	for(auto & v : seam_vertex_edges) {
		count += v.second.size();
	}
	return count / 2;
}

enum SeamAwareDegree
{
	NoUVShapePreserving,
	UVShapePreserving,
	Seamless
};

/*
Decimates a triangle mesh down to a target number of vertices,
preserving the UV parameterization.
TODO Q: Should we do a version that does not preserve the UV parameterization exactly,
        but instead returns a sequence of TC/FTC that can be used to transform a UV
        point between the parameterizations of the decimated and undecimated mesh?
Input parameters:
    V: The 3D positions of the input mesh (3 columns)
    TC: The 2D texture coordinates of the input mesh (2 columns)
    F: Indices into `V` for the three vertices of each triangle.
    FTC: Indices into `TC` for the three vertices of each triangle.
Output parameters:
    Vout: The 3D positions of the decimated mesh (3 columns),
          where #vertices is as close as possible to `target_num_vertices`)
    TCout: The texture coordinates of the decimated mesh (2 columns)
    Fout: Indices into `Vout` for the three vertices of each triangle.
    FTCout: Indices into `TCout` for the three vertices of each triangle.
Returns:
    True if the routine succeeded, false if an error occurred.
Notes:
    The output mesh will a vertex count as close as possible to `target_num_vertices`.
    The decimated mesh should never have fewer vertices than `target_num_vertices`.
*/
template <typename DerivedV, typename DerivedF, typename DerivedT>
bool decimate_down_to(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const Eigen::PlainObjectBase<DerivedT>& TC,
    const Eigen::PlainObjectBase<DerivedF>& FT,
    int target_num_vertices,
    Eigen::MatrixXd& V_out,
    Eigen::MatrixXi& F_out,
    Eigen::MatrixXd& TC_out,
    Eigen::MatrixXi& FT_out,
    int seam_aware_degree
    )
{
#define DEBUG_DECIMATE_DOWN_TO
    assert( target_num_vertices > 0 );
    assert( target_num_vertices < V.rows() );
    
    /// 3D triangle mesh with UVs.
    // 3D
    assert( V.cols() == 3 );
    // triangle mesh
    assert( F.cols() == 3 );
    // UVs
    assert( TC.cols() == 2 );
    assert( FT.cols() == 3 );
    assert( FT.cols() == F.cols() );
    
    // Print information about seams.
    Eigen::MatrixXi seams, boundaries, foldovers;
    igl::seam_edges( V, TC, F, FT, seams, boundaries, foldovers );
#ifdef DEBUG_DECIMATE_DOWN_TO
    std::cout << "seams: " << seams.rows() << "\n";
    std::cout << seams << std::endl;
    std::cout << "boundaries: " << boundaries.rows() << "\n";
    std::cout << boundaries << std::endl;
    std::cout << "foldovers: " << foldovers.rows() << "\n";
    std::cout << foldovers << std::endl;
#endif
    
    // Collect all vertex indices involved in seams.
    std::unordered_set< int > seam_vertex_indices;
    // Also collect the edges in terms of position vertex indices themselves.
    EdgeMap seam_vertex_edges;
    {
		for( int i = 0; i < seams.rows(); ++i ) {
		    const int v1 = F( seams( i, 0 ),   seams( i, 1 ) );
		    const int v2 = F( seams( i, 0 ), ( seams( i, 1 ) + 1 ) % 3 );
			seam_vertex_indices.insert( v1 );
			seam_vertex_indices.insert( v2 );
			insert_edge( seam_vertex_edges, v1, v2 );
			// The vertices on both sides should match:
			assert( seam_vertex_indices.count( F( seams( i, 2 ),   seams( i, 3 ) ) ) );
			assert( seam_vertex_indices.count( F( seams( i, 2 ), ( seams( i, 3 ) + 1 ) % 3 ) ) );
		}
		for( int i = 0; i < boundaries.rows(); ++i ) {
		    const int v1 = F( boundaries( i, 0 ),   boundaries( i, 1 ) );
		    const int v2 = F( boundaries( i, 0 ), ( boundaries( i, 1 ) + 1 ) % 3 );
			seam_vertex_indices.insert( v1 );
			seam_vertex_indices.insert( v2 );
			insert_edge( seam_vertex_edges, v1, v2 );
		}
		for( int i = 0; i < foldovers.rows(); ++i ) {
		    const int v1 = F( foldovers( i, 0 ),   foldovers( i, 1 ) );
		    const int v2 = F( foldovers( i, 0 ), ( foldovers( i, 1 ) + 1 ) % 3 );
			seam_vertex_indices.insert( v1 );
			seam_vertex_indices.insert( v2 );
			insert_edge( seam_vertex_edges, v1, v2 );
			// The vertices on both sides should match:
			assert( seam_vertex_indices.count( F( foldovers( i, 2 ),   foldovers( i, 3 ) ) ) );
			assert( seam_vertex_indices.count( F( foldovers( i, 2 ), ( foldovers( i, 3 ) + 1 ) % 3 ) ) );
		}
	
	    std::cout << "# seam vertices: " << seam_vertex_indices.size() << std::endl;		
		std::cout << "# seam edges: " << count_seam_edge_num(seam_vertex_edges) << std::endl;
    }
  
    // Compute the per-vertex quadric error metric.
    std::vector< Eigen::MatrixXd > Q;
    bool success = false;
    Eigen::VectorXi J;
    
	MapV5d hash_Q;
	half_edge_qslim_5d(V,F,TC,FT,hash_Q);
	std::cout << "computing initial metrics finished\n" << std::endl;
	success = decimate_halfedge_5d(
		V, F,
		TC, FT,
		seam_vertex_edges,
		hash_Q,
		target_num_vertices,
		seam_aware_degree,
		V_out, F_out,
		TC_out, FT_out
		);
	std::cout << "#seams after decimation: " << count_seam_edge_num(seam_vertex_edges) << std::endl;
    std::cout << "#interior foldeover: " << interior_foldovers.size() << std::endl;
    std::cout << "#exterior foldeover: " << exterior_foldovers.size() << std::endl;
    return success;
}
}

int main( int argc, char* argv[] ) {
    std::vector<std::string> args( argv + 1, argv + argc );
    std::string strictness;
    int seam_aware_degree = int( SeamAwareDegree::Seamless );
    const bool found_strictness = pythonlike::get_optional_parameter( args, "--strict", strictness );	
	if ( found_strictness ) {
		seam_aware_degree = atoi(strictness.c_str());
	}
    
    if( args.size() != 3 && args.size() != 4 )	usage( argv[0] );
    std::string input_path, command, command_parameter;
    pythonlike::unpack( args.begin(), input_path, command, command_parameter );
    args.erase( args.begin(), args.begin() + 3 );
    
    // Does the input path exist?
    Eigen::MatrixXd V, TC, CN;
    Eigen::MatrixXi F, FT, FN;
    if( !igl::readOBJ( input_path, V, TC, CN, F, FT, FN ) ) {
        std::cerr << "ERROR: Could not read OBJ: " << input_path << std::endl;
        usage( argv[0] );
    }

    std::cout << "Loaded a mesh with " << V.rows() << " vertices and " << F.rows() << " faces: " << input_path << std::endl;
    
    // Get the target number of vertices.
    int target_num_vertices = 0;
    if( command == "num-vertices" ) {
        // strto<> returns 0 upon failure, which is fine, since that is invalid input for us.
        target_num_vertices = pythonlike::strto< int >( command_parameter );
    }
    else if( command == "percent-vertices" ) {
        const double percent = pythonlike::strto< double >( command_parameter );
        target_num_vertices = lround( ( percent * V.rows() )/100. );
        std::cout << command_parameter << "% of " << std::to_string( V.rows() ) << " input vertices is " << std::to_string( target_num_vertices ) << " output vertices." << std::endl;
        // Ugh, printf() requires me to specify the types of integers versus longs.
        // printf( "%.2f%% of %d input vertices is %d output vertices.", percent, V.rows(), target_num_vertices );
    }
    else {
        std::cerr << "ERROR: Unknown command: " << command << std::endl;
        usage( argv[0] );
    }
    
    // Check that the target number of vertices is positive and fewer than the input number of vertices.
    if( target_num_vertices <= 0 ) {
        std::cerr << "ERROR: Target number of vertices must be a positive integer: " << argv[4] << std::endl;
        usage( argv[0] );
    }
    if( target_num_vertices >= V.rows() ) {
    	std::string output_path = pythonlike::os_path_splitext( input_path ).first + "-decimated_to_" + std::to_string( V.rows() ) + "_vertices.obj";
        if( !igl::writeOBJ( output_path, V, F, CN, FN, TC, FT ) ) {
			std::cerr << "ERROR: Could not write OBJ: " << output_path << std::endl;
			usage( argv[0] );
		}
   		std::cout << "Wrote: " << output_path << std::endl;
        std::cerr << "ERROR: Target number of vertices must be smaller than the input number of vertices: " << argv[4] << std::endl;
        return 0;
    }
    
    // Make the default output path.
    std::string output_path = pythonlike::os_path_splitext( input_path ).first + "-decimated_to_" + std::to_string( target_num_vertices ) + "_vertices.obj";
    if( !args.empty() ) {
        output_path = args.front();
        args.erase( args.begin() );
    }
    
    // We should have consumed all arguments.
    if( !args.empty() ) usage( argv[0] );
    
    // Decimate!
    Eigen::MatrixXd V_out, TC_out, CN_out;
    Eigen::MatrixXi F_out, FT_out, FN_out;
    const bool success = decimate_down_to( V, F, TC, FT, target_num_vertices, V_out, F_out, TC_out, FT_out, seam_aware_degree );
    if( !success ) {
        std::cerr << "WARNING: decimate_down_to() returned false (target number of vertices may have been unachievable)." << std::endl;
    }
    
    if( !igl::writeOBJ( output_path, V_out, F_out, CN_out, FN_out, TC_out, FT_out ) ) {
        std::cerr << "ERROR: Could not write OBJ: " << output_path << std::endl;
        usage( argv[0] );
    }
    std::cout << "Wrote: " << output_path << std::endl;
    
    return 0;
}
