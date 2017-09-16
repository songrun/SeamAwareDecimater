#include "half_edge.h"
#include <Eigen/Geometry>
#include <iostream>
#include <functional>

VertexBundle::VertexBundle(int new_vi, int new_tci) : vi( new_vi ), tci( new_tci ) {}
bool operator==( const VertexBundle& lhs, const VertexBundle& rhs )
{
    return lhs.vi == rhs.vi && lhs.tci == rhs.tci;
}
bool operator!=( const VertexBundle& lhs, const VertexBundle& rhs )
{
    return !( lhs == rhs );
}

HalfEdge::HalfEdge(int a, int b) : fi(a), ki(b) {}

Bundle get_half_edge_bundle(
    int e,
    Eigen::MatrixXi & E,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & FTC
    )
{
    Bundle result;
    result.reserve(2);
    
    for( int side = 0; side < 2; ++side ) {
        
        const int face_index = EF( e, side );
        const int opposite_vertex = EI( e, side );
        
        // Normally F and FTC have the same number of faces.
        // For meshes with boundaries, though, they won't.
        // The outside-of-the-mesh side of the boundary
        // edge will be connected with new faces to a vertex at infinity.
        // These new faces exist in F, but not FT.
        // In case someone accesses one of these edges to infinity,
        // return a bogus texture coordinate (-1).
        // UPDATE: They now do have the same number of faces;
        //         FTC is also augmented with a vertex to infinity.
        
        HalfEdge he( face_index, opposite_vertex );
        int v1 = F  ( face_index, (opposite_vertex + 1)%3);
        int t1 = FTC( face_index, (opposite_vertex + 1)%3);
        /*
        int t1 = face_index < FTC.rows()
               ? FTC( face_index, (opposite_vertex + 1)%3)
               : -1
               ;
        */
        int v2 = F  ( face_index, (opposite_vertex + 2)%3);
        int t2 = FTC( face_index, (opposite_vertex + 2)%3);
        /*
        int t2 = face_index < FTC.rows()
               ? FTC( face_index, (opposite_vertex + 2)%3)
               : -1
               ;
        */
        he.p[0] = VertexBundle(v1, t1);
        he.p[1] = VertexBundle(v2, t2);
        
        result.push_back( he );
    }
    
    return result;
}

/*
BundleKey bundle_key(
	const int e,
	const Eigen::MatrixXi & E)
{
	assert( E.cols() == 2 );
	assert( e < E.rows() );
	int v1 = E(e,0);
	int v2 = E(e,1);
	if( v1 > v2 ) 	std::swap(v1,v2);
	return std::make_pair(v1,v2);
}
*/

void print_bundle( const Bundle & b )
{
	using namespace std;
	for(int i=0; i<b.size(); i++) {
		cout << "Half Edge #: " << i << endl;
		cout << "vi: " << b[i].p[0].vi << " ti: " << b[i].p[0].tci << endl;
		cout << "vi: " << b[i].p[1].vi << " ti: " << b[i].p[1].tci << endl;
		cout << "fi: " << b[i].fi << endl;
		cout << "ki: " << b[i].ki << endl;
	}
}



bool contains_edge( const EdgeMap& edges, int v1, int v2 )
{
    assert( v1 != v2 );
    assert( v1 >= 0 );
    assert( v2 >= 0 );
    
    const bool result1 = edges.count( v1 ) && edges.at( v1 ).count( v2 );
    const bool result2 = edges.count( v2 ) && edges.at( v2 ).count( v1 );
    assert( result1 == result2 );
    return result1;
}

void collapse_edge( EdgeMap& edges, int vertex_to_remove, int vertex_collapsing_into )
{
    assert( contains_edge( edges, vertex_to_remove, vertex_collapsing_into ) );
    
    // This should be true from the way we use it.
    assert( vertex_collapsing_into < vertex_to_remove );
    
    // For every neighbor of `vertex_to_remove`, replace it with `vertex_collapsing_into`.
    const std::unordered_set< int > neighbors( edges[ vertex_to_remove ] );
    for( const auto& n : neighbors ) {
        edges[ n ].erase( vertex_to_remove );
        edges[ n ].insert( vertex_collapsing_into );
    }
    
    // Add every neighbor of `vertex_to_remove` to `vertex_collapsing_into`.
    edges[ vertex_collapsing_into ].insert( neighbors.begin(), neighbors.end() );
    
    // We have now added `vertex_collapsing_into` as a neighbor of itself. Remove it.
    assert( edges[ vertex_collapsing_into ].count( vertex_collapsing_into ) );
    edges[ vertex_collapsing_into ].erase( vertex_collapsing_into );
    
    // Finally, erase `vertex_to_remove` from `edges`.
    edges.erase( vertex_to_remove );
}
void rename_vertex( EdgeMap& edges, int old_vertex_name, int new_vertex_name )
{
    // We must already know about old_vertex_name.
    assert( edges.count( old_vertex_name ) );
    // We must not already know about new_vertex_name.
    assert( !edges.count( new_vertex_name ) );
    
    // For every neighbor of `old_vertex_name`, replace it with `new_vertex_name`.
    const std::unordered_set< int > neighbors( edges[ old_vertex_name ] );
    for( const auto& n : neighbors ) {
        edges[ n ].erase( old_vertex_name );
        edges[ n ].insert( new_vertex_name );
    }
    
    // Move the neighbors of `old_vertex_name` to `new_vertex_name`.
    edges[ new_vertex_name ] = neighbors;
    edges.erase( old_vertex_name );
}
void insert_edge( EdgeMap& edges, int v1, int v2 )
{
    assert( !contains_edge( edges, v1, v2 ) );
    edges[ v1 ].insert( v2 );
    edges[ v2 ].insert( v1 );
}
