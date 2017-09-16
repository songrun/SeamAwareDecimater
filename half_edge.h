#ifndef HALF_EDGE_H
#define HALF_EDGE_H

#include <Eigen/Core>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <string>

struct VertexBundle
{
    int vi = -31337;
    int tci = -31337;
    
    VertexBundle() {};
    VertexBundle(int vi, int tci);
};
bool operator==( const VertexBundle& lhs, const VertexBundle& rhs );
bool operator!=( const VertexBundle& lhs, const VertexBundle& rhs );

struct HalfEdge
{
	int fi;		// face Id
	int ki;		// the index opposite to the half edge, ki = {0,1,2}
	VertexBundle p[2]; 	// endpoints, each endpoint contain vi and ti
	
	HalfEdge(int fi, int ki);
};

// A map from vertex_index to texcoord_index to metrics.
typedef std::unordered_map< int, std::unordered_map< int, Eigen::MatrixXd > > MapV5d;
typedef std::vector<HalfEdge> Bundle;

Bundle get_half_edge_bundle(
    int e,
    Eigen::MatrixXi & E,
    Eigen::MatrixXi & EF,
    Eigen::MatrixXi & EI,
    Eigen::MatrixXi & F,
    Eigen::MatrixXi & FT
    );

void print_bundle( const Bundle & b );


typedef std::unordered_map< int, std::unordered_set< int > > EdgeMap;
// Modifies `edges` in-place to collapse the edge between the two vertices.
void collapse_edge( EdgeMap& edges, int vertex_to_remove, int vertex_collapsing_into );
// Modifies `edges` in-place to remove a vertex and replace it with another one.
// The new_vertex_name should not already be present.
void rename_vertex( EdgeMap& edges, int old_vertex_name, int new_vertex_name );
// Modifies `edges` in-place to add an edge between the two vertices.
void insert_edge( EdgeMap& edges, int v1, int v2 );
// Returns whether the edge v1,v2 is in edges.
bool contains_edge( const EdgeMap& edges, int v1, int v2 );

#endif
