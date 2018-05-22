#ifndef COLLAPSE_EDGE_SEAM_H
#define COLLAPSE_EDGE_SEAM_H
#include <Eigen/Core>
#include <vector>
#include <set>
#include <unordered_set>
#include <utility> // std::swap
#include "decimate.h"

// Assumes (V,F) is a closed manifold mesh (except for previouslly collapsed
// faces which should be set to: 
// [DUV_COLLAPSE_EDGE_NULL DUV_COLLAPSE_EDGE_NULL DUV_COLLAPSE_EDGE_NULL].
// Collapses exactly two faces and exactly 3 edges from E (e and one side of
// each face gets collapsed to the other). This is implemented in a way that
// it can be repeatedly called until satisfaction and then the garbage in F
// can be collected by removing NULL faces.
//
// Inputs:
//   e  index into E of edge to try to collapse. E(e,:) = [s d] or [d s] so
//     that s<d, then d is collapsed to s.
///  p  dim list of vertex position where to place merged vertex
// Inputs/Outputs:
//   V  #V by dim list of vertex positions, lesser index of E(e,:) will be set
//     to midpoint of edge.
//   F  #F by 3 list of face indices into V.
//   E  #E by 2 list of edge indices into V.
//   EMAP #F*3 list of indices into E, mapping each directed edge to unique
//     unique edge in E
//   EF  #E by 2 list of edge flaps, EF(e,0)=f means e=(i-->j) is the edge of
//     F(f,:) opposite the vth corner, where EI(e,0)=v. Similarly EF(e,1) "
//     e=(j->i)
//   EI  #E by 2 list of edge flap corners (see above).
//   e1  index into E of edge collpased on left
//   e2  index into E of edge collpased on left
//   f1  index into E of edge collpased on left
//   f2  index into E of edge collpased on left
// Returns true if edge was collapsed
#define DUV_COLLAPSE_EDGE_NULL -1

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
    int & a_e2);
        
bool collapse_edge_with_uv(
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
    int seam_aware_degree,
    std::set<std::pair<double,int> > & Q,
    std::vector<std::set<std::pair<double,int> >::iterator > & Qit,
    std::vector< placement_info_5d > & C,
    int & e,
    bool test = false);


#endif
