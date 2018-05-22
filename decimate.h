#ifndef DECIMATE_SEAM_H
#define DECIMATE_SEAM_H
#include <igl/igl_inline.h>
#include <igl/point_mesh_squared_distance.h>
#include <Eigen/Core>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include "half_edge.h"

struct placement_info_5d {
	Eigen::RowVectorXd p;
	std::vector<Eigen::RowVectorXd> tcs;
	std::vector<Eigen::MatrixXd>    metrics;
};

typedef std::set<std::pair<double,int> > PriorityQueue;

  // Assumes (V,F) is a manifold mesh (possibly with boundary) Collapses edges
  // until desired number of faces is achieved. This uses default edge cost and
  // merged vertex placement functions {edge length, edge midpoint}.
  //
  // Inputs:
  //   V  #V by dim list of vertex positions.
  //   F  #F by 3 list of face indices into V.
  //   TC #TC by dim list of texture coordinates.
  //   FT #F by 3 list of face indices into TC.
  //   max_m  desired number of output faces
  // Outputs:
  //   V_out  #V_out by dim list of output vertex positions (can be same ref as V)
  //   F_out  #F_out by 3 list of output face indices into V_out (can be same ref as F)
  //   TC_out #TC_out by dim list of output texture coordinates (can be same ref as TC)
  //   FT_out #F_out by 3 list of output face indices into TC_out (can be same ref as FT)
  //   J      #F_out list of indices into F of birth face
  // Returns true if m was reached (otherwise #F_out > m)
  //
  // Assumes a **closed** manifold mesh. See igl::connect_boundary_to_infinity
  // and igl::decimate in decimate.cpp
  // is handling meshes with boundary by connecting all boundary edges with
  // dummy facets to infinity **and** modifying the stopping criteria.
  //
  // Inputs:
  //   cost_and_placement  function computing cost of collapsing an edge and 3d
  //     position where it should be placed:
  //     cost_and_placement(V,F,E,EMAP,EF,EI,cost,placement);
  //   stopping_condition  function returning whether to stop collapsing edges
  //     based on current state. Guaranteed to be called after _successfully_
  //     collapsing edge e removing edges (e,e1,e2) and faces (f1,f2):
  //     bool should_stop =
  //       stopping_condition(V,F,E,EMAP,EF,EI,Q,Qit,C,e,e1,e2,f1,f2);

bool decimate_halfedge_5d(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & TC,
    const Eigen::MatrixXi & FT,
    EdgeMap & seam_edges,
    MapV5d & Vmetrics,
    int target_num_vertices,
    const int seam_aware_degree,
    Eigen::MatrixXd & V_out,
    Eigen::MatrixXi & F_out,
    Eigen::MatrixXd & TC_out,
    Eigen::MatrixXi & FT_out
    );
    
void clean_mesh(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & TC,
	const Eigen::MatrixXi & FT,
	const int nF,
	Eigen::MatrixXd & V_out,
	Eigen::MatrixXi & F_out,
	Eigen::MatrixXd & TC_out,
	Eigen::MatrixXi & FT_out);
	
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
	std::vector< placement_info_5d > & C);
	
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
	int & prev_e);

static std::unordered_set<int> interior_foldovers;
static std::unordered_set<int> exterior_foldovers;

#endif

