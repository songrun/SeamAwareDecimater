#ifndef COST_AND_PLACEMENT_H
#define COST_AND_PLACEMENT_H

#include <Eigen/Core>
#include <unordered_set>
#include <set>
#include "decimate.h"	
  
void cost_and_placement_qslim5d_halfedge (
	const Bundle & e,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & TC,
	const Eigen::MatrixXi & FT,
	const EdgeMap & seam_edges,
	const MapV5d & Vmetrics,
	const int seam_aware_degree,
	double &                cost,
	placement_info_5d &     new_placement
	);

#endif
