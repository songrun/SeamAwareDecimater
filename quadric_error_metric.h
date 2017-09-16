#ifndef QUADRIC_ERROR_METRIC_H
#define QUADRIC_ERROR_METRIC_H

#include <Eigen/Core>
#include <vector>
#include "half_edge.h"

// Inputs:
//   V  #V by dim list of vertex positions, lesser index of E(e,:) will be set
//     to midpoint of edge.
//   F  #F by 3 list of face indices into V.
// Outputs:
//	 Q  array of metrics, one 4-by-4 metric per vertex  
void quadric_error_metric(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F, 
	std::vector< Eigen::MatrixXd >& Q);

void qslim_5d(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& TC, 
	const Eigen::MatrixXi& FT, 
	std::vector< Eigen::MatrixXd >& Q);	
	
void half_edge_qslim_5d(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& TC, 
	const Eigen::MatrixXi& FT, 
	MapV5d & hash_Q);	
#endif
