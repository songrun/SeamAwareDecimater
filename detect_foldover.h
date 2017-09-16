#ifndef DETECT_FOLDOVER_H
#define DETECT_FOLDOVER_H

#include <Eigen/Core>

bool two_points_on_same_side(
	const Eigen::RowVectorXd uv1,
	const Eigen::RowVectorXd uv2,
	const Eigen::RowVectorXd p1,
	const Eigen::RowVectorXd p2);
	
bool try_attach_to_seam(
	const int e, 
	const int vertex_off_seam,
	const int ti,
	const int tj,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F, 	
	const Eigen::MatrixXi & E,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI,
	const Eigen::MatrixXd & TC,
	const Eigen::MatrixXi & FT
);
	
#endif