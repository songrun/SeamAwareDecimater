#ifndef NEIPHBOR_FACES_AND_BOUNDARY_H
#define NEIPHBOR_FACES_AND_BOUNDARY_H

#include <Eigen/Core>
#include <vector>
#include <utility>

// Altered based on igl/circulation.cpp
// find all the neighboring faces and boundary edges given a collapsed edge.
void neighbor_faces_and_boundary (
	const int e,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXi & E,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI,
	std::vector<int> & neigh_faces,
	std::vector<std::pair<int,int>> & boundary
);

#endif