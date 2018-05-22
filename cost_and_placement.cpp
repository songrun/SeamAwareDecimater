#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <cmath>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <igl/active_set.h>
#include <igl/edge_collapse_is_valid.h>
#include <igl/circulation.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/writeDMAT.h>
#include <igl/barycentric_coordinates.h>
#include <igl/point_simplex_squared_distance.h>
#include <chrono>

#include "cost_and_placement.h"
#include "eiquadprog.h"
#include "neighbor_faces_and_boundary.h"

namespace {
	const double eps = 1e-8;
	enum SOLVER_TYPE {
		IGL_SOLVER = 0,
		EIQUADPROG
	} solver = EIQUADPROG;
}
	
const double DINF = std::numeric_limits<double>::infinity();
// #define DEBUG_MODE

// get plane from three points, and express it as ax + by + cz + d = 0. 
// return coefficients a,b,c,d as result
Eigen::Vector4d face_from_three_points (
	const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3)
{
	Eigen::Vector3d n = (v2-v1).cross(v3-v1);
	n = n/n.norm();
	double d = -n.dot(v1);
	
	Eigen::Vector4d res;
	res << n(0), n(1), n(2), d;
	
	assert( res != Eigen::Vector4d::Zero() && "Three points are colinear." );
	return res;
}

void cost_and_placement_qslim5d_halfedge (
	const Bundle & e,
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & TC,
	const Eigen::MatrixXi & FT,
	const EdgeMap & seam_edges,
	const MapV5d  & Vmetrics,
	const int seam_aware_degree,
	double &                cost,
	placement_info_5d &     new_placement
	)
{

	// Assign this to 0 rather than, say, -1 so that deleted elements will get
	// draw as degenerate elements at vertex 0 (which should always exist and
	// never get collapsed to anything else since it is the smallest index)
	using namespace Eigen;
	using namespace std;

	// If one of the endpoints is the special vertex at infinity, don't touch it.
	assert( e.size() == 2 );		//  each edge has two half-edge
	const bool has_infinity_vertex = V.row( V.rows()-1 ).minCoeff() == DINF;
	if( has_infinity_vertex 
	  	&& (/* F(e[0].fi, e[0].ki) == V.rows()-1 
	   		|| F(e[1].fi, e[1].ki) == V.rows()-1
	   		||*/ e[0].p[0].vi == V.rows()-1
	   		|| e[0].p[1].vi == V.rows()-1) ) {
		cost = DINF;
		return;
	}
		
#ifdef DEBUG_MODE
	cout << "No infinite vertex." << endl;
#endif
	
	solver = EIQUADPROG;
	MatrixXd new_metric;
	
	// two vertex indices on one side of e
	const int vi[2] = {e[0].p[0].vi, e[0].p[1].vi};
	// If vi[0] and vi[1] are in seam_edges, but (vi[0], vi[1]) is not, return infinite cost.
	if( seam_edges.count( vi[0] ) && seam_edges.count( vi[1] ) && !contains_edge( seam_edges, vi[0], vi[1] ) ) {
	    cost = DINF;
	    return;
	}
	
	/// case 1: 
	// (vi[0], vi[1]) in seam_edges, compute each half edge
	if( contains_edge( seam_edges, vi[0], vi[1] ) ) {
		VertexBundle e_p0[2];		// two Vertex5d for both sides at one end
		VertexBundle e_p1[2];		// two Vertex5d for both sides at the other end
		MatrixXd m[2];				// two metrics
		for(int side=0; side<2; side++) {
			e_p0[side] = e[side].p[0];
			e_p1[side] = e[side].p[1];
			m[side] = Vmetrics.at(e_p0[side].vi).at(e_p0[side].tci) 
					+ Vmetrics.at(e_p1[side].vi).at(e_p1[side].tci);
		}
		
		// TODO: Check that the seam edges are collinear with
		//       their neighbors along the seam. If they are not, return DINF.
		const auto & is_collinear = [&TC](const int tci_1, const int tci_2, const int tci_3)
		{
			assert( tci_1 >= 0 && tci_1 <TC.rows() );
			assert( tci_2 >= 0 && tci_2 <TC.rows() );
			assert( tci_3 >= 0 && tci_3 <TC.rows() );
			Eigen::RowVector2d uv1 = TC.row(tci_1);
			Eigen::RowVector2d uv2 = TC.row(tci_2);
			Eigen::RowVector2d uv3 = TC.row(tci_3);
			Eigen::RowVector2d n1 = (uv2-uv1)/(uv2-uv1).norm();
			Eigen::RowVector2d n2 = (uv3-uv1)/(uv3-uv1).norm();
			return 1-fabs(n1.dot(n2)) < eps;
		};
		const auto & edge_ratio = [&TC](const int tci_1, const int tci_2, const int tci_3)
		{
			assert( tci_1 >= 0 && tci_1 <TC.rows() );
			assert( tci_2 >= 0 && tci_2 <TC.rows() );
			assert( tci_3 >= 0 && tci_3 <TC.rows() );
			Eigen::RowVector2d uv1 = TC.row(tci_1);
			Eigen::RowVector2d uv2 = TC.row(tci_2);
			Eigen::RowVector2d uv3 = TC.row(tci_3);
			if( (uv3-uv2).norm() == 0 ) 	return DINF;
			return (uv2-uv1).norm()/(uv3-uv2).norm();
		};
		
		bool is_free[2] = {false, false};	// Correspond to vi[2] 	
		// for each end of one side, search its neighbor seam edges to check if anyone is collinear with e
		for(int end=0; end<2; end++) {
			// An end is free only if it has exactly two neighboring seam edges
			if(seam_edges.at(vi[end]).size() != 2)	continue;
			for(auto vj : seam_edges.at(vi[end])) {		// all the neighboring seam vertices.
				// test if exist one vertex which has two uvs collinear with both sides of e's uvs 
				if( vj == vi[1-end] )	continue;
                double ratio[2] = {DINF, DINF};
				for(auto item : Vmetrics.at(vj)) {		// all the tci for one neighboring seam vertex.
					int tcj = item.first;
					if(is_collinear(tcj, e_p0[0].tci, e_p1[0].tci)) {
						ratio[0] = edge_ratio(tcj, e_p0[0].tci, e_p1[0].tci);
					}
					if(is_collinear(tcj, e_p0[1].tci, e_p1[1].tci)) {
						ratio[1] = edge_ratio(tcj, e_p1[1].tci, e_p0[1].tci);
					}
				}
				switch (seam_aware_degree) {
					case 0: is_free[end] = true; 
						break;
					case 1: if (ratio[0]!=DINF && ratio[1]!=DINF) is_free[end] = true;
						break;
					case 2: if (ratio[0]!=DINF && ratio[1]!=DINF && abs(ratio[0]-ratio[1]) <= 1e-3) is_free[end] = true;
						break;
				}
			}
		}

		// neither end has collinear neighbor seam edge, shouldn't touch it.
		if( !is_free[0] && !is_free[1] ) {
			cost = DINF;
			return;
		}
		// Otherwise, if any end is not free, collapse e to it.
		for(int end=0; end<2; end++) {
			if( !is_free[end] ) {
				cost = 0;
				new_placement.tcs.resize(2);
				new_placement.metrics.resize(2);
				for(int side=0; side<2; side++) {
					RowVectorXd v(6);
					v.setOnes();
					assert(e[side].p[end].vi == vi[end] || e[side].p[end].vi == vi[1-end]);
					v.head(3) = V.row(vi[end]);
					v.segment(3,2) = e[side].p[end].vi == vi[end] ? TC.row(e[side].p[end].tci) : TC.row(e[side].p[1-end].tci);
					cost += v*m[side]*v.transpose();
					new_placement.p = v.head(3);
					new_placement.tcs[side] = v.segment(3,2);
					new_placement.metrics[side] = m[side];
				}
				return;
			}
		}
		// Finally, if both ends are free. Collapse the edge to a point on the seam.
		// And the uv on both sides should be proportional.
		VectorXd Z;
		// The unknowns are x,y,z,u0,v0,u1,v1,1
		MatrixXd G(8,8);
		G.setZero();
		// combine both sides' metric of e
		if( solver == EIQUADPROG ) {
			// build new metric
			G.block(0,0,3,3) = m[0].block(0,0,3,3) + m[1].block(0,0,3,3);
			G.block(3,0,2,3) = m[0].block(3,0,2,3);
			G.block(0,3,3,2) = m[0].block(0,3,3,2);
			G.block(3,3,2,2) = m[0].block(3,3,2,2);
			G.block(5,0,2,3) = m[1].block(3,0,2,3);
			G.block(0,5,3,2) = m[1].block(0,3,3,2);
			G.block(5,5,2,2) = m[1].block(3,3,2,2);
			RowVectorXd b(7);	// b is the linear term in combined metric
			b.segment(0,3) = m[0].block(5,0,1,3) + m[1].block(5,0,1,3);
			b.segment(3,2) = m[0].block(5,3,1,2);
			b.segment(5,2) = m[1].block(5,3,1,2);
			G.block(7,0,1,7) = b;
			G.block(0,7,7,1) = b.transpose();
			G(7,7) = m[0](5,5) + m[1](5,5);
			
			// regularizer
			const double w = 1e-6;
			MatrixXd reg(8,8);
			reg.setIdentity();
			G.block(0,0,8,8) = G.block(0,0,8,8) + w*reg;
			VectorXd g0(8);
			g0.setOnes();
			g0.segment(0,3) = (V.row(e_p0[0].vi)+V.row(e_p1[0].vi))/2;
			g0.segment(3,2) = (TC.row(e_p0[0].tci)+TC.row(e_p1[0].tci))/2;
			g0.segment(5,2) = (TC.row(e_p0[1].tci)+TC.row(e_p1[1].tci))/2;
			g0 = -w*g0;
			
			// Add the constraint that uv0 and uv1 stay on the same uv-space line
			// and the new position should have the same parameter along each edge.
			// equality constraints:
			// x(7) - 1 = 0
			// (x(3) - e0_u0) - t(e0_u1 - e0_u0) = 0
			// (x(4) - e0_v0) - t(e0_v1 - e0_v0) = 0
			// (x(5) - e1_u0) - t(e1_u1 - e1_u0) = 0
			// (x(6) - e1_v0) - t(e1_v1 - e1_v0) = 0
			RowVector2d vec[2] = {TC.row(e_p1[0].tci) - TC.row(e_p0[0].tci),
								  TC.row(e_p1[1].tci) - TC.row(e_p0[1].tci)};
			assert( vec[0].norm() != 0 && vec[1].norm() != 0 );
			MatrixXd CE(8,4);
			VectorXd ce0(4);
			CE.setZero();
			CE(7,0) = 1.0;
			ce0(0) = -1;
			if( vec[0](0) != 0 ) {		 // t = (x(3) - e0_u0)/(e0_u1 - e0_u0)
				CE(3,1) = -vec[0](1);
				CE(4,1) = vec[0](0);
				ce0(1) = vec[0](1)*TC.row(e_p0[0].tci)(0) - vec[0](0)*TC.row(e_p0[0].tci)(1);
				CE(3,2) = -vec[1](0);
				CE(5,2) = vec[0](0);
				ce0(2) = vec[1](0)*TC.row(e_p0[0].tci)(0) - vec[0](0)*TC.row(e_p0[1].tci)(0);
				CE(3,3) = -vec[1](1);
				CE(6,3) = vec[0](0);
				ce0(3) = vec[1](1)*TC.row(e_p0[0].tci)(0) - vec[0](0)*TC.row(e_p0[1].tci)(1);
			}
			else {						// t = (x(4) - e0_v0)/(e0_v1 - e0_v0)
				assert( vec[0](1) != 0 );
				CE(4,1) = -vec[0](0);
				CE(3,1) = vec[0](1);
				ce0(1) = vec[0](0)*TC.row(e_p0[0].tci)(1) - vec[0](1)*TC.row(e_p0[0].tci)(0);
				CE(4,2) = -vec[1](0);
				CE(5,2) = vec[0](1);
				ce0(2) = vec[1](0)*TC.row(e_p0[0].tci)(1) - vec[0](1)*TC.row(e_p0[1].tci)(0);
				CE(4,3) = -vec[1](1);
				CE(6,3) = vec[0](1);
				ce0(3) = vec[1](1)*TC.row(e_p0[0].tci)(1) - vec[0](1)*TC.row(e_p0[1].tci)(1);
			}
			
			// inequality constraints:
			// t >= 0 && t <= 1
			MatrixXd CI(8,2);
			VectorXd ci0(2);
			CI.setZero();
			if( vec[0](0) != 0 ) {
				double sign = vec[0](0) > 0 ? 1 : -1;
				CI(3,0) = sign;
				ci0(0) = -sign*TC.row(e_p0[0].tci)(0);
				CI(3,1) = -sign;
				ci0(1) = sign*(TC.row(e_p0[0].tci)(0)+vec[0](0));
			}
			else {
				double sign = vec[0](1) > 0 ? 1 : -1;
				CI(4,0) = sign;
				ci0(0) = -sign*TC.row(e_p0[0].tci)(1);
				CI(4,1) = -sign;
				ci0(1) = sign*(TC.row(e_p0[0].tci)(1)+vec[0](1));
			}
			
			solve_quadprog(G,g0,CE,ce0,CI,ci0,Z);
		} else {
			assert( false && "Unknown solver type" );
		}
			// set UV coordinates to be the middle point
		assert( isfinite(Z(0)) && isfinite(Z(1)) && isfinite(Z(2)) 
		     && isfinite(Z(3)) && isfinite(Z(4)) && isfinite(Z(5)) && isfinite(Z(6)));
		new_placement.p = Z.head(3);
		new_placement.tcs = {Z.segment(3,2), Z.segment(5,2)};
		new_placement.metrics = {m[0], m[1]};
		RowVectorXd v = Z;
		// Multiply by one half because we added two energy terms, shall we?
		cost = v*G*v.transpose();
				
		return;
	}

	/// case 2:
	// new metric is the summation of the collapsed vertices' metrics
	assert( e[0].p[0] == e[1].p[1] && e[0].p[1] == e[1].p[0] );
	const int tci[2] = {e[0].p[0].tci, e[0].p[1].tci};
	new_metric = Vmetrics.at(vi[0]).at(tci[0]) + Vmetrics.at(vi[1]).at(tci[1]);
	assert( new_metric.transpose() == new_metric ); // all the metrics are symmetric
	
	// If one vertex is on a seam, it will stay fixed. Use Q as the cost.
	for(int end=0; end<2; end++) {
		if(seam_edges.count(vi[end]) && !seam_edges.count(vi[1-end])) {
			RowVectorXd v(6);
			v.setOnes();
			v.head(3) = V.row(vi[end]);
			v.segment(3,2) = TC.row(tci[end]);
			cost = v*new_metric*v.transpose();
			new_placement.p = V.row(vi[end]);
			new_placement.tcs = { TC.row(tci[end]) };
			new_placement.metrics = { new_metric };
			return;
		} 
	}
	
	/// case 3: no attachment to seam	
	// solve 
	VectorXd Z;
	if( solver == EIQUADPROG ) {
		const double w = 1e-6;
		MatrixXd G = new_metric;
		// regularizer
		MatrixXd reg(6,6);
		reg.setIdentity();
		G = G + w*reg;
		MatrixXd CE(6,1);
		CE.setZero();
		CE(5,0) = 1.0;
		VectorXd g0(6), ce0(1);
		g0.setOnes();
		g0.segment(0,3) = (V.row(vi[0])+V.row(vi[1]))/2;
		g0.segment(3,2) = (TC.row(tci[0])+TC.row(tci[1]))/2;
		g0 = -w*g0;
		ce0 << -1;
		MatrixXd CI;
		VectorXd ci0;
		solve_quadprog(G,g0,CE,ce0,CI,ci0,Z);
	} else {
		assert( false && "Unknown solver type" );
	}

	// set UV coordinates to be the middle point
	assert( isfinite(Z(0)) && isfinite(Z(1)) && isfinite(Z(2)) && isfinite(Z(3)) && isfinite(Z(4)));
	new_placement.p = Z.head(3);
	new_placement.tcs = {Z.segment(3,2)};
	new_placement.metrics = {new_metric};
	RowVectorXd v = Z;
	cost = v*new_metric*v.transpose();

}