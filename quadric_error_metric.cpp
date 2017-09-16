#include "quadric_error_metric.h"
#include <Eigen/Geometry>
#include <iostream>

namespace {
	const double eps = 1e-7;
	const double uv_weight = 1e6;
}
void quadric_error_metric(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F, 
	std::vector< Eigen::MatrixXd >& Q)
{
	using namespace std;
	using namespace Eigen;
	
	// array of 4x4 matrices
	Q.resize( V.rows() );
	for( int i=0; i<Q.size(); i++ )
		Q[i] = Matrix4d::Zero();
	
	const auto & face_from_three_points = [](const Vector3d& v1, const Vector3d& v2, const Vector3d& v3)
	{
		Vector3d n = (v2-v1).cross(v3-v1);
		n.normalize();
		double d = -n.dot(v1);
		
		Vector4d res;
		res << n(0), n(1), n(2), d;
		
		return res;
	};
	
	// the metric at each vertex equals to the sum of metric of its attached faces
	for( int i=0; i<F.rows(); i++ ) {
		Vector3d v1 = V.row( F(i,0) );
		Vector3d v2 = V.row( F(i,1) );
		Vector3d v3 = V.row( F(i,2) );
		Vector4d p = face_from_three_points(v1, v2, v3);
		Matrix4d metric = p*p.transpose();
		
		Q[ F(i,0) ] += metric;
		Q[ F(i,1) ] += metric;
		Q[ F(i,2) ] += metric;
	}
	
	// the cost v.T*Q*v should equal to zero
	for( int i=0; i<V.rows(); i++ ) {
		// cout << fabs(v * Q[i] * v.transpose()) << endl;
		assert( fabs(V.row(i).homogeneous() * Q[i] * V.row(i).homogeneous().transpose()) <= eps );
	}
}

void qslim_5d(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& TC, 
	const Eigen::MatrixXi& FT, 
	std::vector< Eigen::MatrixXd >& Q)
{
	using namespace std;
	using namespace Eigen;
	
	assert( F.rows() == FT.rows() );
	const int nF = F.rows();
	// array of 6x6 matrices
	Q.resize( V.rows() );
	for( int i=0; i<Q.size(); i++ ) {
		MatrixXd metric(6,6);
		metric.setZero();
		Q[i] = metric;
	}
	
	for(int i=0; i<nF; i++) {
		VectorXd p1(5),p2(5),p3(5);
		p1.head(3) = V.row( F(i,0) );
		p2.head(3) = V.row( F(i,1) );
		p3.head(3) = V.row( F(i,2) );
		p1.tail(2) = TC.row( FT(i,0) );
		p2.tail(2) = TC.row( FT(i,1) );
		p3.tail(2) = TC.row( FT(i,2) );
		// Paper Section 5.1
		VectorXd e1 = (p2-p1)/(p2-p1).norm();
		VectorXd e2 = p3-p1-(e1.dot(p3-p1))*e1;
		e2 /= e2.norm();
		const double eps = 1e-7;
		assert( fabs(e1.norm() - 1) <= eps );
		assert( fabs(e2.norm() - 1) <= eps );
		
		MatrixXd A(5,5);
		A.setIdentity();
		A = A - e1*e1.transpose() - e2*e2.transpose();
		VectorXd b = p1.dot(e1)*e1 + p1.dot(e2)*e2 - p1;
		double c = p1.dot(p1) - p1.dot(e1)*p1.dot(e1) - p1.dot(e2)*p1.dot(e2);
		
		// Paper Section 3.4
		MatrixXd metric(6,6);
		metric.block(0,0,5,5) = A;
		metric.block(0,5,5,1) = b;
		metric.block(5,0,1,5) = b.transpose();
		metric(5,5) = c;
		
		// add metric to each vertex
		Q[ F(i,0) ] += metric;
		Q[ F(i,1) ] += metric;
		Q[ F(i,2) ] += metric;
	}
}
	
void half_edge_qslim_5d(
	const Eigen::MatrixXd& V, 
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& TC, 
	const Eigen::MatrixXi& FT, 
	MapV5d & hash_Q)
{
	using namespace std;
	using namespace Eigen;
	
	// initialize 5d vertex map, key is (vi,ti), value is zero metric
	assert( F.rows() == FT.rows() );
	const int nF = F.rows();
	for(int i=0; i<nF; i++) {
	
		/// A. compute metric for each face
		VectorXd p1(5),p2(5),p3(5);
		p1.head(3) = V.row( F(i,0) );
		p2.head(3) = V.row( F(i,1) );
		p3.head(3) = V.row( F(i,2) );
		p1.tail(2) = TC.row( FT(i,0) );
		p2.tail(2) = TC.row( FT(i,1) );
		p3.tail(2) = TC.row( FT(i,2) );
		// Paper Section 5.1
		VectorXd e1 = (p2-p1)/(p2-p1).norm();
		VectorXd e2 = p3-p1-(e1.dot(p3-p1))*e1;
		e2 /= e2.norm();
		const double eps = 1e-7;
		assert( fabs(e1.norm() - 1) <= eps );
		assert( fabs(e2.norm() - 1) <= eps );
		
		MatrixXd A(5,5);
		A.setIdentity();
		A = A - e1*e1.transpose() - e2*e2.transpose();
		VectorXd b = p1.dot(e1)*e1 + p1.dot(e2)*e2 - p1;
		double c = p1.dot(p1) - p1.dot(e1)*p1.dot(e1) - p1.dot(e2)*p1.dot(e2);
		
		// Paper Section 3.4
		MatrixXd metric(6,6);
		metric.block(0,0,5,5) = A;
		metric.block(0,5,5,1) = b;
		metric.block(5,0,1,5) = b.transpose();
		metric(5,5) = c;	
	
		/// B. assign the face metric to each 5d vertex, if it hasn't appeared, initialize
		/// it with the metric, otherwise, add the metric to its original metric. 
		for(int j=0; j<3; j++) {
			int vi = F(i,j);
			int ti = FT(i,j);
			if( hash_Q[vi].count(ti) == 0 ) {
				hash_Q[vi][ti] = metric;
			} 
			else {
				hash_Q[vi][ti] += metric;
			}
		}
	}
	
}	
