#include "detect_foldover.h"
#include <vector>

// judge if p1, p2 on the same side of the straight line pass through uv1, uv2
bool two_points_on_same_side(
	const Eigen::RowVectorXd uv1,
	const Eigen::RowVectorXd uv2,
	const Eigen::RowVectorXd p1,
	const Eigen::RowVectorXd p2)
{
	assert( uv1.size() == 2 );
	assert( uv2.size() == 2 );
	assert( p1.size() == 2 );
	assert( p2.size() == 2 );
	
	if( uv1 == uv2 )	return true;
	if( uv1(0) == uv2(0) ) {
		if( (uv1(0) - p1(0))*(uv1(0) - p2(0)) <= 0 ) {
			return false;
		}	
	}
	else {
		double k = (uv2(1)-uv1(1))/(uv2(0)-uv1(0));
		double b = uv1(1)-uv1(0)*k;
		if( (p1(0)*k+b-p1(1)) * (p2(0)*k+b-p2(1)) <= 0 ) {
			return false;
		}
	}
	
	return true;
}

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
)
{
	using namespace std;
	using namespace Eigen;
	
	const int m = F.rows();
	const auto & step = [&](
		const int e,
		const int ff,
		const int dir, 
		int & ne,
		int & nf, 
		int & ev)
		{
			assert( dir == 1 || dir == -1 );
			assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
			const int nside = EF(e,0)==ff?1:0;
			const int nv = EI(e,nside);
			// get next face
			nf = EF(e,nside);
			// get next edge 
			ne = EMAP(nf+m*((nv+dir+3)%3));
			ev = (nv-dir+3)%3;
		};	

	int ei = e;
	assert( E(e,0) == vertex_off_seam || E(e,1) == vertex_off_seam );
	const int nside = 0;
	const int f0 = EF(e,nside);
	int fi = f0;
	int fvi; 
	int dir = -1;
	{	// test dir	
		const int nv = EI(e,1-nside);
		// get next face
		const int nf = EF(e,1-nside);
		// get next edge 
		const int ne = EMAP(nf+m*((nv+dir+3)%3));
		if( E(ne,0) != vertex_off_seam && E(ne,1) != vertex_off_seam )	dir = -dir;
	}
	std::vector<std::pair<int,int>> boundary;		// vertices on boundary of vertex_off_seam's neighborhood, except for the other vertex of e.
	while(true)
	{
		step(ei,fi,dir,ei,fi,fvi);
		boundary.push_back(std::make_pair(fi,fvi));
		// back to start?
		if(f0 == fi)
		{
			assert( ei == e );
			break;
		}	
	}
	// if vertex_off_seam and vertex_on_seam are on different sides for any boundary edge 
	// in UV space, return false
	assert( boundary.size() > 2 );
	for( int i=1; i<boundary.size()-1; i++ ) {
		RowVectorXd uv1 = TC.row(FT(boundary[i].first,(boundary[i].second+1)%3));
		RowVectorXd uv2 = TC.row(FT(boundary[i].first,(boundary[i].second+2)%3));
		if( !two_points_on_same_side(uv1, uv2, TC.row(ti), TC.row(tj)) )
			return false;
	}
	
	return true;
}