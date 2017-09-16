#include "neighbor_faces_and_boundary.h"
#include <unordered_set>

void neighbor_faces_and_boundary (
	const int e,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXi & E,
	const Eigen::VectorXi & EMAP,
	const Eigen::MatrixXi & EF,
	const Eigen::MatrixXi & EI,
	std::vector<int> & neigh_faces,
	std::vector<std::pair<int,int>> & boundary
) {
	const int m = F.rows();
	const auto & step = [&](
		const int e, 
		const int ff,
		const bool ccw,
		int & ne, 
		int & nf,
		int & ev)
		{
			assert((EF(e,1) == ff || EF(e,0) == ff) && "e should touch ff");
			//const int fside = EF(e,1)==ff?1:0;
			const int nside = EF(e,0)==ff?1:0;
			const int nv = EI(e,nside);
			// get next face
			nf = EF(e,nside);
			// get next edge 
			const int dir = ccw?-1:1;
			ne = EMAP(nf+m*((nv+dir+3)%3));
			ev = (nv-dir+3)%3;
		};
	// Always start with first face (ccw in step will be sure to turn right
  	// direction)
  	std::unordered_set<int> nfaces;
  	std::vector<std::pair<int,int>> b1, b2;
	const int f0 = EF(e,0);
	int fi = f0;
	int ei = e;
	int fvi;
	while(true)
	{
		step(ei,fi,true,ei,fi,fvi);
		nfaces.insert(fi);
		b1.push_back(std::make_pair(fi,fvi));
		// back to start?
		if(fi == f0)
		{
			assert(ei == e);
			break;
		}
	}
	while(true)
	{
		step(ei,fi,false,ei,fi,fvi);
		nfaces.insert(fi);
		b2.push_back(std::make_pair(fi,fvi));
		// back to start?
		if(fi == f0)
		{
			assert(ei == e);
			break;
		}
	}
	
	neigh_faces.resize(nfaces.size());
	std::copy(nfaces.begin(), nfaces.end(), neigh_faces.begin());
	
	const int k=b1.size(), n=b2.size();
	assert( k > 2 && n > 2 );
	boundary.resize(k+n-4);
	std::copy(b1.begin()+1, b1.begin()+k-1, boundary.begin());
	std::copy(b2.begin()+1, b2.begin()+n-1, boundary.begin()+k-2);
}