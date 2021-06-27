#include <fixed_point_constraints.h>
#include <algorithm>
#include"iostream"
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG::P.cols:" << P.cols() << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:P.rows:" << P.rows() << std::endl;
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	P.setZero();
	for (unsigned int i = 0; i < q_size; i++)
		tripletList.push_back(T(i, i, 1));
	for (unsigned int k = 0; k < indices.size(); k++)
		tripletList.push_back(T(indices[k], indices[k], 0));
	
	//std::cout << "FIXED_POINT_CONSTRAINTS::q_size:" << q_size << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::indices.size():" << indices.size() << std::endl;
	P.setFromTriplets(tripletList.begin(), tripletList.end());
}