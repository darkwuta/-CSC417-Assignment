#include <fixed_point_constraints.h>
#include <algorithm>
#include"iostream"
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG::P.cols:" << P.cols() << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:P.rows:" << P.rows() << std::endl;
	//P.resize(q_size, q_size);
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:m:" << indices.size() << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:n:" << q_size/3 << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:m*n:" << q_size/3*indices.size() << std::endl;
	//std::cout << "FIXED_POINT_CONSTRAINTS::DEBUG:P.size():" << P.size()/9 << std::endl;

	//typedef Eigen::Triplet<double> T;
	//std::vector<T> tripletList;
	////P.setZero();
	////for (unsigned int i = 0; i < q_size/3; i++)
	////{
	////	tripletList.push_back(T(i*3, i*3, 1));
	////	tripletList.push_back(T(i*3+1, i*3+1, 1));
	////	tripletList.push_back(T(i*3+2, i*3+2, 1));
	////}
	//for (unsigned int i = 0; i < q_size/3; i++)
	//{
	//	bool isFixed = 0;
	//	for (unsigned int k = 0; k < indices.size(); k++)
	//	{
	//		if (i == indices[k])
	//		{
	//			tripletList.push_back(T(i * 3, i * 3, 0));
	//			tripletList.push_back(T(i * 3 +1, i * 3 +1, 0));
	//			tripletList.push_back(T(i * 3 +2, i * 3 +2, 0));
	//			isFixed = 1;
	//			break;
	//		}
	//	}
	//	if (!isFixed)
	//	{
	//		tripletList.push_back(T(i * 3, i * 3, 1));
	//		tripletList.push_back(T(i * 3 +1, i * 3 +1, 1));
	//		tripletList.push_back(T(i * 3 +2, i * 3 +2, 1));
	//	}
	//		
	//}	
	//std::cout << "FIXED_POINT_CONSTRAINTS::q_size:" << q_size << std::endl;
	////std::cout << "FIXED_POINT_CONSTRAINTS::indices.size():" << indices.size() << std::endl;
 //   
	//P.setFromTriplets(tripletList.begin(), tripletList.end());
	////for (int k = 0; k < P.outerSize(); ++k)
	////	for (Eigen::SparseMatrix<double>::InnerIterator it(P, k); it; ++it)
	////	{
	////		if (it.value() == 0)
	////		{
	////			std::cout << "FIXED_POINT_CONSTRAINTS::it.value:" << it.value() << std::endl;
	////			std::cout << "FIXED_POINT_CONSTRAINTS::it.row():" << it.row() << std::endl;
	////			std::cout << "FIXED_POINT_CONSTRAINTS::it.col():" << it.col() << std::endl;
	////			std::cout << "FIXED_POINT_CONSTRAINTS::it.index():" << it.index() << std::endl;
	////		}
	////	}
	P.resize(q_size - 3 * indices.size(), q_size);
	P.setZero();
	int count = 0;


	for (int i = 0; i < q_size - 3 * indices.size(); i++) {
		if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
		P.insert(i, i + 3.0 * count) = 1.0;
	}
}