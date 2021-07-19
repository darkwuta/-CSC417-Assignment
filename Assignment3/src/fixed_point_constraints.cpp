#include <fixed_point_constraints.h>
#include <algorithm>
#include <set>
#include <iostream>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    ////--------------------------------------------
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripleList;

    P.resize(q_size - indices.size() * 3, q_size);
    P.setZero();
    int count = 0;
    //TODO 不知道为什么下面的不行
    //for (int i = 0; i < q_size - 3 * indices.size(); i++) {
    //    if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
    //    tripleList.push_back(T(i, i + 3.0 * count, 1.));
    //    //P.insert(i, i + 3.0 * count) = 1.0;
    //}
    for (int i = 0; i < q_size/3; i++)
    {
        if (indices[count] == i)count++;
        else
        {
            tripleList.push_back(T((i - count) * 3, i*3, 1.));
            tripleList.push_back(T((i - count) * 3+1, i * 3+1, 1.));
            tripleList.push_back(T((i - count) * 3+2, i * 3+2, 1.));
        }
    }
    P.setFromTriplets(tripleList.begin(), tripleList.end());
    //--------------------------------------------
    //int l = indices.size();
    //typedef Eigen::Triplet<double> T;
    //std::vector<T> tripleList;
    ////tripleList.reserve(q_size * q_size);

    //int j = 0;
    //for (int i = 0; i < q_size; i += 3) {

    //    if (j < indices.size() && indices[j] == i / 3)
    //    {
    //        j += 1;
    //    }
    //    else
    //    {
    //        tripleList.push_back(T(i - 3 * j + 0, i + 0, 1.));
    //        tripleList.push_back(T(i - 3 * j + 1, i + 1, 1.));
    //        tripleList.push_back(T(i - 3 * j + 2, i + 2, 1.));
    //    }

    //}
    //P.resize(q_size - (3 * l), q_size);
    //P.setFromTriplets(tripleList.begin(), tripleList.end());
}