#include <fixed_point_constraints.h>
#include <algorithm>
#include <set>
#include <iostream>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    ////--------------------------------------------
    //P.resize(q_size - indices.size() * 3, q_size);
    //P.setZero();
    //int count = 0;
    //for (int i = 0; i < q_size - 3 * indices.size(); i++) {
    //    if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
    //    P.insert(i, i + 3.0 * count) = 1.0;
    //}
    ////--------------------------------------------
    P.resize(q_size - 3 * indices.size(), q_size);
    std::vector<Eigen::Triplet<double>> tripletList;
    std::set<int> fixed;
    // std::cout << " ---- " << std::endl;
    for (int i = 0; i < indices.size(); i++) {
        // std::cout << "fixed index = " << indices[i] << std::endl;
        fixed.insert(indices[i]);
    }
    // row pointer of P
    int cnt_NonFixed = 0;
    int n = q_size / 3;
    for (int col = 0; col < n; col++) {
        // index of vertex that is NOT fixed
        if (fixed.count(col) < 1) {
            for (int j = 0; j < 3; j++) {
                tripletList.push_back({ 3 * cnt_NonFixed + j, 3 * col + j, 1 });
            }
            cnt_NonFixed = cnt_NonFixed + 1;
        }
    }
    P.setFromTriplets(tripletList.begin(), tripletList.end());

    // check if cnt_NonFixed == m
    if (cnt_NonFixed != P.rows() / 3) {
        std::cout << "ERROR construct P: num of non fixed != m" << std::endl;
    }
}