#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    P.resize(q_size - indices.size() * 3, q_size);
    P.setZero();
    int count = 0;
    for (int i = 0; i < q_size - 3 * indices.size(); i++) {
        if (count < indices.size() && i == 3.0 * indices[count] - 3.0 * count) { count++; }
        P.insert(i, i + 3.0 * count) = 1.0;
    }

}