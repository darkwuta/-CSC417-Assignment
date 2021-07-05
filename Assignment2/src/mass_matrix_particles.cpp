#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

    M.resize(q.size(), q.size());
    typedef Eigen::Triplet<double> T;//��i��j��massΪ����������
    std::vector<T> tripletList;
    M.setZero();
    for (int i = 0; i < q.rows(); i++)
        tripletList.push_back(T(i, i, mass));
    M.setFromTriplets(tripletList.begin(), tripletList.end());
}
