#include <T_cloth.h>

void T_cloth(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::SparseMatrixd &M) {
    double h = 1.0;
    T = 0.5 * h * qdot.transpose() * M * qdot;
}