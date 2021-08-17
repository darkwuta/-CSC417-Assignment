#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> X) {

    Eigen::Vector3d xbar = X - p;
    Eigen::Matrix3d cross_xbar;

    cross_xbar<< 0, -xbar.z(), xbar.y(),
        xbar.z(), 0, -xbar.x(),
        -xbar.y(), xbar.x(), 0;

    Eigen::Matrix36d mid;

    mid << cross_xbar.transpose(), Eigen::Matrix3d::Identity();

    Eigen::Matrix66d later;
    later.setZero();
    later.block<3, 3>(0, 0) = R.transpose();
    later.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();

    J = R * mid * later;

}

