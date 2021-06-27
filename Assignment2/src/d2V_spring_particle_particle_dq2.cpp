#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::Ref<Eigen::Matrix66d> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    double l = (q0 - q1).norm();
    Eigen::Matrix66d BtB;
    BtB << Eigen::Matrix3d::Identity(), -Eigen::Matrix3d::Identity(),
        -Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Identity();
    Eigen::Vector6d q;
    q << q0, q1;
    H = stiffness * BtB - stiffness * l0 * ((BtB * l - BtB * q * (BtB * q).transpose()) / (l * l));
    
}