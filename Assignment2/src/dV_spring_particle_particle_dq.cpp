#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    Eigen::Vector6d q;
    q << q0 - q1, q1 - q0;
    double l = (q0 - q1).norm();
    f = stiffness * (l - l0) * q / l;
}