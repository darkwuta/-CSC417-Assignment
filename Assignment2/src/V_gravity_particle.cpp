#include <V_gravity_particle.h>

void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    //-------------------------------------
    V = mass * g.y() * q.y();
    //-------------------------------------
    //V = 0.0;

    //V = mass * q.transpose() * -g;
}