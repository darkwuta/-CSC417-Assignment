#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    ////------------------------------------------------------------------
    V = 0.5 * stiffness * pow(((q0 - q1).norm() - l0), 2);
    ////------------------------------------------------------------------
    //V = 0.0;

    //V = (stiffness * pow(l0 - sqrt(pow(q0[0] - q1[0], 2.0) + pow(q0[1] - q1[1], 2.0) + pow(q0[2] - q1[2], 2.0)), 2.0)) / 2.0;
}