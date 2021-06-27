#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    f.resize(q.rows());
    //std::cout << "ASSEMBLE_FORCES::DEBUG::f.rows():" << f.rows() << std::endl;
    //std::cout << "ASSEMBLE_FORCES::DEBUG::q.rows():" << q.rows() << std::endl;
    for (int y = 0; y < E.rows(); y++)
    {
        Eigen::Vector3d ft(0, 0, 0);
        dV_gravity_particle_dq(ft, mass, Eigen::Vector3d(0, 0, 0));
        int i = E(y, 0);
        int j = E(y, 1);
        Eigen::Vector3d q0, q1;
        q0 << q(i * 3);
        q1 << q(j * 3);
        dV_spring_particle_particle_dq(ft, q0, q1, l0(y), k);
        f<< -ft.x(),-ft.y(),-ft.z();
    }
    //std::cout << "ASSEMBLE_FORCES::DEBUG::q.rows():" << q.rows() << std::endl;
};