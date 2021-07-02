#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
    //std::cout << "ASSEMBLE_FORCES::DEBUG::f.rows():" << f.rows() << std::endl;
    //std::cout << "ASSEMBLE_FORCES::DEBUG::q.rows():" << q.rows() << std::endl;
    f.resize(q.rows());
    f.setZero();
    for (int y = 0; y < E.rows(); y++)
    {
        Eigen::Vector6d ft;
        //ft << 0, 0, 0, 0, 0, 0;
        //ft << 0, -9.6, 0, 0, -9.8, 0;
        //dV_gravity_particle_dq(ft, mass, Eigen::Vector3d(0, 0, 0));
        int i = E(y, 0);
        int j = E(y, 1);
        Eigen::Vector3d q0, q1;
        q0 << q(i * 3), q(i * 3 + 1) , q(i * 3 + 2);
        q1 << q(j * 3), q(j * 3 + 1) , q(j * 3 + 2);
        dV_spring_particle_particle_dq(ft, q0, q1, l0(y), k);
        //f<< -ft.x(),-ft.y(),-ft.z();
        f(i*3)+=ft(0);
        f(i*3+1)+=ft(1);
        f(i*3+2)+=ft(2);
        f(j*3)+=ft(3);
        f(j*3+1)+=ft(4);
        f(j*3+2)+=ft(5);
    }
    //std::cout << "ASSEMBLE_FORCES::DEBUG::q.rows():" << q.rows() << std::endl;
};