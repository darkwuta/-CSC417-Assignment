#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::tmp_force.rows():" << tmp_force.rows() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::tmp_stiffness.rows():" << tmp_stiffness.rows() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::q" << q << std::endl;
    //tmp_stiffness = Eigen::SparseMatrixd(q.rows(), q.rows());
    Eigen::SparseMatrixd A;
    Eigen::SparseMatrixd K;
    Eigen::VectorXd b,x;

    force(tmp_force, q, qdot);
    stiffness(tmp_stiffness, q, qdot);

    A = (mass - dt * dt * tmp_stiffness);
    b = mass * qdot + dt * tmp_force;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if(solver.info()!=Eigen::Success){
        std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::solver.compute(A)::ERROR" << std::endl;
        return;
    }
    x = solver.solve(b);
    qdot = x;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::qdot.rows():" << qdot.rows() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::qdot.cols():" << qdot.cols() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::mass.rows():" << mass.rows() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::mass.cols():" << mass.cols() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::tmp_force.rows():" << tmp_force.rows() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::tmp_force.cols():" << tmp_force.cols() << std::endl;
    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::qdot:" << qdot << std::endl;
    //if(solver.info()!=Eigen::Success){
    //    std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::M*qdot+dt*tmp_force::ERROR" << std::endl;
    //}
    //int rows = q.rows();
    //outerSizeÃ»ÎÊÌâ
    //for (int k = 0; k < tmp_stiffness.outerSize(); ++k)
    //{
    //    for (Eigen::SparseMatrix<double>::InnerIterator it(tmp_stiffness, k); it; ++it)
    //    {
    //        std::cout << "FIXED_POINT_CONSTRAINTS::it.value:" << it.value() << std::endl;
    //        std::cout << "FIXED_POINT_CONSTRAINTS::it.row():" << it.row() << std::endl;
    //        std::cout << "FIXED_POINT_CONSTRAINTS::it.col():" << it.col() << std::endl;
    //        std::cout << "FIXED_POINT_CONSTRAINTS::it.index():" << it.index() << std::endl;
    //    }
    //}
    //
    //for (unsigned int i = 0; i < rows; i++)
    //{
    //    //std::cout << "ASSEMBLE_STIFFNESS::DEBUG::K.coeffRef(5,5):" << tmp_stiffness.coeffRef(5, 5) << std::endl;
    //    //Eigen::SparseMatrix<double>::InnerIterator tm(mass, k);
    //    //Eigen::SparseMatrix<double>::InnerIterator tk(tmp_stiffness, k);

    //    if (q(i) == 0)
    //        continue;
    //    double m = 1;
    //    double k = 1000000;
    //    double f = tmp_force(i);
    //    double tqdot = qdot(i);

    //    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::f:" << f << std::endl;
    //    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::m:" << tm.value() << std::endl;
    //    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::k:" << k << std::endl;
    //    
    //    //std::cout << "LINEARLY_IMPLICIT_EULER::DEBUG::tqdot:" << tqdot << std::endl;

    //    qdot(i) = (m * tqdot + dt * f) / (m - dt * dt * k);
    //}
    //qdot = (mass * qdot + 0.0001 * tmp_force);
    q = q + dt * qdot;

    //q += 0.0000001*q;
}