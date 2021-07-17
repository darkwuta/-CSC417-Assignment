#include <Eigen/Dense>
#include <EigenTypes.h>
#include <newtons_method.h>

//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  energy(q, qdot) -  a function that computes the energy of the FEM system. This takes q and qdot as parameters, returns the energy value.
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_qdot - scratch space for storing velocities 
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename ENERGY, typename FORCE, typename STIFFNESS> 
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  ENERGY &energy, FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_qdot, Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    
    auto f_energy_i = [&](Eigen::VectorXd qdot_i) {

        return 0.5 * (qdot_i - qdot).transpose() * mass * (qdot_i - qdot) + energy(qdot_i);
        // return energy(qdot_i);

    };

    auto f_grad_i = [&](Eigen::VectorXd& g_i, Eigen::VectorXd qdot_i) {

        Eigen::VectorXd tmp_f;
        force(tmp_f, q + dt * qdot_i, qdot_i);

        g_i = mass * (qdot_i - qdot) - dt * tmp_f;
    };

    auto f_Hessian_i = [&](Eigen::SparseMatrixd& H_i, Eigen::VectorXd qdot_i) {

        Eigen::SparseMatrixd tmp_K;
        stiffness(tmp_K, q + dt * qdot_i, qdot_i);

        H_i = mass - pow(dt, 2) * tmp_K;
    };

    Eigen::VectorXd qdot_i = qdot;
    unsigned int max_step = 5;
    newtons_method(qdot_i, f_energy_i, f_grad_i, f_Hessian_i, max_step, tmp_force, tmp_stiffness);

    // update q and qdot at next time stamp
    q = q + dt * qdot_i;
    qdot = qdot_i;
}
