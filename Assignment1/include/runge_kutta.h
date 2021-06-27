//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    //force(mass,q,qdot);
    //kappa1
    Eigen::VectorXd k1q=q;
    Eigen::VectorXd k1qdot=qdot;

    Eigen::VectorXd k2qdot=qdot-dt*100.0*q;
    Eigen::VectorXd k2q=qdot * dt;


    qdot =qdot-dt*50.0*(k1q+k2q);
    q = q+dt*0.5*(k1qdot+k2qdot);

    
}