#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    ////-----------------------------------------
    //Eigen::Vector6d q;
    //q << q0[0] - q1[0],
    //    q0[1] - q1[1],
    //    q0[2] - q1[2],
    //    q1[0] - q0[0],
    //    q1[1] - q0[1],
    //    q1[2] - q0[2];
    ////不能用下面注释的方法，会有bug
    ////q << q0 - q1, 
    ////    q1 - q0;
    //double l = (q0 - q1).norm();
    //f = stiffness * (l - l0) * q / l;
    ////-----------------------------------------

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd neg_I = -1 * Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd B(I.rows(), I.cols() + I.cols());
    B << neg_I, I;
    // cat q0 and q1 vertically
    Eigen::Vector6d q;
    q << q0, q1;

    // first order derivative by hand
    double spring_len = (q.transpose() * B.transpose() * B * q).array().sqrt()(0);
    double scalar = stiffness * (spring_len - l0) / spring_len;

    Eigen::Vector3d upper_half = -1 * (q1 - q0);
    Eigen::Vector3d lower_half = (q1 - q0);
    Eigen::Vector6d temp;
    temp << upper_half, lower_half;

    f = scalar * temp;
}