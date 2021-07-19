#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
    

    ////---------------------------------------------------
   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       Eigen::Matrix43d dphi;
       dphi_linear_tetrahedron_dX(dphi, V, element, X);

       Eigen::Vector3d X0 = V.row(element[0]);
       Eigen::Vector3d X1 = V.row(element[1]);
       Eigen::Vector3d X2 = V.row(element[2]);
       Eigen::Vector3d X3 = V.row(element[3]);

       Eigen::Matrix34d x;
       for (int i = 0; i < 4; i++) {
           x.col(i) = q.segment(element(i) * 3, 3);
       }

       Eigen::Vector9d psi;
       Eigen::Matrix3d F = x * dphi;
       dpsi_neo_hookean_dF(psi, F, C, D);
       
       //解方程F_flatten = B * q
       //求得T
       Eigen::Matrix3d T;
       T.col(0) = X1 - X0;
       T.col(1) = X2 - X0;
       T.col(2) = X3 - X0;

       Eigen::Matrix3d inv_T = T.inverse();

        // construct the B
        Eigen::MatrixXd B(9, 12);
        B.setZero();
        B.block(0, 0, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
        B.block(3, 1, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
        B.block(6, 2, 3, 1) = dphi.block(0, 0, 1, 3).transpose();

        B.block(0, 3, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
        B.block(3, 4, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
        B.block(6, 5, 3, 1) = dphi.block(1, 0, 1, 3).transpose();

        B.block(0, 6, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
        B.block(3, 7, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
        B.block(6, 8, 3, 1) = dphi.block(2, 0, 1, 3).transpose();

        B.block(0, 9, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
        B.block(3, 10, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
        B.block(6, 11, 3, 1) = dphi.block(3, 0, 1, 3).transpose();

       dV = B.transpose() * psi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    ////---------------------------------------------------

    //auto neohookean_linear_tet = [&](Eigen::Vector12d& dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    //    // D
    //    Eigen::Matrix43d dphi;
    //    dphi_linear_tetrahedron_dX(dphi, V, element, X);

    //    // construct the B
    //    Eigen::MatrixXd B(9, 12);
    //    B.setZero();
    //    B.block(0, 0, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    //    B.block(3, 1, 3, 1) = dphi.block(0, 0, 1, 3).transpose();
    //    B.block(6, 2, 3, 1) = dphi.block(0, 0, 1, 3).transpose();

    //    B.block(0, 3, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    //    B.block(3, 4, 3, 1) = dphi.block(1, 0, 1, 3).transpose();
    //    B.block(6, 5, 3, 1) = dphi.block(1, 0, 1, 3).transpose();

    //    B.block(0, 6, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    //    B.block(3, 7, 3, 1) = dphi.block(2, 0, 1, 3).transpose();
    //    B.block(6, 8, 3, 1) = dphi.block(2, 0, 1, 3).transpose();

    //    B.block(0, 9, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
    //    B.block(3, 10, 3, 1) = dphi.block(3, 0, 1, 3).transpose();
    //    B.block(6, 11, 3, 1) = dphi.block(3, 0, 1, 3).transpose();


    //    // q element
    //    Eigen::Vector12d q_el;
    //    q_el << q.segment(3 * element(0), 3), q.segment(3 * element(1), 3), q.segment(3 * element(2), 3), q.segment(3 * element(3), 3);

    //    // deformation gradient
    //    Eigen::Vector9d F_flatten = B * q_el;
    //    //std::cout << "F_flatten\n" << F_flatten << std::endl;
    //    Eigen::Map<Eigen::Matrix3d> F(F_flatten.data(), 3, 3);
    //    F.transposeInPlace();
    //    //std::cout << "F\n" << F<< std::endl;

    //    Eigen::Vector9d dF;
    //    dpsi_neo_hookean_dF(dF, F, C, D);
    //    dV = B.transpose() * dF;
    //};

    //quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);
}