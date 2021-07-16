#include <dV_linear_tetrahedron_dq.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <dpsi_neo_hookean_dF.h>
#include <quadrature_single_point.h>
#include <iostream>

void dV_linear_tetrahedron_dq(Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {
    
    Eigen::Vector3d X0 = V.row(element[0]);
    Eigen::Vector3d X1 = V.row(element[1]);
    Eigen::Vector3d X2 = V.row(element[2]);
    Eigen::Vector3d X3 = V.row(element[3]);

    Eigen::Matrix34d x;
    x << X0, X1, X2, X3;

    Eigen::Matrix43d dphi;

   auto neohookean_linear_tet = [&](Eigen::Vector12d &dV, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
       dphi_linear_tetrahedron_dX(dphi, V, element, X);
       Eigen::Vector9d psi;
       Eigen::Matrix3d F = x * dphi;
       dpsi_neo_hookean_dF(psi, F, C, D);
       Eigen::MatrixXd B(9, 12);
       
       //解方程F_flatten = B * q
       Eigen::Vector3d X0 = V.row(element[0]);
       Eigen::Vector3d X1 = V.row(element[1]);
       Eigen::Vector3d X2 = V.row(element[2]);
       Eigen::Vector3d X3 = V.row(element[3]);
       //求得T
       Eigen::Matrix3d T;
       T.col(0) = X1 - X0;
       T.col(1) = X2 - X0;
       T.col(2) = X3 - X0;

       Eigen::Matrix3d inv_T = T.inverse();

       Eigen::Vector9d F_flatten;
       F_flatten << F(0, 0), F(0, 1), F(0, 2), F(1, 1), F(1, 2), F(1, 3), F(2, 0), F(2, 1), F(2, 2);

       B = F_flatten * q.transpose() * (q * q.transpose()).inverse();

       dV = B.transpose() * psi;
    };

    quadrature_single_point(dV, q, element, volume, neohookean_linear_tet);  
    
}