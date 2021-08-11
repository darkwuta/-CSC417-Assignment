#include <dphi_linear_tetrahedron_dX.h>
#include <phi_linear_tetrahedron.h>
#include <iostream>
void dphi_linear_tetrahedron_dX(Eigen::Matrix43d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    ////-------------------------------
    //Eigen::Vector3d X0 = V.row(element[0]);
    //Eigen::Vector3d X1 = V.row(element[1]);
    //Eigen::Vector3d X2 = V.row(element[2]);
    //Eigen::Vector3d X3 = V.row(element[3]);

    //Eigen::Matrix3d T;
    //T.col(0) = X1 - X0;
    //T.col(1) = X2 - X0;
    //T.col(2) = X3 - X0;

    //Eigen::Matrix3d inv_T = T.inverse();

    //dphi.setZero();
    //dphi.block(0, 0, 1, 3) = -inv_T.colwise().sum();//每一列的和
    //dphi.block(1, 0, 3, 3) = inv_T;
    ////dphi.block<1, 3>(0, 0) = -Eigen::Vector3d::Ones() * inv_T;
    ////dphi.block<3, 3>(1, 0) = inv_T;
    ////-------------------------------
    Eigen::Vector3d X_0, X_1, X_2, X_3;
    X_0 = V.row(element(0));
    X_1 = V.row(element(1));
    X_2 = V.row(element(2));
    X_3 = V.row(element(3));

    Eigen::Matrix3d T, T_inv;
    T.setZero();
    T.block(0, 0, 3, 1) = X_1 - X_0;
    T.block(0, 1, 3, 1) = X_2 - X_0;
    T.block(0, 2, 3, 1) = X_3 - X_0;

    T_inv = T.inverse();

    Eigen::RowVector4d I;
    I << 1, -1, -1, -1;



    dphi.setZero();
    dphi.block(0, 0, 1, 3) = T_inv.colwise().sum();
    dphi.block(1, 0, 3, 3) = T_inv;


}