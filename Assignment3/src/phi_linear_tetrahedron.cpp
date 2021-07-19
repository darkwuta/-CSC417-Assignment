#include <phi_linear_tetrahedron.h>

void phi_linear_tetrahedron(Eigen::Vector4d &phi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> x) {
    ////------------------------------------
    Eigen::Vector3d X0 = V.row(element[0]);
    Eigen::Vector3d X1 = V.row(element[1]);
    Eigen::Vector3d X2 = V.row(element[2]);
    Eigen::Vector3d X3 = V.row(element[3]);

    Eigen::Matrix3d T;
    T.col(0) = X1 - X0;
    T.col(1) = X2 - X0;
    T.col(2) = X3 - X0;

    Eigen::Vector3d phi_partial = T.inverse() * (x - X0);
    phi(0) = 1.0 - phi_partial.sum();
    phi(1) = phi_partial[0];
    phi(2) = phi_partial[1];
    phi(3) = phi_partial[2]; 
    ////------------------------------------
    // Get the position of each vertices
    //Eigen::Vector3d X_0, X_1, X_2, X_3;
    //X_0 = V.row(element(0));
    //X_1 = V.row(element(1));
    //X_2 = V.row(element(2));
    //X_3 = V.row(element(3));

    //// Construct the T
    //Eigen::Matrix3d T;
    //T.setZero();
    //T.block(0, 0, 3, 1) = X_1 - X_0;
    //T.block(0, 1, 3, 1) = X_2 - X_0;
    //T.block(0, 2, 3, 1) = X_3 - X_0;

    //// phi_123 = T^-1(X - X_0)
    //Eigen::Vector3d phi_123;
    //phi_123 = T.inverse() * (x - X_0);

    //// phi_0 = 1 - phi_1 - phi_2 - phi_3
    //double phi_0 = 1 - phi_123(0) - phi_123(1) - phi_123(2);

    //phi(0) = phi_0;
    //phi(1) = phi_123(0);
    //phi(2) = phi_123(1);
    //phi(3) = phi_123(2);
}