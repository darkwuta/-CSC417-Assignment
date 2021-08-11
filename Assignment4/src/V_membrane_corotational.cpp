#include <V_membrane_corotational.h>

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    Eigen::Vector3d x0;
    Eigen::Vector3d x1;
    Eigen::Vector3d x2;
    x0 = V.row(element(0));
    x1 = V.row(element(1));
    x2 = V.row(element(2));


    Eigen::Matrix32d T;
    T.row(0) = x1 - x0;
    T.row(1) = x2 - x0;

    Eigen::Vector3d n = (x1 - x0).cross(x2 - x0).normalized();

    Eigen::Matrix3d F;
    Eigen::Matrix34d X;
    X.col(0) = x0;
    X.col(1) = x1;
    X.col(2) = x2;
    X.col(3) = n;

    Eigen::Matrix43d t;
    t.row(0) = Eigen::Vector2d::Ones() * (T.transpose() * T).inverse() * T.transpose();
    t.block<2, 3>(1, 0) = (T.transpose() * T).inverse() * T.transpose();
    t.row(4) = n;
    
    

}
