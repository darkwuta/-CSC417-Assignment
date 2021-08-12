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

    Eigen::Vector3d q0;
    Eigen::Vector3d q1;
    Eigen::Vector3d q2;
    q0 = q.segment(3 * element(0), 3);
    q1 = q.segment(3 * element(1), 3);
    q2 = q.segment(3 * element(2), 3);


    Eigen::Vector3d N = (x1 - x0).cross(x2 - x0).normalized();//未形变时的法向量
    Eigen::Vector3d n = (q1 - q0).cross(q2 - q0).normalized();

    Eigen::Matrix3d F;
    Eigen::Matrix34d X;
    X.col(0) = x0;
    X.col(1) = x1;
    X.col(2) = x2;
    X.col(3) = n;

    Eigen::Matrix43d t;
    t.block<3, 3>(0, 0) = dX;
    t.row(4) = N.transpose();
    
    F = X * t;

    // svd目前不太懂
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d singularVals = svd.singularValues();
    double psi = mu * (pow(singularVals[0] - 1.0, 2) + pow(singularVals[1] - 1.0, 2) + pow(singularVals[2] - 1.0, 2)) +
        lambda / 2.0 * pow(singularVals.sum() - 3.0, 2);
    energy = area * psi;

}
