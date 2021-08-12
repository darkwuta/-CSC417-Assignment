#include <dV_membrane_corotational_dq.h>
#include <iostream>

void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
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

    Eigen::Matrix99d Bj;
    Bj.setZero();
    Bj << dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0, 0,
        dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0, 0,
        dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 1), 0, 0,
        0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0,
        0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0,
        0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0,
        0, 0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0),
        0, 0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1),
        0, 0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2);

    Eigen::Vector9d Qj;
    Qj << q0(0),
        q0(1),
        q0(2),
        q1(0),
        q1(1),
        q1(2),
        q2(0),
        q2(1),
        q2(2);

    Eigen::Matrix93d Nj;
    Nj << N(0), 0, 0,
        N(1), 0, 0,
        N(2), 0, 0,
        0, N(0), 0,
        0, N(1), 0,
        0, N(2), 0,
        0, 0, N(0),
        0, 0, N(1),
        0, 0, N(2);

    Eigen::Matrix39d dn_dq;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d cross_dx1, cross_dx2;
    Eigen::Vector3d dx1, dx2;
    dx1 = q1 - q0;
    dx2 = q2 - q0;
    cross_dx1 << 0, -1.0 * dx1[2], dx1[1],
        dx1[2], 0, -1.0 * dx1[0],
        -1.0 * dx1[1], dx1[0], 0;
    cross_dx2 << 0, -1.0 * dx2[2], dx2[1],
        dx2[2], 0, -1.0 * dx2[0],
        -1.0 * dx2[1], dx2[0], 0;

    Eigen::Matrix39d Identity1;
    Identity1 << -I, Eigen::Matrix3d::Zero(), I;
    Eigen::Matrix39d Identity2;
    Identity2 << -I, I, Eigen::Matrix3d::Zero();

    dn_dq = 1.0 / n.norm() * (I - n * n.transpose()) * (cross_dx1 * Identity1 - cross_dx2 * Identity2);




    // svd目前不太懂
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //Eigen::Vector3d singularVals = svd.singularValues();
    //double psi = mu * (pow(singularVals[0] - 1.0, 2) + pow(singularVals[1] - 1.0, 2) + pow(singularVals[2] - 1.0, 2)) +
    //    lambda / 2.0 * pow(singularVals.sum() - 3.0, 2);

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient 

}
