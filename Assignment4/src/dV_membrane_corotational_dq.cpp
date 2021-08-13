#include <dV_membrane_corotational_dq.h>
#include <iostream>
void compute_dn_dq_j(Eigen::MatrixXd& dn_dq_j, Eigen::Vector3d n, Eigen::Vector3d q0, Eigen::Vector3d q1, Eigen::Vector3d q2) {

    Eigen::Vector3d dx1 = q1 - q0;
    Eigen::Vector3d dx2 = q2 - q0;

    Eigen::Matrix3d skew_x1;
    skew_x1 << 0, -1.0 * dx1[2], dx1[1],
        dx1[2], 0, -1.0 * dx1[0],
        -1.0 * dx1[1], dx1[0], 0;

    Eigen::Matrix3d skew_x2;
    skew_x2 << 0, -1.0 * dx2[2], dx2[1],
        dx2[2], 0, -1.0 * dx2[0],
        -1.0 * dx2[1], dx2[0], 0;

    Eigen::MatrixXd concat_Identity1;
    concat_Identity1.resize(3, 9);
    concat_Identity1 << -1.0 * Eigen::MatrixXd::Identity(3, 3),
        Eigen::MatrixXd::Zero(3, 3),
        Eigen::MatrixXd::Identity(3, 3);

    Eigen::MatrixXd concat_Identity2;
    concat_Identity2.resize(3, 9);
    concat_Identity2 << -1.0 * Eigen::MatrixXd::Identity(3, 3),
        Eigen::MatrixXd::Identity(3, 3),
        Eigen::MatrixXd::Zero(3, 3);

    dn_dq_j = (1.0 / (dx1.cross(dx2)).norm()) * (Eigen::MatrixXd::Identity(3, 3) - n * n.transpose()) * (skew_x1 * concat_Identity1 - skew_x2 * concat_Identity2);

}
void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //// Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    // TODO: SVD Here
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
    X.col(0) = q0;
    X.col(1) = q1;
    X.col(2) = q2;
    X.col(3) = n;

    Eigen::Matrix43d t;
    t.block<3, 3>(0, 0) = dX;
    t.row(3) = N.transpose();//是第四行

    F = X * t;


    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();

    //// Fix for inverted elements (thanks to Danny Kaufman)
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
    
    // TODO: energy model gradient 

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

    Eigen::Matrix99d dF_dq;
    dF_dq = Bj + Nj * dn_dq;

    Eigen::Vector3d ds_diagonal;
    for(int i=0;i<3;i++)
        ds_diagonal[i]=2.0 * mu * (S[i] - 1.0) + lambda * (S.sum() - 3.0) * 1.0;

    Eigen::Matrix3d dphi_dF = U * ds_diagonal.asDiagonal() * W.transpose();

    Eigen::Vector9d vec_dphi_dF;
    vec_dphi_dF.segment(0, 3) = dphi_dF.row(0);
    vec_dphi_dF.segment(3, 3) = dphi_dF.row(1);
    vec_dphi_dF.segment(6, 3) = dphi_dF.row(2);

    dV = area * dF_dq.transpose() * vec_dphi_dF;//这里按照公式应该是负的：-area * dF_dq.transpose() * vec_dphi_dF;
                                                //原因：正负在assemble两个函数中决定

}
