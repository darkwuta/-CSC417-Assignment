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

    ////Deformation Gradient
    //Eigen::Matrix3d dx; //deformed tangent matrix 
    //Eigen::Matrix3d U;
    //Eigen::Vector3d S; 
    //Eigen::Matrix3d W; 

    ////TODO: SVD Here
    //Eigen::Vector3d x0;
    //Eigen::Vector3d x1;
    //Eigen::Vector3d x2;
    //x0 = V.row(element(0));
    //x1 = V.row(element(1));
    //x2 = V.row(element(2));

    //Eigen::Vector3d q0;
    //Eigen::Vector3d q1;
    //Eigen::Vector3d q2;
    //q0 = q.segment(3 * element(0), 3);
    //q1 = q.segment(3 * element(1), 3);
    //q2 = q.segment(3 * element(2), 3);


    //Eigen::Vector3d N = (x1 - x0).cross(x2 - x0).normalized();//δ�α�ʱ�ķ�����
    //Eigen::Vector3d n = (q1 - q0).cross(q2 - q0).normalized();

    //Eigen::Matrix3d F;
    //Eigen::Matrix34d X;
    //X.col(0) = x0;
    //X.col(1) = x1;
    //X.col(2) = x2;
    //X.col(3) = n;

    //Eigen::Matrix43d t;
    //t.block<3, 3>(0, 0) = dX;
    //t.row(4) = N.transpose();

    //F = X * t;


    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    //U = svd.matrixU();
    //S = svd.singularValues();
    //W = svd.matrixV();

    ////Fix for inverted elements (thanks to Danny Kaufman)
    //double det = S[0]*S[1];
    //
    // if(det <= -1e-10)
    //{
    //    if(S[0] < 0) S[0] *= -1;
    //    if(S[1] < 0) S[1] *= -1;
    //    if(S[2] < 0) S[2] *= -1;
    //}
    //
    //if(U.determinant() <= 0)
    //{
    //    U(0, 2) *= -1;
    //    U(1, 2) *= -1;
    //    U(2, 2) *= -1;
    //}
    //
    //if(W.determinant() <= 0)
    //{
    //    W(0, 2) *= -1;
    //    W(1, 2) *= -1;
    //    W(2, 2) *= -1;
    //}
    //
    ////TODO: energy model gradient 

    //Eigen::Matrix99d Bj;
    //Bj.setZero();
    //Bj << dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0, 0,
    //    dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0, 0,
    //    dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 1), 0, 0,
    //    0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0,
    //    0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0,
    //    0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0,
    //    0, 0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0),
    //    0, 0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1),
    //    0, 0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2);

    //Eigen::Vector9d Qj;
    //Qj << q0(0),
    //    q0(1),
    //    q0(2),
    //    q1(0),
    //    q1(1),
    //    q1(2),
    //    q2(0),
    //    q2(1),
    //    q2(2);

    //Eigen::Matrix93d Nj;
    //Nj << N(0), 0, 0,
    //    N(1), 0, 0,
    //    N(2), 0, 0,
    //    0, N(0), 0,
    //    0, N(1), 0,
    //    0, N(2), 0,
    //    0, 0, N(0),
    //    0, 0, N(1),
    //    0, 0, N(2);

    //Eigen::Matrix39d dn_dq;
    //Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    //Eigen::Matrix3d cross_dx1, cross_dx2;
    //Eigen::Vector3d dx1, dx2;
    //dx1 = q1 - q0;
    //dx2 = q2 - q0;
    //cross_dx1 << 0, -1.0 * dx1[2], dx1[1],
    //    dx1[2], 0, -1.0 * dx1[0],
    //    -1.0 * dx1[1], dx1[0], 0;
    //cross_dx2 << 0, -1.0 * dx2[2], dx2[1],
    //    dx2[2], 0, -1.0 * dx2[0],
    //    -1.0 * dx2[1], dx2[0], 0;

    //Eigen::Matrix39d Identity1;
    //Identity1 << -I, Eigen::Matrix3d::Zero(), I;
    //Eigen::Matrix39d Identity2;
    //Identity2 << -I, I, Eigen::Matrix3d::Zero();

    //dn_dq = 1.0 / n.norm() * (I - n * n.transpose()) * (cross_dx1 * Identity1 - cross_dx2 * Identity2);

    //Eigen::Matrix99d dF_dq;
    //dF_dq = Bj + Nj * dn_dq;

    //Eigen::Vector3d ds_diagonal;
    //for(int i=0;i<3;i++)
    //    ds_diagonal[i]=2.0 * mu * (S[i] - 1.0) + lambda * (S.sum() - 3.0) * 1.0;

    //Eigen::Vector9d dpshi_dF = U * ds_diagonal.asDiagonal() * W.transpose();

    //dV = -area * dF_dq.transpose() * dpshi_dF;
///////////////////////////////////////////////////////////////////////


//Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix
    Eigen::Matrix3d U;
    Eigen::Vector3d S;
    Eigen::Matrix3d W;

    //TODO: SVD Here

    Eigen::Vector3d q0, q1, q2;
    q0 = q.segment(3 * element(0), 3);
    q1 = q.segment(3 * element(1), 3);
    q2 = q.segment(3 * element(2), 3);

    Eigen::Vector9d q_j;
    q_j << q0, q1, q2;

    Eigen::Vector3d X0, X1, X2;
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));

    // deformed gradient matrix
    Eigen::Matrix3d F;
    // deformed space unit normal vec
    Eigen::Vector3d n = (q1 - q0).cross(q2 - q0);
    n.normalize();

    // ** UN **deformed space unit normal vec
    Eigen::Vector3d N = (X1 - X0).cross(X2 - X0);
    N.normalize();

    // world space position matrix
    Eigen::Matrix34d x_world;
    x_world.col(0) = q0;
    x_world.col(1) = q1;
    x_world.col(2) = q2;
    x_world.col(3) = n;

    // ** UN **deformed space position matrix
    Eigen::Matrix43d X_undeform;
    X_undeform.block<3, 3>(0, 0) = dX;
    X_undeform.block<1, 3>(3, 0) = N.transpose();
    X_undeform(3,0) = N[0];
    X_undeform(3,1) = N[1];
    X_undeform(3,2) = N[2];


    F = x_world * X_undeform;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();


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



   Eigen::Vector3d ds_diag;
   for (int i=0; i < ds_diag.size(); i++) {
       ds_diag[i] = 2.0 * mu * (S[i] - 1.0) + lambda * (S.sum() - 3.0) * 1.0;
   }
   Eigen::DiagonalMatrix<double, 3> ds = ds_diag.asDiagonal();


    // dpsi_dF matrix
    Eigen::Matrix3d dpsi_dF = U * ds * W.transpose();

    // expand dpsi_dF row by row to a vector of size 9
    Eigen::Vector9d vec_dpsi_dF;
    for (int j=0; j < dpsi_dF.cols(); j++) {
        vec_dpsi_dF.segment(3 * j, 3) = dpsi_dF.row(j);
    }

    // compute dF_dq_j
    Eigen::Matrix99d dF_dq_j;
    // require compute dn_dq_j
    Eigen::MatrixXd dn_dq_j;
    dn_dq_j.resize(3, 9);
    compute_dn_dq_j(dn_dq_j, n, q0, q1, q2);


    Eigen::Matrix99d B_j;
    B_j <<  dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0, 0,
            dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0, 0,
            dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0, 0,
            0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0), 0,
            0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1), 0,
            0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2), 0,
            0, 0, dX(0, 0), 0, 0, dX(1, 0), 0, 0, dX(2, 0),
            0, 0, dX(0, 1), 0, 0, dX(1, 1), 0, 0, dX(2, 1),
            0, 0, dX(0, 2), 0, 0, dX(1, 2), 0, 0, dX(2, 2);

    Eigen::MatrixXd N_j = Eigen::MatrixXd::Zero(9, 3);
    N_j.block<3,1>(0, 0) = N;
    N_j.block<3,1>(3, 1) = N;
    N_j.block<3,1>(6, 2) = N;

    dF_dq_j = B_j + N_j * dn_dq_j;

    dV = area * dF_dq_j.transpose() * vec_dpsi_dF;
}
