#include <d2V_membrane_corotational_dq2.h>
#include <iostream>

void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    //
    ////Compute SVD of F here
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

    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }

    Eigen::Tensor3333d dU;
    Eigen::Tensor333d dS;
    Eigen::Tensor3333d dW;
    dsvd(dU, dS, dW, F);

    //
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

    ////TODO: compute H, the hessian of the corotational energy
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

    dn_dq = 1.0 / dx1.cross(dx2).norm() * (I - n * n.transpose()) * (cross_dx1 * Identity1 - cross_dx2 * Identity2);

    Eigen::Matrix99d dF_dq;
    dF_dq = Bj + Nj * dn_dq;


    Eigen::Matrix99d d2phi_d2F;
    Eigen::Matrix3d d2phi_d2F_ij;
    Eigen::Vector3d ds_diagonal;
    for (int i = 0; i < 3; i++)
        ds_diagonal[i] = 2.0 * mu * (S[i] - 1.0) + lambda * (S.sum() - 3.0) * 1.0;

    Eigen::DiagonalMatrix<double, 3>ds = ds_diagonal.asDiagonal();

    Eigen::Matrix3d d2s;
    d2s << 2.0 * mu + lambda, lambda, lambda,
        lambda, 2.0 * mu + lambda, lambda,
        lambda, lambda, 2.0 * mu + lambda;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Eigen::Matrix3d dU_ij = dU[i][j];
            Eigen::Matrix3d dW_ij = dW[i][j];
            Eigen::Vector3d d2s_ij = d2s * dS[i][j];
            d2phi_d2F_ij = dU_ij * ds * W.transpose() + U * d2s_ij.asDiagonal() * W.transpose() + U * ds * dW_ij.transpose();
            Eigen::Vector9d d2phi_d2F_ij_flatten;
            d2phi_d2F_ij_flatten << d2phi_d2F_ij(0, 0), d2phi_d2F_ij(0, 1), d2phi_d2F_ij(0, 2), d2phi_d2F_ij(1, 0), d2phi_d2F_ij(1, 1), d2phi_d2F_ij(1, 2), d2phi_d2F_ij(2, 0), d2phi_d2F_ij(2, 1), d2phi_d2F_ij(2, 2);
            d2phi_d2F.row(3 * i + j) = d2phi_d2F_ij_flatten;
        }
    }

    H = area * dF_dq.transpose() * d2phi_d2F * dF_dq;

    

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();

}
