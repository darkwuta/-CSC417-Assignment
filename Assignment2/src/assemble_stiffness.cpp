#include <assemble_stiffness.h>
#include"iostream"
//自己的有问题
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    //-----------------------------------------------------------------------------------
    typedef Eigen::Triplet<double> T;
    //std::vector<T> tripletList;
    //tripletList.reserve(q.rows() * q.rows());
    //K.resize(q.rows(), q.rows());
    //K.setZero();
    ////std::cout << "ASSEMBLE_STIFFNESS::DEBUG::K.rows():" << K.rows() << std::endl;
    //for (int y = 0; y < E.rows(); y++)
    //{
    //    Eigen::Matrix66d H;
    //    int i = E(y, 0);
    //    int j = E(y, 1);
    //    Eigen::Vector3d q0, q1;
    //    q0 << q(i * 3), q(i * 3 + 1) , q(i * 3 + 2);
    //    q1 << q(j * 3), q(j * 3 + 1) , q(j * 3 + 2);
    //    d2V_spring_particle_particle_dq2(H, q0, q1, l0(y), k);
    //    Eigen::Matrix3d Kaa, Kab, Kba, Kbb;
    //    Kaa << H(0, 0), H(0, 1), H(0, 2),
    //        H(1, 0), H(1, 1), H(1, 2),
    //        H(2, 0), H(2, 1), H(2, 2);
    //    Kab << H(0, 3), H(0, 4), H(0, 5),
    //        H(1, 3), H(1, 4), H(1, 5),
    //        H(2, 3), H(2, 4), H(2, 5);
    //    Kba << H(3, 0), H(3, 1), H(3, 2),
    //        H(4, 0), H(4, 1), H(4, 2),
    //        H(5, 0), H(5, 1), H(5, 2);
    //    Kbb << H(3, 3), H(3, 4), H(3, 5),
    //        H(4, 3), H(4, 4), H(4, 5),
    //        H(5, 3), H(5, 4), H(5, 5);
    //    //Kaa
    //    for (int a = 0; a < 3; a++)
    //    {
    //        for (int b = 0; b < 3; b++)
    //        {
    //            tripletList.push_back(T(i * 3 + a, i * 3 + b, -Kaa(a, b)));
    //        }
    //    }
    //    //Kab
    //    for (int a = 0; a < 3; a++)
    //    {
    //        for (int b = 0; b < 3; b++)
    //        {
    //            tripletList.push_back(T(i * 3 + a, j * 3 + b, -Kab(a, b)));
    //        }
    //    }
    //    //Kba
    //    for (int a = 0; a < 3; a++)
    //    {
    //        for (int b = 0; b < 3; b++)
    //        {
    //            tripletList.push_back(T(j * 3 + a, i * 3 + b, -Kba(a, b)));
    //        }
    //    }
    //    //Kbb
    //    for (int a = 0; a < 3; a++)
    //    {
    //        for (int b = 0; b < 3; b++)
    //        {
    //            tripletList.push_back(T(j * 3 + a, j * 3 + b, -Kbb(a, b)));
    //        }
    //    }
    //}

    //

    //K.setFromTriplets(tripletList.begin(), tripletList.end());
    //-----------------------------------------------------------------------------------
    Eigen::Matrix66d H;
    Eigen::Vector3d  q0, q1;
    K.resize(q.size(), q.size());
    K.setZero();
    std::vector<T> tripletList;
    tripletList.reserve(9 * V.rows() * V.rows());

    for (int i = 0; i < E.rows(); i++) {
        q0[0] = q[3 * E(i, 0)];
        q0[1] = q[3 * E(i, 0) + 1];
        q0[2] = q[3 * E(i, 0) + 2];
        q1[0] = q[3 * E(i, 1)];
        q1[1] = q[3 * E(i, 1) + 1];
        q1[2] = q[3 * E(i, 1) + 2];
        d2V_spring_particle_particle_dq2(H, q0, q1, l0[i], k);
        for (int m = 0; m < 3; m++) {
            for (int n = 0; n < 3; n++) {
                tripletList.push_back(T(3 * E(i, 0) + m, 3 * E(i, 0) + n, H(m, n)));
                tripletList.push_back(T(3 * E(i, 0) + m, 3 * E(i, 1) + n, H(m, n + 3.0)));
                tripletList.push_back(T(3 * E(i, 1) + m, 3 * E(i, 0) + n, H(m + 3.0, n)));
                tripletList.push_back(T(3 * E(i, 1) + m, 3 * E(i, 1) + n, H(m + 3.0, n + 3.0)));
            }
        }
    }
    K.setFromTriplets(tripletList.begin(), tripletList.end());
};