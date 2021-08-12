#include <assemble_stiffness.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) { 
    K.resize(q.size(), q.size());

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> tripleList;
    tripleList.reserve(q.size() * q.size());

    for (unsigned int i = 0; i < q.size() / 3; i++)
    {
        Eigen::Matrix99d tmp_H;
        Eigen::RowVector3i element;
        element = F.row(i);

        d2V_membrane_corotational_dq2(tmp_H, q, dX, V, element, a0(i), mu, lambda);

        Eigen::Matrix3d H_00 = tmp_H.block<3, 3>(0, 0);
        Eigen::Matrix3d H_01 = tmp_H.block<3, 3>(0, 3);
        Eigen::Matrix3d H_02 = tmp_H.block<3, 3>(0, 6);

        Eigen::Matrix3d H_10 = tmp_H.block<3, 3>(3, 0);
        Eigen::Matrix3d H_11 = tmp_H.block<3, 3>(3, 6);
        Eigen::Matrix3d H_12 = tmp_H.block<3, 3>(3, 9);

        Eigen::Matrix3d H_20 = tmp_H.block<3, 3>(6, 0);
        Eigen::Matrix3d H_21 = tmp_H.block<3, 3>(6, 6);
        Eigen::Matrix3d H_22 = tmp_H.block<3, 3>(6, 9);

        for (int x = 0; x < 3; x++)
        {
            for (int y = 0; y < 3; y++)
            {
                tripleList.push_back(Trip(3 * element[0] + x, 3 * element[0] + y, H_00.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[0] + x, 3 * element[1] + y, H_01.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[0] + x, 3 * element[2] + y, H_02.coeffRef(x, y)));

                tripleList.push_back(Trip(3 * element[1] + x, 3 * element[0] + y, H_10.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[1] + x, 3 * element[1] + y, H_11.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[1] + x, 3 * element[2] + y, H_12.coeffRef(x, y)));

                tripleList.push_back(Trip(3 * element[2] + x, 3 * element[0] + y, H_20.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[2] + x, 3 * element[1] + y, H_21.coeffRef(x, y)));
                tripleList.push_back(Trip(3 * element[2] + x, 3 * element[2] + y, H_22.coeffRef(x, y)));

            }
        }

        
    }
    K.setFromTriplets(tripleList.begin(), tripleList.end());
       
        
    };
