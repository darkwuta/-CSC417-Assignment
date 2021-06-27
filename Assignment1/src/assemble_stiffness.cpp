#include <assemble_stiffness.h>
#include"iostream"

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
    typedef Eigen::Triplet<double> T;//以i，j，mass为索引的容器
    std::vector<T> tripletList;
    std::cout << "ASSEMBLE_STIFFNESS::DEBUG::K.rows():" << K.rows() << std::endl;
    for (int y = 0; y < E.rows(); y++)
    {
        Eigen::Matrix66d H;
        int i = E(y, 0);
        int j = E(y, 1);
        Eigen::Vector3d q0, q1;
        q0 << q(i * 3);
        q1 << q(j * 3);
        //std::cout << "ASSEMBLE_STIFFNESS::DEBUG::q0:" <<q0 << std::endl;
        //std::cout << "ASSEMBLE_STIFFNESS::DEBUG::q1:" << q1 << std::endl;
        //std::cout << "ASSEMBLE_STIFFNESS::DEBUG::num:" << y << std::endl;
        //std::cout << "ASSEMBLE_STIFFNESS::DEBUG::H:" << H << std::endl;
        d2V_spring_particle_particle_dq2(H, q0, q1, l0(y), k);
        for (int k = 0; k < 6; k++)
        {
            for (int s = 0; s < 6; s++)
            {
                tripletList.push_back(T(i * 6 + k, j * 6 + s, H(k, s)));
            }
        }
    }
    std::cout << "ASSEMBLE_STIFFNESS::DEBUG::tripletList.size():" << tripletList.size() << std::endl;
    /*std::cout << "ASSEMBLE_STIFFNESS::DEBUG::K.rows:" << K.rows()<<" K.cols"<<K.cols() << std::endl;
    std::cout << "ASSEMBLE_STIFFNESS::DEBUG::K.IsRowMajor():" << K.IsRowMajor()<< std::endl;
    std::cout << "ASSEMBLE_STIFFNESS::DEBUG::sm2.IsRowMajor():" << K.IsRowMajor() << std::endl;*/
    //sm2.setFromTriplets(1, 5);

    K.setFromTriplets(tripletList.begin(), tripletList.end());
        
};