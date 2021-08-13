#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) { 
    f.resize(q.rows());
    f.setZero();
       for(unsigned int i =0;i<F.rows();i++)
       {
           Eigen::Vector9d tmp_f;
           Eigen::RowVector3i element;
           element = F.row(i);

           Eigen::Matrix<double, 1, 9> tmp_row;
           tmp_row = dX.row(i); //ei is the triangle index. 
           Eigen::Matrix3d dphi_dX = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());
          
           dV_membrane_corotational_dq(tmp_f, q, dphi_dX, V, element, a0(i), mu, lambda);



           f.segment(element(0) * 3, 3) -= tmp_f.segment(0, 3);
           f.segment(element(1) * 3, 3) -= tmp_f.segment(3, 3);
           f.segment(element(2) * 3, 3) -= tmp_f.segment(6, 3);
       }
    };
