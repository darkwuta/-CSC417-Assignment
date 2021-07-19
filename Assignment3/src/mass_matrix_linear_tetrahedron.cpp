 
 #include <mass_matrix_linear_tetrahedron.h>


 void mass_matrix_linear_tetrahedron(Eigen::Matrix1212d &M, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
     ////--------------------------------------------------           
     Eigen::Matrix1212d M_partial;
     M_partial.setZero();
     for (int i = 0; i < 4; i++) {
         for (int j = 0; j < 4; j++) {
             if (i == j) {
                 M_partial.block<3, 3>(i * 3, j * 3) = 1.0 / 60 * Eigen::Matrix3d::Identity();
             }
             else {
                 M_partial.block<3, 3>(i * 3, j * 3) = 1.0 / 120 * Eigen::Matrix3d::Identity();
             }
         }
     }

     M = 6.0 * density * volume * M_partial;
     ////--------------------------------------------------

    /* M.setZero();
     double coeff_10 = (1.0 / 10.0) * density * volume;
     double coeff_20 = (1.0 / 20.0) * density * volume;


     M.block(0, 0, 3, 3) = coeff_10 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(0, 3, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(0, 6, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(0, 9, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();

     M.block(3, 0, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(3, 3, 3, 3) = coeff_10 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(3, 6, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(3, 9, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();

     M.block(6, 0, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(6, 3, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(6, 6, 3, 3) = coeff_10 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(6, 9, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();

     M.block(9, 0, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(9, 3, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(9, 6, 3, 3) = coeff_20 * Eigen::Matrix<double, 3, 3>::Identity();
     M.block(9, 9, 3, 3) = coeff_10 * Eigen::Matrix<double, 3, 3>::Identity();*/
 }