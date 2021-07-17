#include <T_linear_tetrahedron.h>
#include <mass_matrix_linear_tetrahedron.h>
#include <iostream>
void T_linear_tetrahedron(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::RowVectorXi> element, double density, double volume) {
    ////----------------------------------------------------
    //Eigen::Matrix1212d M0;
    //Eigen::Vector12d v;
    //for (unsigned int i = 0; i < 4; i++)
    //    v.segment(i * 3, 3) = qdot.segment(element(i) * 3, 3);
    //mass_matrix_linear_tetrahedron(M0, qdot, element, density, volume);
    //T = 0.5 * v.transpose() * M0 * v;
    ////----------------------------------------------------

    Eigen::Vector12d qdot_tet;
    for (int i = 0; i < element.size(); i++) {
        qdot_tet.segment(i * 3, 3) = qdot.segment(element(i) * 3, 3);
    }
    Eigen::Matrix1212d M_j;
    mass_matrix_linear_tetrahedron(M_j, qdot_tet, element, density, volume);

    T = (0.5 * (qdot_tet.transpose() * M_j * qdot_tet))(0, 0);
}