#include <psi_neo_hookean.h>
#include <dphi_linear_tetrahedron_dX.h>
#include <iostream>
void psi_neo_hookean(double &psi, 
                     Eigen::Ref<const Eigen::Matrix3d> F,
                     double C, double D) {
    double J = F.determinant();
    double trace = (F.transpose() * F).trace();
    psi = C * (pow(J, -2 / 3) * trace - 3) + D * pow(J - 1, 2);
}