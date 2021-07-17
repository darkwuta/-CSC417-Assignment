#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

#include<iostream>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    
    ////------------------------------------------------
    //auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    //    
    //    Eigen::Vector3d X0 = V.row(element[0]);
    //    Eigen::Vector3d X1 = V.row(element[1]);
    //    Eigen::Vector3d X2 = V.row(element[2]);
    //    Eigen::Vector3d X3 = V.row(element[3]);

    //    Eigen::Matrix34d x;
    //    x << X0, X1, X2, X3;

    //    //std::cout << "V_linear_tetrahedron::DEBUG::x:" << x << std::endl;

    //    Eigen::Matrix43d dphi;
    //    dphi_linear_tetrahedron_dX(dphi,V,element,X);
    //    psi_neo_hookean(e, x*dphi, C, D); 
    //};

    //quadrature_single_point(energy, q, element, volume, neohookean_linear_tet); 
    ////------------------------------------------------

    auto neohookean_linear_tet = [&](double& e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        Eigen::MatrixXd x;
        x.resize(3, 4);
        for (int i = 0; i < 4; i++) {
            x.col(i) = q.segment(element(i) * 3, 3);
        }

        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi, V, element, X);

        Eigen::Matrix3d F;
        F = x * dphi;

        psi_neo_hookean(e, F, C, D);

    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);
}