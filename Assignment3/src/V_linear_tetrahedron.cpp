#include <V_linear_tetrahedron.h>

#include <dphi_linear_tetrahedron_dX.h>
#include <psi_neo_hookean.h>
#include <quadrature_single_point.h>

#include<iostream>

void V_linear_tetrahedron(double &energy, Eigen::Ref<const Eigen::VectorXd> q, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                          double C, double D) {

    
    ////------------------------------------------------
    auto neohookean_linear_tet = [&](double &e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
        // TODO 这里用下面这段就是图上不显示，为什么？
        //Eigen::Vector3d X0 = V.row(element[0]);
        //Eigen::Vector3d X1 = V.row(element[1]);
        //Eigen::Vector3d X2 = V.row(element[2]);
        //Eigen::Vector3d X3 = V.row(element[3]);
        //Eigen::Matrix34d x;
        //x << X0, X1, X2, X3;

        Eigen::Matrix34d x;
        for (int i = 0; i < 4; i++) {
            x.col(i) = q.segment(element(i) * 3, 3);
        }

        Eigen::Matrix43d dphi;
        dphi_linear_tetrahedron_dX(dphi,V,element,X);

        psi_neo_hookean(e, x * dphi, C, D);
    };

    quadrature_single_point(energy, q, element, volume, neohookean_linear_tet); 
    ////------------------------------------------------
    //Eigen::Matrix34d x;
    //x.setZero();

    //x.block(0, 0, 3, 1) = q.segment(3 * element(0), 3);
    //x.block(0, 1, 3, 1) = q.segment(3 * element(1), 3);
    //x.block(0, 2, 3, 1) = q.segment(3 * element(2), 3);
    //x.block(0, 3, 3, 1) = q.segment(3 * element(3), 3);


    //Eigen::Matrix43d dphi;


    //auto neohookean_linear_tet = [&](double& e, Eigen::Ref<const Eigen::VectorXd>q, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    //    dphi_linear_tetrahedron_dX(dphi, V, element, X);
    //    psi_neo_hookean(e, x * dphi, C, D); // deformation gradient = F = x*dphi

    //};

    //quadrature_single_point(energy, q, element, volume, neohookean_linear_tet);
}