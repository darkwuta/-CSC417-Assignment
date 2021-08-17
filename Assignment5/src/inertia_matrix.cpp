#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector3d x0, x1, x2;
        Eigen::Vector3i element = F.row(i);
        x0 = V.row(element[0]);
        x1 = V.row(element[1]);
        x2 = V.row(element[2]);
        Eigen::Vector3d dx1, dx2;
        dx1 = x1 - x0;
        dx2 = x2 - x0;
        //double area = dx1.cross(dx2).norm()/2;
        Eigen::Matrix3d n = dx1.cross(dx2).normalized();

        mass += density *1.0/6.0*(x0.x()+x1.x()+x2.x())*n.x();
    }
    

    Eigen::Vector3d integral_center;
    for (int i = 0; i < V.rows(); i++)
    {
        integral_center += density * V.row(i);
    }
    center = integral_center / mass;



}