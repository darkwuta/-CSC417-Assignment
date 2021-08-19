#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {

    I.setZero();
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
        double area = dx1.cross(dx2).norm()/2;
        Eigen::Vector3d n = dx1.cross(dx2).normalized();

        mass += density *1.0/6.0*(x0.x()+x1.x()+x2.x())*n.x();

        // use matlab
        // 不要为了省地方把代码写在一行，容易出错！
        double x_0, x_1, x_2;
        x_0 = x0.x(); 
        x_1 = x1.x(); 
        x_2 = x2.x();
        double y_0, y_1, y_2;
        y_0 = x0.y();
        y_1 = x1.y();
        y_2 = x2.y();
        double z_0, z_1, z_2;
        z_0 = x0.z();
        z_1 = x1.z();
        z_2 = x2.z();

        double xxx = x_0 * (x_1 * x_1) * (1.0 / 2.0E1) + (x_0 * x_0) * x_1 * (1.0 / 2.0E1) + x_0 * (x_2 * x_2) * (1.0 / 2.0E1) + (x_0 * x_0) * x_2 * (1.0 / 2.0E1) + x_1 * (x_2 * x_2) * (1.0 / 2.0E1) + (x_1 * x_1) * x_2 * (1.0 / 2.0E1) + (x_0 * x_0 * x_0) * (1.0 / 2.0E1) + (x_1 * x_1 * x_1) * (1.0 / 2.0E1) + (x_2 * x_2 * x_2) * (1.0 / 2.0E1) + x_0 * x_1 * x_2 * (1.0 / 2.0E1);
        double yyy = y_0 * (y_1 * y_1) * (1.0 / 2.0E1) + (y_0 * y_0) * y_1 * (1.0 / 2.0E1) + y_0 * (y_2 * y_2) * (1.0 / 2.0E1) + (y_0 * y_0) * y_2 * (1.0 / 2.0E1) + y_1 * (y_2 * y_2) * (1.0 / 2.0E1) + (y_1 * y_1) * y_2 * (1.0 / 2.0E1) + (y_0 * y_0 * y_0) * (1.0 / 2.0E1) + (y_1 * y_1 * y_1) * (1.0 / 2.0E1) + (y_2 * y_2 * y_2) * (1.0 / 2.0E1) + y_0 * y_1 * y_2 * (1.0 / 2.0E1);
        double zzz = z_0 * (z_1 * z_1) * (1.0 / 2.0E1) + (z_0 * z_0) * z_1 * (1.0 / 2.0E1) + z_0 * (z_2 * z_2) * (1.0 / 2.0E1) + (z_0 * z_0) * z_2 * (1.0 / 2.0E1) + z_1 * (z_2 * z_2) * (1.0 / 2.0E1) + (z_1 * z_1) * z_2 * (1.0 / 2.0E1) + (z_0 * z_0 * z_0) * (1.0 / 2.0E1) + (z_1 * z_1 * z_1) * (1.0 / 2.0E1) + (z_2 * z_2 * z_2) * (1.0 / 2.0E1) + z_0 * z_1 * z_2 * (1.0 / 2.0E1);
        double xxy = (x_0 * x_0)* y_0* (1.0 / 2.0E1) + (x_0 * x_0) * y_1 * (1.0 / 6.0E1) + (x_1 * x_1) * y_0 * (1.0 / 6.0E1) + (x_0 * x_0) * y_2 * (1.0 / 6.0E1) + (x_1 * x_1) * y_1 * (1.0 / 2.0E1) + (x_2 * x_2) * y_0 * (1.0 / 6.0E1) + (x_1 * x_1) * y_2 * (1.0 / 6.0E1) + (x_2 * x_2) * y_1 * (1.0 / 6.0E1) + (x_2 * x_2) * y_2 * (1.0 / 2.0E1) + x_0 * x_1 * y_0 * (1.0 / 3.0E1) + x_0 * x_1 * y_1 * (1.0 / 3.0E1) + x_0 * x_2 * y_0 * (1.0 / 3.0E1) + x_0 * x_1 * y_2 * (1.0 / 6.0E1) + x_0 * x_2 * y_1 * (1.0 / 6.0E1) + x_1 * x_2 * y_0 * (1.0 / 6.0E1) + x_0 * x_2 * y_2 * (1.0 / 3.0E1) + x_1 * x_2 * y_1 * (1.0 / 3.0E1) + x_1 * x_2 * y_2 * (1.0 / 3.0E1);
        double yyz = (y_0 * y_0) * z_0 * (1.0 / 2.0E1) + (y_0 * y_0) * z_1 * (1.0 / 6.0E1) + (y_1 * y_1) * z_0 * (1.0 / 6.0E1) + (y_0 * y_0) * z_2 * (1.0 / 6.0E1) + (y_1 * y_1) * z_1 * (1.0 / 2.0E1) + (y_2 * y_2) * z_0 * (1.0 / 6.0E1) + (y_1 * y_1) * z_2 * (1.0 / 6.0E1) + (y_2 * y_2) * z_1 * (1.0 / 6.0E1) + (y_2 * y_2) * z_2 * (1.0 / 2.0E1) + y_0 * y_1 * z_0 * (1.0 / 3.0E1) + y_0 * y_1 * z_1 * (1.0 / 3.0E1) + y_0 * y_2 * z_0 * (1.0 / 3.0E1) + y_0 * y_1 * z_2 * (1.0 / 6.0E1) + y_0 * y_2 * z_1 * (1.0 / 6.0E1) + y_1 * y_2 * z_0 * (1.0 / 6.0E1) + y_0 * y_2 * z_2 * (1.0 / 3.0E1) + y_1 * y_2 * z_1 * (1.0 / 3.0E1) + y_1 * y_2 * z_2 * (1.0 / 3.0E1);
        double xxz = (x_0 * x_0) * z_0 * (1.0 / 2.0E1) + (x_0 * x_0) * z_1 * (1.0 / 6.0E1) + (x_1 * x_1) * z_0 * (1.0 / 6.0E1) + (x_0 * x_0) * z_2 * (1.0 / 6.0E1) + (x_1 * x_1) * z_1 * (1.0 / 2.0E1) + (x_2 * x_2) * z_0 * (1.0 / 6.0E1) + (x_1 * x_1) * z_2 * (1.0 / 6.0E1) + (x_2 * x_2) * z_1 * (1.0 / 6.0E1) + (x_2 * x_2) * z_2 * (1.0 / 2.0E1) + x_0 * x_1 * z_0 * (1.0 / 3.0E1) + x_0 * x_1 * z_1 * (1.0 / 3.0E1) + x_0 * x_2 * z_0 * (1.0 / 3.0E1) + x_0 * x_1 * z_2 * (1.0 / 6.0E1) + x_0 * x_2 * z_1 * (1.0 / 6.0E1) + x_1 * x_2 * z_0 * (1.0 / 6.0E1) + x_0 * x_2 * z_2 * (1.0 / 3.0E1) + x_1 * x_2 * z_1 * (1.0 / 3.0E1) + x_1 * x_2 * z_2 * (1.0 / 3.0E1);

        //// compute inertia matrix
        I(0, 0) += 2 * area * 1.0 / 3.0 * density * (n.y() * yyy + n.z() * zzz);
        I(0, 1) += 2 * area * -1.0 / 2.0 * density * n.x() * xxy;
        I(0, 2) += 2 * area * -1.0 / 2.0 * density * n.x() * xxz;
                  
        I(1, 0) += 2 * area * -1.0 / 2.0 * density * n.x() * xxy;
        I(1, 1) += 2 * area * 1.0 / 3.0 * density * (n.x() * xxx + n.z() * zzz);
        I(1, 2) += 2 * area * -1.0 / 2.0 * density * n.y() * yyz;
                 
        I(2, 0) += 2 * area * -1.0 / 2.0 * density * n.x() * xxz;
        I(2, 1) += 2 * area * 1.0 / 3.0 * density * n.y() * yyz;
        I(2, 2) += 2 * area * -1.0 / 2.0 * density * (n.x() * xxx + n.y() * yyy);
    }
    

    Eigen::Vector3d integral_center;
    integral_center.setZero();
    for (int i = 0; i < F.rows(); i++)
    {
        integral_center += density * V.row(i);
    }
    center = integral_center / mass;



}