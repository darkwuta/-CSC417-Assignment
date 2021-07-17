#include <Eigen/Dense>
#include <EigenTypes.h>

//Input:
//  q - generalized coordinates of the FEM system
//  element - vertex indices for the tetrahedron
// volume - volume of the tetrahedron
// integrand(out, q, X) - function to be integrated, returns value in out.
//Output:
//  integrated - the value of the integrated function
template<typename Ret, typename Integrand_Func>
inline void quadrature_single_point(Ret &&integrated, Eigen::Ref<const Eigen::VectorXd> q, 
                                               Eigen::Ref<const Eigen::RowVectorXi> element, double volume,
                                               Integrand_Func integrand) {
    ////--------------------------------------
    //Eigen::Vector12d x;
    //x << q.segment(3 * element(0), 3), q.segment(3 * element(1), 3), q.segment(3 * element(2), 3), q.segment(3 * element(3), 3);

    //Eigen::Vector3d X;

    //X(0) = (x(0) + x(3) + x(6) + x(9)) / 4;
    //X(1) = (x(1) + x(4) + x(7) + x(10)) / 4;
    //X(2) = (x(2) + x(5) + x(8) + x(11)) / 4;

    //integrand(integrated, q, element, X);
    //integrated *= volume;
    ////--------------------------------------

    Eigen::Vector3d X = Eigen::VectorXd::Zero(3);
    integrand(integrated, q, element, X);

    integrated = volume * integrated;
}

