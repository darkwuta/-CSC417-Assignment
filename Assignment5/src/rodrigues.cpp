#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {
	
	Eigen::Vector3d a;

	a = omega / omega.norm();
	
	Eigen::Matrix3d cross_a;

	cross_a << 0, -a.z(), a.y(),
		a.z(), 0, -a.x(),
		-a.y(), a.x(), 0;

	double theta = omega.norm();

	R = Eigen::Matrix3d::Identity() + sin(theta) * cross_a + (1 - cos(theta)) * cross_a * cross_a;
	
}