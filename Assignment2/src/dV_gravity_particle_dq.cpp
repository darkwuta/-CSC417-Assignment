#include <dV_gravity_particle_dq.h>
#include<iostream>
void dV_gravity_particle_dq(Eigen::Ref<Eigen::Vector3d> f,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
	////--------------------------------------------
	//std::cout << "DV_GRAVITY_PARTICLE_DQ::DEBUG::g:" << g << std::endl;
	f += g;
	////--------------------------------------------

	//f = mass * -g;
}