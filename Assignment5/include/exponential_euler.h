#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces) {
	for (int i = 0; i < q.rows() / 12; i++) {
		Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12.0 * i).data());

		Eigen::Matrix3d inertiaM, expR;

		inertiaM.setZero();
		expR.setZero();

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				inertiaM(m, n) = masses.at(i)(m, n);
			}
		}

		Eigen::Vector3d p = Eigen::Map<const Eigen::Vector3d>(q.segment<3>(12 * i + 9.0).data());

		Eigen::Vector3d w = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i).data());
		Eigen::Vector3d pdot = Eigen::Map<const Eigen::Vector3d>(qdot.segment<3>(6.0 * i + 3.0).data());
		Eigen::Vector3d text = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6.0 * i).data());
		Eigen::Vector3d fext = Eigen::Map<const Eigen::Vector3d>(forces.segment<3>(6.0 * i + 3.0).data());

		rodrigues(expR, dt * w);

		Eigen::SparseMatrixd A;
		Eigen::Vector3d b, x, k;
		k = (R * inertiaM * R.transpose()) * w;
		k = (dt * w).cross(k);
		A = (R * inertiaM * R.transpose()).sparseView();
		b = (R * inertiaM * R.transpose()) * w + k + dt * text;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		solver.compute(A);
		x = solver.solve(b);
		w = x;
		pdot = pdot + (1.0 / masses.at(i)(5, 5)) * dt * fext;

		R = expR * R;

		p = p + dt * pdot;


		qdot[6.0 * i] = w[0];
		qdot[6.0 * i + 1.0] = w[1];
		qdot[6.0 * i + 2.0] = w[2];

		qdot[6.0 * i + 3.0] = pdot[0];
		qdot[6.0 * i + 4.0] = pdot[1];
		qdot[6.0 * i + 5.0] = pdot[2];

		for (int m = 0; m < 3; m++) {
			for (int n = 0; n < 3; n++) {
				q[12 * i + 3.0 * m + n] = R(n, m);
			}
		}
		q[12.0 * i + 9.0] = p[0];
		q[12.0 * i + 10.0] = p[1];
		q[12.0 * i + 11.0] = p[2];
	}
    //for (int i = 0; i < q.rows() / 12;i++)
    //{
    //    Eigen::Matrix3d R = Eigen::Map<const Eigen::Matrix3d>(q.segment<9>(12 * i).data());
    //    Eigen::Vector3d omega = qdot.segment<3>(6*i);
    //    Eigen::Vector3d p = q.segment<3>(12 * i + 9);
    //    Eigen::Vector3d pdot = qdot.segment<3>(6 * i + 3);
    //    Eigen::Vector3d T_ext = forces.segment<3>(6 * i);
    //    Eigen::Vector3d F_ext = forces.segment<3>(6 * i + 3);

    //    Eigen::Matrix3d R_next;

    //    Eigen::Matrix3d inertia_matrix = masses[i].block<3,3>(0,0);

    //    rodrigues(R_next, dt * omega);

    //    Eigen::SparseMatrixd A;
    //    Eigen::Vector3d B;
    //    Eigen::Vector3d x;
    //    A = (R * inertia_matrix * R.transpose()).sparseView();
    //    B = R * inertia_matrix * R.transpose() * omega + dt * omega.cross((R * inertia_matrix * R.transpose()) * omega) + dt * T_ext;

    //    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    //    solver.compute(A);
    //    x = solver.solve(B);

    //    R = R_next * R;
    //    omega = x;
    //    pdot = pdot + (1.0 / masses.at(i)(5, 5)) * dt * F_ext;

    //    p = p + dt * pdot;
    //    q.segment<9>(12*i) = Eigen::Map<const Eigen::Vector9d>(R.data());
    //    q.segment<3>(12 * i + 9) = p;
    //    qdot.segment<3>(6 * i) = omega;
    //    qdot.segment<3>(6 * i + 3) = pdot;
    //}
    
}