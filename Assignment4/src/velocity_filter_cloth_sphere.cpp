#include <velocity_filter_cloth_sphere.h>

void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    int count = 0;
    for (unsigned int i : indices)
    {
        Eigen::Vector3d v = qdot.segment(i * 3, 3);
        Eigen::Vector3d n = normals[count++];
        double alpha = 0; 
        if (n.dot(v) >= 0)
            alpha = 0;
        else
            alpha = -n.dot(v);
        Eigen::Vector3d v_filtered = v + alpha * n;
        qdot.segment(i * 3, 3) = v_filtered;
    }
    ///////////////////////////////////////////////////////////
    //for (int v = 0; v < indices.size(); v++) {

    //    unsigned int vert_idx = indices[v];
    //    Eigen::Vector3d vert_dot = qdot.segment(3 * vert_idx, 3);
    //    Eigen::Vector3d n = normals[v];

    //    double alpha;
    //    // check if vertex volocity moves inward object while collision
    //    if (n.transpose() * vert_dot < 0) {
    //        alpha = -1.0 * n.transpose() * vert_dot;
    //    }

    //    Eigen::Vector3d v_filtered = vert_dot + alpha * n;

    //    qdot.segment(3 * vert_idx, 3) = v_filtered;

    //}

}