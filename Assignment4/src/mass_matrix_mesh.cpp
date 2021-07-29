#include <mass_matrix_mesh.h>

void mass_triangle_mesh(Eigen::Matrix99d m, double density, double area)
{
    // TODO ¥Ê“…
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
                m.block<3, 3>(i * 3, j * 3) = (1.0 / 12.0) * Eigen::Matrix3d::Identity();
            else
                m.block<3, 3>(i * 3, j * 3) = (1.0 / 24.0) * Eigen::Matrix3d::Identity();
        }
    }
    m = 2.0 * density * area * m;
}

void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    M.resize(q.rows(), q.rows());
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int tri = 0; tri < F.rows(); tri++)
    {
        Eigen::Matrix99d m;
        mass_triangle_mesh(m, density, areas[tri]);
        Eigen::RowVector3i element= F.row(tri);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int x = element(i);
                int y = element(j);

                for (int k = 0; k < 3; k++)
                {
                    tripletList.push_back({ x * 3 + k, y * 3 + k, m(i * 3 + k, j * 3 + k) });
                }
            }
        }
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());
   
}
 