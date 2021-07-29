#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {

    for (int i = 0; i < fg.rows(); i++)
        fg(i) = g(i % 3);

    //������˵��Ǹ��ļ�-M * fg;����ȷ��
    fg = M * fg;
}
