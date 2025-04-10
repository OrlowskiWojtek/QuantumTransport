#include "utils/include/solver.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/src/Core/Matrix.h>
#include <iostream>
#include <memory>

QDSolver::QDSolver()
    : p(QDParams::getInstance())
    , eigensolver(std::make_unique<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>>())
    , qd(std::make_unique<QDBasis>()){}

void QDSolver::solve() {
    Eigen::MatrixXd H = qd->get_H_matrix();
    Eigen::MatrixXd S = qd->get_S_matrix();

    eigensolver->compute(H, S);
}

void QDSolver::eigenvalues_to_file(){


}
