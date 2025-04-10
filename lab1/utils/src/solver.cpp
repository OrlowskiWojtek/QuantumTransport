#include "utils/include/solver.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/src/Core/Matrix.h>
#include <fstream>
#include <iostream>
#include <memory>

QDSolver::QDSolver()
    : p(QDParams::getInstance())
    , eigensolver(std::make_unique<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>>())
    , qd(std::make_unique<QDBasis>()) {}

void QDSolver::solve() {
    qd->calculate_indexes();

    Eigen::MatrixXd H = qd->get_H_matrix();
    Eigen::MatrixXd S = qd->get_S_matrix();

    eigensolver->compute(H, S);
}

void QDSolver::eigenvalues_to_file(int save_n, std::string filename) {
    std::ofstream file("../data/eigenvalues_" + filename);

    for (int i = 0; i < save_n; i++) {
        if (i >= eigensolver->eigenvalues().size()) {
            file << "Wrong size of requested values" << std::endl;
            break;
        }

        file << eigensolver->eigenvalues()(i) * p->hartree_to_meV << "\n";
    }

    file.close();
}

void QDSolver::wavefunctions_to_file(int save_n, std::string filename) {
    std::ofstream file("../data/wavefunctions_" + filename);
    
    for(int i = 0; i < save_n; i++){
        wavefunction_to_stream(i, file);
        file << "\n";
    }

    file.close();
}

void QDSolver::wavefunction_to_stream(int wf_number, std::ofstream& file){
    double step = (p->xmax - p->xmin) / 100;

    for (double x = p->xmin; x < p->xmax; x += step) {
        for (double y = p->xmin; y < p->xmax; y += step) {
            double psi = 0;
            for(int k = 0; k < p->N; k++){
                psi += eigensolver->eigenvectors().col(wf_number)(k) * qd->base_gaussian(x, y, k);
            }
            file << psi << "\t";
        }
        file << "\n";
    }
}
