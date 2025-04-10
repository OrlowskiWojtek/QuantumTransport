#ifndef QD_SOLVER_HPP
#define QD_SOLVER_HPP

#include "basis_funcs.hpp"
#include "params.hpp"
#include <memory>

class QDSolver{
public:
    void solve();
    QDSolver();

private:
    std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>> eigensolver; 
    std::unique_ptr<QDBasis> qd;
    QDParams* p;

    void eigenvalues_to_file();
};

#endif
