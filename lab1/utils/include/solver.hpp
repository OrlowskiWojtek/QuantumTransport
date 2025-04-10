#ifndef QD_SOLVER_HPP
#define QD_SOLVER_HPP

#include "basis_funcs.hpp"
#include "params.hpp"
#include <memory>

class QDSolver{
public:
    void solve();
    QDSolver();

    // saving 'n' first eigenvalues and wavefunctions to file
    void eigenvalues_to_file(int save_n, std::string filename);
    void wavefunctions_to_file(int save_n, std::string filename);
private:
    std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>> eigensolver; 
    std::unique_ptr<QDBasis> qd;
    QDParams* p;

    void wavefunction_to_stream(int wf_number, std::ofstream& file);
};

#endif
