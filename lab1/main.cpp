#include "utils/include/basis_funcs.hpp"
#include "utils/include/solver.hpp"
#include <iostream>
#include <memory>

int main(){
    // task_one
    QDBasis t1_tests;
    t1_tests.base_gaussian_to_file(0);
    t1_tests.base_gaussian_to_file(8);
    t1_tests.base_gaussian_to_file(9);


    std::unique_ptr<QDSolver> solver = std::make_unique<QDSolver>();
    solver->solve();
}
