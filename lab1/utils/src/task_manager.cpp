#include "utils/include/task_manager.hpp"
#include <memory>

TaskManager::TaskManager()
    : qd_basis(std::make_unique<QDBasis>())
    , qd_solver(std::make_unique<QDSolver>()) {}

void TaskManager::task_1() {
    qd_basis->base_gaussian_to_file(0);
    qd_basis->base_gaussian_to_file(8);
    qd_basis->base_gaussian_to_file(9);
}

void TaskManager::task_2() {
    qd_solver->solve();
}

void TaskManager::task_3() {
    qd_solver->solve();
    qd_solver->eigenvalues_to_file(6, "test");
    qd_solver->wavefunctions_to_file(6, "test");
}
