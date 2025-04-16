#include "utils/include/task_manager.hpp"
#include <memory>

TaskManager::TaskManager()
    : p(QDParams::getInstance())
    , qd_basis(std::make_unique<QDBasis>())
    , qd_solver(std::make_unique<QDSolver>()) {
}

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
    qd_solver->eigenvalues_to_file(6, "task_3");
    qd_solver->wavefunctions_to_file(6, "task_3");
}

void TaskManager::task_4() {
    for (int omx = 20; omx <= 500; omx += 20) {
        p->set_omega_x(omx);
        qd_solver->solve();
        qd_solver->eigenvalues_to_file(10, "task_4_hwx_" + std::to_string(omx));
    }
}


void TaskManager::task_5() {
    // From analytical point of view
    // energies of excited states are equal: hwx(nx + 1/2) + hwy(ny + 1/2)
    // If all states from 0-4 must be excited in x direction then 4th x-excited states energy
    // must be smaller than 1st y-excited state
    // arbitrary chosen hwy = 350 satisfies that criteria, that first 4 states are excited in x direction
    // but 5th state is excited in y direction (first state is ground state)
    // for nx = 4, ny = 0 -> E = 4.5 * 80 + 0.5 * 350 = 535
    // for nx = 5, ny = 0 -> E = 5.5 * 80 + 0.5 * 350 = 615
    // for nx = 0, ny = 1 -> E = 0.5 * 80 + 1.5 * 350 = 565
    p->set_omega_x(80);
    p->set_omega_y(350);
    qd_solver->solve();
    qd_solver->eigenvalues_to_file(6, "task_5");
    qd_solver->wavefunctions_to_file(6, "task_5");
}
