#ifndef TASK_MANAGER_HPP
#define TASK_MANAGER_HPP

#include "utils/include/solver.hpp"
#include "utils/include/params.hpp"
#include <memory>

class TaskManager{
public:
    TaskManager();
    void task_1();
    void task_2();
    void task_3();
    void task_4();

private:
    QDParams* p;
    std::unique_ptr<QDBasis> qd_basis;
    std::unique_ptr<QDSolver> qd_solver;

};


#endif
