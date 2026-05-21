#ifndef SOLVER_H
#define SOLVER_H

#include "poisson_base.h"

// Функция решения основной задачи методом МВР (SOR)
SolverResult solve_main_sor(const TaskParams& p);

#endif // SOLVER_H