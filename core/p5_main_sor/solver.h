#ifndef SOLVER_H
#define SOLVER_H

#include "poisson_base.h"

// Функция решения основной задачи методом верхней релаксации (МВР)
SolverResult solve_main_sor(TaskParams params);

#endif // SOLVER_H