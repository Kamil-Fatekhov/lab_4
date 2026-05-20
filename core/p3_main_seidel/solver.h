#ifndef SOLVER_H
#define SOLVER_H

#include "poisson_base.h"

// Функция для решения основной задачи методом Зейделя
// На вход принимает параметры задачи, возвращает структуру с результатом
SolverResult solve_main_seidel(TaskParams params);

#endif // SOLVER_H