#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "solver.h"
#include "grid.h"
#include "utils.h"

void print_table_fragment(const SolverResult& res, const TaskParams& p) {
    Grid g(p);
    std::cout << "\n--- Таблица 4 (Фрагмент значений численного решения) ---\n";
    std::cout << std::setw(10) << "y\\x";
    
    int step_x = std::max(1, p.n / 5);
    int step_y = std::max(1, p.m / 5);

    // Печать шапки x
    for (int i = 0; i <= p.n; i += step_x) 
        std::cout << std::setw(12) << std::fixed << std::setprecision(3) << g.x(i);
    std::cout << "\n";

    // Печать строк (идем от m вниз к 0 для наглядности декартовой системы)
    for (int j = p.m; j >= 0; j -= step_y) {
        std::cout << std::setw(10) << std::fixed << std::setprecision(3) << g.y(j);
        for (int i = 0; i <= p.n; i += step_x) {
            std::cout << std::setw(12) << std::fixed << std::setprecision(6) << res.values[g.index(i, j)];
        }
        std::cout << "\n";
    }
}

int main() {
    TaskParams p;
    // Геометрия области для 7-го варианта
    p.a = 0.0; p.b = 2.0; p.c = 0.0; p.d = 1.0;
    
    // Параметры основной сетки из английского мейна
    p.n = 50; 
    p.m = 50; 
    
    // Параметр МВР из английского мейна
    p.omega = 1.75; 
    p.eps_met = 1e-8; 
    p.Nmax = 60000;

    std::cout << "=== Основная задача (МВР) | Вариант 7 ===\n";
    std::cout << "Требуемая точность: 0.5e-6\n\n";

    // 1. Решение на основной сетке (n x m)
    SolverResult res1 = solve_main_sor(p);
    
    // 2. Решение на сгущенной сетке (2n x 2m)
    TaskParams p2 = p;
    p2.n *= 2; 
    p2.m *= 2;
    SolverResult res2 = solve_main_sor(p2);

    // 3. Вычисление точности eps2 по общим узлам
    double eps2 = 0.0;
    int max_i = 0, max_j = 0;
    Grid g1(p);
    Grid g2(p2);
    
    // Циклы идут по j, внутри по i — строго как в английском мейне
    for (int j = 0; j <= p.m; ++j) {
        for (int i = 0; i <= p.n; ++i) {
            double diff = std::abs(res1.values[g1.index(i, j)] - res2.values[g2.index(2 * i, 2 * j)]);
            if (diff > eps2) {
                eps2 = diff;
                max_i = i; max_j = j;
            }
        }
    }

    // 4. Вывод Справки для отчета (Порядок строк строго соответствует английскому REPORT)
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "---------------- СПРАВКА ----------------\n";
    std::cout << "Метод: Прегонной верхней релаксации (МВР)\n"; // Совпадает со строкой Method:
    std::cout << "Параметр omega: " << p.omega << "\n";            // Совпадает со строкой Omega parameter:
    std::cout << "Размер сетки: " << p.n << " x " << p.m << "\n";   // Совпадает со строкой Grid size:
    std::cout << "Затрачено итераций (N): " << res1.iterations_done << "\n"; // Совпадает со строкой Iterations spent (N):
    std::cout << "Невязка (макс. норма) ||R||: " << res1.final_residual << "\n"; // Совпадает со строкой Residual (max norm) ||R||:
    std::cout << "Затрачено итераций на измельченной сетке (N2): " << res2.iterations_done << "\n"; // Совпадает со строкой Iterations spent on dense grid (N2):
    std::cout << "Невязка на измельченной сетке ||R2||: " << res2.final_residual << "\n"; // Совпадает со строкой Residual on dense grid ||R2||:
    std::cout << "Рассчитанная точность eps2: " << eps2 << "\n"; // Совпадает со строкой Calculated accuracy eps2:
    std::cout << "Максимальное отклонение в узле: x=" << std::setprecision(3) << g1.x(max_i) << ", y=" << g1.y(max_j) << "\n"; // Совпадает со строкой Max deviation at node:
    
    print_table_fragment(res1, p);

    // Проверка на достижение требуемой точности
    if (eps2 <= 0.5e-6) {
        std::cout << "\nСТАТУС: УСПЕХ (Заданная точность достигнута)\n";
    } else {
        std::cout << "\nСТАТУС: ОШИБКА (Требуется больше узлов)\n";
    }
    std::cout << "-----------------------------------------\n";

    return 0;
}