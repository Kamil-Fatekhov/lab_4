#include <iostream>
#include <iomanip>
#include <vector>
#include "solver.h"
#include "grid.h"

void print_table_fragment(const SolverResult& res, const TaskParams& p) {
    Grid g(p);
    std::cout << "\n--- Table 4 Fragment (Values of v(n,m) in nodes) ---\n";
    std::cout << std::setw(10) << "y\\x";
    
    // Выводим заголовки x (каждый 10-й узел для краткости)
    int step_x = std::max(1, p.n / 5);
    int step_y = std::max(1, p.m / 5);

    for (int i = 0; i <= p.n; i += step_x) 
        std::cout << std::setw(12) << std::fixed << std::setprecision(3) << g.x(i);
    std::cout << "\n";

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
    p.a = 0.0; p.b = 2.0; p.c = 0.0; p.d = 1.0;
    
    // СОВЕТ: Попробуй поставить n=100 и m=100, чтобы увидеть как уменьшится eps2
    p.n = 100; 
    p.m = 100; 
    p.omega = 1.0; 
    p.eps_met = 1e-8; 
    p.Nmax = 100000;

    std::cout << "Starting Main Task (Seidel) - Variant 7\n";
    std::cout << "Target accuracy: 0.5e-6\n";

    // 1. Решение на основной сетке
    SolverResult res1 = solve_main_seidel(p);
    
    // 2. Решение на измельченной сетке
    TaskParams p2 = p;
    p2.n *= 2; p2.m *= 2;
    SolverResult res2 = solve_main_seidel(p2);

    // 3. Расчет точности eps2
    double eps2 = 0.0;
    int max_i = 0, max_j = 0;
    Grid g1(p);
    Grid g2(p2);
    for (int j = 0; j <= p.m; ++j) {
        for (int i = 0; i <= p.n; ++i) {
            double diff = std::abs(res1.values[g1.index(i, j)] - res2.values[g2.index(2 * i, 2 * j)]);
            if (diff > eps2) {
                eps2 = diff;
                max_i = i; max_j = j;
            }
        }
    }

    // ВЫВОД ДЛЯ ОТЧЕТА (Раздел V методички)
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "\n================ REPORT ================\n";
    std::cout << "Method: Seidel (SOR with omega=1)\n";
    std::cout << "Grid size: " << p.n << " x " << p.m << "\n";
    std::cout << "Iterations: " << res1.iterations_done << "\n";
    std::cout << "Residual (max norm): " << res1.final_residual << "\n";
    std::cout << "Accuracy eps2: " << eps2 << "\n";
    std::cout << "Max deviation at node: x=" << g1.x(max_i) << ", y=" << g1.y(max_j) << "\n";
    
    print_table_fragment(res1, p);

    if (eps2 <= 0.5e-6) std::cout << "\nSTATUS: SUCCESS\n";
    else std::cout << "\nSTATUS: FAILED (Need more nodes)\n";
    std::cout << "========================================\n";

    return 0;
}