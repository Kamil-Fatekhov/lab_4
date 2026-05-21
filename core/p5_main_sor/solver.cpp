#include "solver.h"
#include "grid.h"
#include "utils.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace {
    const double PI = std::acos(-1.0);

    // Правая часть f(x, y) для 7-го варианта
    double f_func(double x, double y) {
        return std::abs(x * x - 2.0 * y);
    }

    // Граничные условия для 7-го варианта (Основная задача)
    double mu1(double y) { return std::pow(std::sin(PI * y), 2); }
    double mu2(double y) { return std::pow(std::sin(2.0 * PI * y), 2); }
    double mu3(double x) { return std::pow(std::sin(PI * x), 2); }
    double mu4(double x) { return std::pow(std::sin(2.0 * PI * x), 2); }
}

SolverResult solve_main_sor(TaskParams params) {
    Grid grid(params);
    int N_total = (params.n + 1) * (params.m + 1);
    std::vector<double> v(N_total, 0.0);

    double h2inv = 1.0 / (grid.hx * grid.hx);
    double k2inv = 1.0 / (grid.hy * grid.hy);
    double common_denom = 1.0 / (2.0 * (h2inv + k2inv));
    double omega = params.omega;

    // 1. Установка граничных условий и начального приближения
    // Линейная интерполяция по направлению X (согласно методичке)
    for (int j = 0; j <= params.m; ++j) {
        double y_val = grid.y(j);
        double left = mu1(y_val);
        double right = mu2(y_val);
        
        v[grid.index(0, j)] = left;
        v[grid.index(params.n, j)] = right;

        for (int i = 1; i < params.n; ++i) {
            v[grid.index(i, j)] = Utils::interpolate_linear(grid.x(i), params.a, params.b, left, right);
        }
    }
    
    // Граничные условия по Y (перекрывают углы)
    for (int i = 0; i <= params.n; ++i) {
        v[grid.index(i, 0)] = mu3(grid.x(i));
        v[grid.index(i, params.m)] = mu4(grid.x(i));
    }

    // 2. Итерационный процесс МВР (Successive Over-Relaxation)
    int iter = 0;
    double max_diff = 0.0;

    while (iter < params.Nmax) {
        max_diff = 0.0;
        
        // Проход по внутренним узлам сетки
        for (int j = 1; j < params.m; ++j) {
            for (int i = 1; i < params.n; ++i) {
                int idx = grid.index(i, j);
                double old_val = v[idx];

                // Компоненты разностного оператора Лапласа
                double sum_x = (v[grid.index(i - 1, j)] + v[grid.index(i + 1, j)]) * h2inv;
                double sum_y = (v[grid.index(i, j - 1)] + v[grid.index(i, j + 1)]) * k2inv;
                double f_val = f_func(grid.x(i), grid.y(j));

                // Новое значение по Зейделю (с учетом знака Delta u = -f)
                double v_seidel = (sum_x + sum_y + f_val) * common_denom;
                
                // Взвешенное обновление МВР с параметром omega
                double new_val = (1.0 - omega) * old_val + omega * v_seidel;

                v[idx] = new_val;

                double diff = std::abs(new_val - old_val);
                if (diff > max_diff) max_diff = diff;
            }
        }

        iter++;
        // Критерий остановки по точности изменения решения на шаге
        if (max_diff < params.eps_met) break;
    }

    // 3. Расчет итоговой нормы невязки СЛАУ (максимальная норма)
    double max_residual = 0.0;
    for (int j = 1; j < params.m; ++j) {
        for (int i = 1; i < params.n; ++i) {
            double laplace = (v[grid.index(i-1, j)] - 2*v[grid.index(i,j)] + v[grid.index(i+1, j)]) * h2inv +
                             (v[grid.index(i, j-1)] - 2*v[grid.index(i,j)] + v[grid.index(i, j+1)]) * k2inv;
            double r = std::abs(laplace + f_func(grid.x(i), grid.y(j)));
            if (r > max_residual) max_residual = r;
        }
    }

    SolverResult res;
    res.values = v;
    res.iterations_done = iter;
    res.final_residual = max_residual;
    res.error = max_diff; 
    return res;
}