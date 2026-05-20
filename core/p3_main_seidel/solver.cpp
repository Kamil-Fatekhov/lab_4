#include "solver.h"
#include "grid.h"
#include "utils.h"
#include <cmath>
#include <vector>
#include <algorithm>

namespace {
    // Константа PI
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

SolverResult solve_main_seidel(TaskParams params) {
    Grid grid(params);
    int N_total = (params.n + 1) * (params.m + 1);
    std::vector<double> v(N_total, 0.0);

    double h = grid.hx;
    double k = grid.hy;
    double h2inv = 1.0 / (h * h);
    double k2inv = 1.0 / (k * k);
    double common_denom = 1.0 / (2.0 * (h2inv + k2inv));

    // 1. Установка граничных условий и начального приближения
    // Используем линейную интерполяцию по x (согласно методичке)
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
    // Граничные условия по y (перекрывают углы, что корректно)
    for (int i = 0; i <= params.n; ++i) {
        v[grid.index(i, 0)] = mu3(grid.x(i));
        v[grid.index(i, params.m)] = mu4(grid.x(i));
    }

    // 2. Итерационный процесс (Метод Зейделя)
    int iter = 0;
    double max_diff = 0.0;

    while (iter < params.Nmax) {
        max_diff = 0.0;
        
        // Идем по внутренним узлам
        for (int j = 1; j < params.m; ++j) {
            for (int i = 1; i < params.n; ++i) {
                int idx = grid.index(i, j);
                double old_val = v[idx];

                // Формула Зейделя: используем новые значения i-1 и j-1 сразу
                double sum_x = (v[grid.index(i - 1, j)] + v[grid.index(i + 1, j)]) * h2inv;
                double sum_y = (v[grid.index(i, j - 1)] + v[grid.index(i, j + 1)]) * k2inv;
                double f_val = f_func(grid.x(i), grid.y(j));

                v[idx] = (sum_x + sum_y + f_val) * common_denom;

                double diff = std::abs(v[idx] - old_val);
                if (diff > max_diff) max_diff = diff;
            }
        }

        iter++;
        // Проверка критерия остановки (точность метода eps_met)
        if (max_diff < params.eps_met) break;
    }

    // 3. Расчет невязки (Residual) после выхода из цикла
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
    res.error = max_diff; // В данном контексте это точность выхода метода
    return res;
}