#define _USE_MATH_DEFINES
#include "solver.h"
#include "grid.h"
#include "utils.h"
#include <cmath>
#include <algorithm>

// --- Функции для Варианта 7 ---
static double f_func(double x, double y) {
    return std::abs(x * x - 2.0 * y);
}

static double mu1(double y) {
    return std::pow(std::sin(M_PI * y), 2);
}

static double mu2(double y) {
    return std::pow(std::sin(2.0 * M_PI * y), 2);
}

static double mu3(double x) {
    return std::pow(std::sin(M_PI * x), 2);
}

static double mu4(double x) {
    return std::pow(std::sin(2.0 * M_PI * x), 2);
}

// ПЕРЕИМЕНОВАНО: Добавлен суффикс _p4
SolverResult solve_main_sor_p4(const TaskParams& p) {
    Grid g(p);
    SolverResult res;
    
    res.values.assign((p.n + 1) * (p.m + 1), 0.0);
    
    // 1. Установка граничных условий
    for (int j = 0; j <= p.m; ++j) {
        res.values[g.index(0, j)] = mu1(g.y(j));
        res.values[g.index(p.n, j)] = mu2(g.y(j));
    }
    for (int i = 0; i <= p.n; ++i) {
        res.values[g.index(i, 0)] = mu3(g.x(i));
        res.values[g.index(i, p.m)] = mu4(g.x(i));
    }

    // 2. Начальное приближение
    for (int j = 1; j < p.m; ++j) {
        double u_a = mu1(g.y(j));
        double u_b = mu2(g.y(j));
        for (int i = 1; i < p.n; ++i) {
            res.values[g.index(i, j)] = Utils::interpolate_linear(g.x(i), p.a, p.b, u_a, u_b);
        }
    }

    // 3. Коэффициенты разностной схемы
    double h2 = g.hx * g.hx;
    double k2 = g.hy * g.hy;
    double A = 1.0 / h2;
    double B = 1.0 / k2;
    double C = 2.0 * (A + B);

    // 4. Итерационный процесс МВР (SOR)
    int iter = 0;
    double max_diff = 0.0;
    
    do {
        max_diff = 0.0;
        for (int i = 1; i < p.n; ++i) {
            double x = g.x(i);
            for (int j = 1; j < p.m; ++j) {
                double y = g.y(j);
                int idx = g.index(i, j);
                
                double v_old = res.values[idx];
                
                double sum_x = res.values[g.index(i - 1, j)] + res.values[g.index(i + 1, j)];
                double sum_y = res.values[g.index(i, j - 1)] + res.values[g.index(i, j + 1)];
                double v_gs = (A * sum_x + B * sum_y + f_func(x, y)) / C;
                
                res.values[idx] = (1.0 - p.omega) * v_old + p.omega * v_gs;
                
                double diff = std::abs(res.values[idx] - v_old);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
        }
        iter++;
    } while (max_diff > p.eps_met && iter < p.Nmax);

    res.iterations_done = iter;
    res.error = max_diff;

    // 5. Расчет нормы невязки
    double max_res = 0.0;
    for (int i = 1; i < p.n; ++i) {
        double x = g.x(i);
        for (int j = 1; j < p.m; ++j) {
            double y = g.y(j);
            double sum_x = res.values[g.index(i - 1, j)] + res.values[g.index(i + 1, j)];
            double sum_y = res.values[g.index(i, j - 1)] + res.values[g.index(i, j + 1)];
            double current = res.values[g.index(i, j)];
            
            double r_ij = std::abs((A * sum_x + B * sum_y - C * current) + f_func(x, y));
            if (r_ij > max_res) {
                max_res = r_ij;
            }
        }
    }
    res.final_residual = max_res;

    return res;
}