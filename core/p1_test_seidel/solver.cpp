#include "solver.h"
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range2d.h>
#include <cmath>
#include <algorithm>
#include <iostream>

PoissonSolver::PoissonSolver(const TaskConfig& config,
                             std::function<double(double, double)> f_star,
                             std::function<double(double, double)> u_exact)
    : cfg(config), f_star(f_star), u_exact(u_exact), iterations_done(0), last_diff(0.0) 
{
    hx = (cfg.b - cfg.a) / cfg.n;
    hy = (cfg.d - cfg.c) / cfg.m;
    
    int total_nodes = (cfg.n + 1) * (cfg.m + 1);
    v.assign(total_nodes, 0.0);
    f.assign(total_nodes, 0.0);

    for (int i = 0; i <= cfg.n; ++i) {
        for (int j = 0; j <= cfg.m; ++j) {
            double x = cfg.a + i * hx;
            double y = cfg.c + j * hy;
            f[idx(i, j)] = f_star(x, y);
        }
    }
}

void PoissonSolver::init_grid() {
    for (int i = 0; i <= cfg.n; ++i) {
        double x = cfg.a + i * hx;
        v[idx(i, 0)] = u_exact(x, cfg.c);
        v[idx(i, cfg.m)] = u_exact(x, cfg.d);
    }
    for (int j = 0; j <= cfg.m; ++j) {
        double y = cfg.c + j * hy;
        v[idx(0, j)] = u_exact(cfg.a, y);       // Левая граница
        v[idx(cfg.n, j)] = u_exact(cfg.b, y);   // Правая граница
    }
    
    // Внутри сетки оставляем 0.0 (или можно реализовать линейную интерполяцию)
}

void PoissonSolver::solve() {
    double hx2 = hx * hx;
    double hy2 = hy * hy;
    double denom = 2.0 * (1.0 / hx2 + 1.0 / hy2);

    for (iterations_done = 0; iterations_done < cfg.max_iter; ++iterations_done) {
        double current_max_diff = 0.0;

        // Две фазы: phase=0 (красные узлы), phase=1 (черные узлы)
        for (int phase = 0; phase < 2; ++phase) {
            double phase_diff = tbb::parallel_reduce(
                tbb::blocked_range2d<int>(1, cfg.n, 1, cfg.m),
                0.0,
                [&](const tbb::blocked_range2d<int>& r, double local_max_diff) {
                    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                        for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                            if ((i + j) % 2 == phase) {
                                int current_idx = idx(i, j);
                                double v_old = v[current_idx];
                                
                                // Шаблон крест
                                double sum = (v[idx(i - 1, j)] + v[idx(i + 1, j)]) / hx2 +
                                             (v[idx(i, j - 1)] + v[idx(i, j + 1)]) / hy2;
                                             
                                double v_new = (sum + f[current_idx]) / denom;
                                
                                // Верхняя релаксация (если omega=1, это чистый Зейдель)
                                v[current_idx] = v_old + cfg.omega * (v_new - v_old);
                                
                                local_max_diff = std::max(local_max_diff, std::abs(v[current_idx] - v_old));
                            }
                        }
                    }
                    return local_max_diff;
                },
                [](double x, double y) { return std::max(x, y); }
            );
            current_max_diff = std::max(current_max_diff, phase_diff);
        }

        last_diff = current_max_diff;

        // Проверка сходимости
        if (current_max_diff < cfg.epsilon) {
            iterations_done++; // Учитываем последнюю итерацию
            break;
        }
    }
}

double PoissonSolver::calculate_max_error() const {
    double max_err = 0.0;
    for (int i = 0; i <= cfg.n; ++i) {
        for (int j = 0; j <= cfg.m; ++j) {
            double x = cfg.a + i * hx;
            double y = cfg.c + j * hy;
            double exact = u_exact(x, y);
            max_err = std::max(max_err, std::abs(v[idx(i, j)] - exact));
        }
    }
    return max_err;
}