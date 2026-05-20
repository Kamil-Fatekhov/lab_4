#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <functional>

struct TaskConfig {
    double a, b, c, d;
    int n, m;         // Количество разбиений по x и y
    double epsilon;   // Точность для остановки (eps_Mem)
    int max_iter;     // Максимальное число итераций
    double omega;     // Параметр релаксации (1.0 = Зейдель)
};

class PoissonSolver {
public:
    PoissonSolver(const TaskConfig& config,
                  std::function<double(double, double)> f_star,
                  std::function<double(double, double)> u_exact);

    // Инициализация начального приближения (здесь просто нули + точные границы)
    void init_grid();

    // Решение СЛАУ методом Red-Black SOR с использованием TBB
    void solve();

    // Расчет максимальной погрешности по сравнению с точным решением
    double calculate_max_error() const;

    // Геттеры для вывода
    int get_iterations() const { return iterations_done; }
    double get_last_diff() const { return last_diff; }

private:
    TaskConfig cfg;
    std::vector<double> v; // Одномерный массив для сетки (n+1)*(m+1)
    std::vector<double> f; // Правая часть f*
    
    std::function<double(double, double)> f_star;
    std::function<double(double, double)> u_exact;

    double hx, hy;
    int iterations_done;
    double last_diff;

    // Конвертация 2D индексов в 1D
    inline int idx(int i, int j) const {
        return i * (cfg.m + 1) + j;
    }
};

#endif // SOLVER_H