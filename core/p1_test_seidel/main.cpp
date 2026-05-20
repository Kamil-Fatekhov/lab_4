#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <iomanip>
#include "solver.h"

// Точное решение u*(x, y) для 7 варианта
double u_exact_func(double x, double y) {
    double sin_val = std::sin(M_PI * x * y);
    return std::exp(sin_val * sin_val);
}

// Правая часть f*(x,y) = -Delta u
double f_star_func(double x, double y) {
    double pi_xy = M_PI * x * y;
    double sin2 = std::sin(2.0 * pi_xy);
    double cos2 = std::cos(2.0 * pi_xy);
    double term = sin2 * sin2 + 2.0 * cos2;
    
    double u_val = u_exact_func(x, y);
    double delta_u = M_PI * M_PI * (x * x + y * y) * u_val * term;
    
    return -delta_u;
}

int main() {
    // Настройки задачи согласно методичке для тестовой задачи
    TaskConfig cfg;
    cfg.a = 0.0; cfg.b = 2.0;
    cfg.c = 0.0; cfg.d = 1.0;
    cfg.n = 100; // Разбиения по x
    cfg.m = 50;  // Разбиения по y
    cfg.epsilon = 0.5e-6; // Заданная погрешность по методичке
    cfg.max_iter = 100000;
    cfg.omega = 1.0; // 1.0 = Метод Зейделя. Для МВР нужно подбирать (напр. 1.5 - 1.8)

    std::cout << "--- Решение тестовой задачи (Red-Black TBB) ---" << std::endl;
    std::cout << "Сетка: " << cfg.n << " x " << cfg.m << std::endl;
    std::cout << "Параметр релаксации (omega): " << cfg.omega << std::endl;

    PoissonSolver solver(cfg, f_star_func, u_exact_func);
    
    solver.init_grid();
    
    // Запускаем счет
    solver.solve();

    std::cout << "\nСправка для тестовой задачи:" << std::endl;
    std::cout << "Затрачено итераций N: " << solver.get_iterations() << std::endl;
    std::cout << "Достигнута точность итерационного метода eps^(N): " 
              << std::scientific << solver.get_last_diff() << std::endl;
              
    // Рассчитываем и выводим погрешность (эпсилон 1 из методички)
    double eps1 = solver.calculate_max_error();
    std::cout << "Погрешность решения eps_1 = max|u* - v|: " 
              << std::scientific << eps1 << std::endl;

    if (eps1 <= 0.5e-6) {
        std::cout << "[УСПЕХ] Заданная погрешность достигнута!" << std::endl;
    } else {
        std::cout << "[ИНФО] Погрешность выше требуемой. Для ее уменьшения "
                  << "необходимо увеличить число разбиений сетки (n, m)." << std::endl;
    }

    return 0;
}