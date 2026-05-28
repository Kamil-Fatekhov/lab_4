#define _USE_MATH_DEFINES
#include "httplib.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <algorithm>

#include "shared/poisson_base.h"

// Подключаем файлы всех участников
#include "core/p1_test_seidel/solver.h"
#include "core/p2_test_sor/solver.h"
#include "core/p3_main_seidel/solver.h"
#include "core/p4_main_sor/solver.h"
#include "core/p5_main_sor/solver.h"

using namespace httplib;
namespace fs = std::filesystem;

// ==========================================================
// МАТЕМАТИЧЕСКИЙ БЛОК ВАРИАНТА №7
// ==========================================================

double get_u_star(double x, double y) {
    return std::exp(std::pow(std::sin(M_PI * x * y), 2));
}

double get_f_star(double x, double y) {
    double pi = M_PI;
    double s = std::sin(pi * x * y);
    double c = std::cos(pi * x * y);
    double u = std::exp(s * s);
    double val = pi * pi * (x * x + y * y);
    return - 2.0 * u * val * (c * c - s * s + 2.0 * s * s * c * c);
}

double get_f_main(double x, double y) {
    return std::abs(x * x - 2.0 * y);
}

// Функции границ для создания начального приближения V0
double mu1(double y) { return std::pow(std::sin(M_PI * y), 2); }
double mu2(double y) { return std::pow(std::sin(2.0 * M_PI * y), 2); }
double mu3(double x) { return std::pow(std::sin(M_PI * x), 2); }
double mu4(double x) { return std::pow(std::sin(2.0 * M_PI * x), 2); }

// ==========================================================
// ВСПОМОГАТЕЛЬНЫЕ РАСЧЕТЫ
// ==========================================================

// Функция расчета невязки ||R|| (максимальная норма)
double calculate_residual_manual(const std::vector<double>& v, int n, int m, bool is_test) {
    double h = 2.0 / n; 
    double k = 1.0 / m; 
    double h2inv = 1.0 / (h * h);
    double k2inv = 1.0 / (k * k);
    double max_r = 0.0;

    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            double x = i * h;
            double y = j * k;
            double f = is_test ? get_f_star(x, y) : get_f_main(x, y);
            double laplace = (v[(i - 1) * (m + 1) + j] - 2 * v[i * (m + 1) + j] + v[(i + 1) * (m + 1) + j]) * h2inv +
                             (v[i * (m + 1) + (j - 1)] - 2 * v[i * (m + 1) + j] + v[i * (m + 1) + (j + 1)]) * k2inv;
            double r = std::abs(laplace + f);
            if (r > max_r) max_r = r;
        }
    }
    return max_r;
}

// Генерация начального приближения V0 (для расчета R0)
std::vector<double> generate_v0(int n, int m) {
    std::vector<double> v((n + 1) * (m + 1), 0.0);
    double h = 2.0 / n;
    double k = 1.0 / m;
    for (int j = 0; j <= m; ++j) {
        double y = j * k;
        v[0 * (m + 1) + j] = mu1(y);
        v[n * (m + 1) + j] = mu2(y);
        for (int i = 1; i < n; ++i) 
            v[i * (m + 1) + j] = mu1(y) + (i * h) * (mu2(y) - mu1(y)) / 2.0;
    }
    for (int i = 0; i <= n; ++i) {
        double x = i * h;
        v[i * (m + 1) + 0] = mu3(x);
        v[i * (m + 1) + m] = mu4(x);
    }
    return v;
}

double compute_optimal_omega(int n, int m) {
    double h = 2.0 / n;
    double k = 1.0 / m;
    double mu = (k*k * std::cos(M_PI / n) + h*h * std::cos(M_PI / m)) / (h*h + k*k);
    return 2.0 / (1.0 + std::sqrt(1.0 - mu * mu));
}

std::string build_json_response(
    int n, int m, int it1, int it2, 
    double r0, double rn, double r20, double r2n, 
    double final_accuracy, double omega,
    double eps_used, double eps_last, // Добавили эти два параметра
    const std::vector<double>& v, const std::vector<double>& u, const std::vector<double>& e,
    std::string task_type
) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(10);
    ss << "{\"n\":" << n << ",\"m\":" << m << ",\"it1\":" << it1 << ",\"it2\":" << it2 << ",";
    ss << "\"r0\":" << r0 << ",\"rn\":" << rn << ",\"r20\":" << r20 << ",\"r2n\":" << r2n << ",";
    ss << "\"acc\":" << final_accuracy << ",\"omega\":" << omega << ",";
    ss << "\"eps_used\":" << eps_used << ",\"eps_last\":" << eps_last << ","; // Передаем их в JSON
    ss << "\"type\":\"" << task_type << "\",";
    
    auto write_arr = [&](std::string key, const std::vector<double>& arr) {
        ss << "\"" << key << "\":[";
        for (size_t i = 0; i < arr.size(); ++i) ss << arr[i] << (i == arr.size() - 1 ? "" : ",");
        ss << "]";
    };
    write_arr("v", v); ss << ",";
    write_arr("u", u); ss << ",";
    write_arr("e", e);
    ss << "}";
    return ss.str();
}

// ==========================================================
// ОСНОВНОЙ ЦИКЛ СЕРВЕРА
// ==========================================================

int main() {
    Server svr;
    svr.Get("/", [](const Request&, Response& res) {
        std::string path = "../frontend/index.html";
        if (!fs::exists(path)) path = "frontend/index.html";
        std::ifstream ifs(path);
        if (ifs) res.set_content(std::string((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>()), "text/html; charset=utf-8");
    });

    svr.Get("/solve", [](const Request& req, Response& res) {
        int id = std::stoi(req.get_param_value("id"));
        int n = std::stoi(req.get_param_value("n")), m = std::stoi(req.get_param_value("m"));
        double eps = std::stod(req.get_param_value("eps"));
        
        double eps_actual = eps;
        double omega_val = (id == 1 || id == 3) ? 1.0 : compute_optimal_omega(n, m);
        
        TaskParams p; p.a=0; p.b=2; p.c=0; p.d=1; p.n=n; p.m=m; p.eps_met=eps_actual; p.omega=omega_val; p.Nmax=500000;

        std::vector<double> v_final, u_final, e_final;
        int it1 = 0, it2 = 0;
        double r0 = 0, rn = 0, r20 = 0, r2n = 0, final_acc = 0, last_step_diff = 0;

        if (id <= 2) { // --- ТЕСТОВАЯ ЗАДАЧА ---
            r0 = calculate_residual_manual(generate_v0(n, m), n, m, true);
            if (id == 1) { 
                PoissonSolver_p1 s(p, get_f_star, get_u_star); s.init_grid(); s.solve(); 
                v_final = s.get_v(); it1 = s.get_iterations(); last_step_diff = s.get_last_diff(); 
            } else { 
                PoissonSolver_p2 s(p, get_f_star, get_u_star); s.init_grid(); s.solve(); 
                v_final = s.get_v(); it1 = s.get_iterations(); last_step_diff = s.get_last_diff(); 
            }
            rn = calculate_residual_manual(v_final, n, m, true);

            for (int i = 0; i <= n; ++i) for (int j = 0; j <= m; ++j) {
                double exact = get_u_star(i * (2.0 / n), j * (1.0 / m));
                u_final.push_back(exact);
                double diff = std::abs(v_final[i * (m + 1) + j] - exact);
                e_final.push_back(diff);
                if (diff > final_acc) final_acc = diff;
            }
            res.set_content(build_json_response(n, m, it1, 0, r0, rn, 0, 0, final_acc, omega_val, eps_actual, last_step_diff, v_final, u_final, e_final, "test"), "application/json");
        } 
        else { // --- ОСНОВНАЯ ЗАДАЧА ---
            r0 = calculate_residual_manual(generate_v0(n, m), n, m, false);
            SolverResult resN, res2N;
            TaskParams p2 = p; p2.n *= 2; p2.m *= 2; p2.omega = (id == 3) ? 1.0 : compute_optimal_omega(p2.n, p2.m);

            if (id == 3) { resN = solve_main_seidel_p3(p); res2N = solve_main_seidel_p3(p2); }
            else if (id == 4) { resN = solve_main_sor_p4(p); res2N = solve_main_sor_p4(p2); }
            else { resN = solve_main_sor_p5(p); res2N = solve_main_sor_p5(p2); }

            v_final = resN.values; it1 = resN.iterations_done; it2 = res2N.iterations_done;
            rn = resN.final_residual;
            last_step_diff = resN.error; // ВОТ ТУТ МЫ БЕРЕМ РЕАЛЬНУЮ ТОЧНОСТЬ ШАГА

            r20 = calculate_residual_manual(generate_v0(p2.n, p2.m), p2.n, p2.m, false);
            r2n = res2N.final_residual;

            for (int i = 0; i <= n; ++i) for (int j = 0; j <= m; ++j) {
                double v2_val = res2N.values[(2 * i) * (2 * m + 1) + (2 * j)];
                u_final.push_back(v2_val);
                double diff = std::abs(v_final[i * (m + 1) + j] - v2_val);
                e_final.push_back(diff);
                if (diff > final_acc) final_acc = diff;
            }
            res.set_content(build_json_response(n, m, it1, it2, r0, rn, r20, r2n, final_acc, omega_val, eps_actual, last_step_diff, v_final, u_final, e_final, "main"), "application/json");
        }
    });

    std::cout << "SERVER: http://localhost:8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
    return 0;
}