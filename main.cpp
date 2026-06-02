#include <iostream>
#include <iomanip>
#include <cmath>
#include "mesh.h"
#include "rwg.h"
#include "efie.h"
#include "cg_solver.h"
#include "farfield.h"
#include "fmm.h"

int main() {
    system("chcp 65001 > nul");  // 切换终端到 UTF-8，避免中文乱码
    std::cout << "=== 金属小球 RCS 计算（矩量法 + RWG + FMM）===\n";

    // ---- 参数设置 ----
    double radius = 0.5;          // 球半径 0.5 m
    double freq = 300e6;          // 频率 300 MHz
    double lambda = 3e8 / freq;   // 波长 1.0 m
    double k0 = 2.0 * PI / lambda; // 波数

    int refinement = 2;           // 剖分细化次数（2次后约320个三角形）
    int n_theta = 5;              // RCS 计算的 θ 采样点数（0~180°）
    bool use_fmm = true;          // 使用 FMM 加速
    double fmm_box = 0.3 * lambda;// FMM 盒子尺寸
    bool run_dense_compare = false;// 同时跑稠密法对比 (耗时较长)

    std::cout << "频率: " << freq/1e6 << " MHz\n";
    std::cout << "波长: " << lambda << " m\n";
    std::cout << "球半径: " << radius << " m (ka = " << k0*radius << ")\n";

    // 1. 生成网格
    Mesh mesh;
    mesh.generate_sphere(radius, refinement);

    // 2. 构建 RWG 基函数
    RWGData rwg;
    rwg.build(mesh);

    // 3. 构建 EFIE
    EFIE efie(mesh, rwg, k0);

    // 4. 入射方向（单站）
    double inc_theta = 45.0 * PI / 180.0;
    double inc_phi = 0.0;
    Vertex k_inc = {k0 * sin(inc_theta) * cos(inc_phi),
                    k0 * sin(inc_theta) * sin(inc_phi),
                    k0 * cos(inc_theta)};
    Vertex E0_pol = {0, 1, 0};  // φ 极化

    std::vector<Complex> b;
    efie.fill_vector(k_inc, E0_pol, b);

    // 5. 求解 Z I = b
    std::vector<Complex> I;

    // FMM 对象 (作用域需覆盖对比代码)
    FMM fmm(mesh, rwg, k0, fmm_box);

    if (use_fmm) {
        // ---- FMM 加速路径 ----
        std::cout << "\n--- FMM Setup ---\n";
        fmm.build_boxes();
        fmm.precompute_patterns();
        fmm.fill_nearfield(efie);
        fmm.print_stats();

        std::cout << "\n--- FMM-BiCGSTAB Solve ---\n";
        auto mv = fmm.get_matvec_functor();
        bool converged = CGSolver::solve_bicgstab_functor(mv, b, I, 1e-6, 1000);
        if (!converged) {
            std::cout << "警告: FMM 迭代未收敛\n";
        }
    } else {
        // ---- 稠密矩阵路径 ----
        std::vector<std::vector<Complex>> Z;
        efie.fill_matrix(Z);

        bool converged = CGSolver::solve_bicgstab(Z, b, I, 1e-6, 1000);
        if (!converged) {
            std::cout << "警告: 迭代未收敛\n";
        }
    }

    // 6. 计算双站 RCS
    FarField farfield(mesh, rwg, k0);

    std::cout << "\n双站 RCS (dBsm):\n";
    std::cout << "θ(°) | RCS(dBsm)\n";
    std::cout << "----------------\n";

    for (int i = 0; i <= n_theta; i++) {
        double theta = i * 180.0 / n_theta * PI / 180.0;
        double rcs = farfield.compute_rcs(theta, inc_phi, I);
        std::cout << std::setw(4) << i * 180.0 / n_theta
                  << " | " << std::setw(10) << rcs << std::endl;
    }

    // 可选: 稠密法对比
    if (use_fmm && run_dense_compare) {
        std::cout << "\n--- 稠密法对比 ---\n";
        std::vector<std::vector<Complex>> Z_dense;
        efie.fill_matrix(Z_dense);

        std::vector<Complex> I_dense;
        CGSolver::solve_bicgstab(Z_dense, b, I_dense, 1e-6, 1000);

        // 比较电流系数
        double diff_norm = 0, dense_norm = 0;
        for (int i = 0; i < (int)I.size(); i++) {
            Complex diff = I[i] - I_dense[i];
            diff_norm += std::norm(diff);
            dense_norm += std::norm(I_dense[i]);
        }
        double rel_err = std::sqrt(diff_norm / std::max(dense_norm, 1e-30));
        std::cout << "电流相对误差 ||I_fmm - I_dense|| / ||I_dense|| = "
                  << rel_err << "\n";

        // 直接验证 matvec: 用相同的 x=I_fmm, 比较 y=Z*x
        std::vector<Complex> y_fmm;
        fmm.matvec(I, y_fmm);
        std::vector<Complex> y_dense(I.size(), Complex(0,0));
        for (int i = 0; i < (int)I.size(); i++)
            for (int j = 0; j < (int)I.size(); j++)
                y_dense[i] += Z_dense[i][j] * I[j];
        double mv_diff = 0, mv_norm = 0;
        for (int i = 0; i < (int)I.size(); i++) {
            Complex d = y_fmm[i] - y_dense[i];
            mv_diff += std::norm(d);
            mv_norm += std::norm(y_dense[i]);
        }
        std::cout << "Matvec 相对误差 ||y_fmm - y_dense|| / ||y_dense|| = "
                  << std::sqrt(mv_diff / std::max(mv_norm, 1e-30)) << "\n";

        // 比较 RCS
        std::cout << "\nRCS 对比:\n";
        std::cout << "θ(°) | FMM(dBsm) | Dense(dBsm) | Δ(dB)\n";
        std::cout << "----------------------------------------\n";
        for (int i = 0; i <= n_theta; i++) {
            double theta = i * 180.0 / n_theta * PI / 180.0;
            double rcs_fmm = farfield.compute_rcs(theta, inc_phi, I);
            double rcs_dense = farfield.compute_rcs(theta, inc_phi, I_dense);
            std::cout << std::setw(4) << i * 180.0 / n_theta
                      << " | " << std::setw(10) << rcs_fmm
                      << " | " << std::setw(10) << rcs_dense
                      << " | " << std::setw(8) << (rcs_fmm - rcs_dense) << std::endl;
        }
    }

    system("pause");
    return 0;
}