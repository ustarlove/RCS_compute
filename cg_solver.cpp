#include "cg_solver.h"
#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// 共轭梯度法 — 仅适用于 Hermitian 正定矩阵
//tol：精度控制 max_iter：迭代次数上限
// ---------------------------------------------------------------------------
bool CGSolver::solve(const std::vector<std::vector<Complex>>& Z,
                     const std::vector<Complex>& b,
                     std::vector<Complex>& x,
                     double tol,
                     int max_iter) {
    int n = (int)Z.size();
    x.assign(n, Complex(0,0)); //初始解为0
    std::vector<Complex> r = b;  //初始残差
    std::vector<Complex> p = r;  //初始搜索方向

    double rsold = 0;  //上一步残差模长的平方
    for (int i = 0; i < n; i++) rsold += std::norm(r[i]);
    double bnorm = 0;
    for (int i = 0; i < n; i++) bnorm += std::norm(b[i]);
    double tol2 = tol * tol * bnorm;  //平方形式的收敛阈值
    if (bnorm < 1e-30) { std::cout << "CG: zero RHS\n"; return true; } //如果b为0向量，直接返回0解

    for (int iter = 0; iter < max_iter; iter++) {
        std::vector<Complex> Zp(n, Complex(0,0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Zp[i] += Z[i][j] * p[j];

        Complex p_Zp(0,0);
        for (int i = 0; i < n; i++) p_Zp += std::conj(p[i]) * Zp[i];
        double alpha = rsold / p_Zp.real();

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Zp[i];
        }

        double rsnew = 0;
        for (int i = 0; i < n; i++) rsnew += std::norm(r[i]);
        if (rsnew < tol2) {
            std::cout << "CG converged after " << iter+1 << " iterations\n";
            return true;
        }
        double beta = rsnew / rsold;
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rsold = rsnew;
    }
    std::cout << "CG did not converge after " << max_iter << " iterations\n";
    return false;
}

// ---------------------------------------------------------------------------
// BiCGSTAB — 适用于复非对称矩阵（如 EFIE 阻抗阵） 双共轭梯度稳定法
// 算法基于 van der Vorst (1992)
// ---------------------------------------------------------------------------
bool CGSolver::solve_bicgstab(const std::vector<std::vector<Complex>>& Z,
                               const std::vector<Complex>& b,
                               std::vector<Complex>& x,
                               double tol,
                               int max_iter) {
    int n = (int)Z.size();
    x.assign(n, Complex(0,0));

    // 矩阵-向量乘: y = Z * v
    auto matvec = [&](const std::vector<Complex>& v, std::vector<Complex>& y) {
        y.assign(n, Complex(0,0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                y[i] += Z[i][j] * v[j];
    };

    // 向量点积 (a^H * b)
    auto dotc = [&](const std::vector<Complex>& a, const std::vector<Complex>& b) {
        Complex s(0,0);
        for (int i = 0; i < n; i++) s += std::conj(a[i]) * b[i];
        return s;
    };

    // 计算 b 的范数用于收敛判断
    double bnorm = 0;
    for (int i = 0; i < n; i++) bnorm += std::norm(b[i]);
    bnorm = std::sqrt(bnorm);
    double stop_tol = tol * std::max(bnorm, 1e-30); //防止除0

    // r0 = b - A x0 (= b since x0=0)
    std::vector<Complex> r = b;     // 残差
    std::vector<Complex> r_hat = r; // r̂0（通常取为 r0）
    std::vector<Complex> p(n, Complex(0,0));
    std::vector<Complex> v(n, Complex(0,0));
    std::vector<Complex> s(n, Complex(0,0));
    std::vector<Complex> t(n, Complex(0,0));

    Complex rho_old(1, 0);
    Complex alpha(1, 0);
    Complex omega(1, 0);

    for (int iter = 0; iter < max_iter; iter++) {
        Complex rho = dotc(r_hat, r);
        //防止除0
        if (std::abs(rho) < 1e-30) {
            std::cout << "BiCGSTAB breakdown (rho≈0) at iter " << iter << "\n";
            return false;
        }

        Complex beta(1, 0);
        if (iter > 0) {
            beta = (rho / rho_old) * (alpha / omega);
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]); //更新搜索方向，第一次迭代p=r0    
        }

        matvec(p, v);                               // v = Z p
        Complex rhat_v = dotc(r_hat, v);
        if (std::abs(rhat_v) < 1e-30) {
            std::cout << "BiCGSTAB breakdown (r̂·v≈0) at iter " << iter << "\n";
            return false;
        }
        alpha = rho / rhat_v;  //计算步长

        for (int i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];             // 计算中间残差，s = r - α v
        }

        matvec(s, t);                               // t = Z s
        Complex ts = dotc(t, s);
        Complex tt = dotc(t, t);
        if (std::abs(tt) < 1e-30) {
            std::cout << "BiCGSTAB breakdown (t·t≈0) at iter " << iter << "\n";
            return false;
        }
        omega = ts / tt;    //计算稳定化步长

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i] + omega * s[i];    
            r[i] = s[i] - omega * t[i];              // 更新解和残差
        }

        // 收敛检查
        double r_norm = 0;
        for (int i = 0; i < n; i++) r_norm += std::norm(r[i]);
        r_norm = std::sqrt(r_norm);
        if (r_norm < stop_tol) {
            std::cout << "BiCGSTAB converged after " << iter+1 << " iterations, "
                      << "|r|/|b| = " << r_norm / std::max(bnorm, 1e-30) << "\n";
            return true;
        }

        rho_old = rho;
    }

    // 最终残差
    std::vector<Complex> Ax(n, Complex(0,0));
    matvec(x, Ax);
    double final_res = 0;
    for (int i = 0; i < n; i++) {
        Complex diff = b[i] - Ax[i];
        final_res += std::norm(diff);
    }
    std::cout << "BiCGSTAB did not converge after " << max_iter
              << " iterations, final |r|/|b| = "
              << std::sqrt(final_res) / std::max(bnorm, 1e-30) << "\n";
    return false;
}

// ---------------------------------------------------------------------------
// FMM 加速版 BiCGSTAB — 用 matvec 回调替代稠密 Z 矩阵
// ---------------------------------------------------------------------------
bool CGSolver::solve_bicgstab_functor(
    const std::function<void(const std::vector<Complex>&,
                              std::vector<Complex>&)>& matvec,
    const std::vector<Complex>& b,
    std::vector<Complex>& x,
    double tol,
    int max_iter) {
    int n = (int)b.size();
    x.assign(n, Complex(0, 0));

    auto dotc = [&](const std::vector<Complex>& a, const std::vector<Complex>& b_vec) {
        Complex s(0, 0);
        for (int i = 0; i < n; i++) s += std::conj(a[i]) * b_vec[i];
        return s;
    };

    double bnorm = 0;
    for (int i = 0; i < n; i++) bnorm += std::norm(b[i]);
    bnorm = std::sqrt(bnorm);
    double stop_tol = tol * std::max(bnorm, 1e-30);

    std::vector<Complex> r = b;
    std::vector<Complex> r_hat = r;
    std::vector<Complex> p(n, Complex(0, 0));
    std::vector<Complex> v(n, Complex(0, 0));
    std::vector<Complex> s(n, Complex(0, 0));
    std::vector<Complex> t(n, Complex(0, 0));

    Complex rho_old(1, 0);
    Complex alpha(1, 0);
    Complex omega(1, 0);

    for (int iter = 0; iter < max_iter; iter++) {
        Complex rho = dotc(r_hat, r);
        if (std::abs(rho) < 1e-30) {
            std::cout << "[FMM-BiCGSTAB] breakdown (rho≈0) at iter " << iter << "\n";
            return false;
        }

        Complex beta(1, 0);
        if (iter > 0) {
            beta = (rho / rho_old) * (alpha / omega);
        }
        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        matvec(p, v);  // FMM-accelerated matvec
        Complex rhat_v = dotc(r_hat, v);
        if (std::abs(rhat_v) < 1e-30) {
            std::cout << "[FMM-BiCGSTAB] breakdown (r̂·v≈0) at iter " << iter << "\n";
            return false;
        }
        alpha = rho / rhat_v;

        for (int i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];
        }

        matvec(s, t);  // FMM-accelerated matvec
        Complex ts = dotc(t, s);
        Complex tt = dotc(t, t);
        if (std::abs(tt) < 1e-30) {
            std::cout << "[FMM-BiCGSTAB] breakdown (t·t≈0) at iter " << iter << "\n";
            return false;
        }
        omega = ts / tt;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        double r_norm = 0;
        for (int i = 0; i < n; i++) r_norm += std::norm(r[i]);
        r_norm = std::sqrt(r_norm);
        if (r_norm < stop_tol) {
            std::cout << "[FMM-BiCGSTAB] converged after " << iter + 1
                      << " iterations, |r|/|b| = " << r_norm / std::max(bnorm, 1e-30) << "\n";
            return true;
        }

        rho_old = rho;
    }

    std::cout << "[FMM-BiCGSTAB] did not converge after " << max_iter
              << " iterations\n";
    return false;
}
