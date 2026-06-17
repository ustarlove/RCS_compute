#include "cg_solver.h"
#include <algorithm>
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

// ---------------------------------------------------------------------------
// GMRES(m) — 带重启的广义最小残差法，适用于复非对称矩阵
// 残差范数单调下降（理论最优性），每轮迭代 1 次 matvec
// ---------------------------------------------------------------------------
bool CGSolver::solve_gmres_functor(
    const std::function<void(const std::vector<Complex>&,
                              std::vector<Complex>&)>& matvec,
    const std::vector<Complex>& b,
    std::vector<Complex>& x,
    double tol,
    int max_iter,
    int restart) {
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
    if (bnorm < 1e-30) { std::cout << "[GMRES] zero RHS\n"; return true; }
    double stop_tol = tol * bnorm;

    int m = std::min(restart, max_iter);
    // Hessenberg 矩阵 H (m+1)×m
    std::vector<std::vector<Complex>> H(m + 1, std::vector<Complex>(m, Complex(0, 0)));
    // Givens 旋转参数: (c_real, s_complex) × m
    std::vector<double> c_giv(m, 0.0);
    std::vector<Complex> s_giv(m, Complex(0, 0));
    // 右端项 g (长度 m+1)
    std::vector<Complex> g(m + 1, Complex(0, 0));
    // Krylov 基 V (最多 m+1 个向量)
    std::vector<std::vector<Complex>> V;

    std::vector<Complex> r(n), w(n);
    std::vector<Complex> Av(n);  // 用于最终残差验证

    int total_iter = 0;

    for (int outer = 0; outer < max_iter; outer += m) {
        // 计算当前残差 r = b - A*x
        if (outer == 0) {
            r = b;
        } else {
            matvec(x, Av);
            for (int i = 0; i < n; i++) r[i] = b[i] - Av[i];
        }

        double beta = 0;
        for (int i = 0; i < n; i++) beta += std::norm(r[i]);
        beta = std::sqrt(beta);

        if (beta < stop_tol) {
            std::cout << "[GMRES] converged after " << total_iter
                      << " iterations, |r|/|b| = " << beta / bnorm << "\n";
            return true;
        }

        int inner = std::min(m, max_iter - outer);
        V.assign(inner + 1, std::vector<Complex>());

        // v1 = r / β
        V[0].resize(n);
        for (int i = 0; i < n; i++) V[0][i] = r[i] / beta;

        // 重置 H 和 g
        for (int i = 0; i <= inner; i++) {
            for (int j = 0; j < inner; j++) H[i][j] = Complex(0, 0);
            g[i] = Complex(0, 0);
        }
        g[0] = Complex(beta, 0);

        int j;
        for (j = 0; j < inner; j++) {
            // w = A * v_j
            matvec(V[j], w);

            // Arnoldi: Gram-Schmidt 正交化
            for (int i = 0; i <= j; i++) {
                H[i][j] = dotc(V[i], w);
                for (int k = 0; k < n; k++)
                    w[k] -= H[i][j] * V[i][k];
            }

            double h_next = 0;
            for (int i = 0; i < n; i++) h_next += std::norm(w[i]);
            h_next = std::sqrt(h_next);

            H[j + 1][j] = Complex(h_next, 0);

            if (h_next < 1e-30) {
                // Happy breakdown: 解已在当前子空间内
                inner = j + 1;
                break;
            }

            V[j + 1].resize(n);
            for (int i = 0; i < n; i++) V[j + 1][i] = w[i] / h_next;

            // 对 H 的第 j 列施加之前所有的 Givens 旋转
            for (int i = 0; i < j; i++) {
                Complex a = H[i][j];
                Complex b_val = H[i + 1][j];
                H[i][j]     = c_giv[i] * a + s_giv[i] * b_val;
                H[i + 1][j] = -std::conj(s_giv[i]) * a + c_giv[i] * b_val;
            }

            // 计算新的 Givens 旋转，消去 H[j+1][j]
            {
                Complex a = H[j][j];
                Complex b_val = H[j + 1][j];

                double a_abs = std::abs(a);
                double b_abs = std::abs(b_val);
                double r_raw = std::sqrt(a_abs * a_abs + b_abs * b_abs);

                if (r_raw < 1e-30) {
                    c_giv[j] = 1.0;
                    s_giv[j] = Complex(0, 0);
                } else {
                    c_giv[j] = a_abs / r_raw;
                    if (a_abs > 1e-30) {
                        Complex phase_a = a / a_abs;  // a / |a|
                        s_giv[j] = phase_a * std::conj(b_val) / r_raw;
                    } else {
                        s_giv[j] = std::conj(b_val) / r_raw;
                    }
                }

                // 施加新旋转到 H
                H[j][j]     = c_giv[j] * a + s_giv[j] * b_val;
                H[j + 1][j] = Complex(0, 0);

                // 施加新旋转到 g
                Complex g_old = g[j];
                g[j]     = c_giv[j] * g_old + s_giv[j] * g[j + 1];
                g[j + 1] = -std::conj(s_giv[j]) * g_old + c_giv[j] * g[j + 1];
            }

            double resid = std::abs(g[j + 1]);
            total_iter++;

            if (resid < stop_tol) {
                inner = j + 1;
                break;
            }
        }

        // 回代求解 y: H[0:inner-1][0:inner-1] * y = g[0:inner-1]
        // H 已经通过 Givens 旋转化为上三角阵
        std::vector<Complex> y(inner, Complex(0, 0));
        for (int i = inner - 1; i >= 0; i--) {
            Complex sum = g[i];
            for (int k = i + 1; k < inner; k++)
                sum -= H[i][k] * y[k];
            if (std::abs(H[i][i]) > 1e-30)
                y[i] = sum / H[i][i];
            else
                y[i] = Complex(0, 0);
        }

        // x = x + V[:, 0:inner-1] * y
        for (int k = 0; k < inner; k++) {
            for (int i = 0; i < n; i++) {
                x[i] += V[k][i] * y[k];
            }
        }

        // 检查是否收敛
        if (std::abs(g[inner]) < stop_tol) {
            std::cout << "[GMRES] converged after " << total_iter
                      << " iterations, |r|/|b| = " << std::abs(g[inner]) / bnorm << "\n";
            return true;
        }

        // 准备下一轮重启 (r 会在下次循环开头重算)
        std::cout << "[GMRES] restart at iter " << total_iter
                  << ", residual = " << std::abs(g[inner]) / bnorm << "\n";
    }

    // 最终残差
    matvec(x, Av);
    double final_res = 0;
    for (int i = 0; i < n; i++) {
        Complex diff = b[i] - Av[i];
        final_res += std::norm(diff);
    }
    std::cout << "[GMRES] did not converge after " << total_iter
              << " iterations, final |r|/|b| = "
              << std::sqrt(final_res) / bnorm << "\n";
    return false;
}
