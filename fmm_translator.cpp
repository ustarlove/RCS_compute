#include "fmm_translator.h"
#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// Gauss-Legendre 节点和权重 (Newton 迭代求 P_N 的根)
// ---------------------------------------------------------------------------
void Translator::gauss_legendre(int N, std::vector<double>& nodes,
                                 std::vector<double>& weights) {
    nodes.resize(N);
    weights.resize(N);

    for (int i = 0; i < N; i++) {
        // 初始猜测: Chebyshev 近似
        double x = std::cos(PI * (i + 0.75) / (N + 0.5));

        // Newton 迭代
        for (int iter = 0; iter < 20; iter++) {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 0; j < N; j++) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j + 1.0) * x * p2 - j * p3) / (j + 1.0);
            }
            // p1 = P_N(x), p2 = P_{N-1}(x)
            // P'_N(x) = N * (x * P_N(x) - P_{N-1}(x)) / (x² - 1) (pp)
            double pp = N * (x * p1 - p2) / (x * x - 1.0);
            //dx是本次迭代的修正量
            double dx = p1 / pp;
            x -= dx;
            if (std::abs(dx) < 1e-15) break;
        }

        // 再次计算 P_N(x) 和 P'_N(x) 以求权重
        double p1 = 1.0, p2 = 0.0;
        for (int j = 0; j < N; j++) {
            double p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j + 1.0) * x * p2 - j * p3) / (j + 1.0);
        }
        double pp = N * (x * p1 - p2) / (x * x - 1.0);

        nodes[i] = x;
        weights[i] = 2.0 / ((1.0 - x * x) * pp * pp);
    }
}

// ---------------------------------------------------------------------------
// 球 Hankel 函数 h_l^{(2)}(x) 前向递推
//   h_{-1}^{(2)}(x) = cos(x)/x - j sin(x)/x = e^{-jx}/x
//   h_0^{(2)}(x) = (sin x + j cos x)/x
//   h_{l+1} = (2l+1)/x * h_l - h_{l-1},  l = 0, 1, ..., L-1
// ---------------------------------------------------------------------------
void Translator::hankel_h2(int L, double x, std::vector<Complex>& hl) {
    hl.resize(L + 1);
    if (x < 1e-15) x = 1e-15;

    Complex h_m1(std::cos(x) / x, -std::sin(x) / x);  // h_{-1}
    hl[0] = Complex(std::sin(x), std::cos(x)) / x;    // h_0

    for (int l = 0; l < L; l++) {
        hl[l + 1] = (2.0 * l + 1.0) / x * hl[l] - (l == 0 ? h_m1 : hl[l - 1]);
    }
}

// ---------------------------------------------------------------------------
// Legendre 多项式 P_l(x) 前向递推
//   P_0 = 1, P_1 = x
//   (l+1) P_{l+1} = (2l+1) x P_l - l P_{l-1}
// ---------------------------------------------------------------------------
void Translator::legendre_p(int L, double x, std::vector<double>& pl) {
    pl.resize(L + 1);
    pl[0] = 1.0;
    if (L >= 1) {
        pl[1] = x;
        for (int l = 1; l < L; l++) {
            pl[l + 1] = ((2.0 * l + 1.0) * x * pl[l] - l * pl[l - 1]) / (l + 1.0);
        }
    }
}

// ---------------------------------------------------------------------------
// T_L(kR, cosθ) = Σ_{l=0}^{L} (-j)^l (2l+1) h_l^{(2)}(kR) P_l(cosθ)
// ---------------------------------------------------------------------------
Complex Translator::compute(double kR, double cos_theta, int L) {
    std::vector<Complex> hl;
    hankel_h2(L, kR, hl);
    std::vector<double> pl;
    legendre_p(L, cos_theta, pl);

    Complex T(0, 0);
    // (-j)^l 循环: 1, -j, -1, j, 1, ...
    Complex mj_l(1, 0);  // (-j)^0 = 1
    for (int l = 0; l <= L; l++) {
        T += mj_l * (2.0 * l + 1.0) * hl[l] * pl[l];
        mj_l *= Complex(0, -1);  // 乘 -j
    }
    return T;
}

// ---------------------------------------------------------------------------
// 生成球面角谱积分采样点
//   N_θ = L+1 Gauss-Legendre 在 cosθ ∈ [-1, 1] → θ ∈ [0, π]
//   N_φ = 2(L+1) 均匀采样在 φ ∈ [0, 2π)
//   K = N_θ × N_φ
//   总权重和 ≈ 4π
// ---------------------------------------------------------------------------
std::vector<PWSDirection> Translator::generate_quadrature(int L) {
    int N_theta = L + 1;
    int N_phi = 2 * (L + 1);

    std::vector<double> ct_nodes, ct_weights;
    gauss_legendre(N_theta, ct_nodes, ct_weights);

    double dphi = 2.0 * PI / N_phi;
    int K = N_theta * N_phi;
    std::vector<PWSDirection> dirs(K);

    int idx = 0;
    for (int it = 0; it < N_theta; it++) {
        double ct = ct_nodes[it];          // cosθ
        double st = std::sqrt(1.0 - ct * ct);
        // Gauss-Legendre 权重 × dφ (dφ = 2π/N_φ)
        double w = ct_weights[it] * dphi;

        for (int ip = 0; ip < N_phi; ip++) {
            double phi = ip * dphi;
            dirs[idx].k_hat = {st * std::cos(phi), st * std::sin(phi), ct};
            dirs[idx].weight = w;
            idx++;
        }
    }
    return dirs;
}
