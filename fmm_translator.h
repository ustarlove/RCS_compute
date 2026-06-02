#ifndef FMM_TRANSLATOR_H
#define FMM_TRANSLATOR_H

#include "vec3_ops.h"
#include "constants.h"
#include <vector>

struct PWSDirection {
    Vertex k_hat;   // 单位方向向量
    double weight;  // 积分权重（球面积分 d²k̂ 的离散化）
};

class Translator {
public:
    // 生成球面上的角谱积分采样点，L 为截断数
    // K ≈ 2*(L+1)² 个方向
    static std::vector<PWSDirection> generate_quadrature(int L);

    // 计算对角转移算子 T_L(kR, cosθ)
    static Complex compute(double kR, double cos_theta, int L);

private:
    // h_l^{(2)}(x), l = 0..L, 前向递推
    static void hankel_h2(int L, double x, std::vector<Complex>& hl);

    // P_l(x), l = 0..L, 前向递推
    static void legendre_p(int L, double x, std::vector<double>& pl);

    // 计算 N 点 Gauss-Legendre 节点和权重
    static void gauss_legendre(int N, std::vector<double>& nodes,
                               std::vector<double>& weights);
};

#endif
