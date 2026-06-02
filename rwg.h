#ifndef RWG_H
#define RWG_H

#include "mesh.h"
#include "constants.h"

struct Edge {
    int i;          // 边编号
    int t1, t2;     // 两个相邻三角形的索引（t2 = -1 表示边界边）
    int v1, v2;     // 边两端顶点的索引
    double length;  // 边长
    Vertex apex1;   // t1 中与边相对的顶点坐标（RWG ρ⁺的自由度顶点）
    Vertex apex2;   // t2 中与边相对的顶点坐标（RWG ρ⁻的自由度顶点）
};

struct RWGData {
    std::vector<Edge> edges;        // 所有公共边
    std::vector<int> edge_to_t1;    // 边->三角形1 映射
    std::vector<int> edge_to_t2;    // 边->三角形2 映射
    std::vector<double> rwg_coeff1; // 三角形1上的 RWG 系数 l/(2A⁺)
    std::vector<double> rwg_coeff2; // 三角形2上的 RWG 系数 l/(2A⁻)

    void build(const Mesh& mesh);

    // 获取 RWG 基函数 f_n 在给定三角形 tri_idx 上点 r 处的向量值
    // f⁺(r) = l/(2A⁺) * (r - apex⁺)   on T⁺
    // f⁻(r) = -l/(2A⁻) * (r - apex⁻)   on T⁻
    void get_rwg(int edge_idx, int tri_idx,
                 const Vertex& r, Complex* value) const;

    // 获取 RWG 基函数在给定三角形的散度（常量）
    // ∇·f⁺ = l/A⁺,  ∇·f⁻ = -l/A⁻
    double get_divergence(int edge_idx, int tri_idx) const;
};

#endif // RWG_H