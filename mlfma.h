#ifndef MLFMA_H
#define MLFMA_H

#include "vec3_ops.h"
#include "mesh.h"
#include "rwg.h"
#include "constants.h"
#include "fmm_translator.h"
#include <vector>
#include <functional>

class EFIE;

class MLFMA {
public:
    MLFMA(const Mesh& mesh, const RWGData& rwg, double k0, double finest_box_size);

    int num_basis() const { return n_basis_; }
    int num_levels() const { return n_levels_; }

    void build();                       // 构建八叉树层级
    void precompute();                  // 预计算方向图和转移矩阵
    void fill_nearfield(const EFIE& efie); // 填充最细层近场 CSR

    // FMM 加速矩阵-向量积  y = Z * x
    void matvec(const std::vector<Complex>& x, std::vector<Complex>& y) const;

    std::function<void(const std::vector<Complex>&, std::vector<Complex>&)>
    get_matvec_functor() const;

    void print_stats() const;

private:
    const Mesh& mesh_;
    const RWGData& rwg_;
    double k0_, eta_;
    double finest_box_size_;
    int n_basis_;
    int n_levels_;                      // 0 = 最粗层, n_levels_-1 = 最细层

    static constexpr int N_CH = 4;      // x, y, z, div

    // ---- 单层数据 ----
    struct LevelData {
        double box_size;
        int L;                          // 截断数
        int K;                          // 总方向数 = N_theta * N_phi
        int N_theta;
        int N_phi;

        // 球面采样方向
        std::vector<PWSDirection> dirs;
        std::vector<double> ct_nodes;   // cosθ 节点 (Gauss-Legendre, ∈ [-1,1])

        // 3D 网格
        double xmin, ymin, zmin;
        int nx, ny, nz;
        int n_boxes_total;
        std::vector<int> nonempty;      // 非空盒子编号（全局索引）

        // 盒子中心
        std::vector<Vertex> centers;

        // 最细层：盒子 → 基函数映射
        std::vector<std::vector<int>> box_basis;

        // 父子关系
        std::vector<int> parent_box;    // 本层盒子 → 父层盒子（全局索引）
        std::vector<std::vector<int>> child_boxes; // 父层盒子 → 本层子盒子

        // 近邻列表（±1 邻域，同层）
        std::vector<std::vector<int>> near_boxes;

        // 互作用列表（M2L：本层远场但对父层是近场）
        std::vector<std::vector<int>> interaction_list;

        // M2L 转移
        struct M2LPair { int tgt, src, off; };
        std::vector<M2LPair> m2l_pairs;
        std::vector<Complex> translator;
    };

    std::vector<LevelData> levels_;

    // ---- 最细层方向图（方向优先: [K_finest][n_basis_]）----
    std::vector<std::vector<Complex>> Sx_, Sy_, Sz_, Sdiv_;
    std::vector<std::vector<Complex>> Tx_, Ty_, Tz_, Tdiv_;

    // ---- 近场 CSR（仅最细层）----
    std::vector<int> nf_ptr_, nf_col_;
    std::vector<Complex> nf_val_;

    // ---- 每迭代工作缓冲区 ----
    // agg_[level][box_idx_in_nonempty][ch * K + p]
    mutable std::vector<std::vector<std::vector<Complex>>> agg_;
    // down_[level][box_idx_in_nonempty][ch * K + p]
    mutable std::vector<std::vector<std::vector<Complex>>> down_;

    double C_;   // k0² * eta / (16π²)

    // ---- 辅助方法 ----
    int box_idx_at_level(const LevelData& lev, double x, double y, double z) const;

    // 计算基函数对某层某盒子的辐射/接收方向图
    void compute_rwg_patterns(int edge_idx, const LevelData& lev,
                               const Vertex& box_center,
                               std::vector<Complex>& sx, std::vector<Complex>& sy,
                               std::vector<Complex>& sz, std::vector<Complex>& sdiv,
                               std::vector<Complex>& tx, std::vector<Complex>& ty,
                               std::vector<Complex>& tz, std::vector<Complex>& tdiv) const;

    // 双立方插值: 细层 → 粗层（增加角采样 K_fine → K_coarse）
    // src/dst 布局: [N_CH * K]
    void interp_bicubic(const LevelData& lev_fine, const LevelData& lev_coarse,
                        const std::vector<Complex>& src,
                        std::vector<Complex>& dst) const;

    // 双立方反插值: 粗层 → 细层（减少角采样 K_coarse → K_fine）
    void anterp_bicubic(const LevelData& lev_coarse, const LevelData& lev_fine,
                        const std::vector<Complex>& src,
                        std::vector<Complex>& dst) const;

    // 相位平移：将方向图从一个中心移到另一个中心
    // dst[p] = src[p] * exp(j k0 k̂_p · (r_new - r_old))
    void phase_shift(const LevelData& lev, const Vertex& r_old, const Vertex& r_new,
                     const std::vector<Complex>& src,
                     std::vector<Complex>& dst) const;
};

#endif
