#ifndef FMM_H
#define FMM_H

#include "vec3_ops.h"
#include "mesh.h"
#include "rwg.h"
#include "constants.h"
#include "fmm_translator.h"
#include <vector>
#include <functional>

class EFIE;

struct FMMBox {
    int ix, iy, iz;                      //小盒子索引（各方向上第几个）
    int global_idx;                      //总共来说的小盒子编号
    Vertex center;                       // 盒子中心 (世界坐标)
    std::vector<int> basis_idx;          // 盒内基函数编号
    std::vector<int> near_boxes;         // 近场盒子编号列表
    std::vector<int> far_boxes;          // 远场盒子编号列表
};

class FMM {
public:
    // mesh/rwg: 网格和基函数数据
    // k0: 波数
    // box_size: 盒子物理边长 (如 0.3*lambda)
    FMM(const Mesh& mesh, const RWGData& rwg, double k0, double box_size);

    int num_basis() const { return n_basis_; }
    int num_boxes() const { return n_boxes_; }
    int num_directions() const { return K_; }

    // 构建 3D 网格，分配基函数到盒子
    void build_boxes();

    // 预计算所有基函数的辐射/接收方向图
    void precompute_patterns();

    // 填充近场 CSR 稀疏矩阵 (调用 EFIE::compute_interaction)
    void fill_nearfield(const EFIE& efie);

    // FMM 加速的矩阵-向量积  y = Z * x
    void matvec(const std::vector<Complex>& x, std::vector<Complex>& y) const;

    // 返回 matvec 的 std::function 封装，供求解器使用
    std::function<void(const std::vector<Complex>&, std::vector<Complex>&)>
    get_matvec_functor() const;

    // 打印统计信息
    void print_stats() const;

private:
    const Mesh& mesh_;
    const RWGData& rwg_;
    double k0_;
    double eta_;
    double box_size_;
    int L_;              // 截断数
    int K_;              // 角方向数
    int n_basis_;        //基函数数量
    int n_boxes_;

    // 包围盒和网格尺寸
    double xmin_, xmax_, ymin_, ymax_, zmin_, zmax_;
    int nx_, ny_, nz_; //各方向上的盒子数（1开始）

    // 盒子列表 (扁平的 3D 网格)
    std::vector<FMMBox> boxes_;
    std::vector<int> nempty_;  // 非空盒子编号列表

    // 球面积分方向
    std::vector<PWSDirection> dirs_;

    // 辐射方向图 (方向优先存储: [K_][n_basis_])
    // 源方向图 (聚合用)
    std::vector<std::vector<Complex>> Sx_, Sy_, Sz_, Sdiv_;
    // 测试方向图 (解聚用)
    std::vector<std::vector<Complex>> Tx_, Ty_, Tz_, Tdiv_;

    // 预计算的转移矩阵 [n_pairs * K_]
    // pair_list_: {(target_box, source_box, offset_in_translator)}
    struct FarPair {
        int tgt;    // 目标盒子
        int src;    // 源盒子
        int off;    // translator_ 数组中的偏移（这对盒子对）
    };
    std::vector<FarPair> far_pairs_;
    std::vector<Complex> translator_;

    // 近场稀疏矩阵 (CSR)
    std::vector<int> nf_ptr_;       // 行指针, 长度 n_basis_+1
    std::vector<int> nf_col_;       // 列索引, 长度 nnz
    std::vector<Complex> nf_val_;   // 值, 长度 nnz

    // 每迭代可变缓冲区
    // agg_[box][dir][ch], trans_[box][dir][ch]
    // 其中 ch: 0=x, 1=y, 2=z, 3=div
    static constexpr int N_CH = 4;
    mutable std::vector<std::vector<std::vector<Complex>>> agg_;
    mutable std::vector<std::vector<std::vector<Complex>>> trans_;

    double C_;         // k0² * eta / (16π²)

    // 辅助方法
    int box_idx(double x, double y, double z) const;
    void compute_rwg_patterns(int edge_idx,
                              std::vector<Complex>& sx, std::vector<Complex>& sy,
                              std::vector<Complex>& sz, std::vector<Complex>& sdiv,
                              std::vector<Complex>& tx, std::vector<Complex>& ty,
                              std::vector<Complex>& tz, std::vector<Complex>& tdiv) const;
};

#endif
