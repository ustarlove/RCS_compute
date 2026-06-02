#ifndef EFIE_H
#define EFIE_H

#include "mesh.h"
#include "rwg.h"
#include "constants.h"
#include <vector>

class EFIE {
public:
    EFIE(const Mesh& mesh, const RWGData& rwg, double k0);
    
    // 填充阻抗矩阵 Z
    void fill_matrix(std::vector<std::vector<Complex>>& Z);

    // 填充激励向量 b（平面波入射）
    void fill_vector(const Vertex& k_inc, const Vertex& E0_pol,
                     std::vector<Complex>& b);

    // 计算两条 RWG 边之间的相互作用 (公开给 FMM 近场填充使用)
    Complex compute_interaction(int edge_i, int edge_j) const;

    double get_k0() const { return k0_; }
    int num_edges() const { return n_edges_; }

private:
    const Mesh& mesh_;
    const RWGData& rwg_;
    double k0_;          // 波数
    double lambda_;      // 波长
    int n_edges_;

    // 计算点 r 处的平面波电场（入射）
    Complex plane_wave_e(const Vertex& r, const Vertex& k_inc,
                         const Vertex& E0_pol) const;
    
};

#endif // EFIE_H