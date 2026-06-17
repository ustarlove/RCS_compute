#include "mlfma.h"
#include "efie.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ===== 双立方 Lagrange 插值工具（匿名命名空间）=============================

namespace {

// 4 点 Lagrange 插值节点和权重（非周期、非均匀节点）
void lagrange_4_nonperiodic(const std::vector<double>& src_x,
                            const std::vector<double>& dst_x,
                            std::vector<std::vector<int>>& stencil,
                            std::vector<std::vector<double>>& weight) {
    int N_src = (int)src_x.size();
    int N_dst = (int)dst_x.size();
    stencil.assign(N_dst, std::vector<int>(4, 0));
    weight.assign(N_dst, std::vector<double>(4, 0.0));

    for (int j = 0; j < N_dst; j++) {
        double x = dst_x[j];

        // 二分查找最近的源节点区间
        int lo = 0, hi = N_src - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (src_x[mid] <= x) lo = mid;
            else                  hi = mid;
        }
        // 若 x 在 src_x 范围外，lo 和 hi 会卡在边界
        int i0 = lo - 1;
        if (i0 < 0) i0 = 0;
        if (i0 + 4 > N_src) i0 = N_src - 4;
        if (i0 < 0) i0 = 0;

        for (int s = 0; s < 4; s++)
            stencil[j][s] = i0 + s;

        for (int s = 0; s < 4; s++) {
            double xs = src_x[stencil[j][s]];
            double num = 1.0, den = 1.0;
            for (int t = 0; t < 4; t++) {
                if (t == s) continue;
                double xt = src_x[stencil[j][t]];
                num *= (x - xt);
                den *= (xs - xt);
            }
            weight[j][s] = num / den;
        }
    }
}

// 4 点 Lagrange 插值节点和权重（周期均匀节点，φ 方向）
void lagrange_4_periodic(int N_src, int N_dst,
                         std::vector<std::vector<int>>& stencil,
                         std::vector<std::vector<double>>& weight) {
    stencil.assign(N_dst, std::vector<int>(4, 0));
    weight.assign(N_dst, std::vector<double>(4, 0.0));

    double dphi_src = 2.0 * PI / N_src;
    double dphi_dst = 2.0 * PI / N_dst;

    for (int j = 0; j < N_dst; j++) {
        double phi = j * dphi_dst;
        double t = phi / dphi_src;      // 连续源索引
        int ic = (int)std::floor(t);    // 中心源索引

        for (int s = 0; s < 4; s++) {
            int idx = ic - 1 + s;
            stencil[j][s] = (idx % N_src + N_src) % N_src;
        }

        // 用非包裹坐标计算 Lagrange 权重
        double xs[4];
        for (int s = 0; s < 4; s++)
            xs[s] = (ic - 1.0 + s) * dphi_src;

        for (int s = 0; s < 4; s++) {
            double num = 1.0, den = 1.0;
            for (int t = 0; t < 4; t++) {
                if (t == s) continue;
                num *= (phi - xs[t]);
                den *= (xs[s] - xs[t]);
            }
            weight[j][s] = num / den;
        }
    }
}

} // anonymous namespace

// ===== MLFMA 构造 =======================================================

MLFMA::MLFMA(const Mesh& mesh, const RWGData& rwg, double k0, double finest_box_size)
    : mesh_(mesh), rwg_(rwg), k0_(k0), finest_box_size_(finest_box_size) {
    eta_ = std::sqrt(MU0 / EPS0);
    n_basis_ = (int)rwg_.edges.size();
}

// ===== 盒子索引 ==========================================================

int MLFMA::box_idx_at_level(const LevelData& lev, double x, double y, double z) const {
    int ix = (int)((x - lev.xmin) / lev.box_size);
    int iy = (int)((y - lev.ymin) / lev.box_size);
    int iz = (int)((z - lev.zmin) / lev.box_size);
    ix = std::max(0, std::min(ix, lev.nx - 1));
    iy = std::max(0, std::min(iy, lev.ny - 1));
    iz = std::max(0, std::min(iz, lev.nz - 1));
    return iz * (lev.nx * lev.ny) + iy * lev.nx + ix;
}

// ===== 八叉树构建 ========================================================

void MLFMA::build() {
    // 1. 计算包围盒
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = ymin = zmin = 1e30;
    xmax = ymax = zmax = -1e30;
    for (const auto& v : mesh_.vertices) {
        xmin = std::min(xmin, v.x); xmax = std::max(xmax, v.x);
        ymin = std::min(ymin, v.y); ymax = std::max(ymax, v.y);
        zmin = std::min(zmin, v.z); zmax = std::max(zmax, v.z);
    }
    double eps = 1e-6;
    xmin -= eps; xmax += eps;
    ymin -= eps; ymax += eps;
    zmin -= eps; zmax += eps;

    double Lx = xmax - xmin;
    double Ly = ymax - ymin;
    double Lz = zmax - zmin;
    double Lmax = std::max({Lx, Ly, Lz});

    // 2. 确定层数：finest_box_size * 2^(n-1) >= Lmax
    n_levels_ = 1;
    while (finest_box_size_ * std::pow(2.0, n_levels_ - 1) < Lmax)
        n_levels_++;

    levels_.resize(n_levels_);

    // 3. 逐层初始化网格
    for (int l = 0; l < n_levels_; l++) {
        LevelData& lev = levels_[l];
        // level 0 = coarsest (largest box), level n-1 = finest
        lev.box_size = finest_box_size_ * std::pow(2.0, n_levels_ - 1 - l);
        lev.xmin = xmin; lev.ymin = ymin; lev.zmin = zmin;

        lev.nx = std::max(1, (int)std::ceil(Lx / lev.box_size));
        lev.ny = std::max(1, (int)std::ceil(Ly / lev.box_size));
        lev.nz = std::max(1, (int)std::ceil(Lz / lev.box_size));
        lev.n_boxes_total = lev.nx * lev.ny * lev.nz;

        // 初始化盒子中心
        lev.centers.resize(lev.n_boxes_total);
        for (int iz = 0; iz < lev.nz; iz++) {
            for (int iy = 0; iy < lev.ny; iy++) {
                for (int ix = 0; ix < lev.nx; ix++) {
                    int idx = iz * (lev.nx * lev.ny) + iy * lev.nx + ix;
                    lev.centers[idx] = {
                        lev.xmin + (ix + 0.5) * lev.box_size,
                        lev.ymin + (iy + 0.5) * lev.box_size,
                        lev.zmin + (iz + 0.5) * lev.box_size
                    };
                }
            }
        }
    }

    // 4. 最细层：分配基函数到盒子
    LevelData& finest = levels_[n_levels_ - 1];
    finest.box_basis.resize(finest.n_boxes_total);
    for (int j = 0; j < n_basis_; j++) {
        const Edge& e = rwg_.edges[j];
        double cx = 0.5 * (mesh_.vertices[e.v1].x + mesh_.vertices[e.v2].x);
        double cy = 0.5 * (mesh_.vertices[e.v1].y + mesh_.vertices[e.v2].y);
        double cz = 0.5 * (mesh_.vertices[e.v1].z + mesh_.vertices[e.v2].z);
        int b = box_idx_at_level(finest, cx, cy, cz);
        finest.box_basis[b].push_back(j);
    }

    // 收集最细层非空盒子
    finest.nonempty.clear();
    for (int b = 0; b < finest.n_boxes_total; b++)
        if (!finest.box_basis[b].empty())
            finest.nonempty.push_back(b);

    // 5. 父子关系：自底向上：自底向上
    // parent_box: 本层每个盒子 → 父层盒子索引
    // child_boxes: 父层每个盒子 → 本层子盒子列表
    for (int l = n_levels_ - 1; l >= 1; l--) {
        LevelData& child = levels_[l];
        LevelData& parent = levels_[l - 1];

        child.parent_box.assign(child.n_boxes_total, -1);
        parent.child_boxes.assign(parent.n_boxes_total, std::vector<int>());

        for (int b = 0; b < child.n_boxes_total; b++) {
            const Vertex& c = child.centers[b];
            int p = box_idx_at_level(parent, c.x, c.y, c.z);
            child.parent_box[b] = p;
            parent.child_boxes[p].push_back(b);
        }
    }

    // 6. 确定每层非空盒子（含后代有基函数的盒子）
    for (int l = n_levels_ - 2; l >= 0; l--) {
        LevelData& lev = levels_[l];
        lev.nonempty.clear();
        std::vector<bool> active(lev.n_boxes_total, false);

        // 子层的非空盒子决定了本层的非空盒子
        LevelData& child = levels_[l + 1];
        for (int ci = 0; ci < child.n_boxes_total; ci++) {
            if (l + 1 == n_levels_ - 1) {
                if (!child.box_basis[ci].empty()) {
                    int p = child.parent_box[ci];
                    if (p >= 0) active[p] = true;
                }
            } else {
                // 检查子盒子是否在子层的 nonempty 中
                bool child_active = false;
                for (int nb : child.nonempty) {
                    if (nb == ci) { child_active = true; break; }
                }
                if (child_active) {
                    int p = child.parent_box[ci];
                    if (p >= 0) active[p] = true;
                }
            }
        }
        for (int b = 0; b < lev.n_boxes_total; b++)
            if (active[b]) lev.nonempty.push_back(b);
    }

    // 7. 每层构建近邻列表和互作用列表
    for (int l = 0; l < n_levels_; l++) {
        LevelData& lev = levels_[l];

        // 近邻列表 (±1 邻域同层)
        lev.near_boxes.assign(lev.n_boxes_total, std::vector<int>());
        for (int nb : lev.nonempty) {
            int ix = nb % (lev.nx * lev.ny) % lev.nx;
            int iy = (nb % (lev.nx * lev.ny)) / lev.nx;
            int iz = nb / (lev.nx * lev.ny);

            for (int dz = -1; dz <= 1; dz++) {
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        int ni = ix + dx, nj = iy + dy, nk = iz + dz;
                        if (ni < 0 || ni >= lev.nx ||
                            nj < 0 || nj >= lev.ny ||
                            nk < 0 || nk >= lev.nz) continue;
                        int nbidx = nk * (lev.nx * lev.ny) + nj * lev.nx + ni;
                        if (l == n_levels_ - 1) {
                            if (!lev.box_basis[nbidx].empty())
                                lev.near_boxes[nb].push_back(nbidx);
                        } else {
                            // 粗层：只要盒子活跃就是 "近场"
                            bool active = false;
                            for (int a : lev.nonempty)
                                if (a == nbidx) { active = true; break; }
                            if (active)
                                lev.near_boxes[nb].push_back(nbidx);
                        }
                    }
                }
            }
        }

        // 互作用列表
        lev.interaction_list.assign(lev.n_boxes_total, std::vector<int>());
        if (l == 0) {
            // 最粗层：所有其他非空且非近邻的盒子
            for (int nb : lev.nonempty) {
                for (int fb : lev.nonempty) {
                    if (fb == nb) continue;
                    auto& nr = lev.near_boxes[nb];
                    if (std::find(nr.begin(), nr.end(), fb) == nr.end())
                        lev.interaction_list[nb].push_back(fb);
                }
            }
        } else {
            // 标准定义：父层近邻的所有子盒子，排除本盒子的同层近邻
            LevelData& parent = levels_[l - 1];
            for (int nb : lev.nonempty) {
                int p = lev.parent_box[nb];
                if (p < 0) continue;

                // 收集父层近邻的所有子盒子
                std::vector<int> cousins;
                for (int pn : parent.near_boxes[p]) {
                    for (int c : parent.child_boxes[pn])
                        cousins.push_back(c);
                }

                // 排除本盒子及其近邻
                std::vector<int>& nr = lev.near_boxes[nb];
                for (int c : cousins) {
                    if (c == nb) continue;
                    // 只保留活跃盒子
                    bool active = false;
                    for (int a : lev.nonempty)
                        if (a == c) { active = true; break; }
                    if (!active) continue;
                    if (std::find(nr.begin(), nr.end(), c) == nr.end())
                        lev.interaction_list[nb].push_back(c);
                }
            }
        }
    }

    std::cout << "MLFMA: " << n_levels_ << " levels\n";
    for (int l = 0; l < n_levels_; l++) {
        std::cout << "  Level " << l << ": box=" << levels_[l].box_size
                  << "m, grid=" << levels_[l].nx << "x" << levels_[l].ny
                  << "x" << levels_[l].nz
                  << ", nonempty=" << levels_[l].nonempty.size() << "\n";
    }
}

// ===== 计算单个基函数的方向图 ==============================================

void MLFMA::compute_rwg_patterns(int edge_idx, const LevelData& lev,
                                  const Vertex& box_center,
                                  std::vector<Complex>& sx, std::vector<Complex>& sy,
                                  std::vector<Complex>& sz, std::vector<Complex>& sdiv,
                                  std::vector<Complex>& tx, std::vector<Complex>& ty,
                                  std::vector<Complex>& tz, std::vector<Complex>& tdiv) const {
    const Edge& e = rwg_.edges[edge_idx];
    double l = e.length;
    double inv_k = 1.0 / k0_;

    Vertex cp = mesh_.triangles[e.t1].center(mesh_.vertices);
    Vertex cm = mesh_.triangles[e.t2].center(mesh_.vertices);

    double vpx = 0.5 * l * (cp.x - e.apex1.x);
    double vpy = 0.5 * l * (cp.y - e.apex1.y);
    double vpz = 0.5 * l * (cp.z - e.apex1.z);
    double vmx = -0.5 * l * (cm.x - e.apex2.x);
    double vmy = -0.5 * l * (cm.y - e.apex2.y);
    double vmz = -0.5 * l * (cm.z - e.apex2.z);
    double div_l = l * inv_k;

    int K = lev.K;
    sx.assign(K, Complex(0, 0));
    sy.assign(K, Complex(0, 0));
    sz.assign(K, Complex(0, 0));
    sdiv.assign(K, Complex(0, 0));
    tx.assign(K, Complex(0, 0));
    ty.assign(K, Complex(0, 0));
    tz.assign(K, Complex(0, 0));
    tdiv.assign(K, Complex(0, 0));

    Vertex dcp = cp - box_center;
    Vertex dcm = cm - box_center;

    for (int p = 0; p < K; p++) {
        const Vertex& khat = lev.dirs[p].k_hat;

        double ph_p = k0_ * dot(khat, dcp);
        double ph_m = k0_ * dot(khat, dcm);
        Complex ep(std::cos(ph_p), std::sin(ph_p));
        Complex em(std::cos(ph_m), std::sin(ph_m));

        sx[p]   = Complex(vpx, 0) * ep + Complex(vmx, 0) * em;
        sy[p]   = Complex(vpy, 0) * ep + Complex(vmy, 0) * em;
        sz[p]   = Complex(vpz, 0) * ep + Complex(vmz, 0) * em;
        sdiv[p] = Complex(div_l, 0) * (ep - em);

        Complex epc(std::cos(ph_p), -std::sin(ph_p));
        Complex emc(std::cos(ph_m), -std::sin(ph_m));
        tx[p]   = Complex(vpx, 0) * epc + Complex(vmx, 0) * emc;
        ty[p]   = Complex(vpy, 0) * epc + Complex(vmy, 0) * emc;
        tz[p]   = Complex(vpz, 0) * epc + Complex(vmz, 0) * emc;
        tdiv[p] = Complex(-div_l, 0) * (epc - emc);
    }
}

// ===== 预计算 ===========================================================

void MLFMA::precompute() {
    // 1. 每层计算 L 和 K，生成球面采样
    for (int l = 0; l < n_levels_; l++) {
        LevelData& lev = levels_[l];
        double a = 0.5 * lev.box_size * std::sqrt(3.0);
        double ka = k0_ * a;
        int L_calc = (int)(ka + 1.8 * std::cbrt(ka) * 3.0);
        lev.L = std::max(2, std::min(L_calc, 20));  // MLFMA 粗层 L 更大
        lev.dirs = Translator::generate_quadrature(lev.L);
        lev.K = (int)lev.dirs.size();
        lev.N_theta = lev.L + 1;
        lev.N_phi = 2 * (lev.L + 1);

        // 保存 cosθ 节点
        lev.ct_nodes.resize(lev.N_theta);
        for (int i = 0; i < lev.N_theta; i++)
            lev.ct_nodes[i] = lev.dirs[i * lev.N_phi].k_hat.z;  // cosθ = z component

        double total_w = 0;
        for (const auto& d : lev.dirs) total_w += d.weight;
        std::cout << "  Level " << l << ": L=" << lev.L << ", K=" << lev.K
                  << ", Σw=" << total_w << " (4π=" << 4*PI << ")\n";
    }

    // 2. 最细层：预计算所有基函数的方向图
    LevelData& finest = levels_[n_levels_ - 1];
    int Kf = finest.K;

    Sx_.assign(Kf, std::vector<Complex>(n_basis_));
    Sy_.assign(Kf, std::vector<Complex>(n_basis_));
    Sz_.assign(Kf, std::vector<Complex>(n_basis_));
    Sdiv_.assign(Kf, std::vector<Complex>(n_basis_));
    Tx_.assign(Kf, std::vector<Complex>(n_basis_));
    Ty_.assign(Kf, std::vector<Complex>(n_basis_));
    Tz_.assign(Kf, std::vector<Complex>(n_basis_));
    Tdiv_.assign(Kf, std::vector<Complex>(n_basis_));

    std::vector<Complex> sx(Kf), sy(Kf), sz(Kf), sdiv(Kf);
    std::vector<Complex> tx(Kf), ty(Kf), tz(Kf), tdiv(Kf);

    for (int j = 0; j < n_basis_; j++) {
        // 找到该基函数所在的最细层盒子
        const Edge& e = rwg_.edges[j];
        double cx = 0.5 * (mesh_.vertices[e.v1].x + mesh_.vertices[e.v2].x);
        double cy = 0.5 * (mesh_.vertices[e.v1].y + mesh_.vertices[e.v2].y);
        double cz = 0.5 * (mesh_.vertices[e.v1].z + mesh_.vertices[e.v2].z);
        int b = box_idx_at_level(finest, cx, cy, cz);
        Vertex c_box = finest.centers[b];

        compute_rwg_patterns(j, finest, c_box, sx, sy, sz, sdiv, tx, ty, tz, tdiv);

        for (int p = 0; p < Kf; p++) {
            Sx_[p][j] = sx[p];      Sy_[p][j] = sy[p];
            Sz_[p][j] = sz[p];      Sdiv_[p][j] = sdiv[p];
            Tx_[p][j] = tx[p];      Ty_[p][j] = ty[p];
            Tz_[p][j] = tz[p];      Tdiv_[p][j] = tdiv[p];
        }
    }

    // 3. 每层预计算 M2L 转移矩阵
    for (int l = 0; l < n_levels_; l++) {
        LevelData& lev = levels_[l];
        int K = lev.K;

        // 统计总配对数
        int total_pairs = 0;
        for (int nb : lev.nonempty)
            total_pairs += (int)lev.interaction_list[nb].size();

        lev.translator.resize(total_pairs * K);
        lev.m2l_pairs.clear();

        int pair_idx = 0;
        for (int nb : lev.nonempty) {
            Vertex c_tgt = lev.centers[nb];
            for (int src : lev.interaction_list[nb]) {
                Vertex c_src = lev.centers[src];
                Vertex R = c_tgt - c_src;
                double kR = k0_ * norm(R);
                Vertex R_hat = normalize(R);

                lev.m2l_pairs.push_back({nb, src, pair_idx * K});

                for (int p = 0; p < K; p++) {
                    double ct = dot(lev.dirs[p].k_hat, R_hat);
                    if (ct > 1.0) ct = 1.0;
                    if (ct < -1.0) ct = -1.0;
                    lev.translator[pair_idx * K + p] = Translator::compute(kR, ct, lev.L);
                }
                pair_idx++;
            }
        }
        std::cout << "  Level " << l << ": " << lev.m2l_pairs.size()
                  << " M2L pairs\n";
    }

    // 4. 常数
    C_ = k0_ * k0_ * eta_ / (16.0 * PI * PI);

    // 5. 分配工作缓冲区
    agg_.resize(n_levels_);
    down_.resize(n_levels_);
    for (int l = 0; l < n_levels_; l++) {
        int n_nempty = (int)levels_[l].nonempty.size();
        int K = levels_[l].K;
        agg_[l].assign(n_nempty, std::vector<Complex>(N_CH * K, Complex(0, 0)));
        down_[l].assign(n_nempty, std::vector<Complex>(N_CH * K, Complex(0, 0)));
    }
}

// ===== 近场填充 =========================================================

void MLFMA::fill_nearfield(const EFIE& efie) {
    // 仅在最细层填近场
    LevelData& finest = levels_[n_levels_ - 1];
    int nfn = (int)finest.nonempty.size();

    // 统计每行非零元
    std::vector<int> row_count(n_basis_, 0);
    for (int ni = 0; ni < nfn; ni++) {
        int nb = finest.nonempty[ni];
        for (int i : finest.box_basis[nb]) {
            for (int nbr : finest.near_boxes[nb]) {
                row_count[i] += (int)finest.box_basis[nbr].size();
            }
        }
    }

    nf_ptr_.resize(n_basis_ + 1);
    nf_ptr_[0] = 0;
    for (int i = 0; i < n_basis_; i++)
        nf_ptr_[i + 1] = nf_ptr_[i] + row_count[i];

    int nnz = nf_ptr_[n_basis_];
    nf_col_.resize(nnz);
    nf_val_.resize(nnz);

    std::cout << "MLFMA near-field: " << nnz << " nonzeros ("
              << 100.0 * nnz / (n_basis_ * n_basis_) << "% of dense)\n";

    // 填充
    std::vector<int> cur_ptr(n_basis_, 0);
    for (int ni = 0; ni < nfn; ni++) {
        int nb = finest.nonempty[ni];
        for (int i : finest.box_basis[nb]) {
            int base = nf_ptr_[i];
            for (int nbr : finest.near_boxes[nb]) {
                for (int j : finest.box_basis[nbr]) {
                    int pos = base + cur_ptr[i]++;
                    nf_col_[pos] = j;
                    nf_val_[pos] = efie.compute_interaction(i, j);
                }
            }
        }
    }

    // 排序每行
    for (int i = 0; i < n_basis_; i++) {
        int beg = nf_ptr_[i], end = nf_ptr_[i + 1];
        for (int a = beg; a < end; a++) {
            for (int b = a + 1; b < end; b++) {
                if (nf_col_[a] > nf_col_[b]) {
                    std::swap(nf_col_[a], nf_col_[b]);
                    std::swap(nf_val_[a], nf_val_[b]);
                }
            }
        }
    }
    std::cout << "  near-field matrix complete.\n";
}

// ===== 双立方插值：细 → 粗（增加 K）======================================

void MLFMA::interp_bicubic(const LevelData& lev_fine, const LevelData& lev_coarse,
                            const std::vector<Complex>& src,
                            std::vector<Complex>& dst) const {
    int Kf = lev_fine.K, Kc = lev_coarse.K;
    int Ntf = lev_fine.N_theta, Ntc = lev_coarse.N_theta;
    int Npf = lev_fine.N_phi,   Npc = lev_coarse.N_phi;

    dst.assign(N_CH * Kc, Complex(0, 0));

    // 预计算插值节点和权重（可缓存，这里每次计算，成本相对 matvec 可忽略）
    std::vector<std::vector<int>> phi_st, theta_st;
    std::vector<std::vector<double>> phi_wt, theta_wt;
    lagrange_4_periodic(Npf, Npc, phi_st, phi_wt);
    lagrange_4_nonperiodic(lev_fine.ct_nodes, lev_coarse.ct_nodes, theta_st, theta_wt);

    for (int ch = 0; ch < N_CH; ch++) {
        // Step 1: φ 插值 — 对每个细 θ，插到粗 φ
        // temp[is][jc] for is=0..Ntf-1, jc=0..Npc-1
        std::vector<Complex> temp(Ntf * Npc, Complex(0, 0));
        for (int is = 0; is < Ntf; is++) {
            for (int jc = 0; jc < Npc; jc++) {
                Complex sum(0, 0);
                for (int s = 0; s < 4; s++) {
                    int js = phi_st[jc][s];
                    sum += phi_wt[jc][s] * src[ch * Kf + is * Npf + js];
                }
                temp[is * Npc + jc] = sum;
            }
        }

        // Step 2: θ 插值 — 对每个粗 φ，从细 θ 插到粗 θ
        for (int ic = 0; ic < Ntc; ic++) {
            for (int jc = 0; jc < Npc; jc++) {
                Complex sum(0, 0);
                for (int s = 0; s < 4; s++) {
                    int is = theta_st[ic][s];
                    sum += theta_wt[ic][s] * temp[is * Npc + jc];
                }
                dst[ch * Kc + ic * Npc + jc] = sum;
            }
        }
    }
}

// ===== 双立方反插值：粗 → 细（减少 K，插值的伴随算子）=====================

void MLFMA::anterp_bicubic(const LevelData& lev_coarse, const LevelData& lev_fine,
                            const std::vector<Complex>& src,
                            std::vector<Complex>& dst) const {
    int Kc = lev_coarse.K, Kf = lev_fine.K;
    int Ntc = lev_coarse.N_theta, Ntf = lev_fine.N_theta;
    int Npc = lev_coarse.N_phi,   Npf = lev_fine.N_phi;

    dst.assign(N_CH * Kf, Complex(0, 0));

    // 预计算插值节点和权重（精细→粗糙 方向）
    std::vector<std::vector<int>> phi_st, theta_st;
    std::vector<std::vector<double>> phi_wt, theta_wt;
    lagrange_4_periodic(Npf, Npc, phi_st, phi_wt);
    lagrange_4_nonperiodic(lev_fine.ct_nodes, lev_coarse.ct_nodes, theta_st, theta_wt);

    for (int ch = 0; ch < N_CH; ch++) {
        // Step 1: θ 反插值 — 对每个粗 φ，从粗 θ 分布到细 θ
        // temp[is_fine][jc_coarse]
        std::vector<Complex> temp(Ntf * Npc, Complex(0, 0));
        for (int ic = 0; ic < Ntc; ic++) {
            for (int jc = 0; jc < Npc; jc++) {
                Complex val = src[ch * Kc + ic * Npc + jc];
                for (int s = 0; s < 4; s++) {
                    int is = theta_st[ic][s];
                    double w = theta_wt[ic][s];
                    temp[is * Npc + jc] += w * val;
                }
            }
        }

        // Step 2: φ 反插值 — 对每个细 θ，从粗 φ 分布到细 φ
        for (int is = 0; is < Ntf; is++) {
            for (int jc = 0; jc < Npc; jc++) {
                Complex val = temp[is * Npc + jc];
                for (int s = 0; s < 4; s++) {
                    int jf = phi_st[jc][s];
                    double w = phi_wt[jc][s];
                    dst[ch * Kf + is * Npf + jf] += w * val;
                }
            }
        }
    }
}

// ===== 相位平移 ==========================================================

void MLFMA::phase_shift(const LevelData& lev, const Vertex& r_old, const Vertex& r_new,
                         const std::vector<Complex>& src,
                         std::vector<Complex>& dst) const {
    int K = lev.K;
    dst.resize(N_CH * K);

    Vertex dr = r_old - r_new;   // r_new 相对于 r_old 的位移
    for (int ch = 0; ch < N_CH; ch++) {
        for (int p = 0; p < K; p++) {
            double ph = k0_ * dot(lev.dirs[p].k_hat, dr);
            Complex factor(std::cos(ph), std::sin(ph));
            dst[ch * K + p] = src[ch * K + p] * factor;
        }
    }
}

// ===== MLFMA 矩阵-向量积 =================================================

void MLFMA::matvec(const std::vector<Complex>& x, std::vector<Complex>& y) const {
    y.assign(n_basis_, Complex(0, 0));

    // ---- 阶段 0: 清理工作缓冲区 ----
    for (int l = 0; l < n_levels_; l++) {
        const LevelData& lev = levels_[l];
        int K = lev.K;
        int nne = (int)lev.nonempty.size();
        for (int ni = 0; ni < nne; ni++) {
            for (int q = 0; q < N_CH * K; q++) {
                agg_[l][ni][q] = Complex(0, 0);
                down_[l][ni][q] = Complex(0, 0);
            }
        }
    }

    // ---- 阶段 1: 近场 (最细层 CSR) ----
    for (int i = 0; i < n_basis_; i++) {
        for (int ptr = nf_ptr_[i]; ptr < nf_ptr_[i + 1]; ptr++)
            y[i] += nf_val_[ptr] * x[nf_col_[ptr]];
    }

    // =====================================================
    // ---- 阶段 2: 上行 (聚合) ----
    // =====================================================

    const LevelData& finest = levels_[n_levels_ - 1];
    int Kf = finest.K;

    // 2a: 最细层：基函数 → 盒子聚合
    for (int ni = 0; ni < (int)finest.nonempty.size(); ni++) {
        int nb = finest.nonempty[ni];
        for (int j : finest.box_basis[nb]) {
            Complex xj = x[j];
            if (std::norm(xj) < 1e-30) continue;
            for (int p = 0; p < Kf; p++) {
                agg_[n_levels_ - 1][ni][0 * Kf + p] += xj * Sx_[p][j];
                agg_[n_levels_ - 1][ni][1 * Kf + p] += xj * Sy_[p][j];
                agg_[n_levels_ - 1][ni][2 * Kf + p] += xj * Sz_[p][j];
                agg_[n_levels_ - 1][ni][3 * Kf + p] += xj * Sdiv_[p][j];
            }
        }
    }

    // 每层 nonempty idx → global box 映射（用于后续各级）
    std::vector<std::vector<int>> ne2box_l(n_levels_);
    for (int l = 0; l < n_levels_; l++) {
        const LevelData& lev = levels_[l];
        ne2box_l[l].assign(lev.n_boxes_total, -1);
        for (int ni = 0; ni < (int)lev.nonempty.size(); ni++)
            ne2box_l[l][lev.nonempty[ni]] = ni;
    }

    // 2b: 粗层：子 → 父 插值 + 相位平移 + 累加
    for (int l = n_levels_ - 2; l >= 0; l--) {
        const LevelData& child = levels_[l + 1];
        const LevelData& parent = levels_[l];

        for (int ci = 0; ci < (int)child.nonempty.size(); ci++) {
            int cb = child.nonempty[ci];
            int pb = child.parent_box[cb];
            if (pb < 0) continue;
            int pni = ne2box_l[l][pb];
            if (pni < 0) continue;

            // 取子盒子聚合数据
            std::vector<Complex>& child_agg = agg_[l + 1][ci];

            // 步骤 A: 相位平移（子中心 → 父中心），用子层的 K
            std::vector<Complex> shifted;
            phase_shift(child, child.centers[cb], parent.centers[pb], child_agg, shifted);

            // 步骤 B: 双立方插值（子层 K → 父层 K）
            std::vector<Complex> interp;
            interp_bicubic(child, parent, shifted, interp);

            // 累加到父盒子
            std::vector<Complex>& parent_agg = agg_[l][pni];
            int Kp = parent.K;
            for (int q = 0; q < N_CH * Kp; q++)
                parent_agg[q] += interp[q];
        }
    }

    // =====================================================
    // ---- 阶段 3: M2L 转移 (每层) ----
    // =====================================================

    for (int l = 0; l < n_levels_; l++) {
        const LevelData& lev = levels_[l];
        int K = lev.K;

        for (const auto& mp : lev.m2l_pairs) {
            int tni = ne2box_l[l][mp.tgt];
            int sni = ne2box_l[l][mp.src];
            if (tni < 0 || sni < 0) continue;

            for (int p = 0; p < K; p++) {
                Complex T = lev.translator[mp.off + p];
                for (int ch = 0; ch < N_CH; ch++) {
                    down_[l][tni][ch * K + p] += T * agg_[l][sni][ch * K + p];
                }
            }
        }
    }

    // =====================================================
    // ---- 阶段 4: 下行 (解聚) ----
    // =====================================================

    // 4a: 粗 → 细：反插值 + 相位平移 + 分发到子
    for (int l = 0; l < n_levels_ - 1; l++) {
        const LevelData& parent = levels_[l];
        const LevelData& child = levels_[l + 1];
        int Kc = child.K;

        for (int pni = 0; pni < (int)parent.nonempty.size(); pni++) {
            int pb = parent.nonempty[pni];

            for (int cb : parent.child_boxes[pb]) {
                int cni = ne2box_l[l + 1][cb];
                if (cni < 0) continue;

                // 步骤 A: 双立方反插值（父层 K → 子层 K）
                std::vector<Complex> anterp;
                anterp_bicubic(parent, child, down_[l][pni], anterp);

                // 步骤 B: 相位平移（父中心 → 子中心），用子层的 K
                std::vector<Complex> shifted;
                phase_shift(child, parent.centers[pb], child.centers[cb], anterp, shifted);

                // 累加到子盒子
                std::vector<Complex>& child_down = down_[l + 1][cni];
                for (int q = 0; q < N_CH * Kc; q++)
                    child_down[q] += shifted[q];
            }
        }
    }

    // 4b: 最细层：方向图 → 基函数解聚
    {
        for (int ni = 0; ni < (int)finest.nonempty.size(); ni++) {
            int nb = finest.nonempty[ni];
            for (int i : finest.box_basis[nb]) {
                Complex sum(0, 0);
                for (int p = 0; p < Kf; p++) {
                    double w = finest.dirs[p].weight;
                    Complex vec = Tx_[p][i] * down_[n_levels_ - 1][ni][0 * Kf + p]
                                + Ty_[p][i] * down_[n_levels_ - 1][ni][1 * Kf + p]
                                + Tz_[p][i] * down_[n_levels_ - 1][ni][2 * Kf + p];
                    Complex div = Tdiv_[p][i] * down_[n_levels_ - 1][ni][3 * Kf + p];
                    sum += w * (vec + div);
                }
                y[i] += C_ * sum;
            }
        }
    }
}

// ===== 函数封装 =========================================================

std::function<void(const std::vector<Complex>&, std::vector<Complex>&)>
MLFMA::get_matvec_functor() const {
    return [this](const std::vector<Complex>& x, std::vector<Complex>& y) {
        this->matvec(x, y);
    };
}

// ===== 统计信息 =========================================================

void MLFMA::print_stats() const {
    int nnz = nf_ptr_.empty() ? 0 : nf_ptr_[n_basis_];
    std::cout << "\n=== MLFMA Statistics ===\n"
              << "  Unknowns:        " << n_basis_ << "\n"
              << "  Levels:          " << n_levels_ << "\n";
    for (int l = 0; l < n_levels_; l++) {
        const LevelData& lev = levels_[l];
        int n_m2l = (int)lev.m2l_pairs.size();
        std::cout << "  Level " << l << ": box=" << lev.box_size
                  << "m, L=" << lev.L << ", K=" << lev.K
                  << ", nonempty=" << lev.nonempty.size()
                  << ", M2L pairs=" << n_m2l << "\n";
    }
    std::cout << "  Near-field NZ:   " << nnz << " ("
              << 100.0 * nnz / (n_basis_ * n_basis_) << "%)\n";
}
