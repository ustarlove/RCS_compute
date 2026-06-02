#include "fmm.h"
#include "efie.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

FMM::FMM(const Mesh& mesh, const RWGData& rwg, double k0, double box_size)
    : mesh_(mesh), rwg_(rwg), k0_(k0), box_size_(box_size)
{
    eta_ = std::sqrt(MU0 / EPS0);
    n_basis_ = (int)rwg_.edges.size();
}



// ---------------------------------------------------------------------------
// 盒子索引 (从世界坐标到 3D 网格索引)
// ---------------------------------------------------------------------------
int FMM::box_idx(double x, double y, double z) const {
    int ix = (int)((x - xmin_) / box_size_);
    int iy = (int)((y - ymin_) / box_size_);
    int iz = (int)((z - zmin_) / box_size_);
    //边界裁剪，nx_ - 1 是 x 方向上的最大有效索引（从0开始）
    ix = std::max(0, std::min(ix, nx_ - 1));
    iy = std::max(0, std::min(iy, ny_ - 1));
    iz = std::max(0, std::min(iz, nz_ - 1));
    return iz * (nx_ * ny_) + iy * nx_ + ix;
}

// ---------------------------------------------------------------------------
// 构建 3D 网格并分配基函数到盒子
// ---------------------------------------------------------------------------
void FMM::build_boxes() {
    // 1. 计算包围盒
    xmin_ = ymin_ = zmin_ = 1e30;
    xmax_ = ymax_ = zmax_ = -1e30;
    for (const auto& v : mesh_.vertices) {
        xmin_ = std::min(xmin_, v.x); xmax_ = std::max(xmax_, v.x);
        ymin_ = std::min(ymin_, v.y); ymax_ = std::max(ymax_, v.y);
        zmin_ = std::min(zmin_, v.z); zmax_ = std::max(zmax_, v.z);
    }
    //添加一个很小的余量 eps，防止点正好落在盒子边界上
    double eps = 1e-6;
    xmin_ -= eps; xmax_ += eps;
    ymin_ -= eps; ymax_ += eps;
    zmin_ -= eps; zmax_ += eps;

    // 2. 网格尺寸
    nx_ = std::max(1, (int)std::ceil((xmax_ - xmin_) / box_size_));
    ny_ = std::max(1, (int)std::ceil((ymax_ - ymin_) / box_size_));
    nz_ = std::max(1, (int)std::ceil((zmax_ - zmin_) / box_size_));
    n_boxes_ = nx_ * ny_ * nz_;

    // 3. 初始化盒子
    boxes_.resize(n_boxes_);
    for (int iz = 0; iz < nz_; iz++) {
        for (int iy = 0; iy < ny_; iy++) {
            for (int ix = 0; ix < nx_; ix++) {
                int idx = iz * (nx_ * ny_) + iy * nx_ + ix;
                boxes_[idx].ix = ix;
                boxes_[idx].iy = iy;
                boxes_[idx].iz = iz;
                boxes_[idx].global_idx = idx;
                boxes_[idx].center = {
                    xmin_ + (ix + 0.5) * box_size_,
                    ymin_ + (iy + 0.5) * box_size_,
                    zmin_ + (iz + 0.5) * box_size_
                };
            }
        }
    }

    // 4. 分配基函数到盒子 (用边的中点)
    for (int j = 0; j < n_basis_; j++) {
        const Edge& e = rwg_.edges[j];
        const Vertex& v1 = mesh_.vertices[e.v1];
        const Vertex& v2 = mesh_.vertices[e.v2];
        double cx = 0.5 * (v1.x + v2.x);
        double cy = 0.5 * (v1.y + v2.y);
        double cz = 0.5 * (v1.z + v2.z);
        int b = box_idx(cx, cy, cz);
        boxes_[b].basis_idx.push_back(j);
    }

    // 5. 收集非空盒子
    nempty_.clear();
    for (int b = 0; b < n_boxes_; b++) {
        if (!boxes_[b].basis_idx.empty()) nempty_.push_back(b);
    }

    // 6. 构建近场/远场列表 (±1 邻域)
    for (int nb : nempty_) {
        FMMBox& box = boxes_[nb];
        int ix = box.ix, iy = box.iy, iz = box.iz;
        for (int dz = -1; dz <= 1; dz++) {
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    int ni = ix + dx, nj = iy + dy, nk = iz + dz;
                    if (ni < 0 || ni >= nx_ || nj < 0 || nj >= ny_ || nk < 0 || nk >= nz_)
                        continue;
                    int nbidx = nk * (nx_ * ny_) + nj * nx_ + ni;
                    if (!boxes_[nbidx].basis_idx.empty())
                        box.near_boxes.push_back(nbidx);
                }
            }
        }
        // 远场 = 所有其他非空盒子 - 近场
        for (int fb : nempty_) {
            if (fb == nb) continue;
            if (std::find(box.near_boxes.begin(), box.near_boxes.end(), fb)
                == box.near_boxes.end()) {
                box.far_boxes.push_back(fb);
            }
        }
    }

    std::cout << "FMM boxes: " << nempty_.size() << " non-empty / "
              << n_boxes_ << " total (" << nx_ << "x" << ny_ << "x" << nz_ << ")\n";
    int nonempty_total = 0;
    for (int nb : nempty_) nonempty_total += (int)boxes_[nb].basis_idx.size();
    //平均每个非空盒子中包含的基函数数量
    std::cout << "  avg basis/box = "
              << (double)nonempty_total / nempty_.size() << "\n";
}

// ---------------------------------------------------------------------------
// 计算单个基函数的辐射/接收方向图 (4 通道 × K 方向)
// ----------------------------------------------------------------------------
void FMM::compute_rwg_patterns(int edge_idx,
                               std::vector<Complex>& sx, std::vector<Complex>& sy,
                               std::vector<Complex>& sz, std::vector<Complex>& sdiv,
                               std::vector<Complex>& tx, std::vector<Complex>& ty,
                               std::vector<Complex>& tz, std::vector<Complex>& tdiv) const {
    const Edge& e = rwg_.edges[edge_idx];
    double l = e.length;
    double inv_k = 1.0 / k0_;

    // 确定该基函数所属的盒子
    const Vertex& v1 = mesh_.vertices[e.v1];
    const Vertex& v2 = mesh_.vertices[e.v2];
    double mid_x = 0.5 * (v1.x + v2.x);
    double mid_y = 0.5 * (v1.y + v2.y);
    double mid_z = 0.5 * (v1.z + v2.z);
    int b = box_idx(mid_x, mid_y, mid_z);
    Vertex c_box = boxes_[b].center;

    // 两个三角形的重心
    Vertex cp = mesh_.triangles[e.t1].center(mesh_.vertices);  // T⁺ centroid
    Vertex cm = mesh_.triangles[e.t2].center(mesh_.vertices);  // T⁻ centroid

    // T⁺ 上的 RWG 向量: f⁺ = l/(2A⁺)*(r - apex⁺), area scaling gives l/2 factor
    double vpx = 0.5 * l * (cp.x - e.apex1.x);
    double vpy = 0.5 * l * (cp.y - e.apex1.y);
    double vpz = 0.5 * l * (cp.z - e.apex1.z);

    // T⁻ 上的 RWG 向量: f⁻ = -l/(2A⁻)*(r - apex⁻), area scaling gives -l/2 factor
    double vmx = -0.5 * l * (cm.x - e.apex2.x);
    double vmy = -0.5 * l * (cm.y - e.apex2.y);
    double vmz = -0.5 * l * (cm.z - e.apex2.z);

    double div_l = l * inv_k;  // (l/k0) for normalized divergence

    for (int p = 0; p < K_; p++) {
        const Vertex& khat = dirs_[p].k_hat;

        // ---- 源方向图 (聚合用) ----
        // 相位: +j k k̂·(r - c_box)
        double ph_p = k0_ * dot(khat, cp - c_box);
        double ph_m = k0_ * dot(khat, cm - c_box);
        Complex ep(std::cos(ph_p), std::sin(ph_p));
        Complex em(std::cos(ph_m), std::sin(ph_m));

        sx[p]    = Complex(vpx, 0) * ep + Complex(vmx, 0) * em;
        sy[p]    = Complex(vpy, 0) * ep + Complex(vmy, 0) * em;
        sz[p]    = Complex(vpz, 0) * ep + Complex(vmz, 0) * em;
        sdiv[p]  = Complex(div_l, 0) * (ep - em);

        // ---- 测试方向图 (解聚用) = conj(源方向图) ⨉ channel sign ----
        // 测试相位: -j k k̂·(r - c_box), 即源相位的复共轭
        Complex epc(std::cos(ph_p), -std::sin(ph_p));  // conj(ep)
        Complex emc(std::cos(ph_m), -std::sin(ph_m));  // conj(em)

        tx[p]    = Complex(vpx, 0) * epc + Complex(vmx, 0) * emc;
        ty[p]    = Complex(vpy, 0) * epc + Complex(vmy, 0) * emc;
        tz[p]    = Complex(vpz, 0) * epc + Complex(vmz, 0) * emc;
        tdiv[p]  = Complex(-div_l, 0) * (epc - emc);  // -div/k 符号（这里把系数-1/k^2分给了聚合1/k，这里发散的-1/k)
    }
}

// ---------------------------------------------------------------------------
// 预计算所有辐射/接收方向图和转移矩阵
// ---------------------------------------------------------------------------
void FMM::precompute_patterns() {
    // 1. 确定截断数 L (a = 盒子外接球半径 = 半对角线)
    double a = 0.5 * box_size_ * std::sqrt(3.0);  // 半对角线
    double ka = k0_ * a;
    // 多余带宽公式: L ≈ ka + 1.8 (ka)^{1/3} log10(1/ε), ε ≈ 1e-3
    int L_calc = (int)(ka + 1.8 * std::cbrt(ka) * 3.0);
    L_ = std::max(2, std::min(L_calc, 15));

    // 2. 生成角谱积分
    dirs_ = Translator::generate_quadrature(L_);
    K_ = (int)dirs_.size();

    double total_w = 0;
    for (const auto& d : dirs_) total_w += d.weight;
    std::cout << "FMM: L=" << L_ << ", K=" << K_
              << " directions, Σw=" << total_w << " (4π=" << 4*PI << ")\n";

    // 3. 分配方向图存储 (方向优先: [K_][n_basis_])
    Sx_.assign(K_, std::vector<Complex>(n_basis_));
    Sy_.assign(K_, std::vector<Complex>(n_basis_));
    Sz_.assign(K_, std::vector<Complex>(n_basis_));
    Sdiv_.assign(K_, std::vector<Complex>(n_basis_));
    Tx_.assign(K_, std::vector<Complex>(n_basis_));
    Ty_.assign(K_, std::vector<Complex>(n_basis_));
    Tz_.assign(K_, std::vector<Complex>(n_basis_));
    Tdiv_.assign(K_, std::vector<Complex>(n_basis_));

    // 4. 对每个基函数计算方向图
    std::vector<Complex> sx(K_), sy(K_), sz(K_), sdiv(K_);
    std::vector<Complex> tx(K_), ty(K_), tz(K_), tdiv(K_);

    for (int j = 0; j < n_basis_; j++) {
        compute_rwg_patterns(j, sx, sy, sz, sdiv, tx, ty, tz, tdiv);
        for (int p = 0; p < K_; p++) {
            Sx_[p][j]    = sx[p];
            Sy_[p][j]    = sy[p];
            Sz_[p][j]    = sz[p];
            Sdiv_[p][j]  = sdiv[p];
            Tx_[p][j]    = tx[p];
            Ty_[p][j]    = ty[p];
            Tz_[p][j]    = tz[p];
            Tdiv_[p][j]  = tdiv[p];
        }
    }

    // 5. 预计算所有远场盒子对的转移矩阵（远场盒子对：目标盒子，源盒子）
    far_pairs_.clear();
    int total_pairs = 0;
    for (int nb : nempty_) {
        total_pairs += (int)boxes_[nb].far_boxes.size();
    }
    translator_.resize(total_pairs * K_);
    int pair_idx = 0;
    for (int nb : nempty_) {
        Vertex c_tgt = boxes_[nb].center;
        for (int fb : boxes_[nb].far_boxes) {
            Vertex c_src = boxes_[fb].center;
            Vertex R = c_tgt - c_src;
            double kR = k0_ * norm(R);
            Vertex R_hat = normalize(R);

            far_pairs_.push_back({nb, fb, pair_idx * K_});

            for (int p = 0; p < K_; p++) {
                double ct = dot(dirs_[p].k_hat, R_hat);
                if (ct > 1.0) ct = 1.0;
                if (ct < -1.0) ct = -1.0;
                translator_[pair_idx * K_ + p] = Translator::compute(kR, ct, L_);
            }
            pair_idx++;
        }
    }

    std::cout << "  precomputed " << far_pairs_.size() << " far-field box pairs\n";

    // 6. 常数 C = k²η / (16π²)
    C_ = k0_ * k0_ * eta_ / (16.0 * PI * PI);

    // 7. 分配迭代缓冲区
    agg_.resize(n_boxes_);
    trans_.resize(n_boxes_);
    for (int nb : nempty_) {
        agg_[nb].assign(K_, std::vector<Complex>(N_CH, Complex(0, 0)));
        trans_[nb].assign(K_, std::vector<Complex>(N_CH, Complex(0, 0)));
    }
}

// ---------------------------------------------------------------------------
// 填充近场稀疏矩阵 (CSR)
// ---------------------------------------------------------------------------
void FMM::fill_nearfield(const EFIE& efie) {
    // 第 1 遍: 统计每行非零元数量
    std::vector<int> row_count(n_basis_, 0);

    for (int nb : nempty_) {
        const FMMBox& box = boxes_[nb];
        for (int i_local = 0; i_local < (int)box.basis_idx.size(); i_local++) {
            int i = box.basis_idx[i_local];  // 测试基函数

            for (int nbr : box.near_boxes) {
                const FMMBox& nbr_box = boxes_[nbr];
                // 只计算 i <= j 的上三角 (矩阵对称性)
                for (int j_local = 0; j_local < (int)nbr_box.basis_idx.size(); j_local++) {
                    row_count[i]++;   //near_box的所有基函数数?
                }
            }
        }
    }

    // 构建行指针
    nf_ptr_.resize(n_basis_ + 1);
    nf_ptr_[0] = 0;
    for (int i = 0; i < n_basis_; i++) {
        nf_ptr_[i + 1] = nf_ptr_[i] + row_count[i];
    }
    //nnz:近场稀疏矩阵中非零元素的总数
    int nnz = nf_ptr_[n_basis_];
    nf_col_.resize(nnz);
    nf_val_.resize(nnz);

    std::cout << "Filling near-field sparse matrix: " << nnz << " nonzeros ("
              << 100.0 * nnz / (n_basis_ * n_basis_) << "% of dense)\n";

    // 第 2 遍: 填充列索引和值
    std::vector<int> cur_ptr(n_basis_, 0);  // 每行当前写入位置

    for (int nb : nempty_) {
        const FMMBox& box = boxes_[nb];
        for (int i_local = 0; i_local < (int)box.basis_idx.size(); i_local++) {
            int i = box.basis_idx[i_local];
            int base = nf_ptr_[i];

            for (int nbr : box.near_boxes) {
                const FMMBox& nbr_box = boxes_[nbr];
                for (int j_local = 0; j_local < (int)nbr_box.basis_idx.size(); j_local++) {
                    int j = nbr_box.basis_idx[j_local];
                    int pos = base + cur_ptr[i]++;
                    nf_col_[pos] = j;
                    nf_val_[pos] = efie.compute_interaction(i, j);
                }
            }
        }
    }

    // 按列索引排序每行 (便于缓存访问),列索引升序
    for (int i = 0; i < n_basis_; i++) {
        int beg = nf_ptr_[i];
        int end = nf_ptr_[i + 1];
        // 冒泡排序 (N 小，足够)
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

// ---------------------------------------------------------------------------
// FMM 矩阵-向量积: y = Z * x
// ---------------------------------------------------------------------------
void FMM::matvec(const std::vector<Complex>& x, std::vector<Complex>& y) const {
    y.assign(n_basis_, Complex(0, 0));

    // ---- 第 1 部分: 近场稀疏 matvec ----
    for (int i = 0; i < n_basis_; i++) {
        for (int ptr = nf_ptr_[i]; ptr < nf_ptr_[i + 1]; ptr++) {
            y[i] += nf_val_[ptr] * x[nf_col_[ptr]];
        }
    }

    // ---- 第 2 部分: FMM 远场 ----

    // 清理缓冲区
    for (int nb : nempty_) {
        for (int p = 0; p < K_; p++) {
            agg_[nb][p][0] = agg_[nb][p][1] = Complex(0, 0);
            agg_[nb][p][2] = agg_[nb][p][3] = Complex(0, 0);
            trans_[nb][p][0] = trans_[nb][p][1] = Complex(0, 0);
            trans_[nb][p][2] = trans_[nb][p][3] = Complex(0, 0);
        }
    }

    // 阶段 1: 聚合 — 每个非空盒子对所有方向求和源贡献
    for (int nb : nempty_) {
        for (int j : boxes_[nb].basis_idx) {
            Complex xj = x[j];
            if (std::norm(xj) < 1e-30) continue;
            for (int p = 0; p < K_; p++) {
                agg_[nb][p][0] += xj * Sx_[p][j];
                agg_[nb][p][1] += xj * Sy_[p][j];
                agg_[nb][p][2] += xj * Sz_[p][j];
                agg_[nb][p][3] += xj * Sdiv_[p][j];
            }
        }
    }

    // 阶段 2: 转移 — 远场盒子对之间的平面波传播
    for (size_t fp = 0; fp < far_pairs_.size(); fp++) {
        int tgt = far_pairs_[fp].tgt;
        int src = far_pairs_[fp].src;
        int off = far_pairs_[fp].off;

        for (int p = 0; p < K_; p++) {
            Complex T = translator_[off + p];
            for (int ch = 0; ch < N_CH; ch++) {
                trans_[tgt][p][ch] += T * agg_[src][p][ch];
            }
        }
    }

    // 阶段 3: 解聚 — 在每个测试基函数处叠加远场贡献
    for (int nb : nempty_) {
        for (int i : boxes_[nb].basis_idx) {
            Complex sum(0, 0);
            for (int p = 0; p < K_; p++) {
                double w = dirs_[p].weight;
                // 矢量贡献: Tx*Bx + Ty*By + Tz*Bz
                Complex vec = Tx_[p][i] * trans_[nb][p][0]
                            + Ty_[p][i] * trans_[nb][p][1]
                            + Tz_[p][i] * trans_[nb][p][2];
                // 散度贡献 (Tdiv 已带负号)
                Complex div = Tdiv_[p][i] * trans_[nb][p][3];
                sum += w * (vec + div);
            }
            y[i] += C_ * sum;   //近场和远场两部分贡献累加
        }
    }
}

// ---------------------------------------------------------------------------
std::function<void(const std::vector<Complex>&, std::vector<Complex>&)>
FMM::get_matvec_functor() const {
    return [this](const std::vector<Complex>& x, std::vector<Complex>& y) {
        this->matvec(x, y);
    };
}

// ---------------------------------------------------------------------------
void FMM::print_stats() const {
    int nnz = nf_ptr_.empty() ? 0 : nf_ptr_[n_basis_];
    int nfp = (int)far_pairs_.size();
    std::cout << "\n=== FMM Statistics ===\n"
              << "  Unknowns:        " << n_basis_ << "\n"
              << "  Box size:        " << box_size_ << " m\n"
              << "  Grid:            " << nx_ << "×" << ny_ << "×" << nz_ << "\n"
              << "  Non-empty boxes: " << nempty_.size() << "\n"
              << "  L / K:           " << L_ << " / " << K_ << "\n"
              << "  Near-field NZ:   " << nnz << " ("
              << 100.0 * nnz / (n_basis_ * n_basis_) << "%)\n"
              << "  Far-field pairs: " << nfp << "\n";
}
