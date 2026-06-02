#include "efie.h"
#include <cmath>
#include <iostream>

// ---------------------------------------------------------------------------
// 三角形上的 7 点 Gauss 积分规则（5 阶代数精度，Dunavant 1985）
// ---------------------------------------------------------------------------
namespace {

struct GaussPt {
    double u, v, w;
    double weight;
};

std::vector<GaussPt> getGauss7() {
    const double a1 = 0.05971587178976981;
    const double b1 = 0.47014206410511505;
    const double a2 = 0.7974269853530872;
    const double b2 = 0.10128650732345633;
    const double w1 = 0.13239415278850616;
    const double w2 = 0.12593918054482713;
    return {
        {1.0/3.0, 1.0/3.0, 1.0/3.0, 0.22500000000000000},
        {a1, b1, b1, w1}, {b1, a1, b1, w1}, {b1, b1, a1, w1},
        {a2, b2, b2, w2}, {b2, a2, b2, w2}, {b2, b2, a2, w2},
    };
}

//将重心坐标转换为实际物理坐标
Vertex baryToCart(double u, double v, double w,
                  const Triangle& tri, const std::vector<Vertex>& verts) {
    const Vertex& p0 = verts[tri.v0];
    const Vertex& p1 = verts[tri.v1];
    const Vertex& p2 = verts[tri.v2];
    return {u*p0.x + v*p1.x + w*p2.x,
            u*p0.y + v*p1.y + w*p2.y,
            u*p0.z + v*p1.z + w*p2.z};
}

//计算自由空间标量亥姆霍兹方程的格林函数，r场点,rp源点
Complex green_kernel(double k, const Vertex& r, const Vertex& rp) {
    double dx = r.x - rp.x, dy = r.y - rp.y, dz = r.z - rp.z;
    double R = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (R < 1e-12) R = 1e-12;
    return std::exp(Complex(0, -k * R)) / (4.0 * PI * R);
}

//计算平滑化的格林函数
Complex green_smooth(double k, const Vertex& r, const Vertex& rp) {
    double dx = r.x - rp.x, dy = r.y - rp.y, dz = r.z - rp.z;
    double R = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (R < 1e-12) return Complex(0, -k) / (4.0 * PI);   //-jk/4π
    return (std::exp(Complex(0, -k * R)) - 1.0) / (4.0 * PI * R);
}

// ===== 解析静态奇异积分：∫_T 1/|r-r'| dS' =====
// 恒等式 1/R = ∇'²_s R（平坦三角形） + 表面散度定理 →
//   ∫_T 1/R dS' = ∮_{∂T} n̂_out·(r'-r)/R dl' = Σ_{i=1}^3 P_i^0·f_i
// 其中 P_i^0 = n̂_out·(A-r), f_i = ln((R_A+R_B+L)/(R_A+R_B-L))
// 观测点 r 须在三角形 T 的平面内（Gauss 点满足此条件）。
static double integrate_1_over_R(const Triangle& T, const std::vector<Vertex>& verts,
                                  const Vertex& r) {
    const Vertex& v0 = verts[T.v0];
    const Vertex& v1 = verts[T.v1];
    const Vertex& v2 = verts[T.v2];

    // 三角形外法向（自动检测朝向：外法向与重心同向）
    double e1x = v1.x-v0.x, e1y = v1.y-v0.y, e1z = v1.z-v0.z;
    double e2x = v2.x-v0.x, e2y = v2.y-v0.y, e2z = v2.z-v0.z;
    double nx = e1y*e2z - e1z*e2y;
    double ny = e1z*e2x - e1x*e2z;
    double nz = e1x*e2y - e1y*e2x;
    double nrm = std::sqrt(nx*nx + ny*ny + nz*nz);
    nx /= nrm; ny /= nrm; nz /= nrm;
    double cx = (v0.x+v1.x+v2.x)/3.0;
    double cy = (v0.y+v1.y+v2.y)/3.0;
    double cz = (v0.z+v1.z+v2.z)/3.0;
    if (nx*cx + ny*cy + nz*cz < 0.0) { nx=-nx; ny=-ny; nz=-nz; }

    const Vertex* va[3] = {&v0, &v1, &v2};
    const Vertex* vb[3] = {&v1, &v2, &v0};

    double result = 0.0;
    for (int k = 0; k < 3; k++) {
        const Vertex& A = *va[k];
        const Vertex& B = *vb[k];

        // 边方向及长度
        double tx = B.x-A.x, ty = B.y-A.y, tz = B.z-A.z;
        double L = std::sqrt(tx*tx + ty*ty + tz*tz);
        tx /= L; ty /= L; tz /= L;

        // 边外法向 n_out = t × n_tri
        double nox = ty*nz - tz*ny;
        double noy = tz*nx - tx*nz;
        double noz = tx*ny - ty*nx;

        // P_i^0 = n_out · (A - r)
        double ux = A.x-r.x, uy = A.y-r.y, uz = A.z-r.z;
        double P0 = nox*ux + noy*uy + noz*uz;
        P0 = std::abs(P0);
        if (std::abs(P0) < 1e-15) continue;

        double RA = std::sqrt(ux*ux + uy*uy + uz*uz);
        double vx = B.x-r.x, vy = B.y-r.y, vz = B.z-r.z;
        double RB = std::sqrt(vx*vx + vy*vy + vz*vz);

        // f_i = ln((R_A+R_B+L)/(R_A+R_B-L))
        double sum = RA + RB;
        double denom = sum - L;
        if (denom < 1e-15 * L) continue;
        result += P0 * std::log((sum + L) / denom);
    }
    return result;
}

} // namespace

// ---------------------------------------------------------------------------
// EFIE 构造函数
// ---------------------------------------------------------------------------
EFIE::EFIE(const Mesh& mesh, const RWGData& rwg, double k0)
    : mesh_(mesh), rwg_(rwg), k0_(k0) {
    lambda_ = 2.0 * PI / k0;
    n_edges_ = rwg_.edges.size();
}

// ---------------------------------------------------------------------------
// 计算一对三角形 (T_m, T_n) 上的双重 RWG-EFIE 积分（非奇异情形）
//
// 对 r ∈ T_m, r' ∈ T_n:
//   kernel = jkη [ f_i(r)·f_j(r') - (1/k²)(∇·f_i)(∇·f_j) ] G(R)
//
// f_i(r) = coeff_m * (r - apex_m)     (矢量)
// f_j(r') = coeff_n * (r' - apex_n)
// ∇·f_i = div_m, ∇·f_n = div_n       (常量)
//Am，An是面积，把标准三角形上的积分结果缩放回物理三角形上的积分结果。
//固定对应端点
// ---------------------------------------------------------------------------
static Complex integrateTrianglePair(
    const Triangle& Tm, const Triangle& Tn,
    const std::vector<Vertex>& verts,
    double Am, double An,
    double coeff_m, double div_m, const Vertex& apex_m,
    double coeff_n, double div_n, const Vertex& apex_n,
    double k, double eta, bool is_singular)
{
    const auto& gauss = getGauss7();
    const double jk_eta = k * eta;          // jkη → 复数用到时乘以 j
    const double inv_k2 = 1.0 / (k * k);    // 1/k²

    Complex acc(0, 0);
    for (const auto& gp_m : gauss) {
        Vertex rm = baryToCart(gp_m.u, gp_m.v, gp_m.w, Tm, verts);
        // f_i(rm) = coeff_m * (rm - apex_m)
        double f_ix = coeff_m * (rm.x - apex_m.x);
        double f_iy = coeff_m * (rm.y - apex_m.y);
        double f_iz = coeff_m * (rm.z - apex_m.z);

        for (const auto& gp_n : gauss) {
            Vertex rn = baryToCart(gp_n.u, gp_n.v, gp_n.w, Tn, verts);
            // f_j(rn) = coeff_n * (rn - apex_n)
            double f_jx = coeff_n * (rn.x - apex_n.x);
            double f_jy = coeff_n * (rn.y - apex_n.y);
            double f_jz = coeff_n * (rn.z - apex_n.z);

            double f_dot = f_ix*f_jx + f_iy*f_jy + f_iz*f_jz;  // f_i·f_j
            double div_prod = div_m * div_n;                     // (∇·f_i)(∇·f_j)

            double qw = gp_m.weight * gp_n.weight;

            Complex G;
            if (is_singular) {
                // 奇异情形：分离静态部分与平滑部分
                // G(k,R) = 1/(4πR) + G_smooth(k,R)
                // 平滑部分用 Gauss 积分；静态部分在外部用解析线积分处理
                G = green_smooth(k, rm, rn);
            } else {
                G = green_kernel(k, rm, rn);
            }

            // kernel = jkη [f_i·f_j - (1/k²) div_i·div_j] * G
            Complex kernel = Complex(0, jk_eta) *
                Complex(f_dot - inv_k2 * div_prod, 0) * G;

            acc += qw * kernel;
        }
    }

    return acc * (Am * An);
}

// ---------------------------------------------------------------------------
// 计算两条 RWG 边 (edge_i, edge_j) 之间的完整 EFIE 相互作用 Z_ij
//
// Z_ij = Σ_{p=±} Σ_{q=±} Z_ij^{pq}
// 每一项 Z_ij^{pq} 是在三角形对 (T_i^p, T_j^q) 上的双重积分
// ---------------------------------------------------------------------------
Complex EFIE::compute_interaction(int edge_i, int edge_j) const {
    const Edge& ei = rwg_.edges[edge_i];
    const Edge& ej = rwg_.edges[edge_j];

    const double k = k0_;
    const double eta = std::sqrt(MU0 / EPS0);  // ≈ 376.73 Ω

    Complex Z_ij(0, 0);

    // 每条边有两个三角形: t1 (T⁺) 和 t2 (T⁻)
    struct TriInfo {
        int idx;
        double coeff;   // ± l/(2A)
        double div;     // ± l/A = ± 2*coeff (绝对值)
        Vertex apex;
        double area;
    };

    //每条边对应的两个三角形
    TriInfo ti_list[2] = {
        { ei.t1,  rwg_.rwg_coeff1[edge_i],   2.0 * rwg_.rwg_coeff1[edge_i], ei.apex1,
          mesh_.triangles[ei.t1].area(mesh_.vertices) },
        { ei.t2, -rwg_.rwg_coeff2[edge_i], -2.0 * rwg_.rwg_coeff2[edge_i], ei.apex2,
          mesh_.triangles[ei.t2].area(mesh_.vertices) },
    };

    TriInfo tj_list[2] = {
        { ej.t1,  rwg_.rwg_coeff1[edge_j],   2.0 * rwg_.rwg_coeff1[edge_j], ej.apex1,
          mesh_.triangles[ej.t1].area(mesh_.vertices) },
        { ej.t2, -rwg_.rwg_coeff2[edge_j], -2.0 * rwg_.rwg_coeff2[edge_j], ej.apex2,
          mesh_.triangles[ej.t2].area(mesh_.vertices) },
    };

    for (int pi = 0; pi < 2; pi++) {
        const TriInfo& ti = ti_list[pi];
        const Triangle& Tri = mesh_.triangles[ti.idx];

        for (int pj = 0; pj < 2; pj++) {
            const TriInfo& tj = tj_list[pj];
            const Triangle& Trj = mesh_.triangles[tj.idx];

            bool same_tri = (ti.idx == tj.idx);

            // 双重面积分（非奇异部分，或平滑部分）
            Complex contrib = integrateTrianglePair(
                Tri, Trj, mesh_.vertices, ti.area, tj.area,
                ti.coeff, ti.div, ti.apex,
                tj.coeff, tj.div, tj.apex,
                k, eta, same_tri);

            if (same_tri) {
                // ─── 静态奇异修正：解析线积分 + 中心点近似 ───
                // ∫∫ K/(4πR) dS' dS ≈ K(center)/(4π) · I_double
                // 内层 ∫_T 1/R dS' 用解析线积分（P_i^0·f_i），
                // 外层用 7 点 Gauss 求 I_double = ∫_T ∫_T 1/R dS' dS

                const double jk_eta = k * eta;
                const double inv_k2 = 1.0 / (k * k);

                // 计算 I_double = ∫_T ∫_T 1/|r-r'| dS' dS
                double I_double = 0.0;
                const auto& gauss = getGauss7();
                for (const auto& gp : gauss) {
                    Vertex rm = baryToCart(gp.u, gp.v, gp.w, Tri, mesh_.vertices);
                    I_double += gp.weight * integrate_1_over_R(Tri, mesh_.vertices, rm);
                }
                I_double *= ti.area;  // 外层 Gauss 面积缩放

                // 核函数在三角形中心点求值
                Vertex rc = Tri.center(mesh_.vertices);
                double f_ix = ti.coeff * (rc.x - ti.apex.x);
                double f_iy = ti.coeff * (rc.y - ti.apex.y);
                double f_iz = ti.coeff * (rc.z - ti.apex.z);
                double f_jx = tj.coeff * (rc.x - tj.apex.x);
                double f_jy = tj.coeff * (rc.y - tj.apex.y);
                double f_jz = tj.coeff * (rc.z - tj.apex.z);
                double fdot = f_ix*f_jx + f_iy*f_jy + f_iz*f_jz;
                double ker_val = fdot - inv_k2 * ti.div * tj.div;

                Complex static_contrib = Complex(0, jk_eta) *
                    Complex(ker_val / (4.0 * PI), 0) * Complex(I_double, 0);
                contrib = contrib + static_contrib;
            }

            Z_ij = Z_ij + contrib;
        }
    }

    return Z_ij;
}

// ---------------------------------------------------------------------------
// 填充阻抗矩阵 Z (N × N)
// ---------------------------------------------------------------------------
void EFIE::fill_matrix(std::vector<std::vector<Complex>>& Z) {
    int N = n_edges_;
    Z.assign(N, std::vector<Complex>(N, Complex(0,0)));

    std::cout << "Filling impedance matrix " << N << "x" << N
              << " with 7-pt Gauss quadrature ...\n";

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Z[i][j] = compute_interaction(i, j);
        }
        //每处理完 50 行，打印一次进度信息，并且在处理完最后一行时也打印一次
        if ((i+1) % 50 == 0 || i == N-1) {
            std::cout << "  row " << i+1 << " / " << N << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// 平面波入射电场: E(r) = E0_pol * exp(-j k_inc·r)
// ---------------------------------------------------------------------------
// Complex EFIE::plane_wave_e(const Vertex& r, const Vertex& k_inc,
//                            const Vertex& E0_pol) const {
//     double phase = k_inc.x * r.x + k_inc.y * r.y + k_inc.z * r.z;
//     Complex exp_phase = std::exp(Complex(0, -phase));
//     // 入射电场方向与极化相同（远场平面波近似）
//     return Complex(E0_pol.x, E0_pol.y) * exp_phase;  // 简化：用 y 分量
//     // 实际应为矢量，但对标量激励只需 E·f 的积分
// }

// ---------------------------------------------------------------------------
// 填充激励向量 b: b_i = ∫_S f_i(r) · E_inc(r) dS
// 使用 Gauss 积分逼近
// ---------------------------------------------------------------------------
void EFIE::fill_vector(const Vertex& k_inc, const Vertex& E0_pol,
                       std::vector<Complex>& b) {
    int N = n_edges_;
    b.assign(N, Complex(0,0));
    const auto& gauss = getGauss7();

    for (int i = 0; i < N; i++) {
        const Edge& e = rwg_.edges[i];

        Complex acc(0, 0);

        // 在 T⁺ 上的贡献: f⁺(r) = coeff1 * (r - apex1), 面积分 × A⁺
        {
            const Triangle& Tp = mesh_.triangles[e.t1];
            double Ap = Tp.area(mesh_.vertices);
            double coeff = rwg_.rwg_coeff1[i];  // l/(2A⁺)
            for (const auto& gp : gauss) {
                Vertex r = baryToCart(gp.u, gp.v, gp.w, Tp, mesh_.vertices);
                // f(r) = coeff * (r - apex)
                double fx = coeff * (r.x - e.apex1.x);
                double fy = coeff * (r.y - e.apex1.y);
                double fz = coeff * (r.z - e.apex1.z);

                // E_inc(r) 的 x 分量与 f 的点积
                // E0_pol 定义为 {E0x, E0y, E0z}
                double phase = k_inc.x*r.x + k_inc.y*r.y + k_inc.z*r.z;
                Complex exp_phase = std::exp(Complex(0, -phase));
                double E_dot_f = E0_pol.x*fx + E0_pol.y*fy + E0_pol.z*fz;
                acc = acc + gp.weight * Complex(E_dot_f, 0) * exp_phase;
            }
            acc = acc * Ap;
        }

        // 在 T⁻ 上的贡献: f⁻(r) = -coeff2 * (r - apex2)
        {
            const Triangle& Tm = mesh_.triangles[e.t2];
            double Am = Tm.area(mesh_.vertices);
            double coeff = -rwg_.rwg_coeff2[i];  // -l/(2A⁻)
            for (const auto& gp : gauss) {
                Vertex r = baryToCart(gp.u, gp.v, gp.w, Tm, mesh_.vertices);
                double fx = coeff * (r.x - e.apex2.x);
                double fy = coeff * (r.y - e.apex2.y);
                double fz = coeff * (r.z - e.apex2.z);

                double phase = k_inc.x*r.x + k_inc.y*r.y + k_inc.z*r.z;
                Complex exp_phase = std::exp(Complex(0, -phase));
                double E_dot_f = E0_pol.x*fx + E0_pol.y*fy + E0_pol.z*fz;
                acc = acc + gp.weight * Complex(E_dot_f, 0) * exp_phase;
            }
            acc = acc * Am;
        }

        b[i] = acc;
    }
}
