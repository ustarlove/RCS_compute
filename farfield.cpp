#include "farfield.h"
#include <cmath>

namespace {

struct GaussPt { double u, v, w, weight; };

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

Vertex baryToCart(double u, double v, double w,
                  const Triangle& tri, const std::vector<Vertex>& verts) {
    const Vertex& p0 = verts[tri.v0];
    const Vertex& p1 = verts[tri.v1];
    const Vertex& p2 = verts[tri.v2];
    return {u*p0.x + v*p1.x + w*p2.x,
            u*p0.y + v*p1.y + w*p2.y,
            u*p0.z + v*p1.z + w*p2.z};
}

// 辐射矢量: N(θ,φ) = Σ I_n ∫ f_n(r) exp(j k k̂·r') dS
void compute_N(double theta, double phi, double k0,
               const Mesh& mesh, const RWGData& rwg,
               const std::vector<Complex>& I,
               Complex& Nx, Complex& Ny, Complex& Nz) {
    Nx = Ny = Nz = Complex(0,0);
    const auto& gauss = getGauss7();

    double kx = k0 * sin(theta) * cos(phi);
    double ky = k0 * sin(theta) * sin(phi);
    double kz = k0 * cos(theta);

    for (size_t e = 0; e < rwg.edges.size(); e++) {
        const Edge& ed = rwg.edges[e];
        Complex Ie = I[e];

        // T⁺ 贡献: f⁺(r) = coeff⁺ * (r - apex1),  coeff⁺ = l/(2A⁺)
        {
            const Triangle& Tp = mesh.triangles[ed.t1];
            double Ap = Tp.area(mesh.vertices);
            double coeff = rwg.rwg_coeff1[e];
            Complex sx(0,0), sy(0,0), sz(0,0);

            for (const auto& gp : gauss) {
                Vertex r = baryToCart(gp.u, gp.v, gp.w, Tp, mesh.vertices);
                double fx = coeff * (r.x - ed.apex1.x);
                double fy = coeff * (r.y - ed.apex1.y);
                double fz = coeff * (r.z - ed.apex1.z);
                double phase = kx*r.x + ky*r.y + kz*r.z;
                Complex eph = std::exp(Complex(0, phase));
                sx = sx + gp.weight * Complex(fx, 0) * eph;
                sy = sy + gp.weight * Complex(fy, 0) * eph;
                sz = sz + gp.weight * Complex(fz, 0) * eph;
            }
            Nx = Nx + Ie * sx * Ap;
            Ny = Ny + Ie * sy * Ap;
            Nz = Nz + Ie * sz * Ap;
        }

        // T⁻ 贡献: f⁻(r) = coeff⁻ * (r - apex2),  coeff⁻ = -l/(2A⁻)
        {
            const Triangle& Tm = mesh.triangles[ed.t2];
            double Am = Tm.area(mesh.vertices);
            double coeff = -rwg.rwg_coeff2[e];
            Complex sx(0,0), sy(0,0), sz(0,0);

            for (const auto& gp : gauss) {
                Vertex r = baryToCart(gp.u, gp.v, gp.w, Tm, mesh.vertices);
                double fx = coeff * (r.x - ed.apex2.x);
                double fy = coeff * (r.y - ed.apex2.y);
                double fz = coeff * (r.z - ed.apex2.z);
                double phase = kx*r.x + ky*r.y + kz*r.z;
                Complex eph = std::exp(Complex(0, phase));
                sx = sx + gp.weight * Complex(fx, 0) * eph;
                sy = sy + gp.weight * Complex(fy, 0) * eph;
                sz = sz + gp.weight * Complex(fz, 0) * eph;
            }
            Nx = Nx + Ie * sx * Am;
            Ny = Ny + Ie * sy * Am;
            Nz = Nz + Ie * sz * Am;
        }
    }
}

} // namespace

FarField::FarField(const Mesh& mesh, const RWGData& rwg, double k0)
    : mesh_(mesh), rwg_(rwg), k0_(k0) {}

Complex FarField::compute_Etheta(double theta, double phi,
                                  const std::vector<Complex>& I) const {
    Complex Nx, Ny, Nz;
    compute_N(theta, phi, k0_, mesh_, rwg_, I, Nx, Ny, Nz);

    double ct = cos(theta), st = sin(theta);
    double cp = cos(phi), sp = sin(phi);

    // ê_θ = (cosθ cosφ, cosθ sinφ, -sinθ)
    Complex N_theta = Nx * (ct * cp) + Ny * (ct * sp) + Nz * (-st);

    // E_θ = -j k η / (4π) * N_θ
    return Complex(0, -k0_ * Z0) / (4.0 * PI) * N_theta;
}

Complex FarField::compute_Ephi(double theta, double phi,
                                const std::vector<Complex>& I) const {
    Complex Nx, Ny, Nz;
    compute_N(theta, phi, k0_, mesh_, rwg_, I, Nx, Ny, Nz);

    double cp = cos(phi), sp = sin(phi);

    // ê_φ = (-sinφ, cosφ, 0)
    Complex N_phi = Nx * (-sp) + Ny * cp;

    return Complex(0, -k0_ * Z0) / (4.0 * PI) * N_phi;
}

//先计算单位入射场，对于其他任意强度的入射场，只需简单缩放即可
double FarField::compute_rcs(double theta, double phi,
                                         const std::vector<Complex>& I) const {
    // σ = 4π (|E_θ|² + |E_φ|²) / |E_inc|²   (|E_inc| = 1 V/m)
    Complex E_th = compute_Etheta(theta, phi, I);
    Complex E_ph = compute_Ephi(theta, phi, I);

    double rcs = 4.0 * PI * (std::norm(E_th) + std::norm(E_ph));
    if (rcs < 1e-30) rcs = 1e-30;
    return 10.0 * log10(rcs);
}
