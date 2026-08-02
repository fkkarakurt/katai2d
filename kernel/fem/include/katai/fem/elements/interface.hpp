#pragma once
// Interface (soil-structure interface) — zero-thickness (Goodman 1968) bilinear Coulomb
// joint (P2.4). Sits on a mesh edge (tri6 edge = 3 nodes). Two sides: the soil side (mesh
// nodes, base DOFs) and the structure side (coincident nodes, extra ux,uy via
// DofMap::add_extra_dof). 12 DOFs (6 soil + 6 structure).
//
// Constitutive law (PLAXIS interface; tension-positive): elastic τ=k_s·Δu_s, σ_n=k_n·Δu_n
// (tension cut-off σ_t); Coulomb shear τ_max = c_i − σ_n·tanφ_i; at yield plastic slip
// Δu_s^p (perfectly plastic, ψ_i=0).
// INTEGRATION: Newton-Cotes (nodal, Day & Potts 1994) — prevents the high-k_n oscillation,
// node pairs act like decoupled springs. Math: docs/references/interface-formulation.md.

#include <array>
#include <cmath>

#include <Eigen/Dense>

namespace katai::core::iface {  // 'interface' is a macro on Windows (#define interface struct) — use iface

inline constexpr int kNodeCount = 3;
inline constexpr int kDofCount = 12;     // 6 soil + 6 structure
inline constexpr int kPointCount = 3;    // Newton-Cotes (at the nodes)

using NodeCoords = Eigen::Matrix<double, 3, 2>;  // row i: (x_i, y_i), [A, B, mid]

struct InterfaceProps {
    double kn = 0.0;        // normal stiffness [kN/m³]
    double ks = 0.0;        // shear stiffness [kN/m³]
    double c_i = 0.0;       // interface cohesion (= R_inter·c_soil)
    double phi_i = 0.0;     // interface friction angle [rad] (tanφ_i = R_inter·tanφ_soil)
    double sigma_t = 0.0;   // tensile strength (tension cut-off); default 0
};

// PLAXIS virtual-thickness stiffness: G_i=R²G_soil, ν_i=0.45, E_oed,i=2G_i(1−ν_i)/(1−2ν_i),
// t_i=δ·avg_size. k_s=G_i/t_i, k_n=E_oed,i/t_i. (caller; helper.)
inline void interface_stiffness(double R_inter, double G_soil, double avg_size,
                                double vt_factor, double& kn, double& ks) {
    const double Gi = R_inter * R_inter * G_soil;
    const double nu_i = 0.45;
    const double Eoed_i = 2.0 * Gi * (1.0 - nu_i) / (1.0 - 2.0 * nu_i);
    const double ti = vt_factor * avg_size;
    ks = Gi / ti;
    kn = Eoed_i / ti;
}

inline Eigen::Vector3d shape_deriv_xi(double xi) {  // [A(ξ=−1), B(ξ=+1), mid(ξ=0)]
    return {xi - 0.5, xi + 0.5, -2.0 * xi};
}

// Newton-Cotes point q: local node index (soil_nodes/struct order [A,B,mid]), weight, ξ.
struct NCPoint { int node; double w; double xi; };
inline std::array<NCPoint, 3> nc_points() {
    return {{ {0, 1.0 / 3.0, -1.0},    // A
              {2, 4.0 / 3.0, 0.0},     // mid
              {1, 1.0 / 3.0, 1.0} }};  // B
}

// Edge Jacobian J=ds/dξ and unit tangent (c,s) at a point.
struct EdgeFrame { double J, c, s; };
inline EdgeFrame edge_frame(const NodeCoords& X, double xi) {
    const Eigen::Vector3d d = shape_deriv_xi(xi);
    const double dx = d.dot(X.col(0)), dy = d.dot(X.col(1));
    const double J = std::sqrt(dx * dx + dy * dy);
    return {J, dx / J, dy / J};
}

// ---- 5-NODE interface (tri15 edge, 5 nodes quarter-spaced) — the PLAXIS 15-node soil counterpart ----
inline constexpr int kNodeCount5 = 5;
inline constexpr int kDofCount5 = 20;     // 10 soil + 10 structure
inline constexpr int kPointCount5 = 5;    // Newton-Cotes (at the 5 nodes)
using NodeCoords5 = Eigen::Matrix<double, 5, 2>;

// Quartic Lagrange derivative at ξ=[−1,−0.5,0,0.5,1] (for the edge_frame5 Jacobian).
inline Eigen::Matrix<double, 5, 1> shape_deriv5_xi(double xi) {
    static const double xn[5] = {-1.0, -0.5, 0.0, 0.5, 1.0};
    Eigen::Matrix<double, 5, 1> dN;
    for (int i = 0; i < 5; ++i) {
        double s = 0.0;
        for (int k = 0; k < 5; ++k)
            if (k != i) {
                double p = 1.0 / (xn[i] - xn[k]);
                for (int j = 0; j < 5; ++j)
                    if (j != i && j != k) p *= (xi - xn[j]) / (xn[i] - xn[j]);
                s += p;
            }
        dN(i) = s;
    }
    return dN;
}
// 5-point closed Newton-Cotes (Boole) — at the nodes, weights {7,32,12,32,7}/45 (Σ=2 =
// the [−1,1] length). Node pairs DECOUPLE (at each NC point only that node's N_i=1) → no
// high-kn oscillation (Day&Potts).
inline std::array<NCPoint, 5> nc_points5() {
    return {{ {0, 7.0 / 45.0, -1.0}, {1, 32.0 / 45.0, -0.5}, {2, 12.0 / 45.0, 0.0},
              {3, 32.0 / 45.0, 0.5}, {4, 7.0 / 45.0, 1.0} }};
}
inline EdgeFrame edge_frame5(const NodeCoords5& X, double xi) {
    const Eigen::Matrix<double, 5, 1> d = shape_deriv5_xi(xi);
    const double dx = d.dot(X.col(0)), dy = d.dot(X.col(1));
    const double J = std::sqrt(dx * dx + dy * dy);
    return {J, dx / J, dy / J};
}

// Coulomb return mapping (per node pair). Δu_s,Δu_n local relative displacement; slip_p_c
// committed plastic slip. sigma_n0: the INITIAL normal stress (K0 install / staged
// carry-over; tension-pos, compression<0) → total σ_n = σ_n0 + k_n·Δu_n; carries the K0
// horizontal stress at Δu_n=0 in a discontinuous (split) wall mesh (otherwise a spurious
// "installation" movement; see docs/references/interface-formulation.md §6).
// Returns: τ, σ_n, the tangents D_s/D_n, the new plastic slip.
struct Return { double tau, sigma_n, Ds, Dn, slip_p_new; };
inline Return coulomb_return(const InterfaceProps& p, double du_s, double du_n, double slip_p_c,
                             double sigma_n0 = 0.0) {
    double sigma_n = sigma_n0 + p.kn * du_n, Dn = p.kn;
    if (sigma_n > p.sigma_t) { sigma_n = p.sigma_t; Dn = 0.0; }  // tension cut-off (separation)
    const double tau_tr = p.ks * (du_s - slip_p_c);
    const double tau_max = std::max(0.0, p.c_i - sigma_n * std::tan(p.phi_i));
    if (std::fabs(tau_tr) <= tau_max)
        return {tau_tr, sigma_n, p.ks, Dn, slip_p_c};            // elastic (stick)
    const double tau = std::copysign(tau_max, tau_tr);
    return {tau, sigma_n, 0.0, Dn, du_s - tau / p.ks};           // slip (plastic)
}

}  // namespace katai::core::iface
