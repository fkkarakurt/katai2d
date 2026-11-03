#pragma once
// Interface (soil-structure interface) — zero-thickness (Goodman 1968) bilinear Coulomb
// joint (P2.4). Sits on a mesh edge (tri6 edge = 3 nodes). Two sides: the soil side (mesh
// nodes, base DOFs) and the structure side (coincident nodes, extra ux,uy via
// DofMap::add_extra_dof). 12 DOFs (6 soil + 6 structure).
//
// Constitutive law (tension-positive): elastic τ=k_s·Δu_s, σ_n=k_n·Δu_n
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
    // The cohesion against a shear of negative sense, where the strength has a sense (an inclined joint
    // in NGI-ADP: active one way, passive the other). Below zero: c_i, either way.
    double c_neg = -1.0;
    // The stiffnesses go as (normal stress / k_pref)^k_power; kn and ks are their values at k_pref. 0: constant.
    double k_power = 0.0;
    double k_pref = 100.0;
};

// What the stiffnesses are multiplied by at a point installed with the normal stress sigma_n0 (tension positive), floored
// at a tenth of the reference stress. A new contact, a joint with no ground across it that takes over no stress (a lining
// in the dug-out ring), has no stress to read and stands at the reference stiffness.
inline double stiffness_factor(const InterfaceProps& p, double sigma_n0) {
    if (p.k_power <= 0.0 || sigma_n0 >= 0.0) return 1.0;
    return std::pow(std::max(-sigma_n0, 0.1 * p.k_pref) / p.k_pref, p.k_power);
}

// The cohesion a shear stress of this sign acts against.
inline double cohesion_for(const InterfaceProps& p, double tau) {
    return tau < 0.0 && p.c_neg >= 0.0 ? p.c_neg : p.c_i;
}

// The strength of a joint in NGI-ADP: the traction criterion of the model on a line of tangent (c, s). In the
// plane strain failure ellipse, with no volume change, the shear stress on the line is limited by
//     T = R (-/+ R_A sin 2b + sqrt(R_B^2 sin^2 2b + R_D^2 cos^2 2b)),
// R_A = (suA - suP)/2, R_B = (suA + suP)/2, R_D = suDSS: suDSS on a level or a vertical line, whichever way it
// shears, and at 45 degrees suP where the upper body goes up the slope and suA where it comes down. Which body
// is upper, in the sign of the joint's own shear: the relative slip of the body on the +n side, n = (-s, c), along
// the tangent, is the slip of the upper body along the tangent read towards +x, for either way the line was
// drawn, and sin 2b does not change when the line is reversed. So a positive shear on a line rising to the right is the upper body
// going up the slope.
inline void ngi_adp_joint_strength(double R, double suA, double suP, double suDSS, double c, double s,
                                   double& c_pos, double& c_neg) {
    const double s2 = 2.0 * s * c, c2 = c * c - s * s;
    const double RA = 0.5 * (suA - suP), RB = 0.5 * (suA + suP);
    const double root = std::sqrt(RB * RB * s2 * s2 + suDSS * suDSS * c2 * c2);
    c_pos = R * (root - RA * s2);
    c_neg = R * (root + RA * s2);
}

// Virtual-thickness stiffness: G_i=R²G_soil, ν_i=0.45, E_oed,i=2G_i(1−ν_i)/(1−2ν_i),
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

// ---- 5-NODE interface (tri15 edge, 5 nodes quarter-spaced) — the 15-node soil counterpart ----
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
// `pw` is the excess pore pressure of the ground beside the joint (tension positive, so a
// pressure is negative and a suction positive). The normal stress the joint transmits is a total
// one -- beside undrained ground it carries that ground's excess pore pressure with it -- and
// friction acts on the EFFECTIVE normal stress, sigma_n - pw. With pw left out, a clay unloaded
// undrained behind a wall, whose pore water goes into suction, had the suction counted as a loss
// of contact stress and its joint slipped at a strength the ground beside it did not have.
// The gap opens on the same stress: when the skeleton's contact, sigma_n - pw, reaches the tensile
// strength. The pore water's share stays across it -- a suction beside the joint still holds the two
// faces together -- so an opened joint transmits sigma_t + pw. Cut off on sigma_n alone, a joint
// beside clay in suction opened while its skeleton was still pressed, and lost its friction early.
inline Return coulomb_return(const InterfaceProps& p, double du_s, double du_n, double slip_p_c,
                             double sigma_n0 = 0.0, double pw = 0.0) {
    const double kf = stiffness_factor(p, sigma_n0), kn = p.kn * kf, ks = p.ks * kf;
    double sigma_n = sigma_n0 + kn * du_n, Dn = kn;
    if (sigma_n - pw > p.sigma_t) { sigma_n = p.sigma_t + pw; Dn = 0.0; }  // tension cut-off (separation)
    const double tau_tr = ks * (du_s - slip_p_c);
    const double tau_max = std::max(0.0, cohesion_for(p, tau_tr) - (sigma_n - pw) * std::tan(p.phi_i));
    if (std::fabs(tau_tr) <= tau_max)
        return {tau_tr, sigma_n, ks, Dn, slip_p_c};              // elastic (stick)
    const double tau = std::copysign(tau_max, tau_tr);
    return {tau, sigma_n, 0.0, Dn, du_s - tau / ks};             // slip (plastic)
}

}  // namespace katai::core::iface
