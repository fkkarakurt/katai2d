#pragma once
// Plate (structural wall/beam) — 3-node quadratic Mindlin/Timoshenko beam element (P2.4).
// Shear-deformable (thin and thick plates both correct). Sits on a mesh edge (tri6 edge =
// 3 nodes). 3 DOFs per node: (u_x, u_y, φ) — the two translations are SHARED with the
// soil, the rotation φ is plate-specific.
//
// Constitutive law (PLAXIS 2D Material Models Manual §18.3, Eq 18-6…18-9):
//   N = EA·ε ,   M = EI·κ ,   Q = kGA'·γ ,   kGA' = k·EA/(2(1+ν)) ,  k = 5/6
// Kinematics (local: s = plate axis, n = normal):
//   ε = du_s/ds ,   κ = dφ/ds ,   γ = du_n/ds − φ   (Mindlin: φ independent of the transverse slope)
// SHEAR LOCKING: REDUCED integration on the shear term (bending/axial 3-point, shear
// 2-point selective reduced — Zienkiewicz&Taylor) → no spurious stiffness in thin plates.
//
// Math: docs/references/structural-plate-formulation.md.

#include <array>
#include <cmath>

#include <Eigen/Dense>

namespace katai::core::plate {

inline constexpr int kNodeCount = 3;
inline constexpr int kDofCount = 9;  // 3 nodes × (u_x, u_y, φ)

using NodeCoords = Eigen::Matrix<double, 3, 2>;     // row i: (x_i, y_i)
using ElementMatrix = Eigen::Matrix<double, 9, 9>;  // K
using Dof = Eigen::Matrix<double, 9, 1>;

// Node order: ξ=−1 (end A), ξ=+1 (end B), ξ=0 (middle). Quadratic 1D shape functions.
inline Eigen::Vector3d shape(double xi) {
    return {0.5 * xi * (xi - 1.0), 0.5 * xi * (xi + 1.0), 1.0 - xi * xi};
}
inline Eigen::Vector3d shape_deriv_xi(double xi) {
    return {xi - 0.5, xi + 0.5, -2.0 * xi};
}

// Plate physical parameters. kGA' = shear_factor·EA, shear_factor = k/(2(1+ν)), k=5/6.
// The MASS fields (rho_A, rho_I) are used ONLY in dynamic analysis
// (analysis/structural_dynamics.hpp assemble_structural_mass); the static path never reads
// them → the 0 default leaves statics bit-for-bit.
// WHEN ADDING A FIELD: append at the END — inserting in the middle silently shifts
// positional aggregate-inits.
struct PlateProps {
    double EA = 0.0;   // axial stiffness [kN/m]
    double EI = 0.0;   // bending stiffness [kNm²/m]
    double nu = 0.0;   // plate Poisson
    double k = 5.0 / 6.0;  // shear correction factor
    // Mass per unit length ρA = w/g [Mg/m] (w = plate weight [kN/m/m], g = 9.81 m/s²) and
    // rotary inertia per unit length ρI = ρA·d²/12 [Mg·m] (d = equivalent thickness
    // √(12EI/EA), rectangular section). 0 = massless (no plate inertia in dynamics).
    double rho_A = 0.0;
    double rho_I = 0.0;
    // Elastoplastic capacities (the PLAXIS MMM §18.3 diamond;
    // structural-plate-formulation.md §10): Mp = plastic moment [kNm/m], Np = plastic axial
    // force [kN/m]. ≤0 ⇒ UNBOUNDED in that direction; both ≤0 ⇒ purely elastic (old
    // behaviour bit-for-bit — plastic() false, the assembly stays on the K·u path).
    double Mp = -1.0;
    double Np = -1.0;
    double shear_rigidity() const { return k * EA / (2.0 * (1.0 + nu)); }  // kGA'
    bool plastic() const { return Mp > 0.0 || Np > 0.0; }
};

// Plate inertia for a Dynamic (seismic) phase, from the plate material weight w [kN/m/m]
// and gravity g [m/s^2] (Stage B6: engine-owned; the caller supplies g so this module
// does not reach upward for a constant):
//   rho_A = w/g            mass per unit length            [Mg/m]
//   rho_I = rho_A d^2/12   rotary inertia per unit length  [Mg m],  d = sqrt(12 EI/EA)
// d is PLAXIS's equivalent thickness (Ref. Man. sec. 5.6: a rectangular section with the
// same EA and EI), so rho_I is that section's second moment of mass. Static paths never
// read these fields; a plate with w = 0 stays massless (its stiffness still takes part
// in the dynamic system). Call AFTER EA/EI are set.
inline void set_plate_mass(PlateProps& pp, double w, double g) {
    pp.rho_A = std::fmax(0.0, w) / g;
    const double d = (pp.EA > 0.0 && pp.EI > 0.0) ? std::sqrt(12.0 * pp.EI / pp.EA) : 0.0;
    pp.rho_I = pp.rho_A * d * d / 12.0;
}

namespace detail {
// At a Gauss point: ds/dξ (Jacobian), unit tangent (c,s), dN/ds.
struct EdgeKinematics {
    double J;                 // ds/dξ
    double c, s;              // unit tangent (cosθ, sinθ)
    Eigen::Vector3d dNds;     // dN_i/ds
    Eigen::Vector3d N;        // N_i
};
inline EdgeKinematics edge_kin(const NodeCoords& X, double xi) {
    const Eigen::Vector3d N = shape(xi), dNdxi = shape_deriv_xi(xi);
    const double dx = dNdxi.dot(X.col(0)), dy = dNdxi.dot(X.col(1));
    const double J = std::sqrt(dx * dx + dy * dy);
    EdgeKinematics e;
    e.J = J; e.c = dx / J; e.s = dy / J; e.N = N; e.dNds = dNdxi / J;
    return e;
}
// Axial/bending/shear B rows (1×9), global DOFs [u_x,u_y,φ] in node order.
inline void b_rows(const EdgeKinematics& e, Eigen::Matrix<double, 1, 9>& Be,
                   Eigen::Matrix<double, 1, 9>& Bk, Eigen::Matrix<double, 1, 9>& Bg) {
    Be.setZero(); Bk.setZero(); Bg.setZero();
    for (int i = 0; i < 3; ++i) {
        const int b = 3 * i;
        // ε = du_s/ds = Σ dN_i/ds (c·u_x + s·u_y)
        Be(0, b + 0) = e.c * e.dNds(i);
        Be(0, b + 1) = e.s * e.dNds(i);
        // κ = dφ/ds
        Bk(0, b + 2) = e.dNds(i);
        // γ = du_n/ds − φ = Σ dN_i/ds (−s·u_x + c·u_y) − Σ N_i φ
        Bg(0, b + 0) = -e.s * e.dNds(i);
        Bg(0, b + 1) = e.c * e.dNds(i);
        Bg(0, b + 2) = -e.N(i);
    }
}
}  // namespace detail

// Local→global stiffness (9×9). Selective reduced integration: axial+bending 3-point
// Gauss, shear 2-point Gauss (prevents shear locking). c,s are embedded in B → no separate
// transformation needed.
inline ElementMatrix stiffness(const NodeCoords& X, const PlateProps& p) {
    ElementMatrix K = ElementMatrix::Zero();
    const double kGA = p.shear_rigidity();
    // 3-point Gauss (axial + bending).
    constexpr double g3 = 0.7745966692414834;  // sqrt(3/5)
    const std::array<double, 3> xi3{-g3, 0.0, g3};
    const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (int q = 0; q < 3; ++q) {
        const auto e = detail::edge_kin(X, xi3[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const double wds = w3[q] * e.J;
        K += wds * (p.EA * Be.transpose() * Be + p.EI * Bk.transpose() * Bk);
    }
    // 2-point Gauss (shear, REDUCED).
    constexpr double g2 = 0.5773502691896257;  // 1/sqrt(3)
    const std::array<double, 2> xi2{-g2, g2};
    for (int q = 0; q < 2; ++q) {
        const auto e = detail::edge_kin(X, xi2[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const double wds = 1.0 * e.J;  // 2-point weight = 1
        K += wds * (kGA * Bg.transpose() * Bg);
    }
    return K;
}

// Structural force recovery at one ξ: N, Q, M (per unit width).
struct PlateForces { double N, Q, M; };
inline PlateForces forces(const NodeCoords& X, const PlateProps& p, const Dof& u, double xi) {
    const auto e = detail::edge_kin(X, xi);
    Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
    detail::b_rows(e, Be, Bk, Bg);
    return {p.EA * (Be * u)(0), p.shear_rigidity() * (Bg * u)(0), p.EI * (Bk * u)(0)};
}

// ===========================================================================================
// ELASTOPLASTIC M-N HINGE (the PLAXIS MMM §18.3 diamond; structural-plate-formulation.md §10).
// Generalized-stress perfect plasticity at the bending/axial Gauss (stress) points:
//   f = |N|/Np + |M|/Mp − 1 ≤ 0  (1/cap := 0 if unbounded), associated flow, Koiter corner
// return. Shear Q is NOT on the yield surface (PLAXIS: M and N only) — stays elastic.
// State = [ε_p, κ_p] per Gauss point; the return map is a PURE function of the committed
// state (line-search safe).
// ===========================================================================================
inline constexpr int kPlasticPerPoint = 2;                       // [ε_p, κ_p]
inline constexpr int kBendGaussCount = 3;                        // 3-node: axial+bending Gauss
inline constexpr int kBendGaussCount5 = 5;                       // 5-node
inline constexpr int kPlasticStateSize = kBendGaussCount * kPlasticPerPoint;    // 6 doubles/element
inline constexpr int kPlasticStateSize5 = kBendGaussCount5 * kPlasticPerPoint;  // 10 doubles/element
// Bending/axial Gauss locations (the SAME literals as in stiffness) — the diagram Lagrange
// expansion (PLAXIS: "extrapolation of the values at the stress points") uses these ξ.
inline constexpr std::array<double, 3> kBendGaussXi{-0.7745966692414834, 0.0, 0.7745966692414834};
inline constexpr std::array<double, 5> kBendGaussXi5{
    -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};

// The M-N return map at one stress point. Input: TOTAL (ε, κ) + committed plastic
// (ε_p, κ_p). Output: capped (N, M), the NEW plastic state and the consistent tangent
// (linear surface + perfect plasticity ⇒ consistent ≡ continuum). Validity-cascaded
// selection: quadrant surface → V_M corner → V_N corner.
struct MnReturn {
    double N = 0.0, M = 0.0;                  // return-mapped generalized stresses
    double ep = 0.0, kp = 0.0;                // new plastic state (written to the trial)
    double Dnn = 0.0, Dnm = 0.0, Dmm = 0.0;   // consistent tangent (symmetric 2×2)
    bool yielded = false;
};
inline MnReturn mn_return(const PlateProps& p, double eps, double kap, double ep_c, double kp_c) {
    MnReturn r;
    r.ep = ep_c; r.kp = kp_c;
    const double Ntr = p.EA * (eps - ep_c), Mtr = p.EI * (kap - kp_c);
    r.N = Ntr; r.M = Mtr; r.Dnn = p.EA; r.Dmm = p.EI;
    const double iN = p.Np > 0.0 ? 1.0 / p.Np : 0.0;   // 1/cap; 0 = unbounded
    const double iM = p.Mp > 0.0 ? 1.0 / p.Mp : 0.0;
    const double f = std::fabs(Ntr) * iN + std::fabs(Mtr) * iM - 1.0;
    if (f <= 0.0) return r;                            // elastic (f = −1 when uncapped)
    r.yielded = true;
    const double sN = Ntr >= 0.0 ? 1.0 : -1.0, sM = Mtr >= 0.0 ? 1.0 : -1.0;
    // (1) The quadrant's surface: s = s_tr − Δλ·D·a, a = (sN·iN, sM·iM), Δλ = f/h.
    const double h = p.EA * iN * iN + p.EI * iM * iM;
    const double Nf = Ntr - (f / h) * p.EA * sN * iN;
    const double Mf = Mtr - (f / h) * p.EI * sM * iM;
    if (sN * Nf >= 0.0 && sM * Mf >= 0.0) {
        r.N = Nf; r.M = Mf;
        // Plastic state from the return point; in the uncapped (unchanged) component the
        // committed value is preserved BIT-FOR-BIT — reading it back from the formula would
        // write ±ulp drift into ε_p on every evaluation (the same lesson as the interface
        // resting on the cap; it breaks the nil-phase identity).
        r.ep = iN > 0.0 ? eps - Nf / p.EA : ep_c;
        r.kp = iM > 0.0 ? kap - Mf / p.EI : kp_c;
        // D_ep = D − (D·a)(aᵀD)/h  (rank-1; with one cap the other component stays elastic).
        const double dn = p.EA * sN * iN, dm = p.EI * sM * iM;
        r.Dnn = p.EA - dn * dn / h;
        r.Dnm = -dn * dm / h;
        r.Dmm = p.EI - dm * dm / h;
    } else if ((std::fabs(Mtr) - p.Mp) * p.Mp / p.EI >= std::fabs(Ntr) * p.Np / p.EA) {
        // (2) The V_M = (0, sM·Mp) corner (reachable only when both caps are finite: with
        // one cap the surface return is always valid). Normal-cone derivation §10.1. Both
        // components fixed → tangent 0.
        r.N = 0.0; r.M = sM * p.Mp;
        r.ep = eps;                          // N = 0 → ε_p = ε (exact)
        r.kp = kap - sM * p.Mp / p.EI;
        r.Dnn = r.Dnm = r.Dmm = 0.0;
    } else {
        // (3) The V_N = (sN·Np, 0) corner.
        r.N = sN * p.Np; r.M = 0.0;
        r.ep = eps - sN * p.Np / p.EA;
        r.kp = kap;                          // M = 0 → κ_p = κ (exact)
        r.Dnn = r.Dnm = r.Dmm = 0.0;
    }
    return r;
}

// Elastoplastic element internal force + (if K != null) consistent tangent — 3-node.
// Axial+bending return mapping at the 3-point Gauss rule (state_c → state_t, [ε_p,κ_p]×3),
// shear elastic at the REDUCED 2-point rule (same scheme as stiffness). With
// props.plastic() false the caller uses the OLD K·u path.
inline void internal_force_plastic(const NodeCoords& X, const PlateProps& p, const Dof& u,
                                   const double* state_c, double* state_t,
                                   Dof& f, ElementMatrix* K) {
    f.setZero();
    if (K) K->setZero();
    const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (int q = 0; q < kBendGaussCount; ++q) {
        const auto e = detail::edge_kin(X, kBendGaussXi[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const auto r = mn_return(p, (Be * u)(0), (Bk * u)(0),
                                 state_c[kPlasticPerPoint * q], state_c[kPlasticPerPoint * q + 1]);
        state_t[kPlasticPerPoint * q] = r.ep;
        state_t[kPlasticPerPoint * q + 1] = r.kp;
        const double wds = w3[q] * e.J;
        f.noalias() += wds * (Be.transpose() * r.N + Bk.transpose() * r.M);
        if (K)
            K->noalias() += wds * (r.Dnn * Be.transpose() * Be + r.Dmm * Bk.transpose() * Bk +
                                   r.Dnm * (Be.transpose() * Bk + Bk.transpose() * Be));
    }
    const double kGA = p.shear_rigidity();
    constexpr double g2 = 0.5773502691896257;
    const std::array<double, 2> xi2{-g2, g2};
    for (int q = 0; q < 2; ++q) {
        const auto e = detail::edge_kin(X, xi2[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const double wds = 1.0 * e.J;   // 2-point weight = 1
        f.noalias() += wds * kGA * Bg.transpose() * (Bg * u);
        if (K) K->noalias() += wds * kGA * Bg.transpose() * Bg;
    }
}

// For the diagram: the capped (N, M) at the bending/axial Gauss points with the committed
// state — the stations expand from these by Lagrange (PLAXIS: nodal output is
// stress-point extrapolation; the cap CAN be exceeded at a station, PLAXIS does not check
// either — §10.2 honest declaration).
inline void gauss_forces_plastic(const NodeCoords& X, const PlateProps& p, const Dof& u,
                                 const double* state_c, std::array<double, 3>& Ng,
                                 std::array<double, 3>& Mg) {
    for (int q = 0; q < kBendGaussCount; ++q) {
        const auto e = detail::edge_kin(X, kBendGaussXi[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const auto r = mn_return(p, (Be * u)(0), (Bk * u)(0),
                                 state_c[kPlasticPerPoint * q], state_c[kPlasticPerPoint * q + 1]);
        Ng[q] = r.N; Mg[q] = r.M;
    }
}

// ===========================================================================================
// 5-NODE QUARTIC Timoshenko beam — sits on a tri15 edge (tri15 edge = 5 nodes,
// quarter-spaced). The counterpart of the 5-node plate element PLAXIS 2D uses with the
// 15-node soil. NATURAL node order: ξ = −1, −0.5, 0, +0.5, +1 (nodes 0..4). (u_x,u_y,φ)
// per node → 15 DOFs. Constitutive law + kinematics as the 3-node one; SELECTIVE REDUCED
// INTEGRATION: axial+bending 5-point Gauss, shear 4-point Gauss (one degree lower → no
// shear locking). Math: structural-plate-formulation.md.
// ===========================================================================================
inline constexpr int kNodeCount5 = 5;
inline constexpr int kDofCount5 = 15;
using NodeCoords5 = Eigen::Matrix<double, 5, 2>;
using ElementMatrix5 = Eigen::Matrix<double, 15, 15>;
using Dof5 = Eigen::Matrix<double, 15, 1>;

namespace detail {
inline constexpr std::array<double, 5> kXi5{-1.0, -0.5, 0.0, 0.5, 1.0};  // node locations

// Quartic Lagrange shape functions and derivatives (the general product formula).
inline Eigen::Matrix<double, 5, 1> shape5(double xi) {
    Eigen::Matrix<double, 5, 1> N;
    for (int i = 0; i < 5; ++i) {
        double v = 1.0;
        for (int j = 0; j < 5; ++j)
            if (j != i) v *= (xi - kXi5[j]) / (kXi5[i] - kXi5[j]);
        N(i) = v;
    }
    return N;
}
inline Eigen::Matrix<double, 5, 1> shape5_deriv(double xi) {
    Eigen::Matrix<double, 5, 1> dN;
    for (int i = 0; i < 5; ++i) {
        double s = 0.0;
        for (int k = 0; k < 5; ++k)
            if (k != i) {
                double p = 1.0 / (kXi5[i] - kXi5[k]);
                for (int j = 0; j < 5; ++j)
                    if (j != i && j != k) p *= (xi - kXi5[j]) / (kXi5[i] - kXi5[j]);
                s += p;
            }
        dN(i) = s;
    }
    return dN;
}

struct EdgeKinematics5 {
    double J, c, s;
    Eigen::Matrix<double, 5, 1> dNds, N;
};
inline EdgeKinematics5 edge_kin5(const NodeCoords5& X, double xi) {
    const Eigen::Matrix<double, 5, 1> N = shape5(xi), dNdxi = shape5_deriv(xi);
    const double dx = dNdxi.dot(X.col(0)), dy = dNdxi.dot(X.col(1));
    const double J = std::sqrt(dx * dx + dy * dy);
    EdgeKinematics5 e;
    e.J = J; e.c = dx / J; e.s = dy / J; e.N = N; e.dNds = dNdxi / J;
    return e;
}
inline void b_rows5(const EdgeKinematics5& e, Eigen::Matrix<double, 1, 15>& Be,
                    Eigen::Matrix<double, 1, 15>& Bk, Eigen::Matrix<double, 1, 15>& Bg) {
    Be.setZero(); Bk.setZero(); Bg.setZero();
    for (int i = 0; i < 5; ++i) {
        const int b = 3 * i;
        Be(0, b + 0) = e.c * e.dNds(i);          // ε = du_s/ds
        Be(0, b + 1) = e.s * e.dNds(i);
        Bk(0, b + 2) = e.dNds(i);                // κ = dφ/ds
        Bg(0, b + 0) = -e.s * e.dNds(i);         // γ = du_n/ds − φ
        Bg(0, b + 1) = e.c * e.dNds(i);
        Bg(0, b + 2) = -e.N(i);
    }
}
}  // namespace detail

// 15×15 local→global stiffness. 5-point Gauss (axial+bending), 4-point Gauss (shear, reduced).
inline ElementMatrix5 stiffness5(const NodeCoords5& X, const PlateProps& p) {
    ElementMatrix5 K = ElementMatrix5::Zero();
    const double kGA = p.shear_rigidity();
    // 5-point Gauss (degree-9 exact): axial + bending.
    constexpr double a5 = 0.5384693101056831, b5 = 0.9061798459386640;
    const std::array<double, 5> xi5{-b5, -a5, 0.0, a5, b5};
    const std::array<double, 5> w5{0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                   0.4786286704993665, 0.2369268850561891};
    for (int q = 0; q < 5; ++q) {
        const auto e = detail::edge_kin5(X, xi5[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
        detail::b_rows5(e, Be, Bk, Bg);
        const double wds = w5[q] * e.J;
        K += wds * (p.EA * Be.transpose() * Be + p.EI * Bk.transpose() * Bk);
    }
    // 4-point Gauss (shear, REDUCED → no shear locking).
    constexpr double a4 = 0.3399810435848563, b4 = 0.8611363115940526;
    const std::array<double, 4> xi4{-b4, -a4, a4, b4};
    const std::array<double, 4> w4{0.3478548451374538, 0.6521451548625461,
                                   0.6521451548625461, 0.3478548451374538};
    for (int q = 0; q < 4; ++q) {
        const auto e = detail::edge_kin5(X, xi4[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
        detail::b_rows5(e, Be, Bk, Bg);
        K += (w4[q] * e.J) * (kGA * Bg.transpose() * Bg);
    }
    return K;
}

// 5-node beam force recovery at one ξ: N, Q, M (per unit width).
inline PlateForces forces5(const NodeCoords5& X, const PlateProps& p, const Dof5& u, double xi) {
    const auto e = detail::edge_kin5(X, xi);
    Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
    detail::b_rows5(e, Be, Bk, Bg);
    return {p.EA * (Be * u)(0), p.shear_rigidity() * (Bg * u)(0), p.EI * (Bk * u)(0)};
}

// Elastoplastic element internal force + tangent — 5-node (same logic as 3-node;
// axial+bending 5-point Gauss return mapping [ε_p,κ_p]×5, shear elastic at the REDUCED
// 4-point rule).
inline void internal_force_plastic5(const NodeCoords5& X, const PlateProps& p, const Dof5& u,
                                    const double* state_c, double* state_t,
                                    Dof5& f, ElementMatrix5* K) {
    f.setZero();
    if (K) K->setZero();
    const std::array<double, 5> w5{0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                   0.4786286704993665, 0.2369268850561891};
    for (int q = 0; q < kBendGaussCount5; ++q) {
        const auto e = detail::edge_kin5(X, kBendGaussXi5[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
        detail::b_rows5(e, Be, Bk, Bg);
        const auto r = mn_return(p, (Be * u)(0), (Bk * u)(0),
                                 state_c[kPlasticPerPoint * q], state_c[kPlasticPerPoint * q + 1]);
        state_t[kPlasticPerPoint * q] = r.ep;
        state_t[kPlasticPerPoint * q + 1] = r.kp;
        const double wds = w5[q] * e.J;
        f.noalias() += wds * (Be.transpose() * r.N + Bk.transpose() * r.M);
        if (K)
            K->noalias() += wds * (r.Dnn * Be.transpose() * Be + r.Dmm * Bk.transpose() * Bk +
                                   r.Dnm * (Be.transpose() * Bk + Bk.transpose() * Be));
    }
    const double kGA = p.shear_rigidity();
    constexpr double a4 = 0.3399810435848563, b4 = 0.8611363115940526;
    const std::array<double, 4> xi4{-b4, -a4, a4, b4};
    const std::array<double, 4> w4{0.3478548451374538, 0.6521451548625461,
                                   0.6521451548625461, 0.3478548451374538};
    for (int q = 0; q < 4; ++q) {
        const auto e = detail::edge_kin5(X, xi4[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
        detail::b_rows5(e, Be, Bk, Bg);
        const double wds = w4[q] * e.J;
        f.noalias() += wds * kGA * Bg.transpose() * (Bg * u);
        if (K) K->noalias() += wds * kGA * Bg.transpose() * Bg;
    }
}

// 5-node Gauss (N, M) samples for the diagram (with the committed state; Lagrange expansion at the caller).
inline void gauss_forces_plastic5(const NodeCoords5& X, const PlateProps& p, const Dof5& u,
                                  const double* state_c, std::array<double, 5>& Ng,
                                  std::array<double, 5>& Mg) {
    for (int q = 0; q < kBendGaussCount5; ++q) {
        const auto e = detail::edge_kin5(X, kBendGaussXi5[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
        detail::b_rows5(e, Be, Bk, Bg);
        const auto r = mn_return(p, (Be * u)(0), (Bk * u)(0),
                                 state_c[kPlasticPerPoint * q], state_c[kPlasticPerPoint * q + 1]);
        Ng[q] = r.N; Mg[q] = r.M;
    }
}

}  // namespace katai::core::plate
