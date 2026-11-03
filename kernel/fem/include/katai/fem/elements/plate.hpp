#pragma once
// Plate (structural wall/beam) — 3-node quadratic Mindlin/Timoshenko beam element (P2.4).
// Shear-deformable (thin and thick plates both correct). Sits on a mesh edge (tri6 edge =
// 3 nodes). 3 DOFs per node: (u_x, u_y, φ) — the two translations are SHARED with the
// soil, the rotation φ is plate-specific.
//
// Constitutive law (per unit width; linear elastic in N, M and Q):
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
    // Elastoplastic capacities (the M-N diamond |N|/Np + |M|/Mp <= 1;
    // structural-plate-formulation.md §10): Mp = plastic moment [kNm/m], Np = plastic axial
    // force [kN/m]. ≤0 ⇒ UNBOUNDED in that direction; both ≤0 ⇒ purely elastic (old
    // behaviour bit-for-bit — plastic() false, the assembly stays on the K·u path).
    double Mp = -1.0;
    double Np = -1.0;
    double shear_rigidity() const { return k * EA / (2.0 * (1.0 + nu)); }  // kGA'
    bool plastic() const { return Mp > 0.0 || Np > 0.0; }
    // A plate whose section was changed in a phase carries the strains its forces no longer
    // account for in the same state an elastoplastic section keeps (see rescale_plate_state),
    // so it is evaluated on the path that reads that state even when it has no capacity.
    bool carries_state = false;
    bool stateful() const { return plastic() || carries_state; }
    // AXIAL SYMMETRY: the plate is a shell of revolution (a cone frustum per element, r = x the
    // radius, z = y the axis), per radian like the soil. Besides the meridional strain, curvature
    // and shear of the plane strain plate it has a hoop strain eps_t = u_r / r and a hoop
    // curvature kap_t = cos(beta) phi / r (beta the meridian's angle to the r axis), and
    //     N_s = EA (eps_s + nu eps_t),  N_t = EA (eps_t + nu eps_s),
    //     M_s = EI (kap_s + nu kap_t),  M_t = EI (kap_t + nu kap_s),
    // so that, with no hoop strain at all, it is the plane strain plate whose EA and EI were
    // entered. With a capacity the meridional N_s, M_s keep the M-N diamond and the hoop force its
    // own cap |N_t| <= Np (shell_return).
    bool axisym = false;
};

// Plate inertia for a Dynamic (seismic) phase, from the plate material weight w [kN/m/m]
// and gravity g [m/s^2] (Stage B6: engine-owned; the caller supplies g so this module
// does not reach upward for a constant):
//   rho_A = w/g            mass per unit length            [Mg/m]
//   rho_I = rho_A d^2/12   rotary inertia per unit length  [Mg m],  d = sqrt(12 EI/EA)
// d is the plate's equivalent thickness (the rectangular section with the
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

namespace detail {
// The hoop rows of a shell of revolution (PlateProps::axisym) at a point of radius r: the hoop
// strain u_r / r and the hoop curvature c phi / r, c = dr/ds (the cosine of the meridian's angle
// to the r axis; it changes sign with the element's direction exactly as dphi/ds does, so the
// energy does not).
template <int Nn>
inline void hoop_rows(const Eigen::Matrix<double, Nn, 1>& N, double c, double r,
                      Eigen::Matrix<double, 1, 3 * Nn>& Bt, Eigen::Matrix<double, 1, 3 * Nn>& Bkt) {
    Bt.setZero(); Bkt.setZero();
    for (int i = 0; i < Nn; ++i) {
        Bt(0, 3 * i + 0) = N(i) / r;
        Bkt(0, 3 * i + 2) = c * N(i) / r;
    }
}
}  // namespace detail

// Local→global stiffness (9×9). Selective reduced integration: axial+bending 3-point
// Gauss, shear 2-point Gauss (prevents shear locking). c,s are embedded in B → no separate
// transformation needed. In axial symmetry, per radian, with the hoop terms (PlateProps::axisym).
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
        if (p.axisym) {
            const double r = e.N.dot(X.col(0));
            Eigen::Matrix<double, 1, 9> Bt, Bkt;
            detail::hoop_rows<3>(e.N, e.c, r, Bt, Bkt);
            const double wds = w3[q] * e.J * r;
            K += wds * (p.EA * (Be.transpose() * Be + Bt.transpose() * Bt +
                                p.nu * (Be.transpose() * Bt + Bt.transpose() * Be)) +
                        p.EI * (Bk.transpose() * Bk + Bkt.transpose() * Bkt +
                                p.nu * (Bk.transpose() * Bkt + Bkt.transpose() * Bk)));
            continue;
        }
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
        const double wds = 1.0 * e.J * (p.axisym ? e.N.dot(X.col(0)) : 1.0);  // 2-point weight = 1
        K += wds * (kGA * Bg.transpose() * Bg);
    }
    return K;
}

// Structural force recovery at one ξ: N, Q, M (per unit width).
// In axial symmetry N and M are the meridional ones, and Nt, Mt the hoop force and moment (per
// unit length of the hoop and of the meridian); zero in plane strain.
struct PlateForces { double N, Q, M; double Nt = 0.0, Mt = 0.0; };
inline PlateForces forces(const NodeCoords& X, const PlateProps& p, const Dof& u, double xi) {
    const auto e = detail::edge_kin(X, xi);
    Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
    detail::b_rows(e, Be, Bk, Bg);
    if (p.axisym) {
        Eigen::Matrix<double, 1, 9> Bt, Bkt;
        detail::hoop_rows<3>(e.N, e.c, e.N.dot(X.col(0)), Bt, Bkt);
        const double es = (Be * u)(0), et = (Bt * u)(0), ks = (Bk * u)(0), kt = (Bkt * u)(0);
        return {p.EA * (es + p.nu * et), p.shear_rigidity() * (Bg * u)(0), p.EI * (ks + p.nu * kt),
                p.EA * (et + p.nu * es), p.EI * (kt + p.nu * ks)};
    }
    return {p.EA * (Be * u)(0), p.shear_rigidity() * (Bg * u)(0), p.EI * (Bk * u)(0)};
}

// ===========================================================================================
// ELASTOPLASTIC M-N HINGE (the linear M-N diamond; structural-plate-formulation.md §10).
// Generalized-stress perfect plasticity at the bending/axial Gauss (stress) points:
//   f = |N|/Np + |M|/Mp − 1 ≤ 0  (1/cap := 0 if unbounded), associated flow, Koiter corner
// return. Shear Q is NOT on the yield surface (only M and N are capped) — stays elastic.
// State = [ε_p, κ_p] per Gauss point; the return map is a PURE function of the committed
// state (line-search safe).
// ===========================================================================================
inline constexpr int kPlasticPerPoint = 2;                       // [ε_p, κ_p]
inline constexpr int kBendGaussCount = 3;                        // 3-node: axial+bending Gauss
inline constexpr int kBendGaussCount5 = 5;                       // 5-node
// After the axial/bending points the state holds one SHEAR STRAIN OFFSET per reduced shear
// point: the part of the shear strain the section does not answer for, Q = kGA (gamma - gamma_0).
// It is zero for every plate except one whose section was changed in a phase, which keeps its
// shear force through it as it keeps N and M through eps_p and kap_p. Nothing yields in shear,
// so the offset is only ever carried.
inline constexpr int kShearGaussCount = 2;                       // 3-node, reduced
inline constexpr int kShearGaussCount5 = 4;                      // 5-node
inline constexpr int kShearStateBegin = kBendGaussCount * kPlasticPerPoint;
inline constexpr int kShearStateBegin5 = kBendGaussCount5 * kPlasticPerPoint;
// Last, a shell of revolution (PlateProps::axisym) keeps the HOOP part of the same state per
// axial/bending point, [eps_t_p, kap_t_p]: the hoop strain and curvature its forces do not answer
// for. eps_t_p grows where the hoop force is at its cap; both are carried through a change of
// section. A plane strain plate has no hoop and leaves them zero.
inline constexpr int kHoopPerPoint = 2;                          // [eps_t_p, kap_t_p]
inline constexpr int kHoopStateBegin = kShearStateBegin + kShearGaussCount;
inline constexpr int kHoopStateBegin5 = kShearStateBegin5 + kShearGaussCount5;
inline constexpr int kPlasticStateSize = kHoopStateBegin + kBendGaussCount * kHoopPerPoint;     // 14 doubles/element
inline constexpr int kPlasticStateSize5 = kHoopStateBegin5 + kBendGaussCount5 * kHoopPerPoint;  // 24 doubles/element
// Bending/axial Gauss locations (the SAME literals as in stiffness) — the diagram Lagrange
// expansion (nodal values extrapolated from the stress points) uses these ξ.
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

// ELASTOPLASTIC SHELL OF REVOLUTION (PlateProps::axisym). The meridional force and moment keep
// the M-N diamond of the plane strain plate, |N_s|/Np + |M_s|/Mp <= 1; the hoop force has a cap of
// its own, |N_t| <= Np (an isotropic section has the same capacity in both directions), which
// takes no part in the diamond; the hoop moment is not capped. With the generalized stresses
// s = (N_s, N_t, M_s, M_t) and strains e = (eps_s, eps_t, kap_s, kap_t),
//     s = D (e - e_p),   D = EA [[1, nu], [nu, 1]]  (+)  EI [[1, nu], [nu, 1]],
// the elastic law of the shell (PlateProps::axisym). The admissible set is convex and bounded by
// at most six planes a_i . s <= 1 (four of the diamond, two of the hoop cap), and with associated
// flow the return is the closest admissible point in the energy norm of D: the one set of active
// planes A whose point s = s_tr - D A lambda, (A^T D A) lambda = A^T s_tr - 1, has lambda >= 0 and
// breaks no other plane. It is found by trying the sets of one, two and three planes in turn
// (a corner of the diamond is two, with the hoop cap three). The plastic strain grows by
// A lambda -- a component no active plane points along keeps its committed value bit for bit --
// and the consistent tangent is D - D A (A^T D A)^-1 A^T D. With no capacity (a section changed
// in a phase, which carries a state but cannot yield) the set is everything and s = D (e - e_p).
struct ShellReturn {
    Eigen::Vector4d s = Eigen::Vector4d::Zero();    // (N_s, N_t, M_s, M_t)
    Eigen::Vector4d ep = Eigen::Vector4d::Zero();   // new plastic (eps_s, eps_t, kap_s, kap_t)
    Eigen::Matrix4d D = Eigen::Matrix4d::Zero();    // consistent tangent
    bool yielded = false;
};
inline Eigen::Matrix4d shell_elasticity(const PlateProps& p) {
    Eigen::Matrix4d D = Eigen::Matrix4d::Zero();
    D(0, 0) = D(1, 1) = p.EA;
    D(0, 1) = D(1, 0) = p.nu * p.EA;
    D(2, 2) = D(3, 3) = p.EI;
    D(2, 3) = D(3, 2) = p.nu * p.EI;
    return D;
}
inline ShellReturn shell_return(const PlateProps& p, const Eigen::Vector4d& e,
                                const Eigen::Vector4d& ep_c) {
    ShellReturn r;
    r.D = shell_elasticity(p);
    r.ep = ep_c;
    const Eigen::Vector4d s_tr = r.D * (e - ep_c);
    r.s = s_tr;
    const double iN = p.Np > 0.0 ? 1.0 / p.Np : 0.0;   // 1/cap; 0 = unbounded
    const double iM = p.Mp > 0.0 ? 1.0 / p.Mp : 0.0;
    std::array<Eigen::Vector4d, 6> a;
    int na = 0;
    if (iN > 0.0 || iM > 0.0)
        for (double sN : {1.0, -1.0})
            for (double sM : {1.0, -1.0}) {
                if ((iN == 0.0 && sN < 0.0) || (iM == 0.0 && sM < 0.0)) continue;   // one cap: two planes
                a[(size_t)na++] = Eigen::Vector4d(sN * iN, 0.0, sM * iM, 0.0);
            }
    if (iN > 0.0) {
        a[(size_t)na++] = Eigen::Vector4d(0.0, iN, 0.0, 0.0);
        a[(size_t)na++] = Eigen::Vector4d(0.0, -iN, 0.0, 0.0);
    }
    double tol = 1.0e-12;   // a stress point on its surface stays where it is
    const auto admissible = [&](const Eigen::Vector4d& s) {
        for (int i = 0; i < na; ++i)
            if (a[(size_t)i].dot(s) > 1.0 + tol) return false;
        return true;
    };
    if (admissible(s_tr)) return r;
    r.yielded = true;
    // Candidate active sets, smallest first; each as a bit mask over the planes. Rounding can
    // leave the exact projection a hair outside a plane it does not touch: the second pass
    // accepts it within 1e-9 of the cap.
    for (; tol < 1.0e-8; tol *= 1.0e3)
    for (int size = 1; size <= 3; ++size)
        for (int mask = 1; mask < (1 << na); ++mask) {
            int act[3], k = 0;
            for (int i = 0; i < na && k <= size; ++i)
                if (mask & (1 << i)) { if (k < 3) act[k] = i; ++k; }
            if (k != size) continue;
            Eigen::Matrix<double, 4, Eigen::Dynamic> A(4, size);
            for (int j = 0; j < size; ++j) A.col(j) = a[(size_t)act[j]];
            const Eigen::MatrixXd DA = r.D * A;
            const Eigen::MatrixXd M = A.transpose() * DA;
            const Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
            if (lu.rank() < size) continue;   // dependent planes (opposite faces): not a corner
            const Eigen::VectorXd lam =
                lu.solve(A.transpose() * s_tr - Eigen::VectorXd::Ones(size));
            if (lam.minCoeff() < 0.0) continue;
            const Eigen::Vector4d s = s_tr - DA * lam;
            if (!admissible(s)) continue;
            r.s = s;
            r.ep = ep_c + A * lam;
            for (int c = 0; c < 4; ++c)
                if (A.row(c).isZero(0.0)) r.ep(c) = ep_c(c);   // untouched: bit for bit
            r.D -= DA * lu.inverse() * DA.transpose();
            return r;
        }
    return r;   // unreachable for a convex set of planes; the trial is kept as it is
}

// The generalized strains of a shell of revolution at one point, (eps_s, eps_t, kap_s, kap_t) =
// B u, and the committed plastic ones from the element state at axial/bending point q.
inline Eigen::Vector4d shell_plastic_state(const double* state, int hoop_begin, int q) {
    return Eigen::Vector4d(state[kPlasticPerPoint * q], state[hoop_begin + kHoopPerPoint * q],
                           state[kPlasticPerPoint * q + 1], state[hoop_begin + kHoopPerPoint * q + 1]);
}
inline void write_shell_plastic_state(double* state, int hoop_begin, int q, const Eigen::Vector4d& ep) {
    state[kPlasticPerPoint * q] = ep(0);
    state[kPlasticPerPoint * q + 1] = ep(2);
    state[hoop_begin + kHoopPerPoint * q] = ep(1);
    state[hoop_begin + kHoopPerPoint * q + 1] = ep(3);
}

// The shell of revolution with a state (a capacity, or a section changed in a phase) -- 3-node:
// the axial/bending part by shell_return at the 3-point rule, per radian (weight r), the shear
// elastic at the reduced 2-point rule less its carried offset, as internal_force_plastic.
inline void internal_force_shell(const NodeCoords& X, const PlateProps& p, const Dof& u,
                                 const double* state_c, double* state_t, Dof& f, ElementMatrix* K) {
    f.setZero();
    if (K) K->setZero();
    const std::array<double, 3> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    for (int q = 0; q < kBendGaussCount; ++q) {
        const auto e = detail::edge_kin(X, kBendGaussXi[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg, Bt, Bkt;
        detail::b_rows(e, Be, Bk, Bg);
        const double r = e.N.dot(X.col(0));
        detail::hoop_rows<3>(e.N, e.c, r, Bt, Bkt);
        Eigen::Matrix<double, 4, 9> B;
        B.row(0) = Be; B.row(1) = Bt; B.row(2) = Bk; B.row(3) = Bkt;
        const auto sr = shell_return(p, B * u, shell_plastic_state(state_c, kHoopStateBegin, q));
        write_shell_plastic_state(state_t, kHoopStateBegin, q, sr.ep);
        const double wds = w3[q] * e.J * r;
        f.noalias() += wds * (B.transpose() * sr.s);
        if (K) K->noalias() += wds * (B.transpose() * sr.D * B);
    }
    const double kGA = p.shear_rigidity();
    constexpr double g2 = 0.5773502691896257;
    const std::array<double, 2> xi2{-g2, g2};
    for (int q = 0; q < 2; ++q) {
        const auto e = detail::edge_kin(X, xi2[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
        detail::b_rows(e, Be, Bk, Bg);
        const double wds = 1.0 * e.J * e.N.dot(X.col(0));   // 2-point weight = 1, per radian
        state_t[kShearStateBegin + q] = state_c[kShearStateBegin + q];
        f.noalias() += wds * kGA * Bg.transpose() * ((Bg * u)(0) - state_c[kShearStateBegin + q]);
        if (K) K->noalias() += wds * kGA * Bg.transpose() * Bg;
    }
}

// The four generalized stresses of a shell of revolution at its axial/bending points -- capped by
// the return, or (capped = false) the elastic ones from the state: what a stiffness on top of it
// carries.
inline void gauss_forces_shell(const NodeCoords& X, const PlateProps& p, const Dof& u,
                               const double* state_c, bool capped,
                               std::array<Eigen::Vector4d, 3>& sg) {
    for (int q = 0; q < kBendGaussCount; ++q) {
        const auto e = detail::edge_kin(X, kBendGaussXi[q]);
        Eigen::Matrix<double, 1, 9> Be, Bk, Bg, Bt, Bkt;
        detail::b_rows(e, Be, Bk, Bg);
        detail::hoop_rows<3>(e.N, e.c, e.N.dot(X.col(0)), Bt, Bkt);
        const Eigen::Vector4d strain((Be * u)(0), (Bt * u)(0), (Bk * u)(0), (Bkt * u)(0));
        const Eigen::Vector4d ep = shell_plastic_state(state_c, kHoopStateBegin, q);
        sg[(size_t)q] = capped ? shell_return(p, strain, ep).s : Eigen::Vector4d(shell_elasticity(p) * (strain - ep));
    }
}

// Elastoplastic element internal force + (if K != null) consistent tangent — 3-node.
// Axial+bending return mapping at the 3-point Gauss rule (state_c → state_t, [ε_p,κ_p]×3),
// shear elastic at the REDUCED 2-point rule (same scheme as stiffness). With
// props.plastic() false the caller uses the OLD K·u path.
inline void internal_force_plastic(const NodeCoords& X, const PlateProps& p, const Dof& u,
                                   const double* state_c, double* state_t,
                                   Dof& f, ElementMatrix* K) {
    if (p.axisym) { internal_force_shell(X, p, u, state_c, state_t, f, K); return; }
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
        state_t[kShearStateBegin + q] = state_c[kShearStateBegin + q];
        f.noalias() += wds * kGA * Bg.transpose() * ((Bg * u)(0) - state_c[kShearStateBegin + q]);
        if (K) K->noalias() += wds * kGA * Bg.transpose() * Bg;
    }
}

// For the diagram: the capped (N, M) at the bending/axial Gauss points with the committed
// state — the stations expand from these by Lagrange (nodal output is stress-point
// extrapolation; the cap CAN be exceeded at a station and is not checked there — §10.2
// honest declaration).
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
// quarter-spaced). The 5-node plate element that pairs with the
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
        if (p.axisym) {
            const double r = e.N.dot(X.col(0));
            Eigen::Matrix<double, 1, 15> Bt, Bkt;
            detail::hoop_rows<5>(e.N, e.c, r, Bt, Bkt);
            const double wds = w5[q] * e.J * r;
            K += wds * (p.EA * (Be.transpose() * Be + Bt.transpose() * Bt +
                                p.nu * (Be.transpose() * Bt + Bt.transpose() * Be)) +
                        p.EI * (Bk.transpose() * Bk + Bkt.transpose() * Bkt +
                                p.nu * (Bk.transpose() * Bkt + Bkt.transpose() * Bk)));
            continue;
        }
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
        K += (w4[q] * e.J * (p.axisym ? e.N.dot(X.col(0)) : 1.0)) * (kGA * Bg.transpose() * Bg);
    }
    return K;
}

// 5-node beam force recovery at one ξ: N, Q, M (per unit width).
inline PlateForces forces5(const NodeCoords5& X, const PlateProps& p, const Dof5& u, double xi) {
    const auto e = detail::edge_kin5(X, xi);
    Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
    detail::b_rows5(e, Be, Bk, Bg);
    if (p.axisym) {
        Eigen::Matrix<double, 1, 15> Bt, Bkt;
        detail::hoop_rows<5>(e.N, e.c, e.N.dot(X.col(0)), Bt, Bkt);
        const double es = (Be * u)(0), et = (Bt * u)(0), ks = (Bk * u)(0), kt = (Bkt * u)(0);
        return {p.EA * (es + p.nu * et), p.shear_rigidity() * (Bg * u)(0), p.EI * (ks + p.nu * kt),
                p.EA * (et + p.nu * es), p.EI * (kt + p.nu * ks)};
    }
    return {p.EA * (Be * u)(0), p.shear_rigidity() * (Bg * u)(0), p.EI * (Bk * u)(0)};
}

// Elastoplastic element internal force + tangent — 5-node (same logic as 3-node;
// axial+bending 5-point Gauss return mapping [ε_p,κ_p]×5, shear elastic at the REDUCED
// 4-point rule).
// The shell of revolution with a state -- 5-node (internal_force_shell).
inline void internal_force_shell5(const NodeCoords5& X, const PlateProps& p, const Dof5& u,
                                  const double* state_c, double* state_t, Dof5& f, ElementMatrix5* K) {
    f.setZero();
    if (K) K->setZero();
    const std::array<double, 5> w5{0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                   0.4786286704993665, 0.2369268850561891};
    for (int q = 0; q < kBendGaussCount5; ++q) {
        const auto e = detail::edge_kin5(X, kBendGaussXi5[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg, Bt, Bkt;
        detail::b_rows5(e, Be, Bk, Bg);
        const double r = e.N.dot(X.col(0));
        detail::hoop_rows<5>(e.N, e.c, r, Bt, Bkt);
        Eigen::Matrix<double, 4, 15> B;
        B.row(0) = Be; B.row(1) = Bt; B.row(2) = Bk; B.row(3) = Bkt;
        const auto sr = shell_return(p, B * u, shell_plastic_state(state_c, kHoopStateBegin5, q));
        write_shell_plastic_state(state_t, kHoopStateBegin5, q, sr.ep);
        const double wds = w5[q] * e.J * r;
        f.noalias() += wds * (B.transpose() * sr.s);
        if (K) K->noalias() += wds * (B.transpose() * sr.D * B);
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
        const double wds = w4[q] * e.J * e.N.dot(X.col(0));   // per radian
        state_t[kShearStateBegin5 + q] = state_c[kShearStateBegin5 + q];
        f.noalias() += wds * kGA * Bg.transpose() * ((Bg * u)(0) - state_c[kShearStateBegin5 + q]);
        if (K) K->noalias() += wds * kGA * Bg.transpose() * Bg;
    }
}

inline void gauss_forces_shell5(const NodeCoords5& X, const PlateProps& p, const Dof5& u,
                                const double* state_c, bool capped,
                                std::array<Eigen::Vector4d, 5>& sg) {
    for (int q = 0; q < kBendGaussCount5; ++q) {
        const auto e = detail::edge_kin5(X, kBendGaussXi5[q]);
        Eigen::Matrix<double, 1, 15> Be, Bk, Bg, Bt, Bkt;
        detail::b_rows5(e, Be, Bk, Bg);
        detail::hoop_rows<5>(e.N, e.c, e.N.dot(X.col(0)), Bt, Bkt);
        const Eigen::Vector4d strain((Be * u)(0), (Bt * u)(0), (Bk * u)(0), (Bkt * u)(0));
        const Eigen::Vector4d ep = shell_plastic_state(state_c, kHoopStateBegin5, q);
        sg[(size_t)q] = capped ? shell_return(p, strain, ep).s : Eigen::Vector4d(shell_elasticity(p) * (strain - ep));
    }
}

inline void internal_force_plastic5(const NodeCoords5& X, const PlateProps& p, const Dof5& u,
                                    const double* state_c, double* state_t,
                                    Dof5& f, ElementMatrix5* K) {
    if (p.axisym) { internal_force_shell5(X, p, u, state_c, state_t, f, K); return; }
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
        state_t[kShearStateBegin5 + q] = state_c[kShearStateBegin5 + q];
        f.noalias() += wds * kGA * Bg.transpose() * ((Bg * u)(0) - state_c[kShearStateBegin5 + q]);
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
