#pragma once
// Geogrid (reinforcement/geosynthetic) — 3-node axial membrane element (P2.4). The
// bending-free (EI=0, no shear, no φ-DOF) form of the plate: only the axial force
// N = EA·ε. Sits on a mesh edge (tri6 edge = 3 nodes), 2 DOFs per node (u_x, u_y) —
// SHARED with the soil.
//
// Constitutive law (PLAXIS 2D Material Models Manual §18.2, Eq 18-2):
//   - TENSION-ONLY: tension only; slackens in compression (N=0, REVERSIBLE slack —
//     nonlinear-elastic).
//   - N_p (optional): max tensile force; yields at N=N_p (PERMANENT plastic elongation ε_p).
//
// Math: docs/references/structural-plate-formulation.md §8.

#include <array>
#include <cmath>

#include <Eigen/Dense>

namespace katai::core::geogrid {

inline constexpr int kNodeCount = 3;
inline constexpr int kDofCount = 6;  // 3 nodes × (u_x, u_y)

using NodeCoords = Eigen::Matrix<double, 3, 2>;     // row i: (x_i, y_i)
using ElementMatrix = Eigen::Matrix<double, 6, 6>;
using Dof = Eigen::Matrix<double, 6, 1>;
using AxialB = Eigen::Matrix<double, 1, 6>;

struct GeogridProps {
    double EA = 0.0;    // axial stiffness [kN/m]
    double Np = -1.0;   // max tensile force [kN/m]; ≤0 ⇒ unbounded (pure tension-only)
};

// Quadratic 1D shape functions (same as plate: ξ=−1 end A, ξ=+1 end B, ξ=0 mid).
inline Eigen::Vector3d shape_deriv_xi(double xi) {
    return {xi - 0.5, xi + 0.5, -2.0 * xi};
}

// At a Gauss point: ds/dξ (Jacobian), unit tangent (c,s), the axial B row (1×6, global
// DOFs). ε = du_s/ds = Σ dN_i/ds (c·u_x + s·u_y).
struct AxialKin {
    double J;          // ds/dξ
    AxialB Be;         // axial B (c,s embedded)
};
inline AxialKin axial_kin(const NodeCoords& X, double xi) {
    const Eigen::Vector3d dNdxi = shape_deriv_xi(xi);
    const double dx = dNdxi.dot(X.col(0)), dy = dNdxi.dot(X.col(1));
    const double J = std::sqrt(dx * dx + dy * dy);
    const double c = dx / J, s = dy / J;
    AxialKin k;
    k.J = J;
    k.Be.setZero();
    for (int i = 0; i < 3; ++i) {
        const double dNds = dNdxi(i) / J;
        k.Be(0, 2 * i + 0) = c * dNds;
        k.Be(0, 2 * i + 1) = s * dNds;
    }
    return k;
}

// Tension-only + N_p return mapping (see §8). ε = the current axial strain, ep_c = the
// committed plastic ε. Returns: N (force), ep_new (new plastic ε), Dt (algorithmic tangent).
struct AxialReturn { double N, ep_new, Dt; };
inline AxialReturn axial_return(const GeogridProps& p, double eps, double ep_c) {
    const double N_tr = p.EA * (eps - ep_c);
    if (p.Np > 0.0 && N_tr >= p.Np)         // tensile yield — permanent
        return {p.Np, eps - p.Np / p.EA, 0.0};
    if (N_tr <= 0.0)                         // compression cut-off — reversible slack
        return {0.0, ep_c, 0.0};
    return {N_tr, ep_c, p.EA};               // elastic
}

// 2-point Gauss (ε is linear on the quadratic edge → exact).
inline constexpr int kGaussCount = 2;
inline std::array<double, 2> gauss_xi() {
    constexpr double g2 = 0.5773502691896257;  // 1/sqrt(3)
    return {-g2, g2};
}

}  // namespace katai::core::geogrid
