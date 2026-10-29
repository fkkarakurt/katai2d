#pragma once
// Geogrid (reinforcement/geosynthetic) — 3-node axial membrane element (P2.4). The
// bending-free (EI=0, no shear, no φ-DOF) form of the plate: only the axial force
// N = EA·ε. Sits on a mesh edge (tri6 edge = 3 nodes), 2 DOFs per node (u_x, u_y) —
// SHARED with the soil.
//
// Constitutive law (elastic, tension-only, with an optional tensile capacity):
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
    // A MEMBRANE OF REVOLUTION (axisymmetric model, x = radius): besides the meridional strain the
    // grid stretches round its circumference, eps_t = u_r / r, and carries the hoop force
    // N_t = EA eps_t in the other direction of the grid -- the same stiffness, no Poisson coupling,
    // tension only like the meridional force. Everything is integrated per radian (weight r).
    // The hoop has no plastic state, so a capacity N_p is not taken here (the driver refuses it).
    bool axisym = false;
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

// The hoop law of a membrane of revolution, eps_t = u_r / r (grid_kin below; a stress point lies
// inside the element, so r > 0 there even when an end of the grid is on the axis): tension only
// and elastic -- slack, reversibly, when the ring would shorten.
inline AxialReturn hoop_return(const GeogridProps& p, double eps_t) {
    if (eps_t <= 0.0) return {0.0, 0.0, 0.0};
    return {p.EA * eps_t, 0.0, p.EA};
}

// 2-point Gauss (ε is linear on the quadratic edge → exact).
inline constexpr int kGaussCount = 2;
inline std::array<double, 2> gauss_xi() {
    constexpr double g2 = 0.5773502691896257;  // 1/sqrt(3)
    return {-g2, g2};
}

// ---------------------------------------------------------------------------------------------
// THE ELEMENT OF EITHER ORDER. A geogrid follows the edges of the ground it lies in: three nodes on
// a 6-noded triangle's edge, FIVE on a 15-noded one. Until 2026-10-27 a grid on a 15-noded mesh was
// built of quadratic pieces over half-edges: joined to the ground at every node, but stretched as a
// quadratic between them while the ground and a load on the line are quartic. Where the grid is the
// stiffer of the two the nodes then carry the quartic edge's load shares (7, 32, 12, 32, 7 / 90)
// against the grid's own (1, 4, 2, 4, 1 / 12), and the force oscillates from node to node: measured
// on a cylinder of grid in soft ground, a hoop force 5.8% high at mid-height, and 0.2% of q b on a
// plane-strain bar. The quartic element follows the ground's edge exactly.
//
// Node order of an element: 3 nodes [A, B, M] at xi = -1, +1, 0; 5 nodes [A, B, M, Q1, Q3] at
// xi = -1, +1, 0, -1/2, +1/2 (Q1 between A and M, Q3 between M and B). Its stress points: 2-point
// Gauss for 3 nodes (the historical element, unchanged), 4-point for 5 (exact for the quartic's
// cubic strain squared). The state vectors keep kStateStride slots per element, whatever its order.
inline constexpr int kMaxNodes = 5;
inline constexpr int kStateStride = 4;
using GridCoords = Eigen::Matrix<double, kMaxNodes, 2>;
using GridRow = Eigen::Matrix<double, 1, 2 * kMaxNodes>;
using GridDof = Eigen::Matrix<double, 2 * kMaxNodes, 1>;

inline int point_count(int n) { return n == 5 ? 4 : 2; }
inline void point(int n, int q, double& xi, double& w) {
    if (n == 5) {
        constexpr double a = 0.3399810435848563, wa = 0.6521451548625461;
        constexpr double b = 0.8611363115940526, wb = 0.3478548451374538;
        constexpr double x[4] = {-b, -a, a, b}, ww[4] = {wb, wa, wa, wb};
        xi = x[q]; w = ww[q];
        return;
    }
    xi = gauss_xi()[(size_t)q]; w = 1.0;
}
// Lagrange shape functions and their xi-derivatives over the element's own node positions.
inline void grid_shape(int n, double xi, double* N, double* dN) {
    if (n == 3) {   // the historical element's own expressions, so that it is reproduced to the bit
        const Eigen::Vector3d d = shape_deriv_xi(xi);
        N[0] = 0.5 * xi * (xi - 1.0); N[1] = 0.5 * xi * (xi + 1.0); N[2] = 1.0 - xi * xi;
        dN[0] = d(0); dN[1] = d(1); dN[2] = d(2);
        return;
    }
    constexpr double X[5] = {-1.0, 1.0, 0.0, -0.5, 0.5};
    for (int i = 0; i < n; ++i) {
        double v = 1.0, d = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double p = 1.0 / (X[i] - X[j]);
            for (int k = 0; k < n; ++k)
                if (k != i && k != j) p *= (xi - X[k]) / (X[i] - X[k]);
            d += p;
            v *= (xi - X[j]) / (X[i] - X[j]);
        }
        N[i] = v; dN[i] = d;
    }
}
// At a stress point: the weight w J, the axial row Be and, for a membrane of revolution, the
// radius r and the hoop row Bt (both on the element's 2n DOFs, padded with zeros).
struct GridKin { double wJ = 0.0, xi = 0.0, r = 1.0; GridRow Be, Bt; };
inline GridKin grid_kin(const GridCoords& X, int n, int q, bool axisym) {
    double N[5], dN[5], xi, w;
    point(n, q, xi, w);
    grid_shape(n, xi, N, dN);
    double dx = 0.0, dy = 0.0, r = 0.0;
    for (int i = 0; i < n; ++i) { dx += dN[i] * X(i, 0); dy += dN[i] * X(i, 1); r += N[i] * X(i, 0); }
    const double J = std::sqrt(dx * dx + dy * dy), c = dx / J, s = dy / J;
    GridKin k;
    k.xi = xi;
    k.wJ = w * J;
    k.Be.setZero(); k.Bt.setZero();
    for (int i = 0; i < n; ++i) {
        const double dNds = dN[i] / J;
        k.Be(0, 2 * i + 0) = c * dNds;
        k.Be(0, 2 * i + 1) = s * dNds;
    }
    if (axisym) {
        k.r = r;
        k.wJ *= r;
        for (int i = 0; i < n; ++i) k.Bt(0, 2 * i) = N[i] / r;
    }
    return k;
}
// The position and the interpolated displacement at a point of the element.
inline void grid_point_at(const GridCoords& X, const GridDof& u, int n, double xi,
                          double& x, double& y, double& ux, double& uy) {
    double N[5], dN[5];
    grid_shape(n, xi, N, dN);
    x = y = ux = uy = 0.0;
    for (int i = 0; i < n; ++i) {
        x += N[i] * X(i, 0); y += N[i] * X(i, 1);
        ux += N[i] * u(2 * i); uy += N[i] * u(2 * i + 1);
    }
}

}  // namespace katai::core::geogrid
