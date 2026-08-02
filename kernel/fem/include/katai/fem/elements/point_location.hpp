#pragma once
// Point location + inverse isoparametric mapping (the embedded beam / pile row
// foundation, Phase A.4). Given a physical point (x,y), finds the element containing it
// and its local (ξ,η) coordinates → the soil shape functions N_s can be evaluated there
// (mesh-NONCONFORMING coupling; PLAXIS Sci.Man §7.5 "virtual node"). The inverse mapping
// is a Newton solve of x(ξ,η)=Σ N_i(ξ,η) X_i; for curved tri6/tri15 the Jacobian is
// J=Xᵀ·(dN/dξ,dN/dη). Element-generic (element_traits).
//
// Also generally useful for result-point probing (reading stress/displacement at an
// arbitrary location).

#include <cmath>

#include <Eigen/Dense>

#include <katai/fem/elements/element_traits.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core::ploc {

struct LocalCoord { double xi = 0.0, eta = 0.0; bool converged = false; };

// Ters izoparametrik haritalama: fiziksel (px,py) → yerel (ξ,η), Newton iterasyonu.
// X: the element node coordinates (E::NodeCoords). Start at the centroid (1/3,1/3).
template <class E>
LocalCoord physical_to_local(const typename E::NodeCoords& X, double px, double py,
                             int max_iter = 30, double tol = 1e-12) {
    LocalCoord lc; lc.xi = 1.0 / 3.0; lc.eta = 1.0 / 3.0;
    const Eigen::Vector2d target(px, py);
    for (int it = 0; it < max_iter; ++it) {
        const auto N = E::shape_functions(lc.xi, lc.eta);            // N×1
        const auto dN = E::shape_derivatives_natural(lc.xi, lc.eta); // N×2 [dξ, dη]
        const Eigen::Vector2d xcur(X.col(0).dot(N), X.col(1).dot(N));
        const Eigen::Vector2d r = xcur - target;
        if (r.norm() <= tol) { lc.converged = true; return lc; }
        Eigen::Matrix2d J;  // J = Xᵀ dN : [∂x/∂ξ ∂x/∂η; ∂y/∂ξ ∂y/∂η]
        J(0, 0) = X.col(0).dot(dN.col(0)); J(0, 1) = X.col(0).dot(dN.col(1));
        J(1, 0) = X.col(1).dot(dN.col(0)); J(1, 1) = X.col(1).dot(dN.col(1));
        const double det = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
        if (std::fabs(det) < 1e-300) break;
        const Eigen::Vector2d d = J.inverse() * r;
        lc.xi -= d(0); lc.eta -= d(1);
    }
    // Convergence check after the last increment.
    const auto N = E::shape_functions(lc.xi, lc.eta);
    const Eigen::Vector2d xcur(X.col(0).dot(N), X.col(1).dot(N));
    lc.converged = (xcur - target).norm() <= 1e-9 * (1.0 + target.norm());
    return lc;
}

// Inside the reference triangle (ξ≥0, η≥0, ξ+η≤1)? Tolerance for edges/adjacency.
inline bool inside_ref_triangle(double xi, double eta, double tol = 1e-9) {
    return xi >= -tol && eta >= -tol && (xi + eta) <= 1.0 + tol;
}

struct Located { int element = -1; double xi = 0.0, eta = 0.0; bool found = false; };

namespace detail {
template <class E>
Located locate_impl(const mesh::Mesh& mesh, double px, double py, double tol) {
    constexpr int N = E::kNodeCount;
    for (int e = 0; e < mesh.element_count; ++e) {
        typename E::NodeCoords X;
        for (int k = 0; k < N; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        // Fast rejection: corner bounding box (corners are the first 3 nodes).
        double xmin = X(0,0), xmax = X(0,0), ymin = X(0,1), ymax = X(0,1);
        for (int k = 1; k < 3; ++k) {
            xmin = std::min(xmin, X(k,0)); xmax = std::max(xmax, X(k,0));
            ymin = std::min(ymin, X(k,1)); ymax = std::max(ymax, X(k,1));
        }
        const double pad = 0.05 * (xmax - xmin + ymax - ymin);
        if (px < xmin - pad || px > xmax + pad || py < ymin - pad || py > ymax + pad) continue;
        const LocalCoord lc = physical_to_local<E>(X, px, py);
        if (lc.converged && inside_ref_triangle(lc.xi, lc.eta, tol))
            return {e, lc.xi, lc.eta, true};
    }
    return {};
}
}  // namespace detail

// Finds the element containing (px,py) + its local coordinates. Element type from mesh.nodes_per_element.
inline Located locate_point(const mesh::Mesh& mesh, double px, double py, double tol = 1e-7) {
    if (mesh.nodes_per_element == Tri15Element::kNodeCount)
        return detail::locate_impl<Tri15Element>(mesh, px, py, tol);
    return detail::locate_impl<Tri6Element>(mesh, px, py, tol);
}

}  // namespace katai::core::ploc
