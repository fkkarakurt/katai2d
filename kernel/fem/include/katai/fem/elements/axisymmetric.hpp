#pragma once
// Axisymmetric (r-z) element kernels — P1.6. One of PLAXIS's two analysis modes (the
// other is plane strain). Geometry and shape functions are the SAME as the base element
// (tri6/tri15); the only differences are (i) the strain's 4th component, the hoop
// εθ = u_r / r, (ii) the 4x4 constitutive matrix (the θ direction is a full strain, not
// the plane-strain σ_zz reaction), (iii) the r factor in the integration weight
// (∫·r dr dz, per radian).
//
// The module is therefore templated over the base element E and reuses
// E::shape_functions / shape_derivatives_natural / gauss_points. Coordinates are
// interpreted as (r, z) (column 0 = r, column 1 = z).
//
// Strain/stress order: [ε_r, ε_z, γ_rz, ε_θ] / [σ_r, σ_z, τ_rz, σ_θ].

#include <Eigen/Core>
#include <Eigen/LU>  // inverse

namespace katai::core::axisym {

template <class E>
struct Gradients {
    Eigen::Matrix<double, 4, E::kDofCount> B = Eigen::Matrix<double, 4, E::kDofCount>::Zero();
    double det_jacobian = 0.0;
    double radius = 0.0;  // r at this point = Σ N_i r_i
};

// The axisymmetric B matrix (4 x kDofCount), Jacobian determinant and r at a point.
template <class E>
Gradients<E> strain_displacement(const typename E::NodeCoords& nodes, double xi,
                                 double eta) {
    const auto dn_nat = E::shape_derivatives_natural(xi, eta);   // N x 2
    const typename E::ShapeValues n = E::shape_functions(xi, eta);

    const Eigen::Matrix2d jacobian = dn_nat.transpose() * nodes;
    const Eigen::Matrix2d inv_jacobian = jacobian.inverse();
    const auto dn_phys = dn_nat * inv_jacobian.transpose();      // row i: [dN/dr, dN/dz]

    Gradients<E> out;
    out.det_jacobian = jacobian.determinant();
    for (int i = 0; i < E::kNodeCount; ++i) out.radius += n(i) * nodes(i, 0);

    for (int i = 0; i < E::kNodeCount; ++i) {
        const double dndr = dn_phys(i, 0), dndz = dn_phys(i, 1);
        out.B(0, 2 * i) = dndr;            // ε_r = ∂u_r/∂r
        out.B(1, 2 * i + 1) = dndz;        // ε_z = ∂u_z/∂z
        out.B(2, 2 * i) = dndz;            // γ_rz = ∂u_r/∂z + ∂u_z/∂r
        out.B(2, 2 * i + 1) = dndr;
        out.B(3, 2 * i) = n(i) / out.radius;  // ε_θ = u_r / r
    }
    return out;
}

// Element stiffness matrix Ke = ∫ B^T D B r dA (per radian; the 2π constant is dropped,
// and being dropped identically in both stiffness and loads it does not affect results).
template <class E>
typename E::ElementMatrix element_stiffness(const typename E::NodeCoords& nodes,
                                            const Eigen::Matrix4d& constitutive) {
    typename E::ElementMatrix ke = E::ElementMatrix::Zero();
    for (const auto& gp : E::gauss_points()) {
        const Gradients<E> g = strain_displacement<E>(nodes, gp.xi, gp.eta);
        const double scale = gp.weight * g.det_jacobian * g.radius;
        ke.noalias() += scale * (g.B.transpose() * constitutive * g.B);
    }
    return ke;
}

} // namespace katai::core::axisym
