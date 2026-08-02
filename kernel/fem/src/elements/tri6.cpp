#include <katai/fem/elements/tri6.hpp>

#include <Eigen/Dense>  // determinant + inverse

namespace katai::core::tri6 {

ShapeValues shape_functions(double xi, double eta) {
    const double l1 = 1.0 - xi - eta;  // L1
    const double l2 = xi;              // L2
    const double l3 = eta;             // L3

    ShapeValues n;
    n(0) = l1 * (2.0 * l1 - 1.0);  // corner 1
    n(1) = l2 * (2.0 * l2 - 1.0);  // corner 2
    n(2) = l3 * (2.0 * l3 - 1.0);  // corner 3
    n(3) = 4.0 * l1 * l2;          // orta-kenar 1-2
    n(4) = 4.0 * l2 * l3;          // orta-kenar 2-3
    n(5) = 4.0 * l3 * l1;          // orta-kenar 3-1
    return n;
}

ShapeDerivs shape_derivatives_natural(double xi, double eta) {
    const double l1 = 1.0 - xi - eta;
    const double l2 = xi;
    const double l3 = eta;

    ShapeDerivs d;
    // column 0: dN/dxi
    d(0, 0) = 1.0 - 4.0 * l1;
    d(1, 0) = 4.0 * l2 - 1.0;
    d(2, 0) = 0.0;
    d(3, 0) = 4.0 * (l1 - l2);
    d(4, 0) = 4.0 * l3;
    d(5, 0) = -4.0 * l3;
    // column 1: dN/deta
    d(0, 1) = 1.0 - 4.0 * l1;
    d(1, 1) = 0.0;
    d(2, 1) = 4.0 * l3 - 1.0;
    d(3, 1) = -4.0 * l2;
    d(4, 1) = 4.0 * l2;
    d(5, 1) = 4.0 * (l1 - l3);
    return d;
}

ShapeGradients strain_displacement(const NodeCoords& nodes, double xi,
                                   double eta) {
    const ShapeDerivs dn_natural = shape_derivatives_natural(xi, eta);

    // Jacobian: J = dN_nat^T * coords  (J(a,b) = sum_i dN_i/dnat_a * coord_i_b)
    const Eigen::Matrix2d jacobian = dn_natural.transpose() * nodes;

    ShapeGradients out;
    out.det_jacobian = jacobian.determinant();

    // Physical derivatives: dN_phys = dN_nat * J^{-T}  (row i: [dN_i/dx, dN_i/dy])
    const Eigen::Matrix2d inv_jacobian = jacobian.inverse();
    const ShapeDerivs dn_physical = dn_natural * inv_jacobian.transpose();

    // B (3x12): columns 2i, 2i+1 for each node i.
    for (int i = 0; i < kNodeCount; ++i) {
        const double dndx = dn_physical(i, 0);
        const double dndy = dn_physical(i, 1);
        out.B(0, 2 * i) = dndx;
        out.B(1, 2 * i + 1) = dndy;
        out.B(2, 2 * i) = dndy;
        out.B(2, 2 * i + 1) = dndx;
    }
    return out;
}

std::array<GaussPoint, kGaussCount> gauss_points() {
    // The classic 3-point rule: the interior image of the edge midpoints, equal weights.
    constexpr double a = 1.0 / 6.0;
    constexpr double b = 2.0 / 3.0;
    constexpr double w = 1.0 / 6.0;  // sum = 1/2 (reference triangle area)
    return {{{a, a, w}, {b, a, w}, {a, b, w}}};
}

ElementMatrix element_stiffness(const NodeCoords& nodes,
                                const ConstitutiveMatrix& constitutive,
                                double thickness) {
    ElementMatrix ke = ElementMatrix::Zero();
    for (const GaussPoint& gp : gauss_points()) {
        const ShapeGradients g = strain_displacement(nodes, gp.xi, gp.eta);
        const double scale = gp.weight * g.det_jacobian * thickness;
        ke.noalias() += scale * (g.B.transpose() * constitutive * g.B);
    }
    return ke;
}

} // namespace katai::core::tri6
