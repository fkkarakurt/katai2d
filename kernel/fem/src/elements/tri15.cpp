#include <katai/fem/elements/tri15.hpp>

#include <Eigen/Dense>  // determinant + inverse

namespace katai::core::tri15 {
namespace {

// Barycentric integer indices (a,b,c), a+b+c=4, in node order (see the header).
constexpr int kBary[15][3] = {
    {4, 0, 0}, {0, 4, 0}, {0, 0, 4},              // corners
    {3, 1, 0}, {2, 2, 0}, {1, 3, 0},              // kenar 1-2
    {0, 3, 1}, {0, 2, 2}, {0, 1, 3},              // kenar 2-3
    {1, 0, 3}, {2, 0, 2}, {3, 0, 1},              // kenar 3-1
    {2, 1, 1}, {1, 2, 1}, {1, 1, 2}};             // interior

// Single-variable Lagrange factor A(k,L) = prod_{p<k} (4L-p)/(k-p) and its derivative
// dA/dL. Closed form in the variable u = 4L (du/dL = 4).
struct AVal { double value; double deriv; };  // {A, dA/dL}
AVal lagrange_factor(int k, double L) {
    const double u = 4.0 * L;
    double a, da_du;  // A ve dA/du
    switch (k) {
        case 0: a = 1.0;                                  da_du = 0.0; break;
        case 1: a = u;                                    da_du = 1.0; break;
        case 2: a = u * (u - 1.0) / 2.0;                  da_du = (2.0 * u - 1.0) / 2.0; break;
        case 3: a = u * (u - 1.0) * (u - 2.0) / 6.0;      da_du = (3.0 * u * u - 6.0 * u + 2.0) / 6.0; break;
        default: a = u * (u - 1.0) * (u - 2.0) * (u - 3.0) / 24.0;
                 da_du = (4.0 * u * u * u - 18.0 * u * u + 22.0 * u - 6.0) / 24.0; break;
    }
    return {a, da_du * 4.0};  // dA/dL = dA/du * du/dL
}

} // namespace

ShapeValues shape_functions(double xi, double eta) {
    const double L[3] = {1.0 - xi - eta, xi, eta};
    ShapeValues n;
    for (int i = 0; i < kNodeCount; ++i) {
        n(i) = lagrange_factor(kBary[i][0], L[0]).value *
               lagrange_factor(kBary[i][1], L[1]).value *
               lagrange_factor(kBary[i][2], L[2]).value;
    }
    return n;
}

ShapeDerivs shape_derivatives_natural(double xi, double eta) {
    const double L[3] = {1.0 - xi - eta, xi, eta};
    ShapeDerivs d;
    for (int i = 0; i < kNodeCount; ++i) {
        const AVal a1 = lagrange_factor(kBary[i][0], L[0]);  // A(a, L1)
        const AVal a2 = lagrange_factor(kBary[i][1], L[1]);  // A(b, L2)
        const AVal a3 = lagrange_factor(kBary[i][2], L[2]);  // A(c, L3)
        // L1 = 1-xi-eta (dL1/dxi = dL1/deta = -1), L2 = xi, L3 = eta.
        d(i, 0) = -a1.deriv * a2.value * a3.value + a1.value * a2.deriv * a3.value;
        d(i, 1) = -a1.deriv * a2.value * a3.value + a1.value * a2.value * a3.deriv;
    }
    return d;
}

ShapeGradients strain_displacement(const NodeCoords& nodes, double xi,
                                   double eta) {
    const ShapeDerivs dn_natural = shape_derivatives_natural(xi, eta);

    // Jacobian: J(a,b) = sum_i dN_i/dnat_a * coord_i_b.
    const Eigen::Matrix2d jacobian = dn_natural.transpose() * nodes;

    ShapeGradients out;
    out.det_jacobian = jacobian.determinant();

    // Physical derivatives: dN_phys = dN_nat * J^{-T} (row i: [dN_i/dx, dN_i/dy]).
    const Eigen::Matrix2d inv_jacobian = jacobian.inverse();
    const ShapeDerivs dn_physical = dn_natural * inv_jacobian.transpose();

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
    // Dunavant (1985) degree-6, 12-point symmetric rule. Weights normalized to unit
    // area; scaled by 0.5 for the reference triangle (area 1/2).
    constexpr double a1 = 0.063089014491502, w1 = 0.050844906370207;
    constexpr double a2 = 0.249286745170910, w2 = 0.116786275726379;
    constexpr double a3 = 0.310352451033785, b3 = 0.053145049844816,
                     c3 = 0.636502499121399, w3 = 0.082851075618374;
    const double W1 = 0.5 * w1, W2 = 0.5 * w2, W3 = 0.5 * w3;
    // (xi, eta) = (L2, L3) for each orbit permutation.
    return {{
        // S21 orbit 1 (a1, a1, 1-2a1)
        {a1, a1, W1}, {1.0 - 2.0 * a1, a1, W1}, {a1, 1.0 - 2.0 * a1, W1},
        // S21 orbit 2
        {a2, a2, W2}, {1.0 - 2.0 * a2, a2, W2}, {a2, 1.0 - 2.0 * a2, W2},
        // S111 orbit (6 permutations of barycentric (a3, b3, c3))
        {b3, c3, W3}, {c3, b3, W3}, {a3, c3, W3},
        {c3, a3, W3}, {a3, b3, W3}, {b3, a3, W3},
    }};
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

} // namespace katai::core::tri15
