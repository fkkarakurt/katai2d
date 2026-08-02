#pragma once
// 15-node quartic (4th-order) triangle element (tri15), plane strain — P1.5.
//
// Follows tri6's verified template (same interface, higher-order basis).
// 15 nodes = 3 corners + 3 per edge (9 total) + 3 interior. The nodes sit at the
// (xi/4, eta/4) lattice points of the natural triangle (barycentric integer
// triples with a+b+c=4).
//
// Node ordering (L1 = 1-xi-eta, L2 = xi, L3 = eta; (a,b,c) = 4*(L1,L2,L3)):
//     corners   : 0 @ (4,0,0)=(0,0)   1 @ (0,4,0)=(1,0)   2 @ (0,0,4)=(0,1)
//     edge 1-2  : 3 (3,1,0)  4 (2,2,0)  5 (1,3,0)        [eta=0,  xi=1/4,2/4,3/4]
//     edge 2-3  : 6 (0,3,1)  7 (0,2,2)  8 (0,1,3)        [xi+eta=1 hypotenuse]
//     edge 3-1  : 9 (1,0,3) 10 (2,0,2) 11 (3,0,1)        [xi=0,  eta=3/4,2/4,1/4]
//     interior  : 12 (2,1,1) 13 (1,2,1) 14 (1,1,2)
//
// DOF order: [u0, v0, u1, v1, ..., u14, v14]. Strain [exx, eyy, gxy].
//
// The shape functions come from the general Lagrange product formula:
//   N(a,b,c) = A(a,L1) * A(b,L2) * A(c,L3),
//   A(k,L)   = prod_{p=0}^{k-1} (4L - p) / (k - p),   A(0,L) = 1.
// This removes the error risk of writing 15 quartic polynomials by hand; verified
// by isolated tests (partition of unity, Kronecker delta, zero derivative sum,
// Gauss degree exactness, patch test).

#include <array>

#include <Eigen/Core>

namespace katai::core::tri15 {

inline constexpr int kNodeCount = 15;
inline constexpr int kDofCount = 30;
inline constexpr int kGaussCount = 12;

using NodeCoords = Eigen::Matrix<double, 15, 2>;        // row i: (x_i, y_i)
using ShapeValues = Eigen::Matrix<double, 15, 1>;       // N_i
using ShapeDerivs = Eigen::Matrix<double, 15, 2>;       // [dN/dxi, dN/deta]
using StrainDispMatrix = Eigen::Matrix<double, 3, 30>;  // B
using ElementMatrix = Eigen::Matrix<double, 30, 30>;    // Ke
using ConstitutiveMatrix = Eigen::Matrix3d;             // D

// Shape functions N_i(xi, eta).
ShapeValues shape_functions(double xi, double eta);

// Derivatives w.r.t. natural coordinates: column 0 = dN/dxi, column 1 = dN/deta.
ShapeDerivs shape_derivatives_natural(double xi, double eta);

// The B matrix and Jacobian determinant at a point.
struct ShapeGradients {
    StrainDispMatrix B = StrainDispMatrix::Zero();
    double det_jacobian = 0.0;
};
ShapeGradients strain_displacement(const NodeCoords& nodes, double xi, double eta);

// Element stiffness matrix Ke = ∫ B^T D B t dA (12-point Gauss, degree 6).
ElementMatrix element_stiffness(const NodeCoords& nodes,
                                const ConstitutiveMatrix& constitutive,
                                double thickness = 1.0);

// 12-point triangle Gauss rule (Dunavant 1985, degree 6 — the quartic element's
// B^T D B integrand is degree 6). Weights sum to 1/2 = the reference triangle area.
struct GaussPoint {
    double xi;
    double eta;
    double weight;
};
std::array<GaussPoint, kGaussCount> gauss_points();

} // namespace katai::core::tri15
