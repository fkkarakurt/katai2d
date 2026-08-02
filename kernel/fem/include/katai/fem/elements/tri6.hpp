#pragma once
// 6-node quadratic triangle element (tri6), plane strain — Decision D10.
//
// Node ordering (area/barycentric coordinates L1 = 1-xi-eta, L2 = xi, L3 = eta):
//     corners  : 1 @ (0,0),  2 @ (1,0),  3 @ (0,1)
//     mid-edge : 4 @ 1-2,    5 @ 2-3,    6 @ 3-1
//
// DOF order: [u1, v1, u2, v2, ..., u6, v6]  (ux, uy per node).
// Strain: [exx, eyy, gxy] = [du/dx, dv/dy, du/dy + dv/dx].

#include <array>

#include <Eigen/Core>

namespace katai::core::tri6 {

inline constexpr int kNodeCount = 6;
inline constexpr int kDofCount = 12;
inline constexpr int kGaussCount = 3;

using NodeCoords = Eigen::Matrix<double, 6, 2>;        // row i: (x_i, y_i)
using ShapeValues = Eigen::Matrix<double, 6, 1>;       // N_i
using ShapeDerivs = Eigen::Matrix<double, 6, 2>;       // [dN/dxi, dN/deta] columns
using StrainDispMatrix = Eigen::Matrix<double, 3, 12>; // B
using ElementMatrix = Eigen::Matrix<double, 12, 12>;   // Ke
using ConstitutiveMatrix = Eigen::Matrix3d;            // D

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

// Element stiffness matrix Ke = ∫ B^T D B t dA (3-point Gauss).
// thickness: unit thickness in plane strain (default 1).
ElementMatrix element_stiffness(const NodeCoords& nodes,
                                const ConstitutiveMatrix& constitutive,
                                double thickness = 1.0);

// 3-point triangle Gauss rule (degree 2 — exact for the quadratic element). The
// weights sum to 1/2 = the reference triangle area.
struct GaussPoint {
    double xi;
    double eta;
    double weight;
};
std::array<GaussPoint, kGaussCount> gauss_points();

} // namespace katai::core::tri6
