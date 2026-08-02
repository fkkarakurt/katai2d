#pragma once
// Linear elastic material — plane-strain constitutive matrix (Decision D8: model #1).
//
// Voigt notation (engineering shear strain):
//     stress = [sxx, syy, sxy]^T,  strain = [exx, eyy, gxy]^T
// for plane strain:
//     D = E / ((1+v)(1-2v)) *
//         | 1-v    v       0        |
//         |  v    1-v      0        |
//         |  0     0    (1-2v)/2    |

#include <Eigen/Core>

namespace katai::core {

struct LinearElastic {
    double youngs_modulus = 0.0;  // E  [force/area]
    double poisson_ratio = 0.0;   // v  [-]

    // Plane-strain elastic constitutive matrix (3x3, symmetric positive definite, v<0.5).
    Eigen::Matrix3d plane_strain_matrix() const {
        const double e = youngs_modulus;
        const double v = poisson_ratio;
        const double factor = e / ((1.0 + v) * (1.0 - 2.0 * v));

        Eigen::Matrix3d d = Eigen::Matrix3d::Zero();
        d(0, 0) = factor * (1.0 - v);
        d(0, 1) = factor * v;
        d(1, 0) = factor * v;
        d(1, 1) = factor * (1.0 - v);
        d(2, 2) = factor * (1.0 - 2.0 * v) / 2.0;
        return d;
    }

    // Axisymmetric elastic matrix (4x4) for strain/stress order
    // [r, z, rz(shear), theta(hoop)]. The hoop direction is a full strain (unlike
    // plane strain where eps_zz = 0), so D couples r, z and theta isotropically.
    Eigen::Matrix4d axisymmetric_matrix() const {
        const double e = youngs_modulus, v = poisson_ratio;
        const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
        Eigen::Matrix4d d = Eigen::Matrix4d::Zero();
        d(0, 0) = d(1, 1) = d(3, 3) = f * (1.0 - v);
        d(0, 1) = d(1, 0) = d(0, 3) = d(3, 0) = d(1, 3) = d(3, 1) = f * v;
        d(2, 2) = f * (1.0 - 2.0 * v) / 2.0;  // shear r-z
        return d;
    }
};

} // namespace katai::core
