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
    // AN ANISOTROPIC OPERATOR, where the material has one (a transversely isotropic jointed rock):
    // the full 4 x 4 on [xx, yy, zz, xy] at the modulus anisotropic_E, which a depth profile of
    // youngs_modulus then scales as a whole. Off, the two numbers above are the whole of it. The
    // linear paths (the linear dynamic one above all) assemble from this struct alone, and an
    // anisotropic rock assembled from its E1 / nu1 would quietly be isotropic there.
    bool anisotropic = false;
    Eigen::Matrix4d anisotropic_D = Eigen::Matrix4d::Zero();
    double anisotropic_E = 0.0;

    // Plane-strain elastic constitutive matrix (3x3, symmetric positive definite, v<0.5).
    Eigen::Matrix3d plane_strain_matrix() const {
        if (anisotropic) {
            const double s = anisotropic_E > 0.0 ? youngs_modulus / anisotropic_E : 1.0;
            const int idx[3] = {0, 1, 3};
            Eigen::Matrix3d d;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) d(i, j) = s * anisotropic_D(idx[i], idx[j]);
            return d;
        }
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
        if (anisotropic) {   // [r, z, rz, theta] from [xx, yy, zz, xy]
            const double s = anisotropic_E > 0.0 ? youngs_modulus / anisotropic_E : 1.0;
            const int idx[4] = {0, 1, 3, 2};
            Eigen::Matrix4d d;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) d(i, j) = s * anisotropic_D(idx[i], idx[j]);
            return d;
        }
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
