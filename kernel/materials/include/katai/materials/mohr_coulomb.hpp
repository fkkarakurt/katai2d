#pragma once
// Mohr-Coulomb elasto-plasticity — kinematic and constitutive kernels (P1.2a).
//
// This header provides the *building blocks* on which the principal-stress
// return-mapping algorithm (Clausen, Damkilde & Andersen, 2006/2007) is erected
// in P1.2b: the plane-strain elastic predictor, the spectral decomposition of
// the trial stress, and the Mohr-Coulomb yield function. Each block is a pure
// function and is verified in isolation before the (heavier, branch-rich) return
// mapping is layered on top.
//
// Sign convention: TENSION POSITIVE, consistent with the continuum-mechanics
// convention used throughout the solver (e.g. a surface pressure p produces
// sigma_yy = -p). Consequently the most tensile principal stress is sigma_1 and
// the most compressive is sigma_3.
//
// Plane-strain subtlety: the kinematic constraint eps_zz == 0 does NOT imply
// sigma_zz == 0. The out-of-plane normal stress evolves with loading and may be
// the intermediate OR an extreme principal stress; it must therefore be carried
// explicitly in the material state and included in the yield evaluation. Because
// tau_xz = tau_yz = 0, sigma_zz is itself a principal stress, and the remaining
// two principals follow from the in-plane 2x2 stress block.

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

namespace katai::core {

// Mohr-Coulomb parameters (angles in radians). Non-associated flow is obtained
// whenever the dilatancy angle psi differs from the friction angle phi.
struct MohrCoulombParams {
    double youngs_modulus = 0.0;   // E
    double poisson_ratio = 0.0;    // nu
    double cohesion = 0.0;         // c >= 0
    double friction_angle = 0.0;   // phi  [rad]
    double dilatancy_angle = 0.0;  // psi  [rad], psi <= phi
    // Rankine tension cap (PLAXIS MMM sec 3.2 Eq 3-11 / sec 3.3.10): three extra yield
    // planes f_t,i = sigma_i - sigma_t <= 0 with ASSOCIATED flow; with sorted principals
    // they reduce to sigma_1 <= sigma_t. Internally clamped to sigma_t <= c*cot(phi)
    // (the most tensile principal the MC set admits). Off by default so every existing
    // caller keeps the bit-identical pure-MC path.
    bool tension_cutoff = false;
    double tensile_strength = 0.0;  // sigma_t >= 0 (tension-positive)
};

// Lame constants derived from (E, nu).
struct LameConstants {
    double lambda = 0.0;
    double mu = 0.0;
};

inline LameConstants lame_from(double youngs_modulus, double poisson_ratio) {
    const double e = youngs_modulus, v = poisson_ratio;
    return {e * v / ((1.0 + v) * (1.0 - 2.0 * v)), e / (2.0 * (1.0 + v))};
}

// Plane-strain stress state carried by the constitutive integrator: the three
// in-plane Voigt components plus the out-of-plane normal stress sigma_zz.
struct PlaneStrainStress {
    Eigen::Vector3d in_plane = Eigen::Vector3d::Zero();  // [sxx, syy, sxy]
    double zz = 0.0;                                     // sigma_zz
};

// Elastic predictor: sigma_trial = sigma_n + D_e : d_eps, enforcing d_eps_zz = 0.
// The in-plane response coincides with the standard 3x3 plane-strain operator;
// the out-of-plane increment is d_sigma_zz = lambda (d_eps_xx + d_eps_yy), which
// is the through-thickness reaction that maintains eps_zz = 0.
inline PlaneStrainStress elastic_predictor(const PlaneStrainStress& sigma_n,
                                           const Eigen::Vector3d& dstrain,
                                           const LameConstants& lame) {
    const double lam = lame.lambda, mu = lame.mu;
    const double dexx = dstrain(0), deyy = dstrain(1), dgxy = dstrain(2);
    const double tr = dexx + deyy;  // d_eps_xx + d_eps_yy (d_eps_zz = 0)

    PlaneStrainStress out;
    out.in_plane(0) = sigma_n.in_plane(0) + lam * tr + 2.0 * mu * dexx;
    out.in_plane(1) = sigma_n.in_plane(1) + lam * tr + 2.0 * mu * deyy;
    out.in_plane(2) = sigma_n.in_plane(2) + mu * dgxy;  // tau = mu * gamma
    out.zz = sigma_n.zz + lam * tr;
    return out;
}

// Ordered principal stresses sigma_1 >= sigma_2 >= sigma_3 together with the
// orientation of the in-plane principal frame. The in-plane principal values are
// the eigenvalues of the symmetric 2x2 block [[sxx, sxy], [sxy, syy]]; sigma_zz
// is appended as the third principal. cos2t/sin2t encode the in-plane principal
// axis (theta measured from x), enabling reconstruction of the in-plane tensor
// after the return mapping modifies the principal magnitudes (used in P1.2b).
struct PrincipalStresses {
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;  // descending order
    double cos2t = 1.0, sin2t = 0.0;      // in-plane principal orientation
};

inline PrincipalStresses principal_stresses(const PlaneStrainStress& s) {
    const double sxx = s.in_plane(0), syy = s.in_plane(1), sxy = s.in_plane(2);
    const double mean = 0.5 * (sxx + syy);
    const double half_diff = 0.5 * (sxx - syy);
    const double radius = std::sqrt(half_diff * half_diff + sxy * sxy);

    PrincipalStresses p;
    // In-plane principal orientation (Mohr's circle); guard the isotropic case.
    if (radius > 0.0) {
        p.cos2t = half_diff / radius;
        p.sin2t = sxy / radius;
    }
    const double in_a = mean + radius;  // larger in-plane principal
    const double in_b = mean - radius;  // smaller in-plane principal

    // Merge the two in-plane principals with sigma_zz and sort descending.
    double v[3] = {in_a, in_b, s.zz};
    std::sort(v, v + 3, std::greater<double>());
    p.s1 = v[0];
    p.s2 = v[1];
    p.s3 = v[2];
    return p;
}

// Mohr-Coulomb yield function in principal stresses (tension positive):
//     f = (sigma_1 - sigma_3) + (sigma_1 + sigma_3) sin(phi) - 2 c cos(phi).
// Admissible (elastic) states satisfy f <= 0; f = 0 is the failure surface.
inline double mc_yield(const PrincipalStresses& p, double cohesion,
                       double friction_angle) {
    const double sphi = std::sin(friction_angle);
    const double cphi = std::cos(friction_angle);
    return (p.s1 - p.s3) + (p.s1 + p.s3) * sphi - 2.0 * cohesion * cphi;
}

// Outcome of the closed-form Mohr-Coulomb return mapping.
struct McReturn {
    PlaneStrainStress stress;  // admissible (returned) stress
    bool plastic = false;      // whether a plastic correction was applied
    bool tension = false;      // whether a tension cut-off plane was active (PLAXIS
                               // "tension point" classification; subset of plastic)
    // In-plane consistent (algorithmic) tangent D_T = d(sigma_inplane)/d(eps),
    // 3x3 in Voigt [sxx, syy, sxy] x [exx, eyy, gxy]. Only meaningful when
    // plastic == true; for an elastic step the caller uses the elastic operator.
    // Non-associated flow (psi != phi) makes it unsymmetric.
    Eigen::Matrix3d tangent = Eigen::Matrix3d::Zero();
    // Algorithmic stress-to-stress Jacobian d(sigma^ret)/d(sigma^tr), Voigt order
    // [xx, yy, xy, zz]. Multiplying by the elastic predictor map gives the
    // consistent tangent in either kinematics (plane strain or axisymmetric). It
    // is the identity for an elastic step. Only meaningful when plastic == true.
    Eigen::Matrix4d algo_jacobian = Eigen::Matrix4d::Identity();
};

// Assemble the plane-strain consistent tangent (and the 4x4 stress-to-stress Jacobian
// Psi = d sigma^ret/d sigma^tr in Voigt [xx,yy,xy,zz]) from any coaxial principal-space
// return mapping. Inputs: the trial spectral decomposition (in-plane Mohr-circle cos2t/
// sin2t/radius), the provenance src[i] of each sorted-descending trial principal
// (0 = larger in-plane, 1 = smaller in-plane, 2 = sigma_zz), the returned sorted
// principals ret[i], and the principal Jacobian J(i,j) = d(returned sorted i)/d(trial
// sorted j). Writing the eigenframe rotation (spin) + sigma_zz coupling ONCE lets both
// Mohr-Coulomb (closed-form region J) and Hardening Soil (J = D_pp C_e from the
// substepping continuum tangent) share it. (Sysala 2016 Sec. 5; both conventions are
// tension-positive, so the same De/reconstruction apply.)
struct PrincipalTangent {
    Eigen::Matrix3d tangent = Eigen::Matrix3d::Zero();           // d sigma_inplane / d eps
    Eigen::Matrix4d algo_jacobian = Eigen::Matrix4d::Identity(); // d sigma^ret / d sigma^tr
};

inline PrincipalTangent principal_consistent_tangent(
    double cos2t, double sin2t, double radius, const int src[3],
    const double ret[3], const Eigen::Matrix3d& J, const LameConstants& lame) {
    enum { kInPlaneA = 0, kInPlaneB = 1, kOutOfPlane = 2 };
    int slot[3] = {0, 0, 0};  // slot[source] = sorted position
    for (int i = 0; i < 3; ++i) slot[src[i]] = i;
    const double pa = ret[slot[kInPlaneA]], pb = ret[slot[kInPlaneB]];
    // Reindex J from sorted positions into source order (ia, ib, zz).
    Eigen::Matrix3d Jt;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) Jt(a, b) = J(slot[a], slot[b]);

    // Derivatives of the in-plane trial principals (ia=mean+radius, ib=mean-radius)
    // w.r.t. the in-plane trial Voigt components (sxx, syy, sxy).
    const Eigen::RowVector3d dia(0.5 + 0.5 * cos2t, 0.5 - 0.5 * cos2t, sin2t);
    const Eigen::RowVector3d dib(0.5 - 0.5 * cos2t, 0.5 + 0.5 * cos2t, -sin2t);
    const Eigen::RowVector3d dpa =
        Jt(kInPlaneA, kInPlaneA) * dia + Jt(kInPlaneA, kInPlaneB) * dib;
    const Eigen::RowVector3d dpb =
        Jt(kInPlaneB, kInPlaneA) * dia + Jt(kInPlaneB, kInPlaneB) * dib;
    const double dpa_zz = Jt(kInPlaneA, kOutOfPlane);
    const double dpb_zz = Jt(kInPlaneB, kOutOfPlane);
    const Eigen::Vector3d gpa(0.5 + 0.5 * cos2t, 0.5 - 0.5 * cos2t, 0.5 * sin2t);
    const Eigen::Vector3d gpb(0.5 - 0.5 * cos2t, 0.5 + 0.5 * cos2t, -0.5 * sin2t);
    Eigen::Matrix<double, 3, 4> Phi;
    Phi.leftCols(3) = gpa * dpa + gpb * dpb;
    Phi.col(3) = gpa * dpa_zz + gpb * dpb_zz;

    // Spin part: rotation of the in-plane principal frame. beta = (pa-pb)/(ia-ib) with
    // ia-ib = 2 radius keeps the 1/radius singularity removable; in the in-plane-isotropic
    // limit it tends to the coaxial principal-tangent difference.
    double beta;
    if (2.0 * radius > 1.0e-12 * (std::fabs(pa) + std::fabs(pb) + 1.0)) {
        beta = (pa - pb) / (2.0 * radius);
    } else {
        beta = 0.5 * (Jt(kInPlaneA, kInPlaneA) - Jt(kInPlaneA, kInPlaneB) -
                      Jt(kInPlaneB, kInPlaneA) + Jt(kInPlaneB, kInPlaneB));
    }
    const Eigen::Vector3d uc(1.0, -1.0, 0.0), us(0.0, 0.0, 1.0);
    const Eigen::RowVector3d dcos(0.5 * sin2t * sin2t, -0.5 * sin2t * sin2t,
                                  -cos2t * sin2t);
    const Eigen::RowVector3d dsin(-0.5 * sin2t * cos2t, 0.5 * sin2t * cos2t,
                                  cos2t * cos2t);
    Phi.leftCols(3) += beta * (uc * dcos + us * dsin);

    const Eigen::RowVector3d dpz =
        Jt(kOutOfPlane, kInPlaneA) * dia + Jt(kOutOfPlane, kInPlaneB) * dib;
    const double dpz_zz = Jt(kOutOfPlane, kOutOfPlane);

    PrincipalTangent out;
    out.algo_jacobian.topRows(3) = Phi;
    out.algo_jacobian.row(3) << dpz, dpz_zz;
    const double lam = lame.lambda, G = lame.mu;
    Eigen::Matrix3d De;
    De << lam + 2.0 * G, lam, 0.0,
          lam, lam + 2.0 * G, 0.0,
          0.0, 0.0, G;
    const Eigen::RowVector3d dszz(lam, lam, 0.0);
    out.tangent = Phi.leftCols(3) * De + Phi.col(3) * dszz;
    return out;
}

// Closed-form implicit (backward-Euler) return mapping for perfectly-plastic
// Mohr-Coulomb in principal stress space, after Sysala & Cermak (2016, arXiv:
// 1508.07435) and Clausen, Damkilde & Andersen (2007). Given an elastic trial
// stress it (i) tests admissibility, and otherwise (ii) selects the active
// return region from a-priori criteria computed solely from the trial state --
// the smooth face (sigma_1 > sigma_2 > sigma_3), the left edge (sigma_1 =
// sigma_2), the right edge (sigma_2 = sigma_3), or the apex (sigma_1 = sigma_2 =
// sigma_3) -- and returns the stress in closed form (no iteration: the plastic
// multiplier is affine for perfect plasticity). Coaxiality is exploited to
// reconstruct the in-plane tensor from the returned principal magnitudes. When a
// plastic correction is applied it also returns the closed-form consistent
// (algorithmic) tangent (Sysala Sec. 5): the active region's stress-update
// Jacobian carried through the eigen-decomposition, including the in-plane frame
// rotation (spin) and the out-of-plane (sigma_zz) coupling. The governing
// equations are catalogued in docs/references/mohr-coulomb-formulation.md.
//
// With params.tension_cutoff the admissible set becomes MC intersected with the
// Rankine cap sigma_1 <= sigma_t (three associated planes, PLAXIS MMM Eq 3-11).
// The return then runs a validity-cascaded region selection (formulation doc
// sec 7): pure-MC regions (unchanged, bit-identical when the cap is inactive),
// the tension face / tension-tension edge / tension apex, the MC-tension line
// and its two vertices. Every region is affine, so the region Jacobian stays
// constant-exact and the same principal-space tangent assembly applies.
McReturn mc_return_mapping(const PlaneStrainStress& trial,
                           const MohrCoulombParams& params);

} // namespace katai::core
