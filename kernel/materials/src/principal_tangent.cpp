#include <katai/materials/mohr_coulomb.hpp>

namespace katai::core {

PrincipalTangent principal_consistent_tangent(
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

}  // namespace katai::core
