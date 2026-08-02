#include <katai/materials/mohr_coulomb.hpp>

#include <algorithm>
#include <array>
#include <cmath>

// Closed-form Mohr-Coulomb return mapping in principal stress space (P1.2b).
// Notation and equation numbers refer to Sysala & Cermak (2016, arXiv:1508.07435);
// the catalogue lives in docs/references/mohr-coulomb-formulation.md. Perfectly
// plastic (no hardening, H = 0): every region yields an affine equation for the
// plastic multiplier, so the update is fully closed-form.

namespace katai::core {
namespace {

// Provenance of a principal value: the larger / smaller in-plane principal, or
// the out-of-plane normal stress. Needed to scatter the returned principals back
// onto the (coaxial) trial eigenframe.
enum Source { kInPlaneA = 0, kInPlaneB = 1, kOutOfPlane = 2 };

struct PrincipalValue {
    double value;
    int source;
};

} // namespace

McReturn mc_return_mapping(const PlaneStrainStress& trial,
                           const MohrCoulombParams& mp) {
    const LameConstants lame = lame_from(mp.youngs_modulus, mp.poisson_ratio);
    const double G = lame.mu;
    const double K = lame.lambda + 2.0 * lame.mu / 3.0;
    const double sphi = std::sin(mp.friction_angle);
    const double cphi = std::cos(mp.friction_angle);
    const double spsi = std::sin(mp.dilatancy_angle);
    const double c = mp.cohesion;

    // --- Spectral decomposition of the trial stress (in-plane + sigma_zz) -----
    const double sxx = trial.in_plane(0), syy = trial.in_plane(1),
                 sxy = trial.in_plane(2);
    const double mean = 0.5 * (sxx + syy);
    const double half_diff = 0.5 * (sxx - syy);
    const double radius = std::sqrt(half_diff * half_diff + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;  // in-plane principal orientation
    if (radius > 0.0) {
        cos2t = half_diff / radius;
        sin2t = sxy / radius;
    }
    std::array<PrincipalValue, 3> pv = {{{mean + radius, kInPlaneA},
                                         {mean - radius, kInPlaneB},
                                         {trial.zz, kOutOfPlane}}};
    std::sort(pv.begin(), pv.end(),
              [](const PrincipalValue& a, const PrincipalValue& b) {
                  return a.value > b.value;
              });
    const double s1 = pv[0].value, s2 = pv[1].value, s3 = pv[2].value;

    // --- Admissibility (elastic step) ----------------------------------------
    const double f_trial = (1.0 + sphi) * s1 - (1.0 - sphi) * s3 - 2.0 * c * cphi;
    const double scale =
        (1.0 + sphi) * std::fabs(s1) + (1.0 - sphi) * std::fabs(s3) + 2.0 * c * cphi + 1.0;

    // Rankine tension cap (formulation doc sec 7): with sorted principals the three
    // PLAXIS planes (MMM Eq 3-11, associated flow) reduce to s1 <= st. st is clamped
    // to the MC apex bound c*cot(phi) -- the most tensile principal the MC set admits
    // -- so the cap is never redundant and the corner geometry below stays well-posed
    // (sm <= st). Off => every path below is the bit-identical pure-MC return.
    const bool cutoff = mp.tension_cutoff;
    double st = std::max(0.0, mp.tensile_strength);
    if (cutoff && sphi > 1.0e-14) st = std::min(st, c * cphi / sphi);
    const double tol_t = 1.0e-12 * (std::fabs(s1) + st + 1.0);
    const bool viol_mc = f_trial > 1.0e-12 * scale;
    const bool viol_t = cutoff && (s1 - st > tol_t);

    McReturn out;
    out.stress = trial;
    out.plastic = false;
    if (!viol_mc && !viol_t) return out;
    out.plastic = true;

    // --- Region selection and closed-form return -----------------------------
    // (2/3)(3K-2G) = 2*lambda; coefficients follow Sysala (4.8)-(4.27).
    const double a23 = (2.0 / 3.0) * (3.0 * K - 2.0 * G);
    const double a43 = (4.0 / 3.0) * (3.0 * K - 2.0 * G);

    const double g_sl = (s1 - s2) / (2.0 * G * (1.0 + spsi));  // smooth|left (4.13)
    const double g_sr = (s2 - s3) / (2.0 * G * (1.0 - spsi));  // smooth|right (4.13)

    // Smooth-face multiplier (4.12, H=0): numerator is exactly f_trial.
    const double den_s = a43 * spsi * sphi + 4.0 * G * (1.0 + spsi * sphi);
    const double dl_s = f_trial / den_s;

    double r1 = s1, r2 = s2, r3 = s3;  // returned principals

    // Jacobian of the returned principals w.r.t. the (sorted) trial principals,
    // J(i,j) = d(sigma_i^ret)/d(sigma_j^tr), evaluated for the *active* region.
    // Within a region the closed-form update is affine in (s1,s2,s3) (perfect
    // plasticity), so J is constant and exact -- this is precisely what a forward
    // finite difference fails to recover near a region boundary (the perturbation
    // straddles into a neighbouring branch). It is the kernel of the consistent
    // tangent assembled below.
    Eigen::Matrix3d J = Eigen::Matrix3d::Identity();
    const Eigen::RowVector3d e0(1.0, 0.0, 0.0), e2(0.0, 0.0, 1.0);

    auto return_smooth = [&] {  // sigma_1 > sigma_2 > sigma_3, (4.9)-(4.11)
        const Eigen::Vector3d cf(a23 * spsi + 2.0 * G * (1.0 + spsi),
                                 a23 * spsi,
                                 a23 * spsi - 2.0 * G * (1.0 - spsi));
        r1 = s1 - dl_s * cf(0);
        r2 = s2 - dl_s * cf(1);
        r3 = s3 - dl_s * cf(2);
        // dl_s = f_trial/den_s, with d(f_trial)/ds = (1+sphi, 0, -(1-sphi)).
        const Eigen::RowVector3d df((1.0 + sphi) / den_s, 0.0,
                                    -(1.0 - sphi) / den_s);
        J = Eigen::Matrix3d::Identity() - cf * df;
    };
    auto return_apex = [&] {  // sigma_1 = sigma_2 = sigma_3 = p, (4.24)-(4.27)
        const double sum = s1 + s2 + s3;
        const double den_a = 4.0 * K * spsi * sphi;
        double p;
        if (std::fabs(den_a) > 1.0e-14 * (K + 1.0) && sphi > 1.0e-14) {
            const double dl = ((2.0 / 3.0) * sum * sphi - 2.0 * c * cphi) / den_a;
            p = sum / 3.0 - 2.0 * K * spsi * dl;
            // p collapses to the fixed apex point c*cot(phi), independent of the
            // trial -- the vertex carries no stiffness: J = 0.
            J.setZero();
        } else {
            // Degenerate (psi = 0 or phi = 0): the deviatoric flow cannot reach
            // the tip, so project directly onto the apex point p = c*cot(phi).
            if (sphi > 1.0e-14) {
                p = c * cphi / sphi;  // constant projection: J = 0
                J.setZero();
            } else {
                p = sum / 3.0;  // pure volumetric return: J = (1/3) 11^T
                J = Eigen::Matrix3d::Constant(1.0 / 3.0);
            }
        }
        r1 = r2 = r3 = p;
    };

    // Only two orderings of the region boundaries are possible (Lemma 4.2):
    // g_sl <= ... (left-edge branch) or g_sr <= ... (right-edge branch).
    bool done = false;
    if (viol_mc) {
        if (g_sl <= g_sr) {
            if (dl_s <= g_sl) {
                return_smooth();
            } else {
                // Left edge sigma_1 = sigma_2 (4.14)-(4.18), else apex.
                const double num = 0.5 * (1.0 + sphi) * (s1 + s2) -
                                   (1.0 - sphi) * s3 - 2.0 * c * cphi;
                const double den = a43 * spsi * sphi +
                                   G * (1.0 + spsi) * (1.0 + sphi) +
                                   2.0 * G * (1.0 - spsi) * (1.0 - sphi);
                const double dl = num / den;
                const double g_la = (s1 + s2 - 2.0 * s3) / (2.0 * G * (3.0 - spsi));
                if (dl < g_la) {
                    const double cl = a23 * spsi + G * (1.0 + spsi);    // edge coeff
                    const double c3 = a23 * spsi - 2.0 * G * (1.0 - spsi);
                    const double edge = 0.5 * (s1 + s2) - dl * cl;
                    r1 = r2 = edge;
                    r3 = s3 - dl * c3;
                    // d(num)/ds = (1/2(1+sphi), 1/2(1+sphi), -(1-sphi)).
                    const Eigen::RowVector3d df(0.5 * (1.0 + sphi) / den,
                                                0.5 * (1.0 + sphi) / den,
                                                -(1.0 - sphi) / den);
                    J.row(0) = Eigen::RowVector3d(0.5, 0.5, 0.0) - cl * df;
                    J.row(1) = J.row(0);
                    J.row(2) = e2 - c3 * df;
                } else {
                    return_apex();
                }
            }
        } else {
            if (dl_s <= g_sr) {
                return_smooth();
            } else {
                // Right edge sigma_2 = sigma_3 (4.19)-(4.23), else apex.
                const double num = (1.0 + sphi) * s1 -
                                   0.5 * (1.0 - sphi) * (s2 + s3) - 2.0 * c * cphi;
                const double den = a43 * spsi * sphi +
                                   2.0 * G * (1.0 + spsi) * (1.0 + sphi) +
                                   G * (1.0 - spsi) * (1.0 - sphi);
                const double dl = num / den;
                const double g_ra = (2.0 * s1 - s2 - s3) / (2.0 * G * (3.0 + spsi));
                if (dl < g_ra) {
                    const double c1 = a23 * spsi + 2.0 * G * (1.0 + spsi);
                    const double ce = a23 * spsi - G * (1.0 - spsi);    // edge coeff
                    r1 = s1 - dl * c1;
                    const double edge = 0.5 * (s2 + s3) - dl * ce;
                    r2 = r3 = edge;
                    // d(num)/ds = (1+sphi, -1/2(1-sphi), -1/2(1-sphi)).
                    const Eigen::RowVector3d df((1.0 + sphi) / den,
                                                -0.5 * (1.0 - sphi) / den,
                                                -0.5 * (1.0 - sphi) / den);
                    J.row(0) = e0 - c1 * df;
                    J.row(1) = Eigen::RowVector3d(0.0, 0.5, 0.5) - ce * df;
                    J.row(2) = J.row(1);
                } else {
                    return_apex();
                }
            }
        }
        // The MC-only result stands unless the cap is on and the returned major
        // principal still exceeds it (KKT: the cap is then inactive at the solution
        // even if the TRIAL violated it -- no escalation needed).
        done = !cutoff || (r1 <= st + tol_t);
    }

    if (!done) {
        // --- Tension cut-off regions (formulation doc sec 7) -------------------
        // Reached with s1 > st guaranteed: either the trial violated only the cap,
        // or the MC-only return (which never raises r1 above s1) still exceeds it.
        out.tension = true;
        const double lam = lame.lambda;
        const double tol_mc = 1.0e-12 * scale;
        auto fmc = [&](double sa, double sb) {
            return (1.0 + sphi) * sa - (1.0 - sphi) * sb - 2.0 * c * cphi;
        };

        // T-only candidate: the associated Rankine return to the box s_i <= st
        // (7a face -> 7b edge -> 7c apex). Each failed ordering check is EXACTLY
        // the next region's positive-multiplier condition (duality checked in the
        // formulation doc), so this cascade is the exact box return.
        double q1, q2, q3;
        Eigen::Matrix3d Jt = Eigen::Matrix3d::Zero();
        const double dlt = (s1 - st) / (lam + 2.0 * G);
        const double face2 = s2 - lam * dlt;
        if (face2 <= st + tol_t) {                       // 7a: tension face
            q1 = st; q2 = face2; q3 = s3 - lam * dlt;
            const double k = -lam / (lam + 2.0 * G);
            Jt(1, 0) = k; Jt(1, 1) = 1.0;
            Jt(2, 0) = k; Jt(2, 2) = 1.0;
        } else {
            const double S = (s1 + s2 - 2.0 * st) / (2.0 * lam + 2.0 * G);
            const double edge3 = s3 - lam * S;
            if (edge3 <= st + tol_t) {                   // 7b: tension-tension edge
                q1 = st; q2 = st; q3 = edge3;
                const double k = -lam / (2.0 * lam + 2.0 * G);
                Jt(2, 0) = k; Jt(2, 1) = k; Jt(2, 2) = 1.0;
            } else {                                     // 7c: tension apex (J = 0)
                q1 = st; q2 = st; q3 = st;
            }
        }

        if (fmc(q1, q3) <= tol_mc) {
            // Pure-tension region: MC inactive at the returned point. (The apex
            // candidate always lands here: fmc(st,st) <= 0 by the c*cot(phi) clamp.)
            r1 = q1; r2 = q2; r3 = q3;
            J = Jt;
        } else {
            // Both surfaces active: return to the MC-tension line r = (st, r2, sm)
            // running from V_R = (st, sm, sm) to V_L = (st, st, sm), doc 7d/7e.
            // Multipliers a (MC smooth flow, non-associated) and b (tension,
            // associated) solve an exact 2x2 system; the coefficients reuse the
            // smooth-region constants (A11 = den_s, A21 = cf(0)) so the b = 0 and
            // a = 0 limits coincide with the neighbouring regions to the ulp.
            // det = 4G(lam+G)(1-sphi)(1-spsi) > 0 for phi, psi < 90 deg.
            const double A11 = den_s;
            const double A12 = a23 * sphi + 2.0 * G * (1.0 + sphi);
            const double A21 = a23 * spsi + 2.0 * G * (1.0 + spsi);
            const double A22 = lam + 2.0 * G;
            const double det = A11 * A22 - A12 * A21;
            const double a = (A22 * f_trial - A12 * (s1 - st)) / det;
            const double b = (A11 * (s1 - st) - A21 * f_trial) / det;
            const double sm = ((1.0 + sphi) * st - 2.0 * c * cphi) / (1.0 - sphi);
            const double p2 = s2 - a * a23 * spsi - b * lam;
            if (p2 > st + tol_t) {           // past V_L: s2 hits the cap as well
                r1 = st; r2 = st; r3 = sm;
                J.setZero();
            } else if (p2 < sm - tol_t) {    // past V_R: s2 joins the MC right edge
                r1 = st; r2 = sm; r3 = sm;
                J.setZero();
            } else {                          // line interior: only r2 varies
                r1 = st; r2 = p2; r3 = sm;
                const Eigen::RowVector3d df(1.0 + sphi, 0.0, -(1.0 - sphi));
                const Eigen::RowVector3d e1r(1.0, 0.0, 0.0);
                const Eigen::RowVector3d da = (A22 * df - A12 * e1r) / det;
                const Eigen::RowVector3d db = (A11 * e1r - A21 * df) / det;
                J.setZero();
                J.row(1) =
                    Eigen::RowVector3d(0.0, 1.0, 0.0) - a23 * spsi * da - lam * db;
            }
        }
    }

    // --- Reconstruct the in-plane tensor on the (coaxial) trial eigenframe -----
    const double returned[3] = {r1, r2, r3};
    double pa = 0.0, pb = 0.0, pz = 0.0;
    int slot[3] = {0, 0, 0};  // slot[source] = position in the sorted triplet
    for (int i = 0; i < 3; ++i) {
        slot[pv[i].source] = i;
        switch (pv[i].source) {
            case kInPlaneA: pa = returned[i]; break;
            case kInPlaneB: pb = returned[i]; break;
            case kOutOfPlane: pz = returned[i]; break;
        }
    }
    const double m = 0.5 * (pa + pb);
    const double r = 0.5 * (pa - pb);
    out.stress.in_plane(0) = m + r * cos2t;
    out.stress.in_plane(1) = m - r * cos2t;
    out.stress.in_plane(2) = r * sin2t;
    out.stress.zz = pz;

    // --- Consistent tangent: shared principal-space assembly (spin + szz coupling).
    // J(i,j)=d(returned sorted i)/d(trial sorted j) is the region map (affine, exact);
    // chain eps -> trial principals -> returned principals -> reconstructed sigma is
    // written once in principal_consistent_tangent (also used by Hardening Soil).
    const int src[3] = {pv[0].source, pv[1].source, pv[2].source};
    const double ret[3] = {r1, r2, r3};
    const PrincipalTangent pt =
        principal_consistent_tangent(cos2t, sin2t, radius, src, ret, J, lame);
    out.tangent = pt.tangent;
    out.algo_jacobian = pt.algo_jacobian;
    return out;
}

} // namespace katai::core
