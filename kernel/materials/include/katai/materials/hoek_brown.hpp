#pragma once
// HOEK-BROWN (rock behaviour) -- the material-point core. PLAXIS 2D Material Models Manual §4
// (the 2002 edition of the criterion: Hoek, Carranza-Torres & Corkum 2002; the FE implementation
// including strength factorisation follows Benz, Schwab, Vermeer & Kauther 2007).
//
// WHY A ROCK MODEL IS NOT A SOIL MODEL WITH DIFFERENT NUMBERS. Rock is stiff enough that the
// stress dependency of its stiffness is negligible, so Hooke's law is kept -- but its STRENGTH is
// strongly non-linear in the confining stress, and it has real tensile strength. A Mohr-Coulomb
// fit is a straight line through a curve: it can be made to match over a narrow band of confining
// stress and is wrong outside it, in both directions. That is the whole reason this model exists.
//
// THE CRITERION (manual Eq 4-1, tension POSITIVE and compression negative -- the same convention
// this solver uses throughout):
//     sigma'_1 = sigma'_3 - |sigma_ci| ( m_b sigma'_3 / sigma_ci + s )^a
// with the manual's ordering sigma'_1 <= sigma'_2 <= sigma'_3, i.e. its sigma'_1 is the MOST
// COMPRESSIVE. THIS HEADER USES THE SOLVER'S ORDERING INSTEAD -- s1 >= s2 >= s3, s1 the most
// TENSILE (mohr_coulomb.hpp says so and every caller sorts that way) -- so the manual's
// (sigma'_1, sigma'_3) are this file's (s3, s1). Getting that mapping backwards produces a model
// that runs, converges and is wrong in the direction that matters; it is stated here, asserted in
// the tests against the closed-form envelope, and never re-derived at a call site.
//
// Yield (manual Eq 4-8 in this file's ordering):
//     f = s3 - s1 + fbar(s1),    fbar(x) = |sigma_ci| ( m_b (-x / |sigma_ci|) + s )^a
//
// Flow (manual Eq 4-10/4-12), non-associated, through the transformed stresses
//     S_i = -sigma_i / (m_b |sigma_ci|) + s / m_b^2
//     g = S(s3) - N_psi S(s1),   N_psi = (1 + sin psi_mob) / (1 - sin psi_mob)
// S is LINEAR in sigma, so g is linear in sigma and the plastic corrector direction is CONSTANT:
//     dg/dsigma  proportional to  [ N_psi, 0, -1 ]   in (s1, s2, s3)
// which is the Mohr-Coulomb direction with the mobilised dilatancy. The consequence is worth
// stating because it shapes the algorithm: the return is a STRAIGHT LINE in principal stress
// space and the only non-linearity left is the yield surface itself, so the corrector is a 1-D
// Newton on the plastic multiplier rather than a general closest-point projection.
//
// SCOPE OF THIS VERSION, stated rather than discovered: the return covers the MAIN surface
// (f_13), the EDGE where the manual's second function f_12 becomes active, and the tensile APEX
// where the criterion's own bracket closes. What is NOT here: strength factorisation for Safety
// phases (refused at both seams, with the reason), and the anisotropic Jointed Rock model -- the
// manual keeps that one apart too.
//
// The TANGENT the caller uses is not written here either, but it is no longer the elastic
// operator: material_model.hpp differentiates this return by finite difference when the step is
// plastic, the way the Hardening Soil and Soft Soil branches differentiate theirs. That changed
// after the elastic operator was measured against a boundary value problem rather than a stress
// point -- a tunnel unloaded into a rock mass -- where it stalled the equilibrium iteration as
// soon as the plastic annulus passed a few percent of the opening radius. The return itself is
// cheap enough to differentiate (a bisection on one scalar), which is why the FD is affordable
// here where it is a real cost for the substepped models.

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

#include <katai/materials/mohr_coulomb.hpp>   // LameConstants / PlaneStrainStress / lame_from

namespace katai::core::hoekbrown {

// The eight parameters of the manual's §4.3, in its own units and signs.
struct Params {
    double E = 0.0;        // rock mass Young's modulus E_rm [kPa]
    double nu = 0.0;       // Poisson's ratio [-]
    double sigci = 0.0;    // |sigma_ci|, uni-axial compressive strength of the INTACT rock, > 0
    double mi = 0.0;       // intact rock parameter m_i [-]
    double gsi = 0.0;      // Geological Strength Index [-] (0..100)
    double D = 0.0;        // disturbance factor [-] (0..1)
    double psi = 0.0;      // dilatancy angle at sigma'_3 = 0 [rad]
    double sig_psi = 0.0;  // confining stress at which the dilatancy has died out [kPa, >= 0]
    // OPTIONAL TENSION CUT-OFF (manual sec 4.3.7). The criterion already has a tensile strength of
    // its own -- sigma_t of Eq 4-6, where the bracket closes -- and that one is a CONSEQUENCE of
    // sigma_ci, s and m_b rather than an input. The manual lets the user cap it lower: "users can
    // put a tensile strength value, and if that value is lower than sigma_t, the tensile capacity
    // will be cut-off at that value". Lower only: a value above sigma_t would claim a tensile
    // strength the criterion does not have, and is ignored (the min below).
    bool tension_cutoff = false;
    double sigt_user = 0.0;   // [kPa, >= 0] magnitude; read only when tension_cutoff
};

// The rock-mass constants derived from GSI and D (manual Eq 4-2, 4-3, 4-4).
struct Constants {
    double mb = 0.0, s = 0.0, a = 0.5;
    double sigc = 0.0;   // uni-axial compressive strength of the ROCK MASS (Eq 4-5), negative
    double sigt = 0.0;   // tensile strength of the rock mass (Eq 4-6), positive
};

inline Constants constants_of(const Params& P) {
    Constants C;
    C.mb = P.mi * std::exp((P.gsi - 100.0) / (28.0 - 14.0 * P.D));           // Eq 4-2
    C.s = std::exp((P.gsi - 100.0) / (9.0 - 3.0 * P.D));                      // Eq 4-3
    C.a = 0.5 + (std::exp(-P.gsi / 15.0) - std::exp(-20.0 / 3.0)) / 6.0;      // Eq 4-4
    C.sigc = -std::fabs(P.sigci) * std::pow(C.s, C.a);                        // Eq 4-5 (compression)
    C.sigt = C.mb > 0.0 ? C.s * std::fabs(P.sigci) / C.mb : 0.0;              // Eq 4-6 (tension)
    // The user's cut-off, sec 4.3.7. Everything downstream -- the yield functions' second branch,
    // the apex return region, the tensile branch of the mobilised dilatancy -- is written in terms
    // of C.sigt, so capping it HERE is the whole implementation and no path can miss it.
    if (P.tension_cutoff) C.sigt = std::min(C.sigt, std::fabs(P.sigt_user));
    return C;
}

// fbar(x) of Eq 4-8, with x the MOST TENSILE principal stress (the manual's sigma'_3). The
// bracket is the criterion's own domain: it reaches zero exactly at the isotropic tensile apex
// x = sigma_t, and there is no surface beyond that -- which is why the apex is a return region
// and not a special case bolted on.
inline double fbar(double x, const Params& P, const Constants& C) {
    const double sci = std::fabs(P.sigci);
    const double bracket = C.mb * (-x / sci) + C.s;
    if (bracket <= 0.0) return 0.0;
    return sci * std::pow(bracket, C.a);
}
inline double dfbar_dx(double x, const Params& P, const Constants& C) {
    const double sci = std::fabs(P.sigci);
    const double bracket = C.mb * (-x / sci) + C.s;
    if (bracket <= 0.0) return 0.0;
    return -C.a * C.mb * std::pow(bracket, C.a - 1.0);
}

// The yield function in the SOLVER's ordering (s1 >= s2 >= s3, tension positive) -- AND WITH THE
// SIGN THE REST OF THIS SOLVER USES, which is the opposite of the manual's. Eq 4-7 is written so
// that an unstressed point gives f = +|sigma_ci| s^a > 0: its admissible set is f >= 0, and its
// yielding is f <= 0. Every other yield function in this tree is the other way round (f <= 0
// admissible), so the sign is flipped ONCE, here, rather than at each call site where it would
// eventually be flipped only in most of them:
//     f = s1 - s3 - fbar(s1)      <= 0 admissible, = 0 on the envelope
// A test that only checked "the envelope point gives f = 0" would pass either way; the one that
// caught this asked whether a MORE loaded state is inadmissible, which is the question that has
// a side.
inline double yield(double s1, double s3, const Params& P, const Constants& C) {
    // TWO CONDITIONS, NOT ONE. Past the isotropic tensile apex the criterion's bracket has closed
    // and fbar is zero, so the first term alone would read f = s1 - s3 = 0 for a state in
    // isotropic tension well beyond the rock's tensile strength -- admissible, which it is not.
    // The apex is not an artefact of the algebra: Eq 4-6 says the surface ENDS at s1 = sigma_t,
    // so the admissible set is the criterion INTERSECTED with s1 <= sigma_t, and the yield
    // function of an intersection is the larger of the two.
    return std::max(s1 - s3 - fbar(s1, P, C), s1 - C.sigt);
}

// The second yield function (manual Eq 4-9). In the manual's ordering it pairs sigma'_1 with the
// INTERMEDIATE stress; in this file's that is s3 with s2, so it becomes active exactly where a
// return on the main surface would push s2 past s1. (The other corner, s2 = s3, is not a corner
// of this criterion at all: there Eq 4-9 evaluates to -fbar < 0, comfortably admissible.)
inline double yield12(double s2, double s3, const Params& P, const Constants& C) {
    return std::max(s2 - s3 - fbar(s2, P, C), s2 - C.sigt);
}

// Mobilised dilatancy (manual Eq 4-13/4-14), a function of the most tensile principal stress --
// the manual's sigma'_3, this file's s1. Compression negative, so "more confinement" is s1 more
// negative. Below -sigma_psi the rock has stopped dilating; in the tensile range the manual uses
// an artificially increased value so that plastic expansion remains possible there.
inline double psi_mobilised(double s1, const Params& P, const Constants& C) {
    constexpr double kHalfPi = 1.57079632679489661923;   // 90 degrees, the manual's cap
    if (s1 >= 0.0) {                                      // tensile zone, Eq 4-14
        if (C.sigt <= 0.0) return P.psi;
        const double t = std::min(1.0, s1 / C.sigt);
        return P.psi + t * (kHalfPi - P.psi);
    }
    if (P.sig_psi <= 0.0) return P.psi;                   // no decay declared
    const double sp = P.sig_psi;                          // Eq 4-13: (sigma_psi + sigma'_3)/sigma_psi
    return std::max(0.0, (sp + s1) / sp * P.psi);
}

// The result of a return, in principal space.
struct Return {
    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    bool plastic = false;
    bool apex = false;     // returned to the isotropic tensile apex
    bool edge = false;     // the second yield function was active (ordering would have broken)
};

// Return the trial principal stresses to the criterion. `sorted` must satisfy s1 >= s2 >= s3.
//
// The corrector is sigma = sigma_tr - lambda * D_e * m with m = [N_psi, 0, -1] constant, so
// f(lambda) is a scalar function and the return is a 1-D Newton -- with a bisection fallback,
// because f is smooth but its derivative vanishes at the apex where the bracket closes and a bare
// Newton would step past it.
inline Return return_mapping(double s1t, double s2t, double s3t, const Params& P,
                             const Constants& C) {
    Return R{s1t, s2t, s3t, false, false, false};
    if (yield(s1t, s3t, P, C) <= 0.0) return R;   // elastic

    const LameConstants lame = lame_from(P.E, P.nu);
    const double lam = lame.lambda, mu = lame.mu;
    const double npsi = [&] {
        const double pm = psi_mobilised(s1t, P, C);
        const double sp = std::sin(pm);
        return (1.0 + sp) / (1.0 - sp);
    }();
    // D_e * m for m = [N_psi, 0, -1] (isotropic elasticity in principal space).
    const double d1 = (lam + 2.0 * mu) * npsi - lam;
    const double d2 = lam * (npsi - 1.0);
    const double d3 = lam * npsi - (lam + 2.0 * mu);

    // THE APEX IS A DOMAIN BOUNDARY, NOT A SPECIAL CASE. The criterion's bracket closes exactly at
    // s1 = sigma_t (that is what Eq 4-6 says), so there is no surface on the tensile side of the
    // isotropic tension point. A trial that is already past it cannot be returned to a surface
    // that does not exist, and goes to the apex itself.
    if (s1t >= C.sigt) {
        R.s1 = R.s2 = R.s3 = C.sigt;
        R.plastic = true; R.apex = true;
        return R;
    }

    const auto f_of = [&](double L) {
        return yield(s1t - L * d1, s3t - L * d3, P, C);
    };
    // The corrector drives s1 down and s3 up, so f falls monotonically and without bound: a
    // doubling search finds a bracket on any input, and bisection cannot leave it. (Newton would
    // be faster and would step past the apex where the bracket's derivative vanishes; this is one
    // of the places where the slower method is the one that always answers.)
    double lo = 0.0, hi = 1.0;
    const double scale = std::fabs(s1t) + std::fabs(s3t) + std::fabs(C.sigc) + 1.0;
    hi = scale / std::max(1.0, d1 - d3);
    for (int i = 0; i < 200 && f_of(hi) > 0.0; ++i) hi *= 2.0;
    for (int i = 0; i < 80; ++i) {          // bisection: monotone in lambda, and it cannot escape
        const double mid = 0.5 * (lo + hi);
        (f_of(mid) > 0.0 ? lo : hi) = mid;
    }
    const double L = 0.5 * (lo + hi);
    R.s1 = s1t - L * d1;
    R.s2 = s2t - L * d2;
    R.s3 = s3t - L * d3;
    R.plastic = true;

    // THE ORDER IS PART OF THE ANSWER. The corrector moves s2 as well, so a state that started
    // ordered can leave the region where f_13 alone governs -- which is exactly where the manual's
    // second yield function takes over (Eq 4-9). Clamping s2 back would report a stress the
    // criterion never admitted: measured against this tree's own Mohr-Coulomb return in the Tresca
    // degeneration, the clamp was 25 kPa out on a material whose cohesion is 50. So the edge is
    // returned to properly: BOTH surfaces active, two multipliers, two conditions.
    //
    //     sigma = sigma_tr - la * D m_a - lb * D m_b,   m_a = [N,0,-1],  m_b = [0,N,-1]
    //     f_13(sigma) = 0  and  f_12(sigma) = 0
    //
    // The directions are constant, so this is a 2x2 system in (la, lb) alone. It is solved by a
    // damped Newton on a numerical Jacobian -- two unknowns, smooth residuals, and the main-surface
    // solution as a starting point, which is inside the basin by construction.
    if (R.s2 > R.s1 + 1e-12 || R.s2 < R.s3 - 1e-12) {
        R.edge = true;
        const double e1 = lam * npsi - lam;                    // D m_b, component 1
        const double e2 = (lam + 2.0 * mu) * npsi - lam;       // component 2
        const double e3 = lam * npsi - (lam + 2.0 * mu);       // component 3
        double la = L, lb = 0.0;
        const auto residual = [&](double a, double b, double& r1, double& r2) {
            const double x1 = s1t - a * d1 - b * e1;
            const double x2 = s2t - a * d2 - b * e2;
            const double x3 = s3t - a * d3 - b * e3;
            r1 = yield(x1, x3, P, C);
            r2 = yield12(x2, x3, P, C);
        };
        double r1 = 0.0, r2 = 0.0;
        residual(la, lb, r1, r2);
        const double h = 1e-7 * (std::fabs(la) + std::fabs(lb) + 1e-12);
        for (int it = 0; it < 40; ++it) {
            if (std::fabs(r1) + std::fabs(r2) < 1e-10 * (std::fabs(C.sigc) + 1.0)) break;
            double a1, a2, b1, b2;
            residual(la + h, lb, a1, a2);
            residual(la, lb + h, b1, b2);
            const double J11 = (a1 - r1) / h, J12 = (b1 - r1) / h;
            const double J21 = (a2 - r2) / h, J22 = (b2 - r2) / h;
            const double det = J11 * J22 - J12 * J21;
            if (std::fabs(det) < 1e-300) break;
            const double da = (-r1 * J22 + r2 * J12) / det;
            const double db = (-r2 * J11 + r1 * J21) / det;
            double step = 1.0;
            for (int k = 0; k < 20; ++k) {                     // damping: never leave la, lb >= 0
                const double na = la + step * da, nb = lb + step * db;
                if (na >= 0.0 && nb >= 0.0) {
                    double q1, q2;
                    residual(na, nb, q1, q2);
                    if (std::fabs(q1) + std::fabs(q2) < std::fabs(r1) + std::fabs(r2)) {
                        la = na; lb = nb; r1 = q1; r2 = q2;
                        break;
                    }
                }
                step *= 0.5;
            }
        }
        R.s1 = s1t - la * d1 - lb * e1;
        R.s2 = s2t - la * d2 - lb * e2;
        R.s3 = s3t - la * d3 - lb * e3;
    }
    return R;
}


// ------------------------------------------------------------------------------------------
// The FE-facing entry: a plane-strain stress in, an admissible plane-strain stress out.
//
// The decomposition and the reconstruction are the same ones the Mohr-Coulomb return uses --
// sigma_zz is itself a principal stress because tau_xz = tau_yz = 0 in plane strain, so the
// triplet is {the two in-plane principals, sigma_zz} and the returned triplet goes back onto the
// SAME eigenframe (the return is coaxial: the corrector is a combination of principal directions,
// so it cannot rotate them). What is deliberately NOT here is the algorithmic tangent: the caller
// uses the elastic operator, which is what the Hardening Soil branch does for its own reasons --
// a robust tangent, and one that costs iterations rather than accuracy.
struct PlaneReturn {
    PlaneStrainStress stress;
    bool plastic = false, apex = false, edge = false;
};

inline PlaneReturn plane_return(const PlaneStrainStress& trial, const Params& P,
                                const Constants& C) {
    enum { kA = 0, kB = 1, kZ = 2 };
    const double sxx = trial.in_plane(0), syy = trial.in_plane(1), sxy = trial.in_plane(2);
    const double mean = 0.5 * (sxx + syy), half = 0.5 * (sxx - syy);
    const double radius = std::sqrt(half * half + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = half / radius; sin2t = sxy / radius; }

    struct PV { double v; int src; };
    PV pv[3] = {{mean + radius, kA}, {mean - radius, kB}, {trial.zz, kZ}};
    std::sort(pv, pv + 3, [](const PV& x, const PV& y) { return x.v > y.v; });

    const Return R = return_mapping(pv[0].v, pv[1].v, pv[2].v, P, C);
    PlaneReturn out;
    out.plastic = R.plastic; out.apex = R.apex; out.edge = R.edge;
    if (!R.plastic) { out.stress = trial; return out; }

    const double ret[3] = {R.s1, R.s2, R.s3};
    double pa = 0.0, pb = 0.0, pz = 0.0;
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == kA) pa = ret[i];
        else if (pv[i].src == kB) pb = ret[i];
        else pz = ret[i];
    }
    const double m = 0.5 * (pa + pb), r = 0.5 * (pa - pb);
    out.stress.in_plane(0) = m + r * cos2t;
    out.stress.in_plane(1) = m - r * cos2t;
    out.stress.in_plane(2) = r * sin2t;
    out.stress.zz = pz;
    return out;
}

} // namespace katai::core::hoekbrown
