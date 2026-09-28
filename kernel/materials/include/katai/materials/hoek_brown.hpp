#pragma once
// HOEK-BROWN (rock behaviour) -- the material-point core: the 2002 edition of the criterion
// (Hoek, Carranza-Torres & Corkum 2002). A strength-factorised form of the same yield function
// exists (Benz, Schwab, Kauther & Vermeer 2008, Int. J. Rock Mech. Min. Sci. 45(2), 210-222) and is not implemented in this build.
//
// WHY A ROCK MODEL IS NOT A SOIL MODEL WITH DIFFERENT NUMBERS. Rock is stiff enough that the
// stress dependency of its stiffness is negligible, so Hooke's law is kept -- but its STRENGTH is
// strongly non-linear in the confining stress, and it has real tensile strength. A Mohr-Coulomb
// fit is a straight line through a curve: it can be made to match over a narrow band of confining
// stress and is wrong outside it, in both directions. That is the whole reason this model exists.
//
// THE CRITERION (written with tension POSITIVE and compression negative -- the same convention
// this solver uses throughout):
//     sigma'_1 = sigma'_3 - |sigma_ci| ( m_b sigma'_3 / sigma_ci + s )^a
// with the ordering sigma'_1 <= sigma'_2 <= sigma'_3, i.e. its sigma'_1 is the MOST
// COMPRESSIVE. THIS HEADER USES THE SOLVER'S ORDERING INSTEAD -- s1 >= s2 >= s3, s1 the most
// TENSILE (mohr_coulomb.hpp says so and every caller sorts that way) -- so the criterion's
// (sigma'_1, sigma'_3) are this file's (s3, s1). Getting that mapping backwards produces a model
// that runs, converges and is wrong in the direction that matters; it is stated here, asserted in
// the tests against the closed-form envelope, and never re-derived at a call site.
//
// Yield (in this file's ordering):
//     f = s3 - s1 + fbar(s1),    fbar(x) = |sigma_ci| ( m_b (-x / |sigma_ci|) + s )^a
//
// Flow, non-associated, through the transformed stresses
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
// (f_13), the EDGE where the second yield function f_12 becomes active, and the tensile APEX
// where the criterion's own bracket closes. What is NOT here: strength factorisation for Safety
// phases (refused at both seams, with the reason), and anisotropic jointed-rock behaviour, which
// is a separate model.
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

// The model's eight input parameters, in the units and signs rock mechanics quotes them in.
struct Params {
    double E = 0.0;        // rock mass Young's modulus E_rm [kPa]
    double nu = 0.0;       // Poisson's ratio [-]
    double sigci = 0.0;    // |sigma_ci|, uni-axial compressive strength of the INTACT rock, > 0
    double mi = 0.0;       // intact rock parameter m_i [-]
    double gsi = 0.0;      // Geological Strength Index [-] (0..100)
    double D = 0.0;        // disturbance factor [-] (0..1)
    double psi = 0.0;      // dilatancy angle at sigma'_3 = 0 [rad]
    double sig_psi = 0.0;  // confining stress at which the dilatancy has died out [kPa, >= 0]
    // OPTIONAL TENSION CUT-OFF. The criterion already has a tensile strength of its own --
    // sigma_t = s sigma_ci / m_b, where the bracket closes -- and that one is a CONSEQUENCE of
    // sigma_ci, s and m_b rather than an input. The user may cap it lower: a tensile strength
    // below sigma_t cuts the tensile capacity off at that value. Lower only: a value above
    // sigma_t would claim a tensile strength the criterion does not have, and is ignored (the
    // min below).
    bool tension_cutoff = false;
    double sigt_user = 0.0;   // [kPa, >= 0] magnitude; read only when tension_cutoff
};

// The rock-mass constants derived from GSI and D (Hoek, Carranza-Torres & Corkum 2002).
struct Constants {
    double mb = 0.0, s = 0.0, a = 0.5;
    double sigc = 0.0;   // uni-axial compressive strength of the ROCK MASS, sigma_ci s^a, < 0
    double sigt = 0.0;   // tensile strength of the rock mass, s sigma_ci / m_b, positive
};

inline Constants constants_of(const Params& P) {
    Constants C;
    C.mb = P.mi * std::exp((P.gsi - 100.0) / (28.0 - 14.0 * P.D));           // m_b
    C.s = std::exp((P.gsi - 100.0) / (9.0 - 3.0 * P.D));                      // s
    C.a = 0.5 + (std::exp(-P.gsi / 15.0) - std::exp(-20.0 / 3.0)) / 6.0;      // a
    C.sigc = -std::fabs(P.sigci) * std::pow(C.s, C.a);                        // sigma_c (< 0)
    C.sigt = C.mb > 0.0 ? C.s * std::fabs(P.sigci) / C.mb : 0.0;              // sigma_t (tension)
    // The user's optional cut-off. Everything downstream -- the yield functions' second branch,
    // the apex return region, the tensile branch of the mobilised dilatancy -- is written in terms
    // of C.sigt, so capping it HERE is the whole implementation and no path can miss it.
    if (P.tension_cutoff) C.sigt = std::min(C.sigt, std::fabs(P.sigt_user));
    return C;
}

// fbar(x) of the yield function, with x the MOST TENSILE principal stress (the criterion's
// sigma'_3). The bracket is the criterion's own domain: it reaches zero exactly at the isotropic
// tensile apex x = sigma_t, and there is no surface beyond that -- which is why the apex is a
// return region and not a special case bolted on.
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
// SIGN THE REST OF THIS SOLVER USES, which is the opposite of the criterion's usual yield form.
// Written as f = sigma'_1 - sigma'_3 + fbar(sigma'_3) in the criterion's ordering, an unstressed
// point gives f = +|sigma_ci| s^a > 0: its admissible set is f >= 0, and its yielding is f <= 0.
// Every other yield function in this tree is the other way round (f <= 0 admissible), so the
// sign is flipped ONCE, here, rather than at each call site where it would eventually be
// flipped only in most of them:
//     f = s1 - s3 - fbar(s1)      <= 0 admissible, = 0 on the envelope
// A test that only checked "the envelope point gives f = 0" would pass either way; the one that
// caught this asked whether a MORE loaded state is inadmissible, which is the question that has
// a side.
inline double yield(double s1, double s3, const Params& P, const Constants& C) {
    // TWO CONDITIONS, NOT ONE. Past the isotropic tensile apex the criterion's bracket has closed
    // and fbar is zero, so the first term alone would read f = s1 - s3 = 0 for a state in
    // isotropic tension well beyond the rock's tensile strength -- admissible, which it is not.
    // The apex is not an artefact of the algebra: the bracket closes at s1 = sigma_t, the
    // surface ENDS there, so the admissible set is the criterion INTERSECTED with
    // s1 <= sigma_t, and the yield
    // function of an intersection is the larger of the two.
    return std::max(s1 - s3 - fbar(s1, P, C), s1 - C.sigt);
}

// The second yield function, f_12 = s2 - s3 - fbar(s2). In the criterion's ordering it pairs
// sigma'_1 with the INTERMEDIATE stress; in this file's that is s3 with s2, so it becomes active
// exactly where a return on the main surface would push s2 past s1. (The other corner, s2 = s3,
// is not a corner of this criterion at all: there f_12 evaluates to -fbar < 0, comfortably
// admissible.)
inline double yield12(double s2, double s3, const Params& P, const Constants& C) {
    return std::max(s2 - s3 - fbar(s2, P, C), s2 - C.sigt);
}

// Mobilised dilatancy, a function of the most tensile principal stress -- the criterion's
// sigma'_3, this file's s1. Compression negative, so "more confinement" is s1 more negative.
// In compression psi_mob = psi (sigma_psi + sigma'_3)/sigma_psi, floored at 0: below -sigma_psi
// the rock has stopped dilating. In the tensile range the value is increased artificially, from
// psi at 0 to 90 degrees at sigma_t, so that plastic expansion remains possible there.
inline double psi_mobilised(double s1, const Params& P, const Constants& C) {
    constexpr double kHalfPi = 1.57079632679489661923;   // 90 degrees, the upper limit
    if (s1 >= 0.0) {                                      // tensile zone: psi -> 90 deg at sigma_t
        if (C.sigt <= 0.0) return P.psi;
        const double t = std::min(1.0, s1 / C.sigt);
        return P.psi + t * (kHalfPi - P.psi);
    }
    if (P.sig_psi <= 0.0) return P.psi;                   // no decay declared
    const double sp = P.sig_psi;                          // (sigma_psi + sigma'_3)/sigma_psi
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
Return return_mapping(double s1t, double s2t, double s3t, const Params& P,
                      const Constants& C);


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

PlaneReturn plane_return(const PlaneStrainStress& trial, const Params& P,
                         const Constants& C);

} // namespace katai::core::hoekbrown
