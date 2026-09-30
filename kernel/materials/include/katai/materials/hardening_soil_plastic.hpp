#pragma once
// Hardening Soil — P2.3d (part 1): the multiaxial SHEAR hardening return mapping.
// Strain-driven predictor-corrector (FE constitutive form) in principal stress space.
// Elastic: stress-dependent Eur (frozen at σ3 at the start of the step). Shear yield
// surface (Schanz 1999):
//   f = f̄(q) − γ^p ,   f̄(q) = (2/Ei) q/(1−q/qa) − 2q/Eur ,   q = σ1 − σ3
// Non-associated flow: mobilized dilatancy ψ_m, m_g = (1, R, R), R = ε3^p/ε1^p. ψ_m is Rowe's
// where Rowe is non-negative; below the phase-transformation line HS takes zero and HSsmall
// takes Li & Dafalias instead — one rule, `detail::hs_dilatancy` (see its comment).
//   Rowe flow rule: ε_v^p/ε_q^p = −sinψ_m (DILATION, comp-pos ⇒ ε_v^p<0) ⇒
//   R = −(1+sinψ_m)/(2−sinψ_m) ; hardening γ^p=−(2ε1^p−ε_v^p) ⇒ dγ^p=h·dλ,
//   h = 1−2R = (4+sinψ_m)/(2−sinψ_m). At ψ_m=0, R=−½, h=2 (volumetrically neutral, the
//   hyperbola is preserved).
// Consistent (asymmetric) tangent D^ep = D_e − (D_e m_g)(n_f^T D_e)/(n_f^T D_e m_g + h).
// In the ψ=0 check the elastoplastic machinery PRODUCES the drained triaxial hyperbola
// (not closed-form) → the machinery's correctness (test_hs_shear). Math:
// hardening-soil-formulation.md.
//
// CONVENTION: compression-POSITIVE principal stresses σ1 ≥ σ2 ≥ σ3 (sigma[0..2]). The sign
// flips in FE (σ_HS = −σ_solver) — that integration is wired in a later part.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include <Eigen/Dense>

#include <katai/materials/hardening_soil.hpp>

namespace katai::core {

struct HsShearStep {
    Eigen::Vector3d stress;   // updated principal stresses
    double gamma_p;           // updated hardening parameter
    Eigen::Matrix3d tangent;  // consistent D^ep (principal space 3×3)
    bool plastic;
};

namespace detail {

// Isotropic elastic principal-space stiffness (E, ν).
inline Eigen::Matrix3d hs_elastic(double E, double nu) {
    const double f = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    Eigen::Matrix3d d;
    d << f * (1 - nu), f * nu, f * nu,
         f * nu, f * (1 - nu), f * nu,
         f * nu, f * nu, f * (1 - nu);
    return d;
}

// Mobilized dilatancy ψ_m as a function of the mobilized friction angle φ_m. This rule has
// FOUR callers (the two principal returns and their two cap-coupled siblings) and used to be
// written out four times; a rule that is right in three copies and stale in the fourth is a
// silently wrong answer on whichever stress path reaches the fourth. One definition, four
// bindings.
//
// HS (Rowe 1962; Schanz & Vermeer 1996):  sinψ_m = (sinφ_m − sinφ_cv)/(1 − sinφ_m·sinφ_cv),
// cut off to [0, sinψ], and only above sinφ_m = 3/4 sinφ (ψ_m = 0 below it); a non-positive ψ is
// taken as it is. Where Rowe returns a NEGATIVE value the Hardening Soil model takes ψ_m = 0.
//
// HSsmall (after Li & Dafalias (2000)): that zero cut-off can give too little plastic
// volumetric strain, so wherever Rowe is negative the small-strain model puts a small
// CONTRACTION there instead of nothing:
//   sinψ_m = (1/10)·(−M_c·exp[(1/15)·ln((M_c/M_d)·(q/q_a))] + M_d)          (i)
//   M_c    = 6·sinφ_cv/(3 − sinφ_cv)                                        (ii)
//   M_d    = 6·sinφ_m /(3 − sinφ_m )                                        (iii)
//   q/q_a  = max([(1−sinφ_cv)/sinφ_cv]·[sinφ_m/(1−sinφ_m)], 1e−4)           (iv)
//   sinφ_m ≥ sinφ/(2 − sinφ)                                                (v)
// The q/q_a of (iv) is its own definition, NOT the model's deviatoric hyperbola ratio q/q_a
// (q_a = q_f/R_f); the two are unrelated quantities that share a name.
//
// The two branches MEET EXACTLY, and that is the property worth knowing: at φ_m = φ_cv we get
// M_d = M_c and q/q_a = 1, so the logarithm vanishes and (i) returns exactly 0 — which is
// also precisely where Rowe changes sign. The composite rule is continuous by CONSTRUCTION,
// not to a tolerance. The floor (v) then bounds how contractant it can get, since sinψ_m falls
// monotonically as φ_m drops: the floor is the most negative ψ_m the model can produce.
// (test_hssmall pins both.)
//
// WORKED VALUES (φ=35°, ψ=5°): the upper branch is pure Rowe, the curve vanishes at
// φ_cv = 30.80°, and the plateau starts at the floor (v), φ_m = 23.71°, where (i) gives
// ψ_m = −1.679°. We implement the rule exactly as written above, with the leading 1/10. A
// plotted version of this rule whose amplitude is 1.29× larger (plateau −2.167°, as if the
// leading 1/10 were 1/7.74) is declared as a gap in hssmall-formulation.md; both readings are
// small contractions and the difference between them is far smaller than the difference
// between either and the ψ_m = 0 this replaces.
struct HsDilatancy {
    double sin_cs = 0.0;        // sinφ_cv (Rowe, from the input φ and ψ)
    double sin_psi = 0.0;       // sinψ — the upper cut-off
    double sin_phi = 0.0;       // sinφ — the 3/4 sinφ threshold of the HS rule
    double sin_psi_in = 0.0;    // sinψ as entered, sign kept: the HS rule's ψ <= 0 branch
    double c_cot = 0.0;         // c·cotφ — the cohesion shift of the mobilized-φ definition
    double Mc = 0.0;            // M_c, (ii) (HSsmall branch only)
    double sphi_m_floor = 0.0;  // floor (v) (HSsmall branch only)
    bool li_dafalias = false;   // HSsmall (G0_ref>0) with a usable φ_cv

    double operator()(double sphi_m) const {
        if (!li_dafalias) {
            // THE HARDENING SOIL RULE, in the form the model is implemented and documented in
            // practice -- Rowe's law (Schanz & Vermeer 1996) switched on only above a mobilised
            // friction of 3/4 sin(phi), which neither Rowe nor Benz (2007) states:
            //   sin phi_m <  3/4 sin phi              -> psi_m = 0
            //   sin phi_m >= 3/4 sin phi, psi > 0     -> sin psi_m = max(Rowe, 0)
            //   sin phi_m >= 3/4 sin phi, psi <= 0    -> psi_m = psi
            //   phi = 0                               -> psi_m = 0
            // Until 2026-09-30 the threshold was missing -- Rowe was taken as soon as it turned
            // positive, i.e. from phi_cv, which lies BELOW 3/4 phi whenever psi is large (phi =
            // psi = 41 degrees puts phi_cv at 0: dilation from the first increment of shear) --
            // and a negative psi, which the input accepts, was silently read as 0.
            if (!(sin_phi > 1e-12) || sphi_m < 0.75 * sin_phi) return 0.0;
            if (!(sin_psi_in > 0.0)) return sin_psi_in;
            const double rowe = (sphi_m - sin_cs) / (1.0 - sphi_m * sin_cs);
            return rowe > 0.0 ? std::min(rowe, sin_psi) : 0.0;
        }
        const double rowe = (sphi_m - sin_cs) / (1.0 - sphi_m * sin_cs);
        if (rowe >= 0.0) return std::min(rowe, sin_psi);
        const double s = std::max(sphi_m, sphi_m_floor);                // floor (v)
        const double Md = 6.0 * s / (3.0 - s);                          // M_d, (iii)
        const double qqa =
            std::max((1.0 - sin_cs) / sin_cs * (s / (1.0 - s)), 1e-4);  // q/q_a, (iv)
        return 0.1 * (-Mc * std::exp(std::log((Mc / Md) * qqa) / 15.0) + Md);  // sinψ_m, (i)
    }

    // sinφ_m = q/(σ1 + σ3 + 2c·cotφ) from the mobilized stress state (σ1 = q + σ3, comp-positive).
    double from_q(double q, double s3v) const {
        const double s1 = q + s3v;
        const double denom = s1 + s3v + 2.0 * c_cot;
        return (*this)(denom > 1e-12 ? q / denom : 0.0);
    }
};

inline HsDilatancy hs_dilatancy(const HardeningSoilParams& p) {
    HsDilatancy d;
    const double sphi = std::sin(p.friction), cphi = std::cos(p.friction);
    const double sps = std::sin(p.dilatancy);
    d.sin_cs = (sphi - sps) / (1.0 - sphi * sps);  // critical state
    d.sin_psi = sps > 0.0 ? sps : 0.0;
    d.sin_phi = sphi;
    d.sin_psi_in = p.dilatancy_cut ? 0.0 : sps;
    d.c_cot = (sphi > 1e-12) ? p.cohesion * cphi / sphi : 0.0;
    // φ_cv = 0 (a φ = 0 Tresca soil) would divide by zero in (iv) — but it also makes Rowe
    // non-negative everywhere, so the branch is unreachable there; the guard says so rather
    // than relying on it.
    d.li_dafalias = p.G0_ref > 0.0 && !p.dilatancy_cut && d.sin_cs > 1e-12 && sphi < 1.0;
    if (d.li_dafalias) {
        d.Mc = 6.0 * d.sin_cs / (3.0 - d.sin_cs);  // M_c, (ii)
        d.sphi_m_floor = sphi / (2.0 - sphi);      // floor (v)
    }
    return d;
}

}  // namespace detail

// One shear-hardening step: committed (σ_n, γ_n) + strain increment dε → updated state +
// consistent tangent. Eur/Ei/qa are frozen at the start-of-step σ3 = σ_n[2] (explicit
// stress-dependent stiffness). q is capped at qf at failure (perfectly
// plastic MC).
HsShearStep hs_shear_step(const HardeningSoilParams& p,
                          const Eigen::Vector3d& sigma_n, double gamma_p_n,
                          const Eigen::Vector3d& deps);

// Stress-driven SHEAR return (for the FE wrapping): given the elastic trial principal
// stress (compression-positive, σ1≥σ2≥σ3) + committed γ^p, returns the corrected principal
// stress + γ^p (NO tangent — numerical at the Voigt level). Stiffness frozen at σ3
// (sigma3_stiff = committed σ3). The trial-based sibling of hs_shear_step (the predictor
// lives outside).
struct HsPrincipalReturn {
    Eigen::Vector3d stress;  // corrected principal stress (compression-positive)
    double gamma_p;
};

HsPrincipalReturn hs_shear_correct(const HardeningSoilParams& p,
                                   const Eigen::Vector3d& sig_tr,
                                   double gamma_p_n, double sigma3_stiff);

// --- Cap (volumetric) yield surface -----------------------------------------------
// f_c = q̃²/α² + p² − p_p²,  p = (σ1+σ2+σ3)/3,  q̃ = σ1+(δ−1)σ2−δσ3, δ=(3+sinφ)/(3−sinφ)
// (q̃=q at triaxial σ2=σ3). ASSOCIATED flow g_c=f_c ⇒ dε^p=λ ∂f_c/∂σ, dε_v^pc=λ·2p.
// POWER-LAW hardening ε_v^pc = (β/(1−m))(p_c/p_ref)^(1−m), modulus
// H_cap=(p_ref/β)(p_c/p_ref)^m. In isotropic compression (q̃=0) the tangent is
// K_iso = K_e·H_cap/(K_e+H_cap) (springs in series). The consistent tangent (associated) is
// ~symmetric.
// Schanz, Vermeer & Bonnier (1999); hardening-soil-formulation.md §4.
struct HsCapStep {
    Eigen::Vector3d stress;
    double pp;                // updated preconsolidation pressure
    Eigen::Matrix3d tangent;  // consistent (symmetric) D^ep
    bool plastic;
};

HsCapStep hs_cap_step(const HardeningSoilParams& p,
                      const Eigen::Vector3d& sigma_n, double pp_n,
                      const Eigen::Vector3d& deps);

// --- Two-surface (shear + cap) return mapping -------------------------------------
// Strain-driven, principal-space (fixed principal directions — for the single-element
// triaxial/oedometer drivers and calibration). Active set: shear-only / cap-only / both
// (Koiter). In the both case (oedometer) σ, λ_s, λ_c are solved coupled (inner fixed
// point + outer 2×2 Newton). The σ3 stiffness is frozen at the committed minor. pp
// hardens with the cap volumetric λ_c·2p, γ^p with the shear λ_s.
struct HsState {
    Eigen::Vector3d stress;
    double gamma_p;
    double pp;
};

// sigma3_stiff: the minor principal that freezes the elastic Eur (committed σ3). The FE
// wrapping passes the committed σ3; negative sentinel ⇒ sigma_n(2) is used
// (old call).
HsState hs_return_principal(const HardeningSoilParams& p,
                            const Eigen::Vector3d& sigma_n, double gamma_p_n,
                            double pp_n, const Eigen::Vector3d& dstrain,
                            double sigma3_stiff = -1e300);

// --- ROBUST multi-surface integrator: explicit substepping + Koiter (Sloan/Potts&Gens) ---
// The path the literature takes for complex soil models: the outer strain increment is split
// into small substeps; at each
// substep the active surfaces (shear/cap) are determined, the plastic multipliers are
// solved from the 2×2 system with Koiter multi-surface flow (negative multiplier → that
// surface is dropped), stress+hardening are updated, and a drift correction pulls back to
// the surface. Unlike an implicit nested Newton it does NOT DIVERGE (explicit, small step)
// — robust including high stiffness (E=100+MPa). Strain-driven (the shared path of
// integrate_point + single-element triaxial/oedometer + calibration).
// Source: Sloan, Abbo & Sheng (2001); Potts & Zdravković; Schanz, Vermeer & Bonnier (1999).
struct HsIntegrated {
    Eigen::Vector3d stress;
    double gamma_p;
    double pp;
    Eigen::Matrix3d tangent;  // continuum elastoplastic (last active set)
    bool plastic;
    int nsub;                 // substeps the error control asked for
    int saturated = 0;        // 1 when the guard ceiling clipped the subdivision, i.e. the
                              // integration did NOT meet its tolerance on this call
};

// The subdivision one integration used. A perturbed run (the numerical consistent tangent's
// three forward differences) is REPLAYED along exactly this subdivision: with automatic
// substepping the perturbation would otherwise change the subdivision itself, and that
// difference is larger than the ~E*h signal the difference quotient is trying to read. This is
// the automatic-substepping replacement for the fixed substep count the previous rule pinned.
// The subdivision is UNIFORM (one size, chosen once, with a partial last step), so recording it
// takes two numbers rather than a list -- and the replay can never overflow a buffer.
struct HsSubstepPlan {
    double dT = 0.0;   // pseudo-time size of every substep but the remainder; 0 = nothing recorded
    int n = 0;         // substeps taken
    void clear() { dT = 0.0; n = 0; }
};

// Integrate the Hardening Soil model over `dstrain` from the committed state, with the local
// integration error held under `stol`.
//
// Scheme: explicit substepping with AUTOMATIC ERROR CONTROL (Sloan 1987; Sloan, Abbo & Sheng
// 2001, Engineering Computations 18(1):121-194). Each substep is evaluated twice -- forward
// Euler from the start of the substep, and again from its Euler end point -- and the two are
// averaged (modified Euler, second order). Their difference IS the local truncation error, and
// it is what sizes the subdivision: the size is measured before the walk and then held (see
// "HOW MANY SUBSTEPS" below for why this scheme does not accept/reject as the paper does). A
// drift correction pulls each state back onto the surfaces it was yielding on.
//
// What this replaces, and why. Until 2026-08-20 the substep count was
// `ceil(||De.deps|| / (0.01 (max|sigma| + p_ref)))` -- a fixed fraction of a reference stress,
// with no error estimate, no tolerance and no rejection. The code cited the paper above
// without implementing its contribution. Measured on the corpus oedometer (KV-CST-002),
// refining that fraction from 1% to 0.03% moved the answer by TWO PERCENTAGE POINTS and did
// not converge -- an error axis larger than the model deviation the numerics record attributes
// to the cap calibration, and one no sweep in this project had ever touched. It also made the
// answer depend on the iteration path: a different line search, or a different linear-solver
// backend, lands on a different subdivision, which is how two backends came to disagree on a
// published number (docs/validation/numerical-uncertainty.md, and the 0.9.0 audit).
//
// `stol` is a NUMERICAL CONTROL, not a material property: it belongs to the phase, like the
// tolerated error and the load-step count. 0 means "the default below".
// TEMPORARY (0.9.0 N-1, second half): the tolerance's permanent home is the phase's numerical
// controls, next to the tolerated error and the load-step count, so that a published run carries
// the integration accuracy it was computed with. Until that plumbing lands it is a build-time
// default with an environment override, so the tolerance can be SWEPT -- which is the whole
// point of having one.
// The guard ceiling on one increment's subdivision, and the measurement that set it. On the
// corpus oedometer (KV-CST-002, range 100->200 kPa, STOL 1e-5) the ceiling costs and buys:
//     100 -> -0.9539%, 24 s | 200 -> -0.9861%, 32 s | 500 -> -1.0689%, 67 s | 2000 -> -1.0554%, 67 s
// Nearly all of the extra work goes into TRIAL iterates far from equilibrium, whose stresses are
// discarded; the answer moves by about a tenth of a percentage point across a twentyfold ceiling,
// which is this fixture's path indeterminacy, not integration error. 200 is the cheapest setting
// on that plateau. Saturation is reported (`HsIntegrated::saturated`) so a run that hit it
// is never mistaken for one that met its tolerance.
// TEMPORARY measurement seam (scaffold): the environment override, so the ceiling stays swept.
inline int hs_max_substeps() {
    static const int v = [] {
        const char* e = std::getenv("KATAI_HS_MAXSUB");
        const int d = e ? std::atoi(e) : 200;
        return d > 0 ? d : 200;
    }();
    return v;
}

// The seam that used to stand here (KATAI_HS_MCPROJ, returning the failure bound along the
// consistent direction De.n_s instead of clamping sigma1) is GONE, swept and removed on
// 2026-08-24. It was kept because a measurement said it fixed a 0.24% shortfall on the failure
// plateau at a cost of one boundary-value case. Re-measured, the fix is not there: both returns
// leave the deviator at the SAME 99.686% of q_f, and the cost still is. The shortfall was never
// the return direction -- see the bound below.

// The default integration tolerance, and the measurement that set it. KV-CST-002 (100->200 kPa,
// 40+160 increments, equilibrium tolerance 1e-6); the same material walked at the STRESS POINT
// says -1.003%, so that column is the target:
//     STOL    deviation   Newton iterations (seating/staged)   wall clock
//     1e-3    -1.536%     1887 / 4883                          91 s
//     1e-4    -1.310%     1907 / 1427                          76 s
//     1e-5    -0.986%      432 / 1186                          31 s
//     1e-6    -0.982%      387 / 3379                         117 s
// A LOOSER integration is not cheaper. Below about 1e-5 the integration noise sits above the
// residual the phase is asking for, and the equilibrium iteration grinds against it -- four times
// the seating iterations at 1e-3, for an answer half a percentage point further from the model's
// own. 1e-5 is the fastest setting AND the first one that reproduces the constitutive routine to
// the digits this case is published in.
inline double hs_default_substep_tol() {
    static const double v = [] {
        const char* e = std::getenv("KATAI_HS_STOL");
        const double d = e ? std::atof(e) : 1.0e-5;
        return d > 0.0 ? d : 1.0e-5;
    }();
    return v;
}

// The same value read as a whole-run OVERRIDE rather than as a default: 0 when the variable is
// not set. Since .k2d v15 a phase can name its own integration tolerance, and a study sweeping
// this axis has to be able to re-run a file that does -- otherwise the one thing the seam exists
// for, asking what a published number owes to its integration, is exactly what it cannot do on
// the files that state it. So the environment wins over the file here, which is the same
// precedence KATAI_CONV_NOLOCAL has over the phase's convergence setting.
inline double hs_env_substep_tol() {
    static const double v = [] {
        const char* e = std::getenv("KATAI_HS_STOL");
        const double d = e ? std::atof(e) : 0.0;
        return d > 0.0 ? d : 0.0;
    }();
    return v;
}

HsIntegrated hs_integrate(const HardeningSoilParams& p,
                          const Eigen::Vector3d& sigma_n, double gamma_p_n,
                          double pp_n, const Eigen::Vector3d& dstrain,
                          double stol = 0.0,
                          HsSubstepPlan* plan_out = nullptr,
                          const HsSubstepPlan* plan_in = nullptr);

// Cap preconsolidation initialization (FE initial state): pp = p_eq · OCR,
// p_eq = √(3J2/α² + p²) (the isotropic-equivalent pressure at which the cap passes through
// σ0, the state parameter p_eq). NC (OCR=1) ⇒ the initial state sits on the
// cap (f_c=0); this keeps the cap well-defined in FE (pp_n=0 would make the cap yield
// always). sig_comp_pos = compression-positive principal stress (σ_HS = −σ_solver).
inline double hs_initial_pp(const HardeningSoilParams& p,
                            const Eigen::Vector3d& sig_comp_pos, double OCR = 1.0) {
    const double pm = (sig_comp_pos(0) + sig_comp_pos(1) + sig_comp_pos(2)) / 3.0;
    // The cap's own deviatoric measure q~ = s1 + (delta - 1) s2 - delta s3 (hs_integrate, Benz 2007
    // Eqn 7.21), so the seeded state is ON the cap the integrator reads. It used to be the von
    // Mises q, which is the same number only where two principals are equal; with sigma_xy != 0 (a
    // slope, a POP history under an inclined surface) the seed sat outside the integrator's cap and
    // the first increment spent itself pulling it back.
    Eigen::Vector3d s = sig_comp_pos;
    std::sort(s.data(), s.data() + 3, [](double x, double y) { return x > y; });
    const double sphi = std::sin(p.friction);
    const double delta = (3.0 + sphi) / (3.0 - sphi);
    const double qt = s(0) + (delta - 1.0) * s(1) - delta * s(2);
    const double a = p.cap_alpha;
    return std::sqrt(qt * qt / (a * a) + pm * pm) * OCR;
}

// Shear hardening initialization (FE initial state): γ^p = f̄(q0). The geostatic K0 state
// carries q0>0; starting with γ^p=0 the shear surface (f_s=f̄(q0)−γ^p) is violated FROM THE
// START → the first return makes a large inconsistent correction (the solver diverges).
// γ^p=f̄(q0) seats the state on the surface (admissible pre-stress), the shear counterpart
// of the pp init (cap). sig_comp_pos = compression-positive principals.
inline double hs_initial_gamma_p(const HardeningSoilParams& p,
                                 const Eigen::Vector3d& sig_comp_pos) {
    // The same laws the integrator reads: the stiffness at the floored minor stress, the strength
    // and the hyperbola's asymptote at the stress itself, guarded at a thousandth of the floor
    // (hs_integrate, stiff_at). At a stress-free point the raw stiffness is zero and q/q_a 0/0.
    const double s1 = sig_comp_pos.maxCoeff(), s3 = sig_comp_pos.minCoeff();
    const double s3k = std::max(s3, p.p_limit());
    const double Ei = p.Ei(s3k), Eur = p.Eur(s3k), qf = p.q_failure(s3);
    const double qa = std::max(qf, p.q_failure(1e-3 * p.p_limit())) / p.Rf;
    if (qf <= 0.0) return 0.0;
    const double q = std::min(std::max(s1 - s3, 0.0), qf);
    const double fbar = (2.0 / Ei) * q / (1.0 - q / qa) - 2.0 * q / Eur;
    return std::max(fbar, 0.0);
}

// The initial hardening variables of an in-situ state and its overconsolidation, seeded from
// the PRE-CONSOLIDATION stress state of the K0 procedure, stated in full:
//     sigma'_yy,c = OCR sigma'_yy          (or sigma'_yy + POP)
//     sigma'_xx,c = sigma'_zz,c = K0nc sigma'_yy,c ,   sigma_xy,c = sigma_xy
// i.e. the state the soil was loaded to along the normally consolidated path before it was
// unloaded to where it is. pp is the cap through that state. gamma_p is the shear hardening that
// primary loading to that state leaves behind: the same history, read by the other surface, so
// an overconsolidated soil is elastic in shear up to the deviator it has already carried.
//
// What this replaces, measured. pp was the cap through the CURRENT state scaled by OCR -- the
// in-situ K0, not K0nc -- and a POP was turned into that ratio for every component; gamma_p was
// the current deviator's. On a clay at K0 = 0.61 with POP = 32 (sigma'_v = 126.5) that gave
// pp = 149.8 against 156.3 here, and an undrained compression that left the elastic range at
// q = 71 where a published element test of the same parameter set does at 87.
//
// Either state is taken if it is the larger: a soil given an in-situ K0 above K0nc x OCR would
// otherwise start outside its own surfaces. Compression positive, Cartesian (xx, yy vertical,
// xy, zz); mode 0 = neither, 1 = OCR, 2 = POP.
struct HsSeed { double pp, gamma_p; };
HsSeed hs_seed_from_history(const HardeningSoilParams& p, double sxx, double syy,
                            double sxy, double szz, int mode, double OCR, double POP);

// --- Oedometer probe + cap calibration (on the robust hs_integrate) ------------------
// Oedometer (ε_h=0) NC primary loading (with hs_integrate, robust); returns the tangent
// Eoed at σ1=p_ref and the lateral ratio K0. (α,β) → (K0,Eoed).
// The tolerance the CALIBRATION integrates at. Deliberately independent of, and tighter than,
// whatever a run uses: alpha and beta are properties of the MATERIAL (they are what makes the
// model reproduce Eoed_ref and K0_NC), so they must not move when a phase asks for a looser or
// tighter integration. Before 2026-08-20 the probe inherited the integrator's fixed 1% substep
// rule, which tied every calibrated cap -- and therefore every Hardening Soil answer this
// program has ever produced -- to a numerical constant, in a second and completely invisible
// way. This is the fix for that half of it.
inline constexpr double kHsCalibrationTol = 1.0e-6;

void hs_oedometer_probe(const HardeningSoilParams& p, double& Eoed_pref,
                        double& K0);

// Cap hardening modulus K_p = K1·K2/(K1−K2) — CLOSED FORM from Eoed_ref (elastic and
// plastic bulk springs in series, 1/K2 = 1/K1 + 1/K_p). K1=Eur_ref/(3(1−2ν)) unloading bulk;
// K2=Eoed_ref(1+2K0nc)/3. This ties the cap hardening directly to Eoed (no numerical
// calibration of β needed) ⇒ the K0-Eoed coupling is resolved: β=p_ref/(k·K_p), only α
// (and a small k correction) is calibrated.
inline double hs_cap_Kp(const HardeningSoilParams& p, double K0_NC) {
    const double K1 = p.Eur_ref / (3.0 * (1.0 - 2.0 * p.nu_ur));
    const double K2 = p.Eoed_ref * (1.0 + 2.0 * K0_NC) / 3.0;
    return (K1 > K2) ? K1 * K2 / (K1 - K2) : K1;  // K1>K2 (Eur≫Eoed) typical
}

// Cap parameters (α, β) from the standard HS inputs (Eoed_ref, K0_NC). β = p_ref/(k·K_p) (K_p
// closed form, sets Eoed in CLOSED FORM → the coupling is resolved); α by outer bisection for
// K0_NC; k by inner bisection as an Eoed_ref fine correction (k≈1). Both are fitted by
// simulated oedometer tests. If K0_NC is unreachable, the nearest α.
void hs_calibrate_cap(HardeningSoilParams& p, double K0_NC);

}  // namespace katai::core
