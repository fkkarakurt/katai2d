#pragma once
// Material-point constitutive interface (P1.1) — the nonlinear-analysis core.
//
// Design (see docs/ARCHITECTURE + the DOD decision): NO virtual calls (vtable) on the hot
// path. A material is a tagged flat data structure; integration is a free function
// (`integrate_point`) that switches on the type — it inlines, is branch-predictor
// friendly, and can be vectorized later.
//
// The integrate_point contract (return-mapping ready):
//   input : the converged (committed) state + this step's TOTAL strain increment
//           Δε = [Δexx, Δeyy, Δgxy]
//   output: the trial state (current stress + history) + the consistent tangent D_T
// The caller does committed → integrate(Δε) → trial each Newton iteration; when the step
// converges it commits the trial.

#include <Eigen/Core>

#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/mohr_coulomb.hpp>
#include <katai/materials/soft_soil.hpp>
#include <katai/materials/soft_soil_creep.hpp>

namespace katai::core {

// Depth-VARYING stiffness/strength profile (indexed by material id; PLAXIS "E'_inc /
// c'_inc" + y_ref). The most common real-soil case is stiffness growing with depth — in
// seismics E(y) is directly the Vs profile itself:
//     E(y) = E_ref + E_inc·(y_ref − y)      (decreases ABOVE y_ref)
//     c(y) = c_ref + c_inc·(y_ref − y)      (Mohr-Coulomb strength only)
// Evaluated PER STRESS (Gauss) POINT — PLAXIS does it per stress point too. An element
// average would silently drift on a coarse mesh; the whole meaning of a gradient is that
// it varies within an element. uniform() ⇒ callers use the old constant-E path →
// BIT-FOR-BIT the same result.
struct MaterialProfile {
    double E_inc = 0.0;   // dE/d(depth) [kN/m²/m]; 0 = uniform (default)
    double c_inc = 0.0;   // dc/d(depth) [kN/m²/m]; 0 = uniform
    double y_ref = 0.0;   // reference elevation [m] where E = E_ref, c = c_ref
    bool uniform() const { return E_inc == 0.0 && c_inc == 0.0; }
};

// The profile value at elevation y. CLAMPED AT ZERO: negative stiffness/strength is
// unphysical and breaks the solver (K goes singular/indefinite) — it really happens when
// a user puts y_ref outside the domain with a large gradient, so it is prevented once,
// centrally, here.
inline double profile_at(double ref, double inc, double y_ref, double y) {
    const double v = ref + inc * (y_ref - y);
    return v > 0.0 ? v : 0.0;
}

enum class MaterialType {
    LinearElastic,
    MohrCoulomb,
    HardeningSoil,
    SoftSoil,
    SoftSoilCreep,
};

// State of an integration (Gauss) point. Besides the in-plane Voigt stress we
// carry the out-of-plane normal stress sigma_zz: under the plane-strain
// constraint eps_zz = 0 it is generally non-zero and enters the Mohr-Coulomb
// yield evaluation as a principal stress. (Hardening / plastic-strain history
// will extend this struct in later phases.)
struct GaussState {
    Eigen::Vector3d stress = Eigen::Vector3d::Zero();  // Voigt [sxx, syy, sxy]
    double stress_zz = 0.0;                            // sigma_zz
    // Accumulated (committed) volumetric strain eps_v = eps_xx+eps_yy(+eps_theta).
    // For Undrained (A) the excess pore pressure is u = -(Kw/n) eps_vol (tension-
    // positive). Zero / unused for drained materials. (See effective-stress-
    // formulation.md.)
    double eps_vol = 0.0;
    // Hardening Soil history: shear hardening parameter gamma_p and cap pre-
    // consolidation pressure pp (compression-positive). Zero / unused for other models.
    double gamma_p = 0.0;
    double pp = 0.0;
    // HSsmall: accumulated deviatoric shear strain γ_hist (monotone; drives the
    // small-strain stiffness degradation). Only for HSsmall (G0_ref>0); 0/unused otherwise.
    double gamma_hist = 0.0;
};

// Flat, tagged material descriptor (no vtable -> suitable for the hot loop).
struct MaterialModel {
    MaterialType type = MaterialType::LinearElastic;
    double youngs_modulus = 0.0;   // E
    double poisson_ratio = 0.0;    // v
    double cohesion = 0.0;         // c        (Mohr-Coulomb)
    double friction_angle = 0.0;   // phi [rad] (Mohr-Coulomb)
    double dilatancy_angle = 0.0;  // psi [rad] (Mohr-Coulomb)
    // Rankine tension cap (PLAXIS MMM Eq 3-11; sigma_1 <= sigma_t, associated flow).
    // Read by the MohrCoulomb branch of integrate_point / integrate_point_axisym;
    // HS/SS do not consume it yet (their integrators lack the extra planes -- the
    // GUI states this honestly). Default OFF keeps every direct caller bit-identical.
    bool tension_cutoff = false;
    double tensile_strength = 0.0;  // sigma_t [kN/m2], tension-positive
    // Dilatancy cut-off (Material Models Manual Eq. 5.16b, Fig. 5.6): "After extensive shearing,
    // dilating materials arrive in a state of critical density where dilatancy has come to an
    // end... As soon as the volume change results in a state of maximum void, the mobilised
    // dilatancy angle is automatically set back to zero." Without it a dense sand dilates for
    // ever and its bearing capacity is over-predicted -- an unsafe number, produced quietly.
    //
    // The void ratio follows the volume change: 1 + e = (1 + e_init) exp(eps_v), expansion
    // positive. The manual writes the same statement through ln((1+e)/(1+e_init)) (Eq. 5.17)
    // with a sign convention that is ambiguous as printed; the exponential form says it once.
    // e_min is deliberately NOT here: the manual states it "is not used within the context of
    // the Hardening-Soil model", so storing it would suggest a rule that does not exist.
    bool dilatancy_cutoff = false;
    double e_init = 0.5, e_max = 1.0;

    // Undrained (A): effective parameters above (E', nu', c', phi'); the pore fluid's
    // volumetric stiffness Kw/n is added to the GLOBAL tangent (D_u = D' + (Kw/n)mm^T)
    // and the internal force uses TOTAL stress = sigma' + (Kw/n) eps_v m, while the
    // constitutive model still works on effective stress. The solver tracks the excess
    // pore pressure in GaussState::eps_vol. undrained_poisson = nu_u (PLAXIS default
    // 0.495; exactly 0.5 is singular).
    bool undrained = false;
    double undrained_poisson = 0.495;  // nu_u
    // The EFFECTIVE elastic pair the pore-fluid derivation uses, when it is not (E, nu).
    // Zero = derive from youngs_modulus / poisson_ratio, which is what Linear Elastic and
    // Mohr-Coulomb want -- and it keeps a depth gradient E'_inc flowing into Kw/n, since the
    // caller hands this struct a per-Gauss copy with the profiled modulus already in place.
    // The Hardening Soil family sets it: that model's elasticity is the unload/reload pair
    // (Eur_ref, nu_ur), while E and nu are boxes it never reads. Deriving the pore fluid's
    // stiffness from those boxes gave an undrained HS soil a water stiffness sized by a
    // default the user never entered -- typically a factor of several, in whichever direction
    // the untouched field happened to sit.
    double undrained_E_ref = 0.0, undrained_nu_ref = 0.0;

    // Undrained (C): the material is analysed in TOTAL stress (MMM section 2.7). The stiffness
    // and strength above are the undrained ones, no pore pressure is generated or carried, and
    // the stress this model returns is total stress wearing the effective stress's name. The
    // flag exists because the difference is invisible in the parameters: an undrained Tresca
    // envelope looks exactly like a drained one with phi = 0, and something has to know which
    // partial factor an EC7 run should apply to that cohesion (gamma_cu, not gamma_c').
    bool total_stress = false;

    // Hardening Soil parameters (used when type == HardeningSoil). Compression-positive
    // internally; the FE wrapper converts (sigma_HS = -sigma_solver). See
    // hardening_soil_plastic.hpp / docs/references/hardening-soil-formulation.md.
    HardeningSoilParams hs;

    // Soft Soil parameters (used when type == SoftSoil; PLAXIS MMM §10, soft_soil.hpp /
    // docs/references/soft-soil-formulation.md). youngs_modulus/poisson_ratio are NOT USED —
    // the stiffness comes from the ln-law (K = p'/κ*). CAUTION (honest Stage-2 limit):
    // SS + undrained (A) is NOT wired yet — kw_over_n() derives from E and would silently
    // be 0 with E=0; Stage 3 will add undrained support via the stress-dependent K_ur. The
    // GUI/project file does not produce this model until Stage 3; the kernel tests use it
    // drained.
    softsoil::Params ssoil;

    // Soft Soil Creep parameters (type == SoftSoilCreep; PLAXIS MMM §11,
    // soft_soil_creep.hpp / docs/references/soft-soil-creep-formulation.md). Time enters
    // via integrate_point's trailing parameter dt_day (0 = no creep — in PLAXIS too, SSC in
    // a phase without a time interval gives only elastic+MC). Undrained/Safety limits as SS.
    softsoilcreep::Params ssc;

    // Plane-strain elastic constitutive matrix (3x3 SPD, v < 0.5).
    Eigen::Matrix3d elastic_plane_strain() const {
        const double e = youngs_modulus, v = poisson_ratio;
        const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
        Eigen::Matrix3d d = Eigen::Matrix3d::Zero();
        d(0, 0) = f * (1.0 - v);
        d(0, 1) = f * v;
        d(1, 0) = f * v;
        d(1, 1) = f * (1.0 - v);
        d(2, 2) = f * (1.0 - 2.0 * v) / 2.0;
        return d;
    }

    // Undrained (A) pore-fluid bulk stiffness Kw/n from the EFFECTIVE parameters and
    // an assumed undrained Poisson ratio (PLAXIS default nu_u = 0.495; exactly 0.5
    // makes the stiffness singular). MMM Eq. 2-50. Derived so that K' + Kw/n = the
    // correct undrained bulk modulus Ku (see effective-stress-formulation.md).
    double kw_over_n(double nu_u) const {
        const bool own = undrained_E_ref > 0.0;   // the model carries its own elastic pair
        const double e = own ? undrained_E_ref : youngs_modulus;
        const double v = own ? undrained_nu_ref : poisson_ratio;
        const double k_eff = e / (3.0 * (1.0 - 2.0 * v));  // K'
        return 3.0 * (nu_u - v) / ((1.0 - 2.0 * nu_u) * (1.0 + v)) * k_eff;
    }

    // K' as the pore-fluid derivation sees it -- the denominator of Skempton's B and the
    // number a GUI or a report needs to show Kw/n as a multiple of the skeleton stiffness.
    double undrained_k_eff() const {
        const bool own = undrained_E_ref > 0.0;
        const double e = own ? undrained_E_ref : youngs_modulus;
        const double v = own ? undrained_nu_ref : poisson_ratio;
        return e / (3.0 * (1.0 - 2.0 * v));
    }

    // Undrained (A) plane-strain stiffness D_u = D' + (Kw/n) m m^T, m = [1,1,0]:
    // the pore fluid adds volumetric (bulk) stiffness, leaving the shear part (G')
    // unchanged. Used in the global stiffness/tangent; the constitutive model still
    // works on EFFECTIVE stress.
    Eigen::Matrix3d undrained_plane_strain(double nu_u) const {
        Eigen::Matrix3d d = elastic_plane_strain();
        const double kwn = kw_over_n(nu_u);
        d(0, 0) += kwn; d(0, 1) += kwn; d(1, 0) += kwn; d(1, 1) += kwn;
        return d;
    }

    // Axisymmetric elastic matrix (4x4), strain/stress order [r, z, rz, theta].
    Eigen::Matrix4d elastic_axisym() const {
        const double e = youngs_modulus, v = poisson_ratio;
        const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
        Eigen::Matrix4d d = Eigen::Matrix4d::Zero();
        d(0, 0) = d(1, 1) = d(3, 3) = f * (1.0 - v);
        d(0, 1) = d(1, 0) = d(0, 3) = d(3, 0) = d(1, 3) = d(3, 1) = f * v;
        d(2, 2) = f * (1.0 - 2.0 * v) / 2.0;
        return d;
    }
};

// --- The undrained stiffness trio (PLAXIS MMM section 2.4, alpha_Biot = 1) -------------------
// Three quantities describe the same pore fluid and the manual gives one equation for each:
//   Kw/n  from nu_u        Eq. 2-50, MaterialModel::kw_over_n above
//   nu_u  from Skempton B  Eq. 2-55, undrained_poisson_from_skempton
//   B     from Kw/n        Eq. 2-57, skempton_from_kw_over_n
// They are mutually consistent -- going round the ring returns the number it started from, which
// is what test_undrained_stiffness measures rather than assumes. Which one the USER supplies is
// the choice PLAXIS offers (nu-undrained definition: Direct or Skempton-B based); the other two
// are then derived and are worth showing, because a B of 0.98 and a Kw/n of 45 K' are the same
// statement and an engineer recognises one of them.

// Eq. 2-55 with alpha_Biot = 1: nu_u = (3 nu' + B (1 - 2 nu')) / (3 - B (1 - 2 nu')).
// B -> 1 gives exactly 0.5 (incompressible, singular), so the caller must keep B < 1.
inline double undrained_poisson_from_skempton(double B, double nu_eff) {
    const double t = B * (1.0 - 2.0 * nu_eff);
    return (3.0 * nu_eff + t) / (3.0 - t);
}

// Eq. 2-57 with alpha_Biot = 1: B = 1 / (1 + n K'/Kw), written on the quantity the engine
// actually carries -- B = (Kw/n) / (K' + Kw/n) = (Kw/n) / Ku. The porosity cancels because
// Kw and n only ever enter as the ratio; PLAXIS's own input box is Kw,ref / n for that reason.
inline double skempton_from_kw_over_n(double kw_over_n, double k_eff) {
    return kw_over_n / (k_eff + kw_over_n);
}

// Frozen unload/reload modulus Eur for the HS forward update: stress-dependent Eur(sigma3)
// evaluated at the committed minor principal (compression-positive), clamped at hs_integrate's
// stiffness floor (0.1 p_ref) so the elastic predictor, the principal compliance and the
// substepping integrator all freeze the SAME Eur -> elastic steps reconstruct the trial exactly.
// Caller (kinematics-specific) builds the elastic trial with this Eur. pe: (HSsmall-adjusted).
inline double hs_frozen_Eur(const HardeningSoilParams& pe,
                            const Eigen::Vector3d& comm_in_plane, double comm_zz) {
    const double cxx = comm_in_plane(0), cyy = comm_in_plane(1), cxy = comm_in_plane(2);
    const double cmean = 0.5 * (cxx + cyy);
    const double cR = std::sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);
    const double s1c = std::max(std::max(cmean + cR, cmean - cR), comm_zz);  // major tension-pos
    const double s3_stiff = std::max(-s1c, 0.1 * pe.p_ref);                  // comp-pos minor
    return pe.Eur(s3_stiff);
}

// HSsmall (G0_ref>0): scale Eur_ref by the small-strain over-stiffness Et/Eur(gamma_hist).
// G0_ref=0 -> returns p unchanged (plain HS, byte-identical). (Material Models Manual sec 7.)
//
// The threshold is the RELOADING one, 2*gamma07 (Masing, Eq 7-11): what this function sets is
// the model's quasi-elastic (unload/reload) stiffness, and the manual keeps that factor
// constant at 2 throughout loading rather than switching it on at a reversal. Riding the
// virgin backbone here instead would degrade the stiffness twice as fast as the manual's
// model: measured 12.9% more heave on the unloading case KV-CST-008 (0.7996 mm against the
// manual's curve, 0.7132 mm) -- and softer is not the safe side when the number being read is
// a wall deflection or a heave.
inline HardeningSoilParams hs_small_strain_params(const HardeningSoilParams& p,
                                                  double gamma_hist) {
    HardeningSoilParams pe = p;
    if (p.G0_ref > 0.0) {
        const double d = 1.0 + HardeningSoilParams::kHDa * gamma_hist / p.gamma07_reload();
        const double E0_ref = 2.0 * (1.0 + p.nu_ur) * p.G0_ref;
        pe.Eur_ref = std::max(E0_ref / (d * d), p.Eur_ref);
    }
    return pe;
}

// KINEMATICS-AGNOSTIC Hardening Soil principal return. The in-plane block (xx,yy,xy) is the
// 2x2 Mohr circle and the out-of-plane normal (sigma_zz plane strain / hoop sigma_theta
// axisymmetric) is the third principal -- the same structure in both modes, so the principal-
// space return is identical and only the caller's elastic predictor / elastic operator differ.
// Inputs: committed stress (comm_*), the elastic TRIAL stress (trial_*, from the caller's
// predictor with the SAME frozen Eur), the (HSsmall-adjusted) params and that Eur. Returns the
// returned (admissible) stress (in-plane Voigt + out-of-plane), the updated history (gamma_p,
// pp), and a PrincipalTangent (3x3 plane-strain D_T + 4x4 algo_jacobian d sigma^ret/d sigma^tr,
// for axisymmetry via algo_jacobian * D_e_axisym). sigma_HS = -sigma_solver.
struct HsReturnCore {
    Eigen::Vector3d in_plane;   // returned [sxx,syy,sxy] / [srr,szz,srz]
    double zz;                  // returned sigma_zz / sigma_theta
    double gamma_p, pp;
    PrincipalTangent tan;
    bool plastic;               // did any surface activate this increment
    int nsub;                   // substep count used by hs_integrate
};

// The material's own tension cut-off, as the cap the principal returns take. PLAXIS switches
// the cut-off ON by default, and so does this schema, so `off` here is a deliberate choice by
// the engineer and not a default nobody looked at.
inline double tension_cap_of(const MaterialModel& m) {
    return m.tension_cutoff ? m.tensile_strength : kNoTensionCap;
}

// nsub_fixed: >0 pins the hs_integrate substep count (so the numerical consistent
// tangent's perturbed runs follow the SAME substep sequence as the base run).
inline HsReturnCore hs_return_core(const HardeningSoilParams& pe, double Eur,
                                   const Eigen::Vector3d& comm_in_plane, double comm_zz,
                                   const Eigen::Vector3d& trial_in_plane, double trial_zz,
                                   double gamma_p_n, double pp_n, int nsub_fixed = 0,
                                   double sigma_t_cap = kNoTensionCap) {
    const double cxx = comm_in_plane(0), cyy = comm_in_plane(1), cxy = comm_in_plane(2);
    const double cmean = 0.5 * (cxx + cyy);
    const double cR = std::sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);

    // Spectral decomposition of the trial stress (tension-positive), tracking sources.
    const double sxx = trial_in_plane(0), syy = trial_in_plane(1), sxy = trial_in_plane(2);
    const double mean = 0.5 * (sxx + syy), hd = 0.5 * (sxx - syy);
    const double radius = std::sqrt(hd * hd + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = hd / radius; sin2t = sxy / radius; }
    struct PV { double v; int src; };
    PV pv[3] = {{mean + radius, 0}, {mean - radius, 1}, {trial_zz, 2}};
    std::sort(pv, pv + 3, [](const PV& a, const PV& b) { return a.v > b.v; });

    // Robust strain-driven substepping return (hs_integrate). Feed the committed principal
    // stress + the principal elastic strain increment that reproduces the trial under the same
    // frozen Eur (deps_p = C_e (sigma_tr - sigma_n)) -> substep 1 reconstructs the trial exactly
    // (elastic steps byte-identical), then substepping traces the loading path (unconditionally
    // stable). The earlier nested-Newton projection diverged on shear-dominated, low-confinement
    // BVPs (strip footings, free surfaces).
    const Eigen::Vector3d sigHS(-pv[2].v, -pv[1].v, -pv[0].v);   // trial, comp-pos desc
    double cp[3] = {cmean + cR, cmean - cR, comm_zz};            // committed, tension-pos
    std::sort(cp, cp + 3, [](double a, double b) { return a > b; });
    const Eigen::Vector3d sig_n_cp(-cp[2], -cp[1], -cp[0]);      // committed, comp-pos desc
    const double nu = pe.nu_ur;
    Eigen::Matrix3d Ce;  // principal elastic compliance at the frozen Eur
    Ce << 1.0, -nu, -nu, -nu, 1.0, -nu, -nu, -nu, 1.0;
    Ce /= Eur;
    const Eigen::Vector3d deps_p = Ce * (sigHS - sig_n_cp);
    const HsIntegrated ret = hs_integrate(pe, sig_n_cp, gamma_p_n, pp_n, deps_p, nsub_fixed);
    double r[3] = {-ret.stress(2), -ret.stress(1), -ret.stress(0)};  // tension, desc
    // Tension cut-off (MMM Eq. 3-11), applied to the principals the model's own return produced.
    // Off by default, and the branch is not taken when the largest principal is admissible, so
    // every run without a cut-off is bit-identical to the one before this existed.
    const bool capped = apply_rankine_cap(r, sigma_t_cap, lame_from(Eur, pe.nu_ur).lambda,
                                          lame_from(Eur, pe.nu_ur).mu);

    // Analytic consistent (continuum) tangent: the shared principal-space assembly (spin +
    // out-of-plane coupling) used by Mohr-Coulomb, fed J(i,j)=J_comp(2-i,2-j) where J_comp =
    // D_pp C_e is the comp-positive stress-to-stress Jacobian (sign flip leaves J invariant,
    // the sort order reverses). Returns both the 3x3 plane-strain D_T and the 4x4 algo_jacobian.
    const Eigen::Matrix3d Jcomp = ret.tangent * Ce;
    Eigen::Matrix3d Jten;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Jten(i, j) = Jcomp(2 - i, 2 - j);
    const int src[3] = {pv[0].src, pv[1].src, pv[2].src};
    const LameConstants lame_ur = lame_from(Eur, nu);

    HsReturnCore out;
    out.tan = principal_consistent_tangent(cos2t, sin2t, radius, src, r, Jten, lame_ur);

    double pa = 0.0, pb = 0.0, pz = 0.0;  // coaxial reconstruction on the trial eigenframe
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == 0) pa = r[i];
        else if (pv[i].src == 1) pb = r[i];
        else pz = r[i];
    }
    const double m2 = 0.5 * (pa + pb), r2 = 0.5 * (pa - pb);
    out.in_plane(0) = m2 + r2 * cos2t;
    out.in_plane(1) = m2 - r2 * cos2t;
    out.in_plane(2) = r2 * sin2t;
    out.zz = pz;
    out.gamma_p = ret.gamma_p;
    out.pp = ret.pp;
    out.plastic = ret.plastic || capped;
    out.nsub = ret.nsub;
    return out;
}

// Hardening Soil FE forward stress update (PLANE STRAIN). Frozen Eur -> plane-strain elastic
// predictor (deps_zz=0) -> kinematics-agnostic principal return -> analytic 3x3 tangent.
// nsub_fixed / plastic_out / nsub_out: the numerical consistent-tangent plumbing (see
// hs_consistent_tangent) — the base run reports its nsub, the perturbed runs pin it.
// Has the dilatancy cut-off been reached at this stress point? (Material Models Manual
// Eq. 5.16b.) The state carries the accumulated volumetric strain, and the void ratio follows
// the volume change: 1 + e = (1 + e_init) exp(eps_v), expansion positive. The question is asked
// of the COMMITTED state -- the void ratio at the start of the increment -- so the answer is
// the same for every iteration of that increment and the return mapping stays a pure function
// of its inputs, which is what makes the tangent consistent and the line search safe.
inline double void_ratio_of(const MaterialModel& m, const GaussState& s) {
    return (1.0 + m.e_init) * std::exp(s.eps_vol) - 1.0;
}
inline bool dilatancy_cut(const MaterialModel& m, const GaussState& s) {
    return m.dilatancy_cutoff && void_ratio_of(m, s) >= m.e_max;
}
// The dilatancy angle the return mapping must use: zero once the soil has dilated to its
// critical void ratio, the material's own value before that.
inline double effective_dilatancy(const MaterialModel& m, const GaussState& s) {
    return dilatancy_cut(m, s) ? 0.0 : m.dilatancy_angle;
}

inline void hs_forward(const MaterialModel& m, const GaussState& committed,
                       const Eigen::Vector3d& de, GaussState& trial,
                       Eigen::Matrix3d* tangent_out = nullptr, int nsub_fixed = 0,
                       bool* plastic_out = nullptr, int* nsub_out = nullptr) {
    HardeningSoilParams pe = hs_small_strain_params(m.hs, committed.gamma_hist);
    // Dilatancy cut-off: psi = 0 clamps the mobilised dilatancy sin(psi_m) to [0, 0] inside the
    // return core, which IS Eq. 5.16b -- the rule enters where the manual puts it, and nothing
    // else in the HS machinery has to know about void ratios. The flag carries the SAME fact a
    // second way because for HSsmall the two stopped being equivalent: sec. 7.9.1 reads psi only
    // through phi_cv, so zeroing psi alone would move phi_cv up to phi and switch the Li &
    // Dafalias contraction on everywhere below failure -- a "stop dilating" option that starts
    // producing volume loss. The cut-off's own words are that psi_m "is automatically set back
    // to zero", so it is set back to zero.
    if (dilatancy_cut(m, committed)) { pe.dilatancy = 0.0; pe.dilatancy_cut = true; }
    const double Eur = hs_frozen_Eur(pe, committed.stress, committed.stress_zz);
    const LameConstants lame_ur = lame_from(Eur, pe.nu_ur);
    const PlaneStrainStress comm{committed.stress, committed.stress_zz};
    const PlaneStrainStress pred = elastic_predictor(comm, de, lame_ur);

    const HsReturnCore c = hs_return_core(pe, Eur, committed.stress, committed.stress_zz,
                                          pred.in_plane, pred.zz, committed.gamma_p,
                                          committed.pp, nsub_fixed, tension_cap_of(m));
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.gamma_p = c.gamma_p;
    trial.pp = c.pp;
    trial.eps_vol = committed.eps_vol;
    if (tangent_out) *tangent_out = c.tan.tangent;
    if (plastic_out) *plastic_out = c.plastic;
    if (nsub_out) *nsub_out = c.nsub;

    // HSsmall: γ_hist += Δγ (monotone accumulation; Eq 7-5 γ=√(1.5 e:e), e=deviatoric; plane strain εzz=0).
    if (m.hs.G0_ref > 0.0) {
        const double tr = de(0) + de(1), em = tr / 3.0;
        const double exx = de(0) - em, eyy = de(1) - em, ezz = -em, exy = 0.5 * de(2);
        const double ee = exx * exx + eyy * eyy + ezz * ezz + 2.0 * exy * exy;
        trial.gamma_hist = committed.gamma_hist + std::sqrt(1.5 * ee);
    } else {
        trial.gamma_hist = committed.gamma_hist;
    }
}

// FD step (relative): the sqrt(machine-epsilon) scale — the Pérez-Foguet &
// Rodríguez-Ferran & Huerta (CMAME 2000) recommendation. Signal/cancellation balance:
// Δσ ≈ E·h ~ 1e-3 kPa, double noise ~1e-14·|σ| → ~1e-8 relative error in the columns
// (ample for Newton).
inline double hs_fd_step(double comp) {
    constexpr double kSqrtEps = 1.4901161193847656e-08;  // sqrt(2^-52)
    return kSqrtEps * (1.0 + std::fabs(comp));
}

// ============================ SOFT SOIL FE wrapper ============================
// KINEMATICS-AGNOSTIC Soft Soil principal return (the SS counterpart of hs_return_core):
// the in-plane 2×2 Mohr block + the out-of-plane normal as the third principal — the SAME
// structure in both kinematics. The elastic law is SS's own exponential ln-law
// (p_tr = p_n·exp(Δε_v/κ*), s_tr = s_n + 2G(p_tr)·Δe) and, being isotropic, is built
// FRAME-INDEPENDENTLY in Voigt; after the principal decomposition, the principal strain
// increment fed to ss_step is the EXACT INVERSE of this law between committed→trial
// (Δε_v = κ*·ln(p_tr/p_n), Δe = Δs_dev/2G, rank-matched coaxial — the adaptation of HS's
// deps_p = C_e(σ_tr−σ_n) trick to the nonlinear law). The core thus reconstructs the trial
// to round-off: elastic steps are the identity map and the constitutive law runs from ONE
// SOURCE (ss_step). de4 = [Δεxx, Δεyy, Δγxy, Δε_out] (solver tension-positive; Δε_out = 0
// in plane strain, the REAL hoop strain in axisymmetry).
struct SsReturnCore {
    Eigen::Vector3d in_plane;   // returned [sxx,syy,sxy] / [srr,szz,srz] (tension-positive)
    double zz;                  // returned σ_zz / σ_θ
    double pp;                  // current preconsolidation pressure
    double K, G;                // elastic moduli frozen at the trial mean (tangent assembly)
    bool plastic;               // did cap/MC activate (SSC: MC or meaningful creep accumulated)
    int nsub;                   // substep count the core used (FD pinning)
};
// ssc != nullptr → the same spectral/inversion skeleton wraps the Soft Soil CREEP core
// (the elastic law is IDENTICAL to SS, so the inversion carries over exactly); dt_day is
// then the time increment. Spectral decomposition + rank matching + coaxial reconstruction
// are ONE source — so the corner/ulp-class errors measured in SS Stage 2 never live in two
// places again.
inline SsReturnCore ss_return_core(const softsoil::Params& P,
                                   const Eigen::Vector3d& comm_in_plane, double comm_zz,
                                   const Eigen::Vector4d& de4, double pp_n,
                                   int nsub_fixed = 0,
                                   const softsoilcreep::Params* ssc = nullptr,
                                   double dt_day = 0.0,
                                   double sigma_t_cap = kNoTensionCap) {
    const double cxx = comm_in_plane(0), cyy = comm_in_plane(1), cxy = comm_in_plane(2);

    // SS elastic trial in Voigt (tension-positive): exponentially EXACT mean + linear
    // deviatoric with G at the trial p — the SAME formula pair as ss_step's predictor
    // (guarantees the inversion below matches to round-off).
    const double mean_c = (cxx + cyy + comm_zz) / 3.0;
    const double p_c = std::max(-mean_c, softsoil::kPmin);            // compression-positive
    const double dv = -(de4(0) + de4(1) + de4(3));                    // compression-positive Δε_v
    const double p_tr = std::max(p_c * std::exp(dv / P.kap_star), softsoil::kPmin);
    const double K_tr = p_tr / P.kap_star;
    const double G = 1.5 * K_tr * (1.0 - 2.0 * P.nu_ur) / (1.0 + P.nu_ur);
    const double em = (de4(0) + de4(1) + de4(3)) / 3.0;
    const double sxx = -p_tr + (cxx - mean_c) + 2.0 * G * (de4(0) - em);
    const double syy = -p_tr + (cyy - mean_c) + 2.0 * G * (de4(1) - em);
    const double szz = -p_tr + (comm_zz - mean_c) + 2.0 * G * (de4(3) - em);
    const double sxy = cxy + G * de4(2);

    // Spectral decomposition of the trial (tension-positive) — the hs_return_core machinery verbatim.
    const double mean2 = 0.5 * (sxx + syy), hd = 0.5 * (sxx - syy);
    const double radius = std::sqrt(hd * hd + sxy * sxy);
    double cos2t = 1.0, sin2t = 0.0;
    if (radius > 0.0) { cos2t = hd / radius; sin2t = sxy / radius; }
    struct PV { double v; int src; };
    PV pv[3] = {{mean2 + radius, 0}, {mean2 - radius, 1}, {szz, 2}};
    std::sort(pv, pv + 3, [](const PV& a, const PV& b) { return a.v > b.v; });

    const double cmean = 0.5 * (cxx + cyy);
    const double cR = std::sqrt(0.25 * (cxx - cyy) * (cxx - cyy) + cxy * cxy);
    double cp[3] = {cmean + cR, cmean - cR, comm_zz};
    std::sort(cp, cp + 3, [](double a, double b) { return a > b; });
    const Eigen::Vector3d sig_n_cp(-cp[2], -cp[1], -cp[0]);      // committed, comp-pos descending
    const Eigen::Vector3d sig_tr_cp(-pv[2].v, -pv[1].v, -pv[0].v);  // trial, comp-pos descending

    // Principal Δε = the exact inverse of the elastic law between the two states (the
    // principal mean = the tensor mean, so the p's match the ones above to round-off).
    const double p_n_p = std::max(sig_n_cp.mean(), softsoil::kPmin);
    const double p_tr_p = std::max(sig_tr_cp.mean(), softsoil::kPmin);
    const double dv_p = P.kap_star * std::log(p_tr_p / p_n_p);
    const Eigen::Vector3d ddev =
        ((sig_tr_cp - Eigen::Vector3d::Constant(sig_tr_cp.mean())) -
         (sig_n_cp - Eigen::Vector3d::Constant(sig_n_cp.mean()))) / (2.0 * G);
    const Eigen::Vector3d deps_p = ddev + Eigen::Vector3d::Constant(dv_p / 3.0);

    Eigen::Vector3d ret_sig;
    double ret_pp;
    bool ret_plastic;
    int ret_nsub;
    if (ssc) {
        const softsoilcreep::StepResult rc =
            softsoilcreep::ssc_step(*ssc, sig_n_cp, pp_n, deps_p, dt_day, nsub_fixed);
        ret_sig = rc.sig;
        ret_pp = rc.pp;
        ret_plastic = rc.mc_active || rc.devc > 1e-12;
        ret_nsub = rc.nsub;
    } else {
        const softsoil::StepResult ret = softsoil::ss_step(P, sig_n_cp, pp_n, deps_p, nsub_fixed);
        ret_sig = ret.sig;
        ret_pp = ret.pp;
        ret_plastic = ret.cap_active || ret.mc_active;
        ret_nsub = ret.nsub;
    }
    double r[3] = {-ret_sig(2), -ret_sig(1), -ret_sig(0)};   // tension-positive, descending rank
    // Tension cut-off (MMM Eq. 3-11); see the Hardening Soil branch. K_tr and G are this step's
    // secant moduli, so the cap is returned on the same elasticity the step was taken with.
    const bool ss_capped = apply_rankine_cap(r, sigma_t_cap, K_tr - 2.0 * G / 3.0, G);

    SsReturnCore out;
    double pa = 0.0, pb = 0.0, pz = 0.0;   // coaxial reconstruction (in the trial eigenframe)
    for (int i = 0; i < 3; ++i) {
        if (pv[i].src == 0) pa = r[i];
        else if (pv[i].src == 1) pb = r[i];
        else pz = r[i];
    }
    const double m2 = 0.5 * (pa + pb), r2 = 0.5 * (pa - pb);
    out.in_plane(0) = m2 + r2 * cos2t;
    out.in_plane(1) = m2 - r2 * cos2t;
    out.in_plane(2) = r2 * sin2t;
    out.zz = pz;
    out.pp = ret_pp;
    out.K = K_tr;
    out.G = G;
    out.plastic = ret_plastic || ss_capped;
    out.nsub = ret_nsub;
    return out;
}

// Soft Soil FE forward update (PLANE STRAIN, Δε_out = 0). Elastic tangent = the (K_tr, G)
// isotropic operator — the EXACT continuum tangent of the exponential mean at the trial
// point is K = p_tr/κ* (the dG/dp cross term is dropped; frozen-modulus, the HS
// tradition). The caller builds the plastic tangent by FD (ss_fd_step).
inline void ss_forward(const MaterialModel& m, const GaussState& committed,
                       const Eigen::Vector3d& de, GaussState& trial,
                       bool* plastic_out = nullptr, double* K_out = nullptr,
                       double* G_out = nullptr, int nsub_fixed = 0, int* nsub_out = nullptr) {
    Eigen::Vector4d de4;
    de4 << de(0), de(1), de(2), 0.0;
    const SsReturnCore c = ss_return_core(m.ssoil, committed.stress, committed.stress_zz, de4,
                                          committed.pp, nsub_fixed, nullptr, 0.0,
                                          tension_cap_of(m));
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.pp = c.pp;
    trial.gamma_p = committed.gamma_p;
    trial.gamma_hist = committed.gamma_hist;
    trial.eps_vol = committed.eps_vol;
    if (plastic_out) *plastic_out = c.plastic;
    if (K_out) *K_out = c.K;
    if (G_out) *G_out = c.G;
    if (nsub_out) *nsub_out = c.nsub;
}

// Soft Soil CREEP FE forward update (PLANE STRAIN): the time-dependent twin of ss_forward —
// the same frame-independent Voigt trial + spectral skeleton (ss_return_core, ssc branch),
// core ssc_step (backward-Euler creep + single-source MC). dt_day = this increment's time
// share [days] (0 → no creep: elastic + MC).
inline void ssc_forward(const MaterialModel& m, const GaussState& committed,
                        const Eigen::Vector3d& de, double dt_day, GaussState& trial,
                        bool* plastic_out = nullptr, double* K_out = nullptr,
                        double* G_out = nullptr, int nsub_fixed = 0, int* nsub_out = nullptr) {
    Eigen::Vector4d de4;
    de4 << de(0), de(1), de(2), 0.0;
    const softsoil::Params S = m.ssc.ss();
    const SsReturnCore c = ss_return_core(S, committed.stress, committed.stress_zz, de4,
                                          committed.pp, nsub_fixed, &m.ssc, dt_day,
                                          tension_cap_of(m));
    trial.stress = c.in_plane;
    trial.stress_zz = c.zz;
    trial.pp = c.pp;
    trial.gamma_p = committed.gamma_p;
    trial.gamma_hist = committed.gamma_hist;
    trial.eps_vol = committed.eps_vol;
    if (plastic_out) *plastic_out = c.plastic;
    if (K_out) *K_out = c.K;
    if (G_out) *G_out = c.G;
    if (nsub_out) *nsub_out = c.nsub;
}

// SS FD step: the noise floor in the core's stress output is the scalar-Newton tolerance
// (1e-11·pp ≈ 1e-9 kPa; the surface gradients are ANALYTIC — soft_soil.hpp). h = 1e-6
// strain carries the signal to ~K·1e-6 ≈ 2.5e-3..2e-2 kPa (≫ the floor), truncation error
// relative ~(h/2)·(1/κ*) ≈ 2.5e-5 — ample for Newton; this larger step (instead of
// sqrt(eps)) is a deliberate choice because of the exponential law's strong curvature
// (K's derivative scales with 1/κ*).
inline double ss_fd_step(double comp) { return 1e-6 * (1.0 + std::fabs(comp)); }

// Tangent request mode (affects Hardening Soil only; LE/MC are always analytic):
//   kNone       — no tangent wanted (the line-search residual path): the costly tangent
//                 work is skipped.
//   kContinuum  — fast path: hs_integrate's last-active-set continuum tangent (linear
//                 global convergence; ~3-4 iterations per step suffice in the easy/
//                 confined regime — measured: performance-baseline.md B4).
//   kConsistent — NUMERICAL CONSISTENT (algorithmic) tangent: the forward-difference
//                 derivative of the fully substepped update (active set + plateau + drift
//                 included) w.r.t. Δε, D(:,j) = [σ(Δε+h·e_j) − σ(Δε)]/h. In the hard
//                 regime (footing edge, low confining pressure) it converges increments
//                 continuum cannot (GUI HS footing load_factor 0.844→1.000 — measured).
//                 The perturbed runs are pinned to the base run's SUBSTEP COUNT
//                 (nsub_fixed) — left free, the ceil() jump pollutes the FD column by the
//                 integration error. On an elastic step the map is linear → FD is
//                 unnecessary and not run.
//                 Source: Pérez-Foguet, Rodríguez-Ferran & Huerta, CMAME 189 (2000)
//                 277-296 + CMAME 190 (2001) 4627-4647; hardening-soil-formulation.md §9.
// The solver uses a hybrid: an increment starts with continuum; if it fails to converge,
// the same increment is retried with consistent (nonlinear_solver.cpp).
enum class TangentMode { kNone, kContinuum, kConsistent };

// Material-point integration. committed: converged state at the start of the
// step; strain_increment: the step's total strain increment Delta-eps; trial:
// updated (output) state; tangent: the 3x3 in-plane operator D_T fed back to the
// element (meaningless while mode==kNone, do not use it).
// dt_day: this increment's TIME share [days] — read only by time-dependent constitutive
// models (SoftSoilCreep); 0 (default) is bit-for-bit the old behaviour for all old callers.
inline void integrate_point(const MaterialModel& m, const GaussState& committed,
                            const Eigen::Vector3d& strain_increment,
                            GaussState& trial, Eigen::Matrix3d& tangent,
                            TangentMode mode = TangentMode::kConsistent, double dt_day = 0.0) {
    const LameConstants lame = lame_from(m.youngs_modulus, m.poisson_ratio);
    const PlaneStrainStress previous{committed.stress, committed.stress_zz};
    const PlaneStrainStress predictor =
        elastic_predictor(previous, strain_increment, lame);

    switch (m.type) {
        case MaterialType::LinearElastic: {
            trial.stress = predictor.in_plane;
            trial.stress_zz = predictor.zz;
            tangent = m.elastic_plane_strain();
            break;
        }
        case MaterialType::MohrCoulomb: {
            const MohrCoulombParams params{m.youngs_modulus, m.poisson_ratio,
                                           m.cohesion, m.friction_angle,
                                           effective_dilatancy(m, committed), m.tension_cutoff,
                                           m.tensile_strength};
            const McReturn base = mc_return_mapping(predictor, params);
            trial.stress = base.stress.in_plane;
            trial.stress_zz = base.stress.zz;
            // Elastic step -> elastic operator; plastic step -> the closed-form
            // consistent (algorithmic) tangent assembled inside the return mapping
            // (the active region's Jacobian carried through the eigen-decomposition,
            // incl. in-plane frame rotation and the sigma_zz coupling). A region-
            // consistent analytic tangent is essential for phi > 0: a forward
            // finite difference straddles region boundaries near the edges
            // (sigma_2 = sigma_3, ubiquitous under gravity) and breaks global
            // Newton convergence.
            tangent = base.plastic ? base.tangent : m.elastic_plane_strain();
            break;
        }
        case MaterialType::HardeningSoil: {
            bool plastic = false;
            int nsub = 0;
            hs_forward(m, committed, strain_increment, trial, &tangent, 0, &plastic, &nsub);
            // kConsistent + plastic step → numerical consistent tangent (TangentMode
            // block): 3 perturbed forward runs, substep count pinned to the base run.
            if (mode == TangentMode::kConsistent && plastic) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j));
                    dep(j) += h;
                    hs_forward(m, committed, dep, pert, nullptr, nsub);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            break;
        }
        case MaterialType::SoftSoil: {
            bool plastic = false;
            double K = 0.0, G = 0.0;
            int nsub = 0;
            ss_forward(m, committed, strain_increment, trial, &plastic, &K, &G, 0, &nsub);
            const double lam = K - 2.0 * G / 3.0;
            tangent << lam + 2.0 * G, lam, 0.0,
                       lam, lam + 2.0 * G, 0.0,
                       0.0, 0.0, G;
            // Plastic step + tangent wanted → NUMERICAL (forward-difference) tangent —
            // unlike HS, in BOTH modes (kContinuum/kConsistent): SS has no "cheap
            // continuum", the elastic operator is λ*/κ* times too stiff on the cap and
            // breaks Newton. kNone: the tangent is unused, skip (the elastic operator
            // stays, harmless). The perturbed runs are pinned to the base run's SUBSTEP
            // count (nsub — the HS lesson).
            if (plastic && mode != TangentMode::kNone) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    ss_forward(m, committed, dep, pert, nullptr, nullptr, nullptr, nsub);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            break;
        }
        case MaterialType::SoftSoilCreep: {
            // Same skeleton as the SS branch; time enters via dt_day and the FD columns run
            // with the SAME dt + SAME nsub (the tangent ∂σ/∂ε is taken at fixed time —
            // which is exactly what Newton needs).
            bool plastic = false;
            double K = 0.0, G = 0.0;
            int nsub = 0;
            ssc_forward(m, committed, strain_increment, dt_day, trial, &plastic, &K, &G, 0, &nsub);
            const double lam = K - 2.0 * G / 3.0;
            tangent << lam + 2.0 * G, lam, 0.0,
                       lam, lam + 2.0 * G, 0.0,
                       0.0, 0.0, G;
            if (plastic && mode != TangentMode::kNone) {
                GaussState pert;
                for (int j = 0; j < 3; ++j) {
                    Eigen::Vector3d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    ssc_forward(m, committed, dep, dt_day, pert, nullptr, nullptr, nullptr, nsub);
                    tangent.col(j) = (pert.stress - trial.stress) / h;
                }
            }
            break;
        }
    }
}

// Axisymmetric material-point integration. Strain/stress order [r, z, rz, theta];
// the hoop is a real strain (unlike plane strain). The Mohr-Coulomb return mapping
// is reused verbatim -- the axisymmetric stress (r-z block + hoop sigma_theta) has
// the same principal structure as plane strain (in-plane block + sigma_zz). The
// consistent 4x4 tangent is D_T = Psi * D_e, where Psi = d(sigma^ret)/d(sigma^tr)
// (McReturn::algo_jacobian) and D_e is the axisymmetric elastic matrix.
inline void integrate_point_axisym(const MaterialModel& m,
                                   const GaussState& committed,
                                   const Eigen::Vector4d& strain_increment,
                                   GaussState& trial, Eigen::Matrix4d& tangent,
                                   TangentMode mode = TangentMode::kConsistent,
                                   double dt_day = 0.0) {
    const Eigen::Matrix4d De = m.elastic_axisym();
    Eigen::Vector4d s_n;
    s_n << committed.stress, committed.stress_zz;       // [r, z, rz, theta]
    const Eigen::Vector4d s_tr = s_n + De * strain_increment;  // elastic predictor

    switch (m.type) {
        case MaterialType::LinearElastic: {
            trial.stress = s_tr.head<3>();
            trial.stress_zz = s_tr(3);
            tangent = De;
            break;
        }
        case MaterialType::MohrCoulomb: {
            PlaneStrainStress predictor;
            predictor.in_plane = s_tr.head<3>();
            predictor.zz = s_tr(3);
            const MohrCoulombParams params{m.youngs_modulus, m.poisson_ratio,
                                           m.cohesion, m.friction_angle,
                                           effective_dilatancy(m, committed), m.tension_cutoff,
                                           m.tensile_strength};
            const McReturn base = mc_return_mapping(predictor, params);
            trial.stress = base.stress.in_plane;
            trial.stress_zz = base.stress.zz;
            tangent = base.plastic ? (base.algo_jacobian * De) : De;
            break;
        }
        case MaterialType::HardeningSoil: {
            // Axisymmetric HS: the (r,z) block is the in-plane Mohr circle and the hoop
            // sigma_theta is the third principal -- the SAME principal structure as plane
            // strain, so hs_return_core is reused verbatim. Only the elastic predictor /
            // operator differ: the trial uses the AXISYMMETRIC elastic matrix at the frozen
            // Eur (NOT m.youngs_modulus; the hoop is a real strain), and the consistent 4x4
            // tangent is D_T = algo_jacobian * D_e_axisym (Psi from the shared principal
            // assembly), exactly mirroring Mohr-Coulomb.
            HardeningSoilParams pe = hs_small_strain_params(m.hs, committed.gamma_hist);
    // Dilatancy cut-off: psi = 0 clamps the mobilised dilatancy sin(psi_m) to [0, 0] inside the
    // return core, which IS Eq. 5.16b -- the rule enters where the manual puts it, and nothing
    // else in the HS machinery has to know about void ratios. The flag carries the SAME fact a
    // second way because for HSsmall the two stopped being equivalent: sec. 7.9.1 reads psi only
    // through phi_cv, so zeroing psi alone would move phi_cv up to phi and switch the Li &
    // Dafalias contraction on everywhere below failure -- a "stop dilating" option that starts
    // producing volume loss. The cut-off's own words are that psi_m "is automatically set back
    // to zero", so it is set back to zero.
    if (dilatancy_cut(m, committed)) { pe.dilatancy = 0.0; pe.dilatancy_cut = true; }
            const double Eur = hs_frozen_Eur(pe, committed.stress, committed.stress_zz);
            const double nu = pe.nu_ur;
            const double f = Eur / ((1.0 + nu) * (1.0 - 2.0 * nu));
            Eigen::Matrix4d De_ur = Eigen::Matrix4d::Zero();
            De_ur(0, 0) = De_ur(1, 1) = De_ur(3, 3) = f * (1.0 - nu);
            De_ur(0, 1) = De_ur(1, 0) = De_ur(0, 3) = De_ur(3, 0) = De_ur(1, 3) =
                De_ur(3, 1) = f * nu;
            De_ur(2, 2) = f * (1.0 - 2.0 * nu) / 2.0;
            const Eigen::Vector4d s_tr_ur = s_n + De_ur * strain_increment;
            const HsReturnCore c = hs_return_core(pe, Eur, committed.stress, committed.stress_zz,
                                                  s_tr_ur.head<3>(), s_tr_ur(3),
                                                  committed.gamma_p, committed.pp, 0,
                                                  tension_cap_of(m));
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.gamma_p = c.gamma_p;
            trial.pp = c.pp;
            trial.eps_vol = committed.eps_vol;
            tangent = c.tan.algo_jacobian * De_ur;  // continuum 4x4 (Psi * D_e)
            // kConsistent + plastic step → numerical consistent 4x4 tangent (same rationale
            // as the plane-strain block; the perturbed runs are pinned to the base nsub).
            if (mode == TangentMode::kConsistent && c.plastic) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = hs_fd_step(strain_increment(j));
                    dep(j) += h;
                    const Eigen::Vector4d s_tr_p = s_n + De_ur * dep;
                    const HsReturnCore cp = hs_return_core(
                        pe, Eur, committed.stress, committed.stress_zz,
                        s_tr_p.head<3>(), s_tr_p(3), committed.gamma_p, committed.pp,
                        c.nsub, tension_cap_of(m));
                    tangent(0, j) = (cp.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cp.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cp.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cp.zz - c.zz) / h;
                }
            }
            // HSsmall: gamma_hist += dgamma (axisym deviator includes the real hoop strain).
            if (m.hs.G0_ref > 0.0) {
                const double tr = strain_increment(0) + strain_increment(1) + strain_increment(3);
                const double em = tr / 3.0;
                const double err = strain_increment(0) - em, ez = strain_increment(1) - em,
                             eth = strain_increment(3) - em, erz = 0.5 * strain_increment(2);
                const double ee = err * err + ez * ez + eth * eth + 2.0 * erz * erz;
                trial.gamma_hist = committed.gamma_hist + std::sqrt(1.5 * ee);
            } else {
                trial.gamma_hist = committed.gamma_hist;
            }
            break;
        }
        case MaterialType::SoftSoil: {
            // Axisymmetric SS: the (r,z) block is the in-plane Mohr circle, the hoop σ_θ is
            // the third principal — ss_return_core takes de4 DIRECTLY as [Δεr, Δεz, Δγrz,
            // Δεθ] (the hoop is a real strain; the volumetric term contains all three).
            // Tangent: elastic (K_tr,G) isotropic 4×4; plastic by FD, same rationale as the
            // plane-strain branch.
            const SsReturnCore c = ss_return_core(m.ssoil, committed.stress, committed.stress_zz,
                                                  strain_increment, committed.pp, 0, nullptr,
                                                  0.0, tension_cap_of(m));
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.pp = c.pp;
            trial.gamma_p = committed.gamma_p;
            trial.gamma_hist = committed.gamma_hist;
            trial.eps_vol = committed.eps_vol;
            const double lam = c.K - 2.0 * c.G / 3.0;
            tangent = Eigen::Matrix4d::Zero();
            tangent(0, 0) = tangent(1, 1) = tangent(3, 3) = lam + 2.0 * c.G;
            tangent(0, 1) = tangent(1, 0) = tangent(0, 3) = tangent(3, 0) = tangent(1, 3) =
                tangent(3, 1) = lam;
            tangent(2, 2) = c.G;
            if (c.plastic && mode != TangentMode::kNone) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    const SsReturnCore cq =
                        ss_return_core(m.ssoil, committed.stress, committed.stress_zz, dep,
                                       committed.pp, c.nsub, nullptr, 0.0, tension_cap_of(m));
                    tangent(0, j) = (cq.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cq.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cq.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cq.zz - c.zz) / h;
                }
            }
            break;
        }
        case MaterialType::SoftSoilCreep: {
            // Same skeleton as the SS axisymmetric branch + time (dt_day); FD columns with
            // the same dt + the base nsub.
            const softsoil::Params S = m.ssc.ss();
            const SsReturnCore c = ss_return_core(S, committed.stress, committed.stress_zz,
                                                  strain_increment, committed.pp, 0, &m.ssc,
                                                  dt_day, tension_cap_of(m));
            trial.stress = c.in_plane;
            trial.stress_zz = c.zz;
            trial.pp = c.pp;
            trial.gamma_p = committed.gamma_p;
            trial.gamma_hist = committed.gamma_hist;
            trial.eps_vol = committed.eps_vol;
            const double lam = c.K - 2.0 * c.G / 3.0;
            tangent = Eigen::Matrix4d::Zero();
            tangent(0, 0) = tangent(1, 1) = tangent(3, 3) = lam + 2.0 * c.G;
            tangent(0, 1) = tangent(1, 0) = tangent(0, 3) = tangent(3, 0) = tangent(1, 3) =
                tangent(3, 1) = lam;
            tangent(2, 2) = c.G;
            if (c.plastic && mode != TangentMode::kNone) {
                for (int j = 0; j < 4; ++j) {
                    Eigen::Vector4d dep = strain_increment;
                    const double h = ss_fd_step(strain_increment(j));
                    dep(j) += h;
                    const SsReturnCore cq =
                        ss_return_core(S, committed.stress, committed.stress_zz, dep,
                                       committed.pp, c.nsub, &m.ssc, dt_day,
                                       tension_cap_of(m));
                    tangent(0, j) = (cq.in_plane(0) - c.in_plane(0)) / h;
                    tangent(1, j) = (cq.in_plane(1) - c.in_plane(1)) / h;
                    tangent(2, j) = (cq.in_plane(2) - c.in_plane(2)) / h;
                    tangent(3, j) = (cq.zz - c.zz) / h;
                }
            }
            break;
        }
    }
}

} // namespace katai::core
