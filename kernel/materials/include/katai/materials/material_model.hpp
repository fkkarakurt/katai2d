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
#include <katai/materials/hoek_brown.hpp>
#include <katai/materials/mohr_coulomb.hpp>
#include <katai/materials/soft_soil.hpp>
#include <katai/materials/soft_soil_creep.hpp>

namespace katai::core {

// Depth-VARYING stiffness/strength profile (indexed by material id; the depth increments
// "E'_inc / c'_inc" + y_ref). The most common real-soil case is stiffness growing with depth — in
// seismics E(y) is directly the Vs profile itself:
//     E(y) = E_ref + E_inc·(y_ref − y)      (decreases ABOVE y_ref)
//     c(y) = c_ref + c_inc·(y_ref − y)      (Mohr-Coulomb strength only)
// Evaluated PER STRESS (Gauss) POINT. An element
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
    HoekBrown,
};

// State of an integration (Gauss) point. Besides the in-plane Voigt stress we
// carry the out-of-plane normal stress sigma_zz: under the plane-strain
// constraint eps_zz = 0 it is generally non-zero and enters the Mohr-Coulomb
// yield evaluation as a principal stress. (Hardening / plastic-strain history
// will extend this struct in later phases.)
struct GaussState {
    Eigen::Vector3d stress = Eigen::Vector3d::Zero();  // Voigt [sxx, syy, sxy]
    double stress_zz = 0.0;                            // sigma_zz
    // Accumulated (committed) volumetric strain eps_v = eps_xx+eps_yy(+eps_theta), in every
    // phase and for every material: the void ratio the dilatancy cut-off reads comes from it.
    double eps_vol = 0.0;
    // Hardening Soil history: shear hardening parameter gamma_p and cap pre-
    // consolidation pressure pp (compression-positive). Zero / unused for other models.
    double gamma_p = 0.0;
    double pp = 0.0;
    // HSsmall: accumulated deviatoric shear strain γ_hist (monotone; drives the
    // small-strain stiffness degradation). Only for HSsmall (G0_ref>0); 0/unused otherwise.
    double gamma_hist = 0.0;
    // THE EXCESS PORE PRESSURE IS A STATE OF THE PHASE CHAIN, not a function of the volume
    // change. Its stress (tension positive, the sign it enters the total stress with) at this
    // point is
    //     pw_carried + (Kw/n) eps_vol_und
    // eps_vol_und is the volumetric strain accumulated while the point GENERATED excess pore
    // pressure -- an Undrained (A)/(B) material in a phase that does not ignore undrained
    // behaviour -- and pw_carried is a pressure handed over by a phase that computes the pressure
    // itself (consolidation, fully coupled). It used to be (Kw/n) eps_vol, the volume change since
    // the initial state whatever produced it, so a drained stage or a consolidation came back in
    // the next undrained phase as a pore pressure nobody had generated: measured on a confined
    // column, a phase that changed nothing heaved 28.7 mm. eps_vol_und accumulates by the same
    // expression eps_vol does, so a chain in which every phase generates is bit-for-bit what it
    // was. (docs/references/effective-stress-formulation.md)
    double eps_vol_und = 0.0;
    double pw_carried = 0.0;
};

// Isotropic elastic constitutive matrix from an (E, nu) PAIR, rather than from the pair a
// material happens to keep in its youngs_modulus/poisson_ratio boxes. The Hardening Soil
// family's elasticity is (Eur, nu_ur) and it never reads those boxes, so every place that
// needs "the elastic operator of THIS point" has to be told which pair to use. Written once
// here: the members below and the axisymmetric HS branch used to carry their own copy of the
// same six lines, and a stiffness formula that exists in three copies is a formula that can
// be right in two of them.
inline Eigen::Matrix3d elastic_plane_strain_of(double e, double v) {
    const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
    Eigen::Matrix3d d = Eigen::Matrix3d::Zero();
    d(0, 0) = f * (1.0 - v);
    d(0, 1) = f * v;
    d(1, 0) = f * v;
    d(1, 1) = f * (1.0 - v);
    d(2, 2) = f * (1.0 - 2.0 * v) / 2.0;
    return d;
}

// Axisymmetric twin, strain/stress order [r, z, rz, theta].
inline Eigen::Matrix4d elastic_axisym_of(double e, double v) {
    const double f = e / ((1.0 + v) * (1.0 - 2.0 * v));
    Eigen::Matrix4d d = Eigen::Matrix4d::Zero();
    d(0, 0) = d(1, 1) = d(3, 3) = f * (1.0 - v);
    d(0, 1) = d(1, 0) = d(0, 3) = d(3, 0) = d(1, 3) = d(3, 1) = f * v;
    d(2, 2) = f * (1.0 - 2.0 * v) / 2.0;
    return d;
}

// Isotropic elastic matrix from a bulk/shear PAIR (the soft-soil integrators report K and G
// rather than E and nu, because their stiffness is the ln-law's and there is no constant E to
// report). Same operator, entered through the other door.
inline Eigen::Matrix3d elastic_plane_strain_kg(double K, double G) {
    const double lam = K - 2.0 * G / 3.0;
    Eigen::Matrix3d d = Eigen::Matrix3d::Zero();
    d(0, 0) = d(1, 1) = lam + 2.0 * G;
    d(0, 1) = d(1, 0) = lam;
    d(2, 2) = G;
    return d;
}

inline Eigen::Matrix4d elastic_axisym_kg(double K, double G) {
    const double lam = K - 2.0 * G / 3.0;
    Eigen::Matrix4d d = Eigen::Matrix4d::Zero();
    d(0, 0) = d(1, 1) = d(3, 3) = lam + 2.0 * G;
    d(0, 1) = d(1, 0) = d(0, 3) = d(3, 0) = d(1, 3) = d(3, 1) = lam;
    d(2, 2) = G;
    return d;
}

// Flat, tagged material descriptor (no vtable -> suitable for the hot loop).
struct MaterialModel {
    MaterialType type = MaterialType::LinearElastic;
    double youngs_modulus = 0.0;   // E
    double poisson_ratio = 0.0;    // v
    double cohesion = 0.0;         // c        (Mohr-Coulomb)
    double friction_angle = 0.0;   // phi [rad] (Mohr-Coulomb)
    double dilatancy_angle = 0.0;  // psi [rad] (Mohr-Coulomb)
    // Rankine tension cap (sigma_1 <= sigma_t, associated flow).
    // Read by the MohrCoulomb branch of integrate_point / integrate_point_axisym;
    // HS/SS do not consume it yet (their integrators lack the extra planes -- the
    // GUI states this honestly). Default OFF keeps every direct caller bit-identical.
    bool tension_cutoff = false;
    double tensile_strength = 0.0;  // sigma_t [kN/m2], tension-positive
    // Dilatancy cut-off: after extensive shearing a dilating material arrives at a state of
    // critical density where dilatancy has come to an end. As soon as the volume change takes the
    // soil to its maximum void ratio e_max, the mobilised dilatancy angle is set back to zero.
    // Without it a dense sand dilates for ever and its bearing capacity is over-predicted -- an
    // unsafe number, produced quietly.
    //
    // The void ratio follows the volume change: 1 + e = (1 + e_init) exp(eps_v), expansion
    // positive. The same statement can be written eps_v = ln((1+e)/(1+e_init)), where the sign
    // convention is easy to get wrong; the exponential form says it once.
    // e_min is deliberately NOT here: the cut-off does not use a minimum void ratio, so storing
    // one would suggest a rule that does not exist.
    bool dilatancy_cutoff = false;
    double e_init = 0.5, e_max = 1.0;

    // Undrained (A): effective parameters above (E', nu', c', phi'); the pore fluid's
    // volumetric stiffness Kw/n is added to the GLOBAL tangent (D_u = D' + (Kw/n)mm^T)
    // and the internal force uses TOTAL stress = sigma' + (Kw/n) eps_v m, while the
    // constitutive model still works on effective stress. The solver tracks the excess
    // pore pressure in GaussState::eps_vol_und. undrained_poisson = nu_u (0.495 by default, a
    // nearly incompressible undrained Poisson's ratio; exactly 0.5 is singular).
    bool undrained = false;
    // The phase ignores undrained behaviour (`ignoreund`): the water's stiffness is not in the
    // tangent and no NEW excess pore pressure is generated, but the pressure generated before
    // stays in the total stress: a phase with no time and no flow has nothing that could dissipate
    // it. It is a flag of its own because clearing `undrained` for the phase used to
    // do both jobs at once and a third by accident: the existing pressure vanished, and an
    // Undrained (B) soil lost the undrained partial factor its cohesion is owed (design_code.hpp).
    bool ignore_undrained = false;
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

    // Undrained (C): the material is analysed in TOTAL stress. The stiffness
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

    // Soft Soil parameters (used when type == SoftSoil; soft_soil.hpp /
    // docs/references/soft-soil-formulation.md). youngs_modulus/poisson_ratio are NOT USED —
    // the stiffness comes from the ln-law (K = p'/κ*). CAUTION (honest Stage-2 limit):
    // SS + undrained (A) is NOT wired yet — kw_over_n() derives from E and would silently
    // be 0 with E=0; Stage 3 will add undrained support via the stress-dependent K_ur. The
    // GUI/project file does not produce this model until Stage 3; the kernel tests use it
    // drained.
    softsoil::Params ssoil;

    // Soft Soil Creep parameters (type == SoftSoilCreep; Vermeer & Neher 1999,
    // soft_soil_creep.hpp / docs/references/soft-soil-creep-formulation.md). Time enters
    // via integrate_point's trailing parameter dt_day (0 = no creep: SSC in a phase without
    // a time interval gives only elastic+MC). Undrained/Safety limits as SS.
    softsoilcreep::Params ssc;

    // Hoek-Brown parameters (type == HoekBrown; hoek_brown.hpp). The elastic part
    // IS youngs_modulus / poisson_ratio -- rock keeps Hooke's law, which is the whole reason the
    // model is only a strength criterion. hb.E / hb.nu are filled from those two at the seam so
    // the core stays self-contained.
    hoekbrown::Params hb;

    // Plane-strain elastic constitutive matrix (3x3 SPD, v < 0.5).
    Eigen::Matrix3d elastic_plane_strain() const {
        return elastic_plane_strain_of(youngs_modulus, poisson_ratio);
    }

    // Undrained (A) pore-fluid bulk stiffness Kw/n from the EFFECTIVE parameters and
    // an assumed undrained Poisson ratio (nu_u = 0.495 by default; exactly 0.5
    // makes the stiffness singular): Kw/n = 3 (nu_u - nu') / ((1 - 2 nu_u)(1 + nu')) K'.
    // Derived so that K' + Kw/n = the correct undrained bulk modulus Ku (see
    // effective-stress-formulation.md).
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
        return elastic_axisym_of(youngs_modulus, poisson_ratio);
    }
};

// --- The undrained stiffness trio (alpha_Biot = 1) ---------------------------------------------
// Three quantities describe the same pore fluid, and each has a closed form:
//   Kw/n  from nu_u        MaterialModel::kw_over_n above
//   nu_u  from Skempton B  undrained_poisson_from_skempton
//   B     from Kw/n        skempton_from_kw_over_n
// They are mutually consistent -- going round the ring returns the number it started from, which
// is what test_undrained_stiffness measures rather than assumes. Which one the USER supplies is
// an input choice (docs/k2d-format.md, `und_mode`: nu_u entered, or Skempton's B); the other two
// are then derived and are worth showing, because a B of 0.98 and a Kw/n of 45 K' are the same
// statement and an engineer recognises one of them.

// With alpha_Biot = 1: nu_u = (3 nu' + B (1 - 2 nu')) / (3 - B (1 - 2 nu')).
// B -> 1 gives exactly 0.5 (incompressible, singular), so the caller must keep B < 1.
inline double undrained_poisson_from_skempton(double B, double nu_eff) {
    const double t = B * (1.0 - 2.0 * nu_eff);
    return (3.0 * nu_eff + t) / (3.0 - t);
}

// With alpha_Biot = 1: B = 1 / (1 + n K'/Kw) (Skempton 1954), written on the quantity the
// engine actually carries -- B = (Kw/n) / (K' + Kw/n) = (Kw/n) / Ku. The porosity cancels
// because Kw and n only ever enter as the ratio, which is why Kw/n is the natural quantity.
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
// G0_ref=0 -> returns p unchanged (plain HS, byte-identical). (Benz 2007.)
//
// The threshold is the RELOADING one, 2*gamma07 (Masing 1926): what this function sets is
// the model's quasi-elastic (unload/reload) stiffness, and the model keeps that factor
// constant at 2 throughout loading rather than switching it on at a reversal. Riding the
// virgin backbone here instead would degrade the stiffness twice as fast as the model
// intends, and measurably overstate the heave on the unloading case KV-CST-008 -- and softer
// is not the safe side when the number being read is a wall deflection or a heave.
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
    // 1 when the substep guard clipped the subdivision, i.e. the integration did NOT meet its
    // tolerance on this call. Carried out of the integrator because a run that hit the ceiling
    // must never be reported as one that met its tolerance -- and until 2026-08-25 the finite
    // element path dropped this flag on the floor, so exactly that was happening.
    int saturated = 0;
};

// The material's own tension cut-off, as the cap the principal returns take. The schema
// switches the cut-off ON by default, so `off` here is a deliberate choice by the engineer and
// not a default nobody looked at.
inline double tension_cap_of(const MaterialModel& m) {
    return m.tension_cutoff ? m.tensile_strength : kNoTensionCap;
}

// plan_in: replay a previous run's accepted substep subdivision, so the numerical consistent
// tangent's perturbed runs walk EXACTLY the base run's substeps (see HsSubstepPlan). plan_out:
// record this run's subdivision for that purpose.
HsReturnCore hs_return_core(const HardeningSoilParams& pe, double Eur,
                            const Eigen::Vector3d& comm_in_plane, double comm_zz,
                            const Eigen::Vector3d& trial_in_plane, double trial_zz,
                            double gamma_p_n, double pp_n,
                            const HsSubstepPlan* plan_in = nullptr,
                            double sigma_t_cap = kNoTensionCap,
                            HsSubstepPlan* plan_out = nullptr,
                            // The constitutive integration error tolerance (STOL), from
                            // the phase (.k2d v15) or 0 for the material class's default.
                            double substep_tol = 0.0);

// Hardening Soil FE forward stress update (PLANE STRAIN). Frozen Eur -> plane-strain elastic
// predictor (deps_zz=0) -> kinematics-agnostic principal return -> analytic 3x3 tangent.
// nsub_fixed / plastic_out / nsub_out: the numerical consistent-tangent plumbing (see
// hs_consistent_tangent) — the base run reports its nsub, the perturbed runs pin it.
// Has the dilatancy cut-off been reached at this stress point (e >= e_max)? The state
// carries the accumulated volumetric strain, and the void ratio follows
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

void hs_forward(const MaterialModel& m, const GaussState& committed,
                const Eigen::Vector3d& de, GaussState& trial,
                Eigen::Matrix3d* tangent_out = nullptr,
                const HsSubstepPlan* plan_in = nullptr,
                bool* plastic_out = nullptr, HsSubstepPlan* plan_out = nullptr,
                double* Eur_out = nullptr, double substep_tol = 0.0,
                int* saturated_out = nullptr);

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
SsReturnCore ss_return_core(const softsoil::Params& P,
                            const Eigen::Vector3d& comm_in_plane, double comm_zz,
                            const Eigen::Vector4d& de4, double pp_n,
                            int nsub_fixed = 0,
                            const softsoilcreep::Params* ssc = nullptr,
                            double dt_day = 0.0,
                            double sigma_t_cap = kNoTensionCap);

// Soft Soil FE forward update (PLANE STRAIN, Δε_out = 0). Elastic tangent = the (K_tr, G)
// isotropic operator — the EXACT continuum tangent of the exponential mean at the trial
// point is K = p_tr/κ* (the dG/dp cross term is dropped; frozen-modulus, the HS
// tradition). The caller builds the plastic tangent by FD (ss_fd_step).
void ss_forward(const MaterialModel& m, const GaussState& committed,
                const Eigen::Vector3d& de, GaussState& trial,
                bool* plastic_out = nullptr, double* K_out = nullptr,
                double* G_out = nullptr, int nsub_fixed = 0, int* nsub_out = nullptr);

// Soft Soil CREEP FE forward update (PLANE STRAIN): the time-dependent twin of ss_forward —
// the same frame-independent Voigt trial + spectral skeleton (ss_return_core, ssc branch),
// core ssc_step (backward-Euler creep + single-source MC). dt_day = this increment's time
// share [days] (0 → no creep: elastic + MC).
void ssc_forward(const MaterialModel& m, const GaussState& committed,
                 const Eigen::Vector3d& de, double dt_day, GaussState& trial,
                 bool* plastic_out = nullptr, double* K_out = nullptr,
                 double* G_out = nullptr, int nsub_fixed = 0, int* nsub_out = nullptr);

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
// The cohesion this material's STRENGTH actually uses, whichever model carries it. The local
// convergence criteria normalise a stress difference by max(tau_max, c, ...), and reading
// MaterialModel::cohesion for every model would report zero for the two soft-soil models,
// which keep theirs in their own parameter blocks. A normaliser that is quietly a factor too
// small makes every point look accurate, which is the one failure mode a convergence check
// must not have.
inline double cohesion_of(const MaterialModel& m) {
    switch (m.type) {
        case MaterialType::HardeningSoil: return m.hs.cohesion;
        case MaterialType::SoftSoil:      return m.ssoil.c;
        case MaterialType::SoftSoilCreep: return m.ssc.c;
        // Rock keeps no cohesion at all -- its strength IS the curve -- so the default below
        // would hand the normaliser a zero, which is precisely the failure the paragraph above
        // describes, on a material whose strengths are megapascals. What the criterion has at
        // zero confinement is the rock mass's uni-axial compressive strength
        // sigma_c = sigma_ci s^a, and a compressive strength sits on a cohesion scale at half of
        // it -- the Tresca relation c = sigma_c / 2, which is the same degeneration KV-CST-014
        // checks this model against. It is a SCALE for a convergence check, not a parameter:
        // nothing else reads it.
        case MaterialType::HoekBrown: {
            hoekbrown::Params hp = m.hb;
            hp.E = m.youngs_modulus; hp.nu = m.poisson_ratio;
            return 0.5 * std::fabs(hoekbrown::constants_of(hp).sigc);
        }
        default:                          return m.cohesion;
    }
}

// The reference pressure the stress-dependent stiffness is defined at. Only the models with a
// stress-dependent elastic stiffness have one, and only those are counted by the non-linear
// elastic criterion; for anything else the value is never read.
inline double p_ref_of(const MaterialModel& m) {
    switch (m.type) {
        case MaterialType::HardeningSoil: return m.hs.p_ref;
        case MaterialType::SoftSoil:
        case MaterialType::SoftSoilCreep: return 100.0;  // the ln-law has no p_ref; 100 kPa
        default:                          return 100.0;
    }
}

// Maximum shear stress at a point: (sigma_1 - sigma_3)/2 over ALL THREE principals, not over
// the in-plane pair. The out-of-plane component is a principal stress in both kinematics
// (sigma_zz in plane strain, sigma_theta in axisymmetry) and under K0 conditions it is
// routinely the largest or the smallest of the three -- taking the in-plane circle alone
// would under-report tau_max exactly where the soil is closest to failure.
inline double tau_max_of(double sxx, double syy, double sxy, double szz) {
    const double mean = 0.5 * (sxx + syy);
    const double r = std::sqrt(0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy);
    const double s1 = std::max(mean + r, szz);
    const double s3 = std::min(mean - r, szz);
    return 0.5 * (s1 - s3);
}

enum class TangentMode { kNone, kContinuum, kConsistent };

// What the material law DID at one stress point, handed back to the caller. Two facts, both
// of them already computed inside every branch below and, until now, thrown away at the
// closing brace:
//
//   plastic          a yield or cap surface was active in this increment. The convergence
//                    machinery needs it because the local error criteria are counted over
//                    PLASTIC points -- an elastic point is accurate by construction.
//   elastic          D^e at this point, in the solver's Voigt order. Not the tangent: the
//                    ELASTIC operator, which is what an equilibrium stress is built with
//                    (sigma_eq = sigma_c,prev + D^e * delta-eps). For the Hardening Soil
//                    family this is (Eur, nu_ur) evaluated where the model evaluates it, not
//                    the (E, nu) boxes those models never read.
//   stress_dependent the elastic stiffness itself is a function of stress. Such a point is
//                    counted separately: it can be inaccurate while carrying no plasticity
//                    at all, which is the whole reason a second local count exists.
//
// Passing nullptr (the default) costs nothing and keeps every existing caller bit-for-bit.
template <class MatT>
struct PointReportT {
    bool plastic = false;
    bool stress_dependent = false;
    // The constitutive integration at this point hit its substep GUARD, so it did not meet the
    // error tolerance it was asked for. This is not a convergence question -- no equilibrium
    // residual can see it -- and it is the one state in which a "converged" answer rests on an
    // integration that was cut short.
    bool integration_saturated = false;
    MatT elastic = MatT::Zero();
};
using PointReport = PointReportT<Eigen::Matrix3d>;    // plane strain [xx, yy, xy]
using PointReport4 = PointReportT<Eigen::Matrix4d>;   // axisymmetric [r, z, rz, theta]

// Material-point integration. committed: converged state at the start of the
// step; strain_increment: the step's total strain increment Delta-eps; trial:
// updated (output) state; tangent: the 3x3 in-plane operator D_T fed back to the
// element (meaningless while mode==kNone, do not use it).
// dt_day: this increment's TIME share [days] — read only by time-dependent constitutive
// models (SoftSoilCreep); 0 (default) is bit-for-bit the old behaviour for all old callers.
void integrate_point(const MaterialModel& m, const GaussState& committed,
                     const Eigen::Vector3d& strain_increment,
                     GaussState& trial, Eigen::Matrix3d& tangent,
                     TangentMode mode = TangentMode::kConsistent, double dt_day = 0.0,
                     PointReport* report = nullptr, double substep_tol = 0.0);

// Axisymmetric material-point integration. Strain/stress order [r, z, rz, theta];
// the hoop is a real strain (unlike plane strain). The Mohr-Coulomb return mapping
// is reused verbatim -- the axisymmetric stress (r-z block + hoop sigma_theta) has
// the same principal structure as plane strain (in-plane block + sigma_zz). The
// consistent 4x4 tangent is D_T = Psi * D_e, where Psi = d(sigma^ret)/d(sigma^tr)
// (McReturn::algo_jacobian) and D_e is the axisymmetric elastic matrix.
void integrate_point_axisym(const MaterialModel& m,
                            const GaussState& committed,
                            const Eigen::Vector4d& strain_increment,
                            GaussState& trial, Eigen::Matrix4d& tangent,
                            TangentMode mode = TangentMode::kConsistent,
                            double dt_day = 0.0,
                            PointReport4* report = nullptr,
                            double substep_tol = 0.0);

} // namespace katai::core
