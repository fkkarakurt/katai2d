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
#include <katai/materials/ngi_adp.hpp>
#include <katai/materials/jointed_rock.hpp>

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
//
// ABOVE THE REFERENCE LEVEL THE VALUE IS THE REFERENCE VALUE. The increment is per unit of
// DEPTH below y_ref: for any elevation above it the stiffness (or the cohesion) is the reference
// one, and below it grows with the depth. It used to be extrapolated upwards as well, so a
// layer whose reference level lay inside it was softer and weaker above that level than the
// value entered for it -- and, with a large increment, had none at all at its top.
inline double profile_at(double ref, double inc, double y_ref, double y) {
    const double v = ref + inc * std::max(0.0, y_ref - y);
    return v > 0.0 ? v : 0.0;
}

enum class MaterialType {
    LinearElastic,
    MohrCoulomb,
    HardeningSoil,
    SoftSoil,
    SoftSoilCreep,
    HoekBrown,
    NgiAdp,
    JointedRock,
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
    // HSsmall: the scalar shear strain gamma_hist the small-strain stiffness is read at, as it
    // stood at the end of the last increment. Only for HSsmall (G0_ref>0); 0/unused otherwise.
    // It is derived from the history tensor below (hs_small_strain_step) and kept for output.
    double gamma_hist = 0.0;
    // HSsmall: the deviatoric strain HISTORY, a symmetric tensor [hxx, hyy, hzz, hxy] (tensor
    // shear; zz is the out-of-plane or hoop direction). A strain increment that reverses the
    // history in one of its own principal directions resets it there, which is what gives the
    // stiffness back on unloading and reloading (Benz 2007, Sec. 6.1 and Appendix A).
    Eigen::Vector4d hs_hist = Eigen::Vector4d::Zero();
    // HSsmall: the smallest stiffness multiplier G/G_ur the point has had on the PRIMARY
    // loading curve (0 = not loaded yet). The hardening laws are scaled with it
    // (HardeningSoilParams::hardening_scale); once it has come down to 1 the model hardens as
    // Hardening Soil does.
    double hs_g_mono = 0.0;
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
    // The point lies in a zone where soil plasticity is excluded -- the ground around the foot of
    // a plate that is given end bearing. It answers with the elastic part of its model and its
    // hardening state does not move. Set by the phase that owns the zone, phase by phase.
    bool plasticity_excluded = false;
    // NGI-ADP: the initial state the yield surface is translated from (ngi_adp.hpp), compression
    // positive, [xx, yy, zz, xy] -- the stress the point carried when the model first acted on
    // it -- and whether it has been taken. The model's plastic shear strain is gamma_p above.
    Eigen::Vector4d adp_s0 = Eigen::Vector4d::Zero();
    bool adp_set = false;
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
    // NGI-ADP (ngi_adp.hpp): the plane strain active strength is `cohesion` and the shear modulus
    // follows from youngs_modulus / poisson_ratio -- so the depth profile of both applies as it
    // does to any model -- and these are the ratios and the total failure strains [-].
    double adp_suP = 0.0, adp_suDSS = 0.0;     // suP / suA, suDSS / suA
    double adp_gfC = 0.0, adp_gfE = 0.0, adp_gfDSS = 0.0;
    // Jointed rock (jointed_rock.hpp). The FIRST joint set's strength is the shared cohesion,
    // friction, dilatancy and tension cut-off above -- so the depth profile of the cohesion, the
    // strength reduction of a Safety phase and the interface's Rinter reach it as they reach any
    // material's -- and its dip orients the anisotropy; E / nu are E1 / nu1, the plane of
    // isotropy's. jr holds the rest: the second and third sets whole, and E2, nu2, G2, which a
    // depth profile of E scales with E1 (jr_E1_ref is the E1 they were entered beside).
    jr::Params jr;
    double jr_E1_ref = 0.0;

    // The jointed rock's parameters at this point: E1 / nu1 and the first set from the shared
    // fields, E2 and G2 scaled with E1.
    jr::Params jr_at_point() const {
        jr::Params P = jr;
        P.E1 = youngs_modulus;
        P.nu1 = poisson_ratio;
        const double s = jr_E1_ref > 0.0 ? youngs_modulus / jr_E1_ref : 1.0;
        P.E2 = jr.E2 * s;
        P.G2 = jr.G2 * s;
        P.p[0].c = cohesion;
        P.p[0].phi = friction_angle;
        P.p[0].psi = dilatancy_angle;
        P.p[0].tcut = tension_cutoff;
        P.p[0].sigt = tensile_strength;
        return P;
    }

    // Plane-strain elastic constitutive matrix (3x3 SPD, v < 0.5). A transversely isotropic rock
    // answers with its own operator: every caller that asks a material for its elasticity
    // (consolidation, the coupled phases, the recovery of pore pressures) then sees the anisotropy.
    Eigen::Matrix3d elastic_plane_strain() const {
        if (type == MaterialType::JointedRock && jr.aniso) {
            const jr::Mat4 D = jr::elastic(jr_at_point());
            const int idx[3] = {0, 1, 3};
            Eigen::Matrix3d d;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) d(i, j) = D(idx[i], idx[j]);
            return d;
        }
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

    // THE PORE FLUID OF A SOIL WHOSE STIFFNESS FOLLOWS THE STRESS. Kw/n is a multiple of the
    // skeleton's bulk modulus K', and the soft-soil models have no constant one: K' =
    // (p' + c cot phi)/kappa*. kw_over_n() is then the value at the reference pressure below, and
    // this is the ratio of the current K' to that -- the factor by which the pore fluid at this
    // stress point is stiffer or softer than the reference one. 1 for every other model. The
    // excess pore pressure is accumulated as an EQUIVALENT strain (GaussState::eps_vol_und
    // advances by this factor times the volumetric strain increment), so that
    // u = kw_over_n() * eps_vol_und stays the one expression every reader of the state uses.
    static constexpr double kSoftSoilUndrainedRefPressure = 100.0;   // [kPa]
    double undrained_stiffness_scale(double p_eff_compression) const {
        if (type != MaterialType::SoftSoil && type != MaterialType::SoftSoilCreep) return 1.0;
        const softsoil::Params sp = type == MaterialType::SoftSoil ? ssoil : ssc.ss();
        const double shift = softsoil::cohesion_shift(sp);
        return (std::max(p_eff_compression, softsoil::kPmin) + shift) /
               (kSoftSoilUndrainedRefPressure + shift);
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
        if (type == MaterialType::JointedRock && jr.aniso) {
            // [r, z, rz, theta] from the core's [xx, yy, zz, xy]: the hoop direction is the
            // out-of-plane one, and it lies in every joint (a cone about the axis)
            const jr::Mat4 D = jr::elastic(jr_at_point());
            const int idx[4] = {0, 1, 3, 2};
            Eigen::Matrix4d d;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) d(i, j) = D(idx[i], idx[j]);
            return d;
        }
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
    return pe.Eur(s3_stiff) * pe.elastic_scale;   // the ELASTIC stiffness of the step
}

// THE PORE FLUID FOLLOWS THE SKELETON, at every stress point. Kw/n is a multiple of the
// skeleton's elastic bulk modulus -- that multiple is what makes the undrained Poisson's ratio
// 0.495 -- and in the hardening and soft-soil models that modulus follows the stress. This is the
// ratio of the current modulus to the reference one MaterialModel::kw_over_n() is built from:
// E_ur(sigma3)/E_ur,ref for Hardening Soil, (p' + c cot phi)/(p_ref + c cot phi) for Soft Soil,
// 1 for a model whose stiffness is a constant.
//
// Until 2026-10 a Hardening Soil material kept the reference value everywhere. Under a few
// metres of soft clay (m = 1) the skeleton is ten to a hundred times softer than at p_ref, so
// the pore fluid was that much stiffer against it than the model intends: an undrained Poisson's
// ratio of 0.4999 instead of 0.495. Measured on an embankment built undrained on soft clay, the
// factor of safety is 1.4020 with the fluid following the skeleton and 1.4019 without; it is
// the definition that changed, not an answer. (It was also tried as the cure for that phase
// stalling at a third of its load on a finer mesh. It is not: the stall is still there.)
inline double undrained_stiffness_scale(const MaterialModel& m, const Eigen::Vector3d& comm_in_plane,
                                        double comm_zz) {
    if (m.type == MaterialType::HardeningSoil && m.hs.Eur_ref > 0.0)
        return hs_frozen_Eur(m.hs, comm_in_plane, comm_zz) / m.hs.Eur_ref;
    const double p_eff = -(comm_in_plane(0) + comm_in_plane(1) + comm_zz) / 3.0;
    return m.undrained_stiffness_scale(p_eff);
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
        pe.elastic_scale = std::max(E0_ref / (d * d), p.Eur_ref) / p.Eur_ref;
    }
    return pe;
}

// THE SMALL-STRAIN STEP OF HSsmall: what one strain increment does to the strain history, and
// the stiffness it is taken with. `exx, eyy, ezz, exy` is the DEVIATORIC strain increment (tensor
// shear; zz the out-of-plane or hoop direction). Follows Benz (2007), Appendix A, statement by
// statement:
//   * the history tensor is rotated into the principal directions of the increment;
//   * in a direction where the history and the increment have opposite signs the history is
//     reset (the scaling 1/sqrt(1 + chi) of that row and column, which takes the elongation chi
//     there to zero and scales the shear components with it);
//   * the increment is added, and the scalar history strain is gamma = sqrt(3) |H de| / |de|
//     before (gamma_1) and after (gamma_2) -- in a triaxial test eps_axial - eps_lateral;
//   * the stiffness of the increment is the SECANT of the Hardin-Drnevich curve between the two,
//     G = G0 [gamma_2/(1 + gamma_2/F) - gamma_1/(1 + gamma_1/F)] / (gamma_2 - gamma_1) with
//     F = gamma_0.7 / 0.385 on the primary curve and twice that on the reloading curve (Masing),
//     never below G_ur (the cut-off at gamma = (sqrt(G0/Gur) - 1) F);
//   * the elastic stiffness is always the RELOADING one -- primary loading is softer through the
//     plasticity, not through the elasticity -- and the smallest value the PRIMARY curve has
//     given is remembered as the multiplier the hardening laws are scaled with.
// Until this was written the history was one scalar that only ever grew: a point unloaded after
// loading stayed on the stiffness the loading had degraded it to.
struct HsSmallStrainStep {
    HardeningSoilParams pe;       // the parameters this increment is integrated with
    Eigen::Vector4d hist = Eigen::Vector4d::Zero();
    double g_mono = 0.0;
    double gamma_hist = 0.0;
};
inline HsSmallStrainStep hs_small_strain_step(const HardeningSoilParams& p, const GaussState& c,
                                              double exx, double eyy, double ezz, double exy) {
    HsSmallStrainStep out;
    out.pe = p;
    out.hist = c.hs_hist;
    out.g_mono = c.hs_g_mono;
    out.gamma_hist = c.gamma_hist;
    if (!(p.G0_ref > 0.0)) return out;
    const double gfac = 2.0 * (1.0 + p.nu_ur) * p.G0_ref / p.Eur_ref;   // G0 / G_ur
    if (!(gfac > 1.0)) return out;
    constexpr double kOff = 100.0;   // the listing's offset of the diagonal; any value >> strains
    const double hd = 0.5 * (exx - eyy), rr = std::sqrt(hd * hd + exy * exy);
    const double d[3] = {0.5 * (exx + eyy) + rr, 0.5 * (exx + eyy) - rr, ezz};
    const double de2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    const double f_load = p.gamma07 / HardeningSoilParams::kHDa;
    const auto tangent_ratio = [&](double gam, double f) {
        const double q = 1.0 + gam / f;
        return std::max(gfac / (q * q), 1.0);
    };
    if (!(de2 > 0.0)) {
        // No deviatoric increment: nothing moves in the history and the stress does not depend
        // on the stiffness, so what is set here is only the operator an iteration STARTS with.
        // With no increment there is no direction, and the two directions answer differently:
        // carrying on, with the reloading curve's tangent where the history stands; turning
        // back, with G0. The start is made with G0, the stiffer of the two: an operator that
        // is too stiff makes a first step that falls short and is corrected by the next, one
        // that is too soft overshoots by G0 / G_t. It was the tangent at the history, and the
        // first time step of a consolidation phase that unloads the ground -- where the search
        // takes its least bad step rather than none -- went from 0.44 of its scale to 1e8.
        out.pe.elastic_scale = gfac;
    } else {
        double cs = 1.0, sn = 0.0;
        if (rr > 0.0) {
            const double t = 0.5 * std::atan2(exy, hd);
            cs = std::cos(t); sn = std::sin(t);
        }
        const double hxx = c.hs_hist(0) + kOff, hyy = c.hs_hist(1) + kOff, hzz = c.hs_hist(2) + kOff,
                     hxy = c.hs_hist(3);
        double h11 = cs * cs * hxx + sn * sn * hyy + 2.0 * cs * sn * hxy;
        double h22 = sn * sn * hxx + cs * cs * hyy - 2.0 * cs * sn * hxy;
        double h12 = cs * sn * (hyy - hxx) + (cs * cs - sn * sn) * hxy;
        double h33 = hzz;
        const double chi1[3] = {h11 - kOff, h22 - kOff, h33 - kOff};
        double xt[3];
        for (int i = 0; i < 3; ++i)
            xt[i] = chi1[i] * d[i] < 0.0 ? std::sqrt(kOff / (chi1[i] + kOff)) : 1.0;
        h11 *= xt[0] * xt[0]; h22 *= xt[1] * xt[1]; h33 *= xt[2] * xt[2]; h12 *= xt[0] * xt[1];
        const double chi2[3] = {h11 - kOff, h22 - kOff, h33 - kOff};
        h11 += d[0]; h22 += d[1]; h33 += d[2];
        const double chi3[3] = {h11 - kOff, h22 - kOff, h33 - kOff};
        double n1 = 0.0, n2 = 0.0;
        for (int i = 0; i < 3; ++i) {
            n1 += d[i] * d[i] * chi2[i] * chi2[i];
            n2 += d[i] * d[i] * chi3[i] * chi3[i];
        }
        const double g1 = std::sqrt(3.0 * n1 / de2), g2 = std::sqrt(3.0 * n2 / de2);
        const auto secant_ratio = [&](double f) {
            const double off = (std::sqrt(gfac) - 1.0) * f, span = g2 - g1;
            const auto area = [&](double g) { return g / (1.0 + g / f); };
            if (span <= 0.0) return tangent_ratio(g1, f);
            if (g1 < off && g2 < off) return gfac * (area(g2) - area(g1)) / span;
            if (g1 < off) return (gfac * (area(off) - area(g1)) + (g2 - off)) / span;
            return 1.0;
        };
        const double g_mono0 = c.hs_g_mono > 0.0 ? c.hs_g_mono : std::min(gfac, 10.0);
        out.g_mono = g_mono0 > 1.0 ? std::max(1.0, std::min(secant_ratio(f_load), g_mono0)) : 1.0;
        // No threshold below which the multiplier is simply G0 / G_ur: the secant goes there by
        // itself as the history empties, and a threshold is a step in the stress (0.2% of the
        // stiffness at a history of 1e-6, which an iteration then cannot get below).
        const double act = secant_ratio(2.0 * f_load);
        out.pe.elastic_scale = std::max(act, 1.0);
        out.hist(0) = cs * cs * h11 + sn * sn * h22 - 2.0 * cs * sn * h12 - kOff;
        out.hist(1) = sn * sn * h11 + cs * cs * h22 + 2.0 * cs * sn * h12 - kOff;
        out.hist(2) = h33 - kOff;
        out.hist(3) = cs * sn * (h11 - h22) + (cs * cs - sn * sn) * h12;
        out.gamma_hist = g2;
    }
    const double gm = out.g_mono > 0.0 ? out.g_mono : std::min(gfac, 10.0);
    const double Ei_ref = 2.0 * p.E50_ref / (2.0 - p.Rf) * p.ei_scale;
    if (gm > 1.0 && Ei_ref > 0.0) out.pe.hardening_scale = std::pow(gm, 1.0 + p.Eur_ref / Ei_ref);
    return out;
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
                int* saturated_out = nullptr,
                const Eigen::Vector3d* de_small = nullptr);

// FD step (relative): the sqrt(machine-epsilon) scale — the Pérez-Foguet &
// Rodríguez-Ferran & Huerta (CMAME 2000) recommendation. Signal/cancellation balance:
// Δσ ≈ E·h ~ 1e-3 kPa, double noise ~1e-14·|σ| → ~1e-8 relative error in the columns
// (ample for Newton).
inline double hs_fd_step(double comp) {
    constexpr double kSqrtEps = 1.4901161193847656e-08;  // sqrt(2^-52)
    return kSqrtEps * (1.0 + std::fabs(comp));
}
// The same step, kept below a thousandth of the increment it differences. The fixed step is
// 1.5e-8 whatever the increment, and an increment of that size is differenced with a step as
// large as itself: the perturbed walk then meets surfaces the increment never does. Measured on
// Hardening Soil at random stress points: for increments below 1e-7 the tangent was more than
// 1% off the difference of the stress at one point in five.
inline double hs_fd_step(double comp, double increment_norm) {
    // ... and not below 1e-11: the stress difference over a smaller step is rounding (a stress of
    // 100 kPa carries 1e-14, and E h has to stand well clear of it).
    const double h = hs_fd_step(comp);
    return increment_norm > 0.0 ? std::max(std::min(h, 1.0e-3 * increment_norm), 1.0e-11) : h;
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
    bool unstable = false;      // Cam-Clay: no stress on the softening ellipse (StepResult::unstable)
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
                double* G_out = nullptr, int nsub_fixed = 0, int* nsub_out = nullptr,
                bool* unstable_out = nullptr);

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
    // The STIFFEST response the point can give to any increment from where it stands, where
    // that is more than `elastic`: a soil with small-strain stiffness answers an increment that
    // turns against its strain history with G0, whatever stiffness the increment at hand was
    // taken with. An iteration on a fixed operator contracts only if the operator is at least
    // as stiff as the response, so the elastic-operator pass assembles this one.
    bool has_stiffest = false;
    MatT stiffest = MatT::Zero();
};
using PointReport = PointReportT<Eigen::Matrix3d>;    // plane strain [xx, yy, xy]
using PointReport4 = PointReportT<Eigen::Matrix4d>;   // axisymmetric [r, z, rz, theta]

// Material-point integration. committed: converged state at the start of the
// step; strain_increment: the step's total strain increment Delta-eps; trial:
// updated (output) state; tangent: the 3x3 in-plane operator D_T fed back to the
// element (meaningless while mode==kNone, do not use it).
// small_strain_increment: where given, the strain increment the SMALL-STRAIN STIFFNESS of the
// step is read from (Hardening Soil with small-strain stiffness; every other model ignores it).
// The stiffness of a step depends on the direction of the step's strain increment, which makes
// the stress a function with slopes of order 1 / |d eps| and steps where the increment's
// principal directions are not unique; an iteration that cannot close on it is given the
// stiffness of a fixed increment close to the answer instead (see small_freeze in the solver).
// The stress is integrated over strain_increment and the strain history is advanced by it.
// dt_day: this increment's TIME share [days] — read only by time-dependent constitutive
// models (SoftSoilCreep); 0 (default) is bit-for-bit the old behaviour for all old callers.
void integrate_point(const MaterialModel& m, const GaussState& committed,
                     const Eigen::Vector3d& strain_increment,
                     GaussState& trial, Eigen::Matrix3d& tangent,
                     TangentMode mode = TangentMode::kConsistent, double dt_day = 0.0,
                     PointReport* report = nullptr, double substep_tol = 0.0,
                     const Eigen::Vector3d* small_strain_increment = nullptr);

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
                            double substep_tol = 0.0,
                            const Eigen::Vector4d* small_strain_increment = nullptr);

} // namespace katai::core
