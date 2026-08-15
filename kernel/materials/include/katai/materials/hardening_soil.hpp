#pragma once
// Hardening Soil (HS) model — P2.3a: stress-dependent stiffness laws + drained triaxial
// hyperbolic primary-loading response (the Duncan-Chang core). This header contains pure
// (closed-form) functions; the multiaxial elastoplastic integration + consistent tangent
// come in later steps (see docs/references/hardening-soil-formulation.md §6–7).
//
// CONVENTION: compression-POSITIVE (classic triaxial), σ3 = cell (minor principal)
// pressure ≥ 0, q = σ1 − σ3 ≥ 0 deviator, −ε1 = axial compression strain. Since the
// solver is tension-positive, the sign flips in the FE integration (σ_HS = −σ_solver).
//
// Sources: Schanz, Vermeer & Bonnier (1999); Duncan & Chang (1970); PLAXIS Material
// Models Manual. Math: docs/references/hardening-soil-formulation.md.

#include <algorithm>
#include <cmath>

namespace katai::core {

// HS parameters (effective). Stiffnesses are reference values at p_ref; m = power exponent.
struct HardeningSoilParams {
    double E50_ref = 0.0;    // secant stiffness (reference, at 50% strength)
    double Eur_ref = 0.0;    // unloading-reloading stiffness (reference)
    double Eoed_ref = 0.0;   // oedometer tangent stiffness (reference)
    double m = 0.5;          // stress-dependency exponent (sand ~0.5, clay ~1)
    double p_ref = 100.0;    // reference pressure (usually 100 kPa)
    double cohesion = 0.0;   // c' (effective)
    double friction = 0.0;   // φ' [rad]
    double dilatancy = 0.0;  // ψ' [rad]
    // The dilatancy cut-off has ENGAGED for this Gauss point (e ≥ e_max). The caller signals it
    // here as well as by zeroing ψ, because for HSsmall those stopped being the same thing when
    // sec. 7.9.1 landed: that rule reads ψ only through φ_cv, so a zeroed ψ moves φ_cv up to φ
    // and switches the Li & Dafalias branch ON everywhere below failure — which would turn a
    // "stop dilating" option into a source of contraction. The manual's cut-off says the
    // MOBILISED angle itself "is automatically set back to zero"; this flag is what keeps
    // ψ_m = 0 meaning ψ_m = 0.
    bool dilatancy_cut = false;
    double Rf = 0.9;         // failure ratio (qa = qf/Rf)
    double nu_ur = 0.2;      // unloading-reloading Poisson

    // Cap (volumetric) yield surface parameters (Rocscience/PLAXIS HS, Eq 15.11-15.13).
    // PLAXIS derives these from K0^NC and Eoed_ref by an oedometer simulation; cap_beta=0 ⇒
    // cap OFF. cap_alpha = α (shape), cap_beta = β (POWER-LAW hardening parameter).
    double cap_alpha = 1.0;  // α  (= PLAXIS cap aspect M; ellipse height M·p_c on the q axis)
    double cap_beta = 0.0;   // β (0 = cap off)

    // Elastic bulk modulus K_e = Eur/(3(1−2ν)) (at σ3).
    double bulk(double sigma3) const {
        return Eur(sigma3) / (3.0 * (1.0 - 2.0 * nu_ur));
    }

    // Cap hardening power law (Eq 15.13): p_c ↔ ε_v^(p-cap). m=1 is the logarithmic special case.
    double cap_ev_from_pc(double pc) const {  // ε_v^(p-cap)(p_c)
        if (std::fabs(1.0 - m) < 1e-12) return cap_beta * std::log(pc / p_ref);
        return cap_beta / (1.0 - m) * std::pow(pc / p_ref, 1.0 - m);
    }
    double cap_pc_from_ev(double ev) const {  // p_c(ε_v^(p-cap))
        if (std::fabs(1.0 - m) < 1e-12) return p_ref * std::exp(ev / cap_beta);
        return p_ref * std::pow((1.0 - m) * ev / cap_beta, 1.0 / (1.0 - m));
    }
    double cap_hardening_modulus(double pc) const {  // H_cap = dp_c/dε_v^pc = (p_ref/β)(p_c/p_ref)^m
        return (p_ref / cap_beta) * std::pow(pc / p_ref, m);
    }

    // Stress-dependent scale factor ((c·cosφ + σ·sinφ)/(c·cosφ + p_ref·sinφ))^m.
    // sigma = the relevant principal stress (σ3 for E50/Eur, σ1 for Eoed).
    double stiffness_factor(double sigma) const {
        const double cc = cohesion * std::cos(friction);
        const double s = std::sin(friction);
        const double num = cc + sigma * s;
        const double den = cc + p_ref * s;
        return std::pow(num / den, m);
    }

    double E50(double sigma3) const { return E50_ref * stiffness_factor(sigma3); }
    double Eur(double sigma3) const { return Eur_ref * stiffness_factor(sigma3); }
    double Eoed(double sigma1) const { return Eoed_ref * stiffness_factor(sigma1); }

    // --- HSsmall (Material Models Manual §7): small-strain stiffness. G0_ref=0 ⇒ OFF (plain HS).
    double G0_ref = 0.0;     // very-small-strain shear modulus (reference); 0 = no HSsmall
    double gamma07 = 1.0e-4; // threshold shear strain where G_s drops to 0.722·G0 (virgin)
    static constexpr double kHDa = 0.385;  // Hardin-Drnevich (Santos&Correia): G_s=0.722G0 @ γ=γ07
    // Masing's rule, Eq 7-11: gamma_0.7,re-loading = 2 gamma_0.7,virgin-loading. The INPUT
    // gamma07 is the virgin-loading threshold (sec. 7.4: "gamma_0.7 is to be supplied for
    // virgin loading"); the curve the model's QUASI-ELASTIC stiffness rides on is that
    // backbone scaled by 2. The manual is explicit that this factor is not switched on at a
    // detected reversal: "the scaling factor for the threshold shear strain is assumed to be
    // constant and equal to 2 throughout loading", because in HS-small the virgin response is
    // elasto-plastic from the start of shearing and the hardening plasticity already supplies
    // the faster virgin decay. The functions below are the BACKBONE (Eq 7-3/7-7/7-8/7-10, at
    // the input gamma07); the FE overlay in hs_small_strain_params rides the reloading curve.
    static constexpr double kMasing = 2.0;
    double gamma07_reload() const { return kMasing * gamma07; }
    // MMM sec. 7.5: "Although Alpan suggests that the ratio E0/Eur can exceed 10 for very soft
    // clays, the maximum ratio E0/Eur or G0/Gur permitted in the HSsmall model is limited to
    // 20." Both moduli follow the same power law, so the ratio is stress-independent and the
    // cap can be applied once, at the reference values.
    static constexpr double kMaxG0Ratio = 20.0;
    double Gur_ref() const { return Eur_ref / (2.0 * (1.0 + nu_ur)); }
    double G0_ref_cap() const { return kMaxG0Ratio * Gur_ref(); }

    double G0(double sigma3) const { return G0_ref * stiffness_factor(sigma3); }     // stress-dependent
    double Gur(double sigma3) const { return Eur(sigma3) / (2.0 * (1.0 + nu_ur)); }  // lower cut-off modulus
    double g_secant(double gamma, double sigma3) const {                            // Eq 7-3
        return G0(sigma3) / (1.0 + kHDa * std::fabs(gamma) / gamma07);
    }
    double g_tangent(double gamma, double sigma3) const {                           // Eq 7-8 (+ cut-off 7-9)
        const double d = 1.0 + kHDa * std::fabs(gamma) / gamma07;
        return std::max(G0(sigma3) / (d * d), Gur(sigma3));
    }
    double gamma_cutoff(double sigma3) const {                                       // Eq 7-10
        return (1.0 / kHDa) * (std::sqrt(G0(sigma3) / Gur(sigma3)) - 1.0) * gamma07;
    }
    // Strain-dependent elastic (unload/reload) Young modulus: replaces HS's constant Eur.
    double Et_small(double gamma, double sigma3) const {
        return 2.0 * (1.0 + nu_ur) * g_tangent(gamma, sigma3);
    }

    // MC failure deviator qf = 2(c·cosφ + σ3·sinφ)/(1 − sinφ)  [= σ3(Kp−1)+2c√Kp].
    double q_failure(double sigma3) const {
        const double s = std::sin(friction);
        return 2.0 * (cohesion * std::cos(friction) + sigma3 * s) / (1.0 - s);
    }

    // Initial stiffness Ei = 2·E50/(2−Rf) and asymptote qa = qf/Rf.
    double Ei(double sigma3) const { return 2.0 * E50(sigma3) / (2.0 - Rf); }
    double q_asymptote(double sigma3) const { return q_failure(sigma3) / Rf; }

    // Drained triaxial PRIMARY loading: deviator q for a given axial compression −ε1 (≥0).
    // The hyperbola −ε1 = (1/Ei) q/(1−q/qa) is inverted: q = Ei·(−ε1)·qa/(qa + Ei·(−ε1));
    // cut at qf at MC failure (perfectly plastic plateau).
    double triaxial_q(double sigma3, double minus_eps1) const {
        if (minus_eps1 <= 0.0) return 0.0;
        const double ei = Ei(sigma3), qa = q_asymptote(sigma3);
        const double a = ei * minus_eps1;
        const double q_hyper = a * qa / (qa + a);
        return std::min(q_hyper, q_failure(sigma3));
    }

    // Primary tangent modulus Et = dq/d(−ε1) = Ei(1−q/qa)² (0 past failure = perfectly plastic).
    double triaxial_tangent(double sigma3, double minus_eps1) const {
        const double q = triaxial_q(sigma3, minus_eps1);
        if (q >= q_failure(sigma3)) return 0.0;
        const double r = 1.0 - q / q_asymptote(sigma3);
        return Ei(sigma3) * r * r;
    }

    // --- Oedometer (1D compression) — the uniaxial response of cap hardening ----------
    // In confined (ε_h=0) primary loading the tangent oedometer modulus Eoed depends on
    // the vertical stress σ1: dε1 = dσ1/Eoed(σ1). This is the uniaxial manifestation of
    // the HS cap yield surface (volumetric hardening) (FE return mapping P2.3d).
    // Closed-form integral (compression-positive axial strain −ε1 ≥ 0; cc=c·cosφ, s=sinφ):
    //   m≠1: −ε1 = D/(Eoed_ref·s) · [(cc+σ1·s)^(1−m) − (cc+σ0·s)^(1−m)]/(1−m)
    //   m=1: −ε1 = D/(Eoed_ref·s) · ln((cc+σ1·s)/(cc+σ0·s))    (classic logarithmic clay)
    // D = (cc+p_ref·s)^m. c=0, m=1 → −ε1 = (p_ref/Eoed_ref)·ln(σ1/σ0).
    double oedometer_strain(double sigma1_0, double sigma1) const {
        const double cc = cohesion * std::cos(friction);
        const double s = std::sin(friction);
        const double D = std::pow(cc + p_ref * s, m);
        const double u0 = cc + sigma1_0 * s, u1 = cc + sigma1 * s;
        const double coef = D / (Eoed_ref * s);
        if (std::fabs(1.0 - m) < 1e-12)
            return coef * std::log(u1 / u0);
        return coef * (std::pow(u1, 1.0 - m) - std::pow(u0, 1.0 - m)) / (1.0 - m);
    }
};

// Stateful drained-triaxial path for the HS SHEAR mechanism (constant cell pressure
// sigma3). Demonstrates the defining HS feature vs Duncan-Chang nonlinear elasticity:
// primary loading is plastic (hyperbola), unloading/reloading is elastic with the stiff
// Eur, leaving permanent plastic strain. Axial strain decomposes as
//   -eps1 = eps1_p + q/Eur ,   eps1_p = hyperbola(q_max) - q_max/Eur   (monotone, history)
// with q_max the hardening memory (largest deviator reached, capped at qf). This is the
// 1D form of the shear yield surface + hardening; the multiaxial FE integration (return
// mapping + consistent tangent) builds on this logic (P2.3d).
// (See docs/references/hardening-soil-formulation.md sec 2-3.)
struct HsTriaxialPath {
    HardeningSoilParams params;
    double sigma3 = 0.0;
    double q_max = 0.0;  // hardening memory (max deviator reached)

    // Primary-loading hyperbola -eps1(q) = (1/Ei) q/(1-q/qa).
    double hyperbola(double q) const {
        return (1.0 / params.Ei(sigma3)) * q / (1.0 - q / params.q_asymptote(sigma3));
    }
    // Permanent (plastic) axial strain accumulated up to the current history q_max.
    double plastic_strain() const {
        return hyperbola(q_max) - q_max / params.Eur(sigma3);
    }
    // Axial compression strain -eps1 at deviator q. Loading beyond q_max grows the
    // hardening memory (plastic); within q_max it is elastic unload/reload (slope 1/Eur).
    // q is capped at qf (perfectly plastic failure).
    double axial_strain(double q) {
        const double qf = params.q_failure(sigma3);
        const double qc = std::min(std::max(q, 0.0), qf);
        if (qc > q_max) q_max = qc;  // primary loading extends the memory
        return plastic_strain() + qc / params.Eur(sigma3);
    }
};

}  // namespace katai::core
