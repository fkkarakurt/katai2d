#pragma once
// van Genuchten (1980) soil-water characteristic curve (SWCC) + Mualem (1976) relative
// permeability — unsaturated flow. The PLAXIS 2D Material Models Manual "van Genuchten"
// model; parameters g_a [1/length], g_n, g_l.
// Suction ψ = −p_w/γ_w  (p_w<0 ⇒ ψ>0 unsaturated; p_w≥0 ⇒ ψ≤0 saturated → S_e=1).
//
// Effective saturation (van Genuchten 1980, Eq):
//   S_e = [ 1 + (g_a·ψ)^{g_n} ]^{−g_m},   g_m = 1 − 1/g_n            (ψ>0; ψ≤0 ⇒ S_e=1)
// Degree of saturation:  S = S_res + (S_sat − S_res)·S_e.
// Relative permeability (Mualem 1976):
//   k_rel = S_e^{g_l} · [ 1 − (1 − S_e^{1/g_m})^{g_m} ]²              (g_l = pore connectivity, def 0.5)
// Moisture capacity (for the transient storage matrix S, Sci.Man Eq 3-33):
//   dS/dψ = (S_sat − S_res)·dS_e/dψ,  dS_e/dψ = −g_m·g_n·u·(1+u)^{−g_m−1}/ψ,  u=(g_a·ψ)^{g_n}  (ψ>0; ≤0 ⇒ 0)
// Sources: van Genuchten (1980) SSSAJ 44:892; Mualem (1976) WRR 12:513; PLAXIS MMM.
// Saturated limit (ψ≤0): S_e=1, k_rel=1, dS/dψ=0 → reduces bit-close to the existing
// saturated seepage/consolidation.

#include <algorithm>
#include <cmath>

namespace katai::core {

// van Genuchten/Mualem water-retention parameters. For saturated soil use g_a→0 (S_e=1
// always) or a call that bypasses retention (saturated zone). k_rel_min = numerical floor
// (PLAXIS ~1e-4; prevents system singularity in the dry zone — Sci.Man Eq 3-15 footnote).
struct WaterRetention {
    double g_a = 2.0;         // [1/length] van Genuchten α (inversely related to air entry)
    double g_n = 2.0;         // van Genuchten n (>1)
    double g_l = 0.5;         // Mualem pore connectivity l
    double S_res = 0.0;       // residual saturation
    double S_sat = 1.0;       // saturated saturation
    double k_rel_min = 1.0e-4;  // relative-permeability floor (numerical)
};

// Effective saturation S_e(ψ) ∈ (0,1].  ψ = suction (≥0 unsaturated); ψ≤0 ⇒ 1 (saturated).
inline double effective_saturation(const WaterRetention& w, double psi) {
    if (psi <= 0.0) return 1.0;
    const double gm = 1.0 - 1.0 / w.g_n;
    const double u = std::pow(w.g_a * psi, w.g_n);
    return std::pow(1.0 + u, -gm);
}

// Degree of saturation S(ψ) = S_res + (S_sat − S_res)·S_e.
inline double saturation(const WaterRetention& w, double psi) {
    return w.S_res + (w.S_sat - w.S_res) * effective_saturation(w, psi);
}

// Moisture capacity dS/dψ ≤ 0 (saturation falls as suction grows).  ψ≤0 ⇒ 0 (saturated, flat).
inline double moisture_capacity(const WaterRetention& w, double psi) {
    if (psi <= 0.0) return 0.0;
    const double gm = 1.0 - 1.0 / w.g_n;
    const double u = std::pow(w.g_a * psi, w.g_n);
    const double dSe_dpsi = -gm * w.g_n * u * std::pow(1.0 + u, -gm - 1.0) / psi;
    return (w.S_sat - w.S_res) * dSe_dpsi;
}

// Relative permeability k_rel(S_e) ∈ [k_rel_min, 1] (Mualem). S_e saturated→1 ⇒ k_rel=1.
inline double relative_permeability(const WaterRetention& w, double Se) {
    const double s = std::clamp(Se, 1.0e-12, 1.0);
    const double gm = 1.0 - 1.0 / w.g_n;
    const double t = 1.0 - std::pow(1.0 - std::pow(s, 1.0 / gm), gm);
    const double kr = std::pow(s, w.g_l) * t * t;
    return std::clamp(kr, w.k_rel_min, 1.0);
}

// Convenience: relative permeability directly from suction, k_rel(S_e(ψ)).
inline double relative_permeability_psi(const WaterRetention& w, double psi) {
    return relative_permeability(w, effective_saturation(w, psi));
}

}  // namespace katai::core
