#pragma once
// Eurocode 7 (EN 1997-1) + TBDY 2018 design-code partial factors -- v0.3 workstream B1.
// Values LOCKED (reference-checked) in docs/references/design-codes-ec7-tbdy.md.
//
// Two factoring paradigms map DIFFERENTLY onto a deformation FEM (reference doc S0):
//   * MATERIAL factoring (EC7 DA1-C2, DA3): reduce shear strength
//       c'_d = c'/gamma_c',  phi'_d = atan(tan phi'/gamma_phi'),  c_u,d = c_u/gamma_cu (Tresca),
//     then re-solve. This MIRRORS strength_reduction.cpp::factor_strength but with the code's
//     FIXED partial factors instead of a bisected strength-reduction factor.
//   * RESISTANCE / ACTION factoring (EC7 DA2, TBDY 2018): factor actions (A-set, multiply loads)
//     and the computed resistance (R-set / gamma_Rv, E_d <= R_d). That lives in the check/report
//     layer, NOT in the constitutive strength -- this header only supplies the numbers.

#include <cmath>

#include <katai/materials/material_model.hpp>

namespace katai::core {

// The design situation. Enum values are file-stable (append only) so serialized projects survive.
enum class DesignApproach {
    None = 0,           // characteristic values -- degenerates to the plain (unfactored) analysis
    EC7_DA1_C1 = 1,     // A1 + M1 + R1  (structural-action governed; material unreduced)
    EC7_DA1_C2 = 2,     // A2 + M2 + R1  (material-factored -> re-solve)
    EC7_DA2 = 3,        // A1 + M1 + R2  (resistance-factored)
    EC7_DA3 = 4,        // A2(geotechnical) + M2 + R3  (material-factored -> re-solve)
    TBDY2018_Static = 5,   // E_d <= R_d = R_k/gamma_Rv, static combination (resistance-factored)
    TBDY2018_Seismic = 6,  // ditto, seismic combination
};

// Partial factors for one design combination. Material factors DIVIDE strength; action factors
// MULTIPLY loads; resistance factors DIVIDE the computed resistance. A value of 1.0 = "no effect".
struct PartialFactors {
    // Material set (M) -- divide strength.
    double gamma_phi = 1.0;    // on tan(phi')  (effective angle of shearing resistance)
    double gamma_c = 1.0;      // on c'         (effective cohesion)
    double gamma_cu = 1.0;     // on c_u        (undrained shear strength, Tresca)
    double gamma_gamma = 1.0;  // on unit weight
    // Action set (A) -- multiply loads.
    double gamma_G_unfav = 1.0;
    double gamma_G_fav = 1.0;
    double gamma_Q_unfav = 1.0;
    double gamma_Q_fav = 0.0;
    // Resistance set (R) -- divide the computed resistance (report/check layer).
    double gamma_Rv = 1.0;  // bearing
    double gamma_Rh = 1.0;  // sliding
    double gamma_Re = 1.0;  // passive / earth resistance (EN 1997-1 Table A.13; slope stability is DA3)
    // True when the resistance factors are NOT yet pinned from the primary standard
    // (TBDY 2018 Table 16.1). The check layer must refuse a pass/fail verdict until they are
    // filled -- we never emit a fabricated safety number.
    bool resistance_pending = false;
};

// EN 1997-1 Annex A recommended values (reference doc S1) + TBDY 2018 framework (S2).
inline PartialFactors design_factors(DesignApproach da) {
    PartialFactors f;
    // Named EN 1997-1 recommended sets.
    auto A1 = [&] { f.gamma_G_unfav = 1.35; f.gamma_G_fav = 1.0; f.gamma_Q_unfav = 1.5; f.gamma_Q_fav = 0.0; };
    auto A2 = [&] { f.gamma_G_unfav = 1.0;  f.gamma_G_fav = 1.0; f.gamma_Q_unfav = 1.3; f.gamma_Q_fav = 0.0; };
    auto M1 = [&] { f.gamma_phi = 1.0;  f.gamma_c = 1.0;  f.gamma_cu = 1.0; f.gamma_gamma = 1.0; };
    auto M2 = [&] { f.gamma_phi = 1.25; f.gamma_c = 1.25; f.gamma_cu = 1.4; f.gamma_gamma = 1.0; };
    auto R1 = [&] { f.gamma_Rv = 1.0; f.gamma_Rh = 1.0; f.gamma_Re = 1.0; };
    auto R2 = [&] { f.gamma_Rv = 1.4; f.gamma_Rh = 1.1; f.gamma_Re = 1.4; };  // Re = passive (Table A.13)
    auto R3 = [&] { f.gamma_Rv = 1.0; f.gamma_Rh = 1.0; f.gamma_Re = 1.0; };
    switch (da) {
        case DesignApproach::None:                          break;
        case DesignApproach::EC7_DA1_C1: A1(); M1(); R1();  break;
        case DesignApproach::EC7_DA1_C2: A2(); M2(); R1();  break;
        case DesignApproach::EC7_DA2:    A1(); M1(); R2();  break;
        case DesignApproach::EC7_DA3:    A2(); M2(); R3();  break;  // A2 = geotechnical actions
        case DesignApproach::TBDY2018_Static:
        case DesignApproach::TBDY2018_Seismic:
            // Resistance-factored (no material reduction), like EC7 DA2. TBDY 2018 Table 16.2
            // (Section 16.8.2): bearing gamma_Rv = 1.40, sliding gamma_Rh = 1.10, passive
            // gamma_Re = 1.40 -- NUMERICALLY IDENTICAL to EC7 DA2 (R2). The same factors apply to
            // static and seismic; the difference is in the load combination and the seismic bearing-
            // capacity formula, not the factor. Values: official Resmi Gazete 18.03.2018 (TBDY 2018),
            // cross-verified against two sources + the EC7 DA2 identity (docs/references/design-codes-ec7-tbdy.md).
            f.gamma_Rv = 1.4; f.gamma_Rh = 1.1; f.gamma_Re = 1.4;
            break;
    }
    return f;
}

// Whether this approach reduces MATERIAL strength (i.e. the solve re-runs with factored c/phi).
inline bool factors_material(DesignApproach da) {
    return da == DesignApproach::EC7_DA1_C2 || da == DesignApproach::EC7_DA3;
}

// Apply the M-set to a material's shear strength, MIRRORING strength_reduction.cpp::factor_strength
// but with distinct gamma_c'/gamma_phi' (and gamma_cu for an undrained Tresca material, phi=0). The
// Hardening Soil failure-surface sub-struct carries its OWN c/phi, so it is reduced too -- otherwise
// an HS material keeps full strength. Dilatancy is clamped to the reduced friction (psi <= phi).
inline void factor_material_strength(MaterialModel& m, const PartialFactors& f) {
    const bool tresca = m.undrained && m.friction_angle == 0.0;  // Undrained (B): c = c_u, phi = 0
    m.cohesion /= (tresca ? f.gamma_cu : f.gamma_c);
    m.friction_angle = std::atan(std::tan(m.friction_angle) / f.gamma_phi);
    if (m.dilatancy_angle > m.friction_angle) m.dilatancy_angle = m.friction_angle;
    // MC tension cap: factored like the cohesion-class strength it is (EC7 has no
    // dedicated tensile-strength factor; mirroring Safety's reduction is the safe
    // direction and keeps the two factoring paths consistent).
    m.tensile_strength /= (tresca ? f.gamma_cu : f.gamma_c);
    const bool hs_tresca = m.undrained && m.hs.friction == 0.0;
    m.hs.cohesion /= (hs_tresca ? f.gamma_cu : f.gamma_c);
    m.hs.friction = std::atan(std::tan(m.hs.friction) / f.gamma_phi);
    if (m.hs.dilatancy > m.hs.friction) m.hs.dilatancy = m.hs.friction;
}

} // namespace katai::core
