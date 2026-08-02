// Eurocode 7 + TBDY 2018 design-code partial factors -- value + material-factoring verification (B1).
// The factor VALUES are locked in docs/references/design-codes-ec7-tbdy.md (EN 1997-1 Annex A
// recommended). Material factoring must MIRROR strength_reduction.cpp::factor_strength (same atan-of-
// tan reduction), and must degenerate to identity when every partial factor is 1.0.
#include <katai/analysis/design_code.hpp>
#include <katai/materials/material_model.hpp>

#include <cmath>
#include <cstdio>

using katai::core::DesignApproach;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::PartialFactors;
using katai::core::design_factors;
using katai::core::factor_material_strength;
using katai::core::factors_material;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol = 1e-12) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}
const double kDeg = std::acos(-1.0) / 180.0;

// (1) Each Design Approach returns the exact EN 1997-1 recommended set combination.
void test_ec7_tables() {
    // DA1-C1 = A1 + M1 + R1
    const PartialFactors c1 = design_factors(DesignApproach::EC7_DA1_C1);
    check(close(c1.gamma_G_unfav, 1.35) && close(c1.gamma_Q_unfav, 1.5), "DA1-C1 actions A1 (1.35 / 1.5)");
    check(close(c1.gamma_phi, 1.0) && close(c1.gamma_c, 1.0), "DA1-C1 material M1 (1.0)");
    check(close(c1.gamma_Rv, 1.0), "DA1-C1 resistance R1 (1.0)");

    // DA1-C2 = A2 + M2 + R1  (the material-factored combination)
    const PartialFactors c2 = design_factors(DesignApproach::EC7_DA1_C2);
    check(close(c2.gamma_G_unfav, 1.0) && close(c2.gamma_Q_unfav, 1.3), "DA1-C2 actions A2 (1.0 / 1.3)");
    check(close(c2.gamma_phi, 1.25) && close(c2.gamma_c, 1.25) && close(c2.gamma_cu, 1.4),
          "DA1-C2 material M2 (1.25 / 1.25 / 1.4)");
    check(close(c2.gamma_Rv, 1.0), "DA1-C2 resistance R1 (1.0)");

    // DA2 = A1 + M1 + R2
    const PartialFactors d2 = design_factors(DesignApproach::EC7_DA2);
    check(close(d2.gamma_G_unfav, 1.35), "DA2 actions A1 (1.35)");
    check(close(d2.gamma_phi, 1.0), "DA2 material M1 (1.0)");
    check(close(d2.gamma_Rv, 1.4) && close(d2.gamma_Rh, 1.1) && close(d2.gamma_Re, 1.4),
          "DA2 resistance R2 (1.4 / 1.1 / 1.4)");

    // DA3 = A2(geotechnical) + M2 + R3
    const PartialFactors d3 = design_factors(DesignApproach::EC7_DA3);
    check(close(d3.gamma_G_unfav, 1.0), "DA3 geotechnical actions A2 (1.0)");
    check(close(d3.gamma_phi, 1.25) && close(d3.gamma_c, 1.25), "DA3 material M2 (1.25)");
    check(close(d3.gamma_Rv, 1.0), "DA3 resistance R3 (1.0)");

    // Which approaches re-solve with factored material.
    check(factors_material(DesignApproach::EC7_DA1_C2) && factors_material(DesignApproach::EC7_DA3),
          "DA1-C2 and DA3 factor material");
    check(!factors_material(DesignApproach::EC7_DA1_C1) && !factors_material(DesignApproach::EC7_DA2),
          "DA1-C1 and DA2 do not factor material");

    // TBDY 2018 Table 16.2 (Section 16.8.2) resistance factors, pinned to the official values:
    // bearing 1.4, sliding 1.1, passive 1.4 -- numerically identical to EC7 DA2.
    const PartialFactors tbdy = design_factors(DesignApproach::TBDY2018_Static);
    check(close(tbdy.gamma_Rv, 1.4) && close(tbdy.gamma_Rh, 1.1) && close(tbdy.gamma_Re, 1.4),
          "TBDY 2018 resistance factors (Table 16.2): 1.4 / 1.1 / 1.4");
    check(!tbdy.resistance_pending, "TBDY 2018 resistance factors are pinned (not pending)");
    check(close(tbdy.gamma_Rv, d2.gamma_Rv) && close(tbdy.gamma_Rh, d2.gamma_Rh) &&
              close(tbdy.gamma_Re, d2.gamma_Re),
          "TBDY 2018 resistance factors == EC7 DA2 (R2)");
    check(design_factors(DesignApproach::TBDY2018_Seismic).gamma_Rv == 1.4,
          "TBDY 2018 seismic uses the same gamma_Rv as static");

    // None = characteristic values (identity).
    const PartialFactors none = design_factors(DesignApproach::None);
    check(close(none.gamma_phi, 1.0) && close(none.gamma_c, 1.0) && close(none.gamma_Rv, 1.0),
          "None = characteristic (all factors 1.0)");
}

// (2) Material factoring reduces strength exactly as EN 1997-1 prescribes.
void test_material_factoring() {
    // Drained Mohr-Coulomb: c' = 10, phi' = 30 deg. DA1-C2 (M2) -> c'/1.25, atan(tan30/1.25).
    MaterialModel m;
    m.type = MaterialType::MohrCoulomb;
    m.cohesion = 10.0;
    m.friction_angle = 30.0 * kDeg;
    m.dilatancy_angle = 5.0 * kDeg;
    m.hs.cohesion = 10.0;
    m.hs.friction = 30.0 * kDeg;
    factor_material_strength(m, design_factors(DesignApproach::EC7_DA1_C2));
    check(close(m.cohesion, 10.0 / 1.25), "drained c' / gamma_c' = 8.0");
    check(close(m.friction_angle, std::atan(std::tan(30.0 * kDeg) / 1.25)),
          "drained tan(phi') / gamma_phi'");
    check(m.dilatancy_angle <= m.friction_angle + 1e-15, "dilatancy clamped to reduced friction");
    check(close(m.hs.cohesion, 8.0), "HS failure-surface c reduced consistently");

    // Undrained Tresca (B): c = c_u = 50, phi = 0. DA1-C2 -> c_u / gamma_cu (1.4); phi stays 0.
    MaterialModel u;
    u.type = MaterialType::MohrCoulomb;
    u.undrained = true;
    u.cohesion = 50.0;
    u.friction_angle = 0.0;
    factor_material_strength(u, design_factors(DesignApproach::EC7_DA1_C2));
    check(close(u.cohesion, 50.0 / 1.4), "undrained Tresca c_u / gamma_cu = 35.714");
    check(close(u.friction_angle, 0.0), "undrained phi stays 0");

    // Identity: an all-1.0 set (M1) leaves strength unchanged (within atan(tan) round-off).
    MaterialModel id;
    id.type = MaterialType::MohrCoulomb;
    id.cohesion = 12.0;
    id.friction_angle = 28.0 * kDeg;
    factor_material_strength(id, design_factors(DesignApproach::EC7_DA1_C1));  // M1: all gamma = 1
    check(close(id.cohesion, 12.0) && close(id.friction_angle, 28.0 * kDeg),
          "M1 (all factors 1.0) is identity");
}

}  // namespace

int main() {
    test_ec7_tables();
    test_material_factoring();
    if (g_failures) {
        std::fprintf(stderr, "%d design-code check(s) FAILED\n", g_failures);
        return 1;
    }
    std::printf("design-code (EC7 tables + material factoring + TBDY framework): all checks passed\n");
    return 0;
}
