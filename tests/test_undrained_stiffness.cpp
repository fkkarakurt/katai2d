// The pore fluid's stiffness, per material.
//
// Undrained (A) and (B) work by adding a bulk stiffness for the water to the soil skeleton:
// Kw/n. It is not a property of water -- PLAXIS says so plainly ("The bulk stiffness of water
// Kw, calculated in this way, is a numerical value related to the soil stiffness") -- it is a
// number derived from the material's own effective stiffness and from one choice the engineer
// makes: how compressible the undrained soil is allowed to be. PLAXIS offers that choice two
// ways, the equivalent undrained Poisson ratio entered directly or Skempton's B, and it is a
// per-material choice because it describes a material.
//
// KATAI used the constant 0.495 for every undrained material in every project, which is PLAXIS's
// DEFAULT presented as if it were the rule, and -- for the Hardening Soil family -- derived K'
// from the E and nu boxes that model never reads. This case pins both halves of the repair.
//
// verify: KV-CST-006
//   oracle:   closed_form
//   source:   PLAXIS 2D 2025.1 Material Models Manual section 2.4 (Undrained effective stress analysis), read from the manual: Eq. 2-50 Kw/n = 3(nu_u - nu')/((1 - 2 nu_u)(1 + nu')) K' = ((0.495 - nu')/(1 + nu')) 300 K' >= 30 K' for alpha_Biot = 1; Eq. 2-54 Ku = 2G(1 + nu_u)/(3(1 - 2 nu_u)); Eq. 2-55 nu_u = (3 nu' + alpha B (1 - 2 nu'))/(3 - alpha B (1 - 2 nu')); Eq. 2-57 B = alpha/(alpha + n(K'/Kw + alpha - 1)), which at alpha = 1 is B = (Kw/n)/(K' + Kw/n). With the 1D confined closed form u = -q H / M_u, M_u = M' + Kw/n (docs/references/effective-stress-formulation.md)
//   locator:  a weightless laterally confined two-layer column loaded undrained (A) through a full-width surcharge, the lower layer given nu_u = 0.497 directly and the upper one Skempton's B = 0.90 -- two materials that could not have differed at all before, since every undrained material was solved at nu_u = 0.495
//   quantity: settlement of the column top and of the layer interface, and the effective vertical stress carried by each layer [m; kPa]
//   expected: u(interface) = -q H1 / M_u1 and u(top) = -q (H1/M_u1 + H2/M_u2), each M_u,i built from that layer's OWN Kw/n; sigma'_yy = -M' q / M_u per layer; the manual's three equations agree with each other (the Kw/n the engine derives through nu_u equals the one Eq. 2-57 gives from B alone); and the interface settlement is unchanged when only the UPPER layer's B is changed, since in a column the layer below cannot know what was done above it
//   band:     2% on both settlements and 3% on the effective stresses, as asserted below (measured +0.00% on all four); 1e-14 relative on the manual's equation ring; 1e-12 relative on the interface under a changed upper layer (measured 2.1e-14, against a 40% change at the surface -- not bit-identity, because the two runs solve different global systems)
//
// The second half -- which K' an advanced model's pore fluid is sized by -- is checked at the
// registry rather than through a BVP, because the quantity IS a construction: what the catalogue
// builds from a Hardening Soil data set whose Linear-elastic boxes were never filled.

#include <katai/analysis/results.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/materials/registry.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::core::MaterialModel;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kE = 1.0e4, kNu = 0.3;      // effective stiffness, both layers
constexpr double kH1 = 6.0, kH2 = 4.0;       // lower / upper layer thickness [m]
constexpr double kW = 2.0, kQ = 50.0;        // column width [m], surcharge [kPa]
constexpr double kNuU1 = 0.497;              // lower layer: nu_u entered directly
constexpr double kB2 = 0.90;                 // upper layer: Skempton's B entered

double k_eff(double E, double nu) { return E / (3.0 * (1.0 - 2.0 * nu)); }
double m_eff(double E, double nu) { return E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); }
// MMM Eq. 2-50, written out here rather than called from the header: the comparison is meant to
// be law against implementation, not implementation against itself.
double kwn_from_nu_u(double nu_u, double E, double nu) {
    return 3.0 * (nu_u - nu) / ((1.0 - 2.0 * nu_u) * (1.0 + nu)) * k_eff(E, nu);
}
// MMM Eq. 2-57 at alpha_Biot = 1, solved for Kw/n: B = (Kw/n)/(K' + Kw/n)  =>  Kw/n = B K'/(1 - B).
// This route never touches nu_u, so agreeing with the engine (which goes B -> nu_u -> Kw/n) is a
// statement about the manual's own equations as much as about the code.
double kwn_from_skempton(double B, double E, double nu) {
    return B * k_eff(E, nu) / (1.0 - B);
}

// A 2 x 10 m weightless column in two layers, laterally confined, loaded by a full-width
// surcharge: the setting in which the undrained 1D closed form is exact.
m::Project two_layer_column(double B_upper) {
    m::Project pr;
    pr.name = "KV-CST-006 undrained stiffness";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH1 + kH2;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;

    m::Material lower;
    lower.name = "Lower clay (nu_u entered)";
    lower.model = m::SoilModel::LinearElastic;
    lower.drainage = m::Drainage::Undrained;
    lower.E = kE; lower.nu = kNu;
    lower.gamma_unsat = 0.0; lower.gamma_sat = 0.0;   // only the surcharge acts
    lower.und_mode = 0;
    lower.nu_u = kNuU1;
    pr.materials.push_back(lower);

    m::Material upper = lower;
    upper.name = "Upper clay (Skempton B entered)";
    upper.und_mode = 1;
    upper.nu_u = 0.495;          // present and deliberately NOT the one in force
    upper.skempton_B = B_upper;
    pr.materials.push_back(upper);

    m::SoilPolygon L;
    L.name = "Lower";
    L.material = 0;
    L.x = {0, kW, kW, 0};
    L.y = {0, 0, kH1, kH1};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(L);

    m::SoilPolygon U;
    U.name = "Upper";
    U.material = 1;
    U.x = {0, kW, kW, 0};
    U.y = {kH1, kH1, kH1 + kH2, kH1 + kH2};
    U.edge_bc = {(int)m::BCType::Free, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(U);

    m::Load q;
    q.kind = m::LoadKind::Distributed;
    q.name = "Surcharge";
    q.x1 = 0; q.y1 = kH1 + kH2; q.x2 = kW; q.y2 = kH1 + kH2;
    q.qx1 = q.qx2 = 0; q.qy1 = q.qy2 = -kQ;
    pr.loads.push_back(q);
    return pr;
}

// Lowest (most negative) vertical displacement on a horizontal line, and the mean effective
// vertical stress of the nodes strictly inside a layer.
double settlement_at(const katai::app::SolveResult& R, double y) {
    double u = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.y[n] - y) < 1e-6) u = std::fmin(u, R.disp[n * 2 + 1]);
    return u;
}
double stress_between(const katai::app::SolveResult& R, double ylo, double yhi) {
    double sum = 0.0;
    int count = 0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > ylo && R.mesh.y[n] < yhi) { sum += R.stress.stress[n](1); ++count; }
    return count ? sum / count : 0.0;
}

// ---------------------------------------------------------------- the manual's equation ring --
void test_equation_ring() {
    std::printf("== the three equations of MMM section 2.4 agree with each other ==\n");
    const double nus[] = {0.0, 0.15, 0.2, 0.3, 0.35};
    const double Bs[] = {0.5, 0.8, 0.9, 0.95, 0.98, 0.999};
    double worst_ring = 0.0, worst_ku = 0.0, min_ratio = 1e300;
    for (double nu : nus) {
        for (double B : Bs) {
            // B -> nu_u (Eq. 2-55) -> Kw/n (Eq. 2-50) -> B (Eq. 2-57): the ring must close.
            const double nu_u = katai::core::undrained_poisson_from_skempton(B, nu);
            const double kwn = kwn_from_nu_u(nu_u, kE, nu);
            const double B_back = katai::core::skempton_from_kw_over_n(kwn, k_eff(kE, nu));
            worst_ring = std::fmax(worst_ring, std::fabs(B_back - B) / B);
            // And Eq. 2-54: K' + Kw/n IS the undrained bulk modulus 2G(1 + nu_u)/(3(1 - 2 nu_u)).
            const double G = kE / (2.0 * (1.0 + nu));
            const double Ku = 2.0 * G * (1.0 + nu_u) / (3.0 * (1.0 - 2.0 * nu_u));
            worst_ku = std::fmax(worst_ku,
                                 std::fabs(k_eff(kE, nu) + kwn - Ku) / Ku);
        }
        // The manual's printed numeric claim at its own default: Kw/n = ((0.495 - nu')/(1 + nu'))
        // 300 K', "larger than 30 K', at least for nu' <= 0.35".
        const double kwn_default = kwn_from_nu_u(0.495, kE, nu);
        const double printed = (0.495 - nu) / (1.0 + nu) * 300.0 * k_eff(kE, nu);
        check(std::fabs(kwn_default - printed) < 1e-9 * printed,
              "Eq. 2-50's two printed forms agree at nu' = " + std::to_string(nu));
        min_ratio = std::fmin(min_ratio, kwn_default / k_eff(kE, nu));
    }
    std::printf("  ring B -> nu_u -> Kw/n -> B closes to %.2e relative\n", worst_ring);
    std::printf("  K' + Kw/n = Ku (Eq. 2-54) to %.2e relative\n", worst_ku);
    std::printf("  smallest Kw/n over nu' in [0, 0.35] at nu_u = 0.495: %.1f K' (manual: >= 30 K')\n",
                min_ratio);
    check(worst_ring < 1e-14, "Eq. 2-55, 2-50 and 2-57 are mutually consistent");
    check(worst_ku < 1e-14, "Eq. 2-50 reproduces the undrained bulk modulus of Eq. 2-54");
    check(min_ratio >= 30.0, "Kw/n >= 30 K' for nu' <= 0.35, as the manual states");
}

// ------------------------------------------------------- which K' an advanced model is sized by --
void test_hardening_soil_reference_stiffness() {
    std::printf("\n== the pore fluid of a Hardening Soil material follows THAT model's elasticity ==\n");
    katai::core::MaterialParams p;
    p.E = 1.3e4; p.nu = 0.3;          // the Linear-elastic boxes: a default nobody entered
    p.E50_ref = 3.0e4; p.Eoed_ref = 3.0e4; p.Eur_ref = 9.0e4;
    p.m = 0.5; p.p_ref = 100.0; p.Rf = 0.9; p.nu_ur = 0.2;
    p.c = 5.0; p.phi_rad = 25.0 * std::acos(-1.0) / 180.0;
    p.drainage = katai::core::DrainageClass::UndrainedA;

    const katai::core::ModelEntry* entry = katai::core::find_model("HardeningSoil");
    check(entry != nullptr, "the Hardening Soil model is registered");
    if (!entry) return;
    const MaterialModel mm = entry->build(p);
    const double got = mm.kw_over_n(mm.undrained_poisson);
    const double want = kwn_from_nu_u(0.495, p.Eur_ref, p.nu_ur);
    const double from_unread_boxes = kwn_from_nu_u(0.495, p.E, p.nu);
    std::printf("  Kw/n from (Eur_ref = %.0f, nu_ur = %.2f) : %.0f kPa\n", p.Eur_ref, p.nu_ur, want);
    std::printf("  Kw/n from the unread (E, nu) boxes        : %.0f kPa  (%.2fx smaller)\n",
                from_unread_boxes, want / from_unread_boxes);
    check(std::fabs(got - want) < 1e-9 * want,
          "Kw/n is built from the unload/reload pair the HS model actually reads");
    check(std::fabs(got - from_unread_boxes) > 0.5 * from_unread_boxes,
          "and is not the number the untouched Linear-elastic boxes would have given");

    // Skempton's B on such a material must be resolved against nu_ur too, or the same
    // inconsistency comes back through the other input.
    p.skempton_mode = true;
    p.skempton_B = 0.90;
    const MaterialModel ms = entry->build(p);
    const double want_s = kwn_from_skempton(0.90, p.Eur_ref, p.nu_ur);
    std::printf("  Skempton B = 0.90 -> nu_u = %.6f, Kw/n = %.0f kPa (B route: %.0f)\n",
                ms.undrained_poisson, ms.kw_over_n(ms.undrained_poisson), want_s);
    check(std::fabs(ms.kw_over_n(ms.undrained_poisson) - want_s) < 1e-9 * want_s,
          "Skempton's B on a Hardening Soil material resolves against nu_ur");
}

}  // namespace

int main() {
    test_equation_ring();
    test_hardening_soil_reference_stiffness();

    std::printf("\n== two undrained layers with different pore-fluid stiffness (KV-CST-006) ==\n");
    const m::Project pr = two_layer_column(kB2);
    const auto report = katai::io::validate_project(pr);
    check(report.ok(), "the two-layer column validates");
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "the column meshed");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }

    const auto run = [&M](const m::Project& p) {
        return katai::app::solve_phases(p, M.mesh,
                                        katai::app::initial_phase_from(p.initial_procedure),
                                        nullptr, nullptr, {});
    };
    const auto res = run(pr);
    check(!res.empty() && res.back().ok, "the undrained column solves");
    if (res.empty() || !res.back().ok) return 1;
    const auto& R = res.back();

    // The closed forms, each layer with its own pore fluid. The upper layer's Kw/n is computed
    // here from B alone (Eq. 2-57), so this compares the engine's nu_u route against the
    // manual's other equation rather than against itself.
    const double Mp = m_eff(kE, kNu);
    const double kwn1 = kwn_from_nu_u(kNuU1, kE, kNu);
    const double kwn2 = kwn_from_skempton(kB2, kE, kNu);
    const double Mu1 = Mp + kwn1, Mu2 = Mp + kwn2;
    const double u_int_ex = -kQ * kH1 / Mu1;
    const double u_top_ex = u_int_ex - kQ * kH2 / Mu2;
    const double s1_ex = -Mp * kQ / Mu1, s2_ex = -Mp * kQ / Mu2;

    const double u_int = settlement_at(R, kH1);
    const double u_top = settlement_at(R, kH1 + kH2);
    const double s1 = stress_between(R, 0.2 * kH1, 0.8 * kH1);
    const double s2 = stress_between(R, kH1 + 0.2 * kH2, kH1 + 0.8 * kH2);

    std::printf("  lower layer: nu_u = %.3f  -> Kw/n = %8.0f kPa = %5.1f K'   M_u = %8.0f kPa\n",
                kNuU1, kwn1, kwn1 / k_eff(kE, kNu), Mu1);
    std::printf("  upper layer: B    = %.3f  -> Kw/n = %8.0f kPa = %5.1f K'   M_u = %8.0f kPa\n",
                kB2, kwn2, kwn2 / k_eff(kE, kNu), Mu2);
    std::printf("  u(interface) = %.6e m (exact %.6e, %+.2f%%)\n",
                u_int, u_int_ex, 100.0 * (u_int - u_int_ex) / u_int_ex);
    std::printf("  u(top)       = %.6e m (exact %.6e, %+.2f%%)\n",
                u_top, u_top_ex, 100.0 * (u_top - u_top_ex) / u_top_ex);
    std::printf("  sigma'_yy lower = %.3f kPa (exact %.3f)   upper = %.3f kPa (exact %.3f)\n",
                s1, s1_ex, s2, s2_ex);

    check(std::fabs(u_int - u_int_ex) < 0.02 * std::fabs(u_int_ex),
          "the interface settles by -q H1 / M_u1, the LOWER layer's undrained modulus");
    check(std::fabs(u_top - u_top_ex) < 0.02 * std::fabs(u_top_ex),
          "the top settles by the sum of two layers with DIFFERENT pore fluids");
    check(std::fabs(s1 - s1_ex) < 0.03 * std::fabs(s1_ex),
          "the lower layer carries -M' q / M_u1 in effective stress");
    check(std::fabs(s2 - s2_ex) < 0.03 * std::fabs(s2_ex),
          "the upper layer carries -M' q / M_u2 -- a different share, because its water differs");
    check(std::fabs(s2 - s1) > 0.5 * std::fabs(s1),
          "the two layers do not carry the same effective stress (they would if Kw/n were global)");

    // In a column the load reaching the lower layer is the surcharge, whatever happens above it.
    // Changing only the upper layer's water must therefore leave the interface EXACTLY where it
    // was -- the sharpest available statement that the parameter is per material.
    const auto res_b = run(two_layer_column(0.95));
    check(!res_b.empty() && res_b.back().ok, "the column solves with a stiffer upper pore fluid");
    if (res_b.empty() || !res_b.back().ok) return 1;
    const double u_int_b = settlement_at(res_b.back(), kH1);
    const double u_top_b = settlement_at(res_b.back(), kH1 + kH2);
    const double rel_int = std::fabs(u_int_b - u_int) / std::fabs(u_int);
    std::printf("  B(upper) 0.90 -> 0.95 : u(interface) %.9e -> %.9e (difference %.1e m, %.1e "
                "relative)\n", u_int, u_int_b, std::fabs(u_int_b - u_int), rel_int);
    std::printf("                          u(top)       %.9e -> %.9e (%.1f%% less settlement)\n",
                u_top, u_top_b, 100.0 * (1.0 - u_top_b / u_top));
    // Not bit-identical, and it should not be claimed: the two runs solve different global
    // systems, so the direct solver's arithmetic differs. What is asserted is the physics --
    // the difference is round-off next to a 40% change one layer higher up.
    check(rel_int < 1e-12,
          "changing the UPPER layer's B leaves the interface where it was, to solver round-off");
    check(std::fabs(u_top_b - u_top) > 0.05 * std::fabs(u_top),
          "and changes the surface settlement, so the parameter is doing something");

    // With both layers left at the schema default the column is the single-material one, which
    // is what every project written before this field existed asked for.
    m::Project def = two_layer_column(kB2);
    def.materials[0].und_mode = 0; def.materials[0].nu_u = 0.495;
    def.materials[1].und_mode = 0; def.materials[1].nu_u = 0.495;
    const auto res_d = run(def);
    check(!res_d.empty() && res_d.back().ok, "the default-stiffness column solves");
    if (!res_d.empty() && res_d.back().ok) {
        const double u_def = settlement_at(res_d.back(), kH1 + kH2);
        const double Mu_def = Mp + kwn_from_nu_u(0.495, kE, kNu);
        const double u_def_ex = -kQ * (kH1 + kH2) / Mu_def;
        std::printf("  both layers at the default nu_u = 0.495: u(top) = %.6e m (exact %.6e, %+.2f%%)\n",
                    u_def, u_def_ex, 100.0 * (u_def - u_def_ex) / u_def_ex);
        check(std::fabs(u_def - u_def_ex) < 0.02 * std::fabs(u_def_ex),
              "at the default the two layers are one material again");
    }

    std::printf(g_failures ? "\nFAILED (%d)\n" : "\nAll checks passed.\n", g_failures);
    return g_failures ? 1 : 0;
}
