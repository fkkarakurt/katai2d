// Undrained (C): the analysis that declines to separate the water from the skeleton.
//
// (A) and (B) are effective-stress analyses -- the pore fluid gets a stiffness, the excess pore
// pressure is computed, and the stresses that come out are effective. (C) is the older practice:
// enter the UNDRAINED stiffness and the UNDRAINED strength, compute in total stress, and report
// no pore pressure at all. PLAXIS keeps it because that is what a site investigation often hands
// an engineer -- Eu and su, no effective parameters -- and because the NGI-ADP family is written
// for it. KATAI had no way to say it: a total-stress model had to be entered as "Drained" with
// undrained numbers, which then had buoyancy subtracted from it and a K0 taken on the wrong
// stress. That is a wrong answer with nothing in the run to mark it.
//
// verify: KV-CST-007
//   oracle:   independent_path
//   source:   PLAXIS 2D 2025.1 Material Models Manual section 2.7 (Undrained total stress analysis) and section 2.4: an undrained Young's modulus converts to an effective one by Eq. 2-59, E' = 2(1 + nu')Eu/3, because the shear modulus is the same in both descriptions of the same soil. The two analyses of that soil must therefore agree, and the manual says exactly how far apart they are allowed to be: the total-stress route is exactly incompressible only at nu_u = 0.5, which "is not possible, since this would lead to singularity", so the difference is the compressibility that nu_u < 0.5 leaves behind. Reference Manual section 6.2.3 for the rest of the type's behaviour ("Pore pressures are not generated"; K0 refers to total stress; a consolidation calculation does not affect such a material)
//   locator:  a weightless laterally confined column loaded by a surcharge, solved twice -- (a) Undrained (A) with the effective pair (E' from Eq. 2-59, nu' = 0.3) and (b) Undrained (C) with the undrained pair (Eu, nu_u) -- for nu_u = 0.495, 0.499 and 0.4999; and a submerged K0 column of the same material read against the total-stress overburden
//   quantity: settlement of the column top under each route [m], and the initial vertical stress of a submerged column [kPa]
//   expected: u(C)/u(A) = 2(1 + nu_u)/3 exactly (the closed-form ratio of the two constrained moduli), so the two routes converge as nu_u -> 0.5: 0.33% apart at 0.495, 0.067% at 0.499, 0.0067% at 0.4999; and the Undrained (C) column's initial vertical stress is the TOTAL overburden -gamma_sat (H - y) with sigma_h = K0 sigma_v, where the same column declared Drained carries the buoyant -gamma' (H - y)
//   band:     0.2% on the ratio against its closed form, as asserted below (measured +0.00% at all three nu_u); 1e-6 relative on the total-stress overburden; the refusals are exact (message match)
//
// The identity is the point. A total-stress analysis and an effective-stress analysis of the same
// soil are two descriptions of one experiment; if the implementation of (C) were wrong -- buoyancy
// subtracted, K0 taken on effective stress, the undrained pair read as an effective one -- the two
// would not meet, and they would not meet in a way that shrinks like (0.5 - nu_u).

#include <katai/analysis/results.hpp>
#include <katai/io/validate.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kEu = 1.5e4;     // undrained Young's modulus [kPa]
constexpr double kNuEff = 0.3;    // effective Poisson's ratio of the same soil
constexpr double kH = 10.0, kW = 2.0, kQ = 50.0;

// MMM Eq. 2-59: E' = 2(1 + nu') Eu / 3. The shear modulus is common to both descriptions.
constexpr double kEeff = 2.0 * (1.0 + kNuEff) * kEu / 3.0;

double constrained(double E, double nu) {
    return E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

// A weightless confined column: the setting where the 1D closed form is exact and the only
// difference between the two runs is the description of the soil.
m::Project column(m::Drainage drainage, double E, double nu, double nu_u) {
    m::Project pr;
    pr.name = "KV-CST-007 undrained C";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;

    m::Material s;
    s.name = drainage == m::Drainage::UndrainedC ? "Clay (total stress)" : "Clay (effective)";
    s.model = m::SoilModel::LinearElastic;
    s.drainage = drainage;
    s.E = E; s.nu = nu;
    s.nu_u = nu_u;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;   // only the surcharge acts
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);

    m::Load q;
    q.kind = m::LoadKind::Distributed;
    q.name = "Surcharge";
    q.x1 = 0; q.y1 = kH; q.x2 = kW; q.y2 = kH;
    q.qx1 = q.qx2 = 0; q.qy1 = q.qy2 = -kQ;
    pr.loads.push_back(q);
    return pr;
}

// A submerged column under its own weight, K0 procedure: what the initial stress IS.
m::Project submerged(m::Drainage drainage) {
    m::Project pr;
    pr.name = "KV-CST-007 submerged";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, kW};
    pr.wy = {kH, kH};                 // water table at the ground surface
    pr.initial_procedure = m::InitialProcedure::K0Procedure;

    m::Material s;
    s.name = "Clay";
    s.model = m::SoilModel::MohrCoulomb;
    s.drainage = drainage;
    s.E = kEu; s.nu = 0.495;
    s.c = 60.0; s.phi = 0.0; s.psi = 0.0;   // su with phi = 0, as both (B) and (C) ask
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.k0_auto = false; s.k0 = 0.8;           // pinned, so the check is about the stress it multiplies
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

katai::app::SolveResult run(const m::Project& pr, bool& ok) {
    const auto M = katai::app::mesh_from_project(pr);
    ok = M.ok;
    if (!M.ok) { katai::app::SolveResult R; R.message = M.message; return R; }
    auto res = katai::app::solve_phases(pr, M.mesh,
                                        katai::app::initial_phase_from(pr.initial_procedure),
                                        nullptr, nullptr, {});
    if (res.empty()) { ok = false; katai::app::SolveResult R; R.message = "no phase result"; return R; }
    ok = res.back().ok;
    return res.back();
}

double top_settlement(const katai::app::SolveResult& R) {
    double u = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > kH - 1e-6) u = std::fmin(u, R.disp[n * 2 + 1]);
    return u;
}

}  // namespace

int main() {
    std::printf("== the same soil, described two ways (KV-CST-007) ==\n");
    std::printf("  Eu = %.0f kPa, nu' = %.2f  ->  E' = 2(1+nu')Eu/3 = %.0f kPa (MMM Eq. 2-59)\n",
                kEu, kNuEff, kEeff);

    const double nu_us[] = {0.495, 0.499, 0.4999};
    for (double nu_u : nu_us) {
        bool ok_a = false, ok_c = false;
        const auto A = run(column(m::Drainage::Undrained, kEeff, kNuEff, nu_u), ok_a);
        const auto C = run(column(m::Drainage::UndrainedC, kEu, nu_u, nu_u), ok_c);
        check(ok_a && ok_c, "both routes solve at nu_u = " + std::to_string(nu_u));
        if (!ok_a || !ok_c) { std::printf("      (%s / %s)\n", A.message.c_str(), C.message.c_str()); continue; }
        const double uA = top_settlement(A), uC = top_settlement(C);
        // The closed-form ratio of the two constrained moduli. M_u(A) = E'(1 - nu_u)/((1 + nu')
        // (1 - 2 nu_u)) once Kw/n is added to M'; M_u(C) = Eu(1 - nu_u)/((1 + nu_u)(1 - 2 nu_u)).
        // With E' from Eq. 2-59 everything cancels but 2(1 + nu_u)/3.
        const double want = 2.0 * (1.0 + nu_u) / 3.0;
        const double got = uC / uA;
        std::printf("  nu_u = %-7.4f  u(A) = %.6e m   u(C) = %.6e m   ratio %.7f "
                    "(closed form %.7f, %+.4f%%)\n",
                    nu_u, uA, uC, got, want, 100.0 * (got - want) / want);
        check(std::fabs(got - want) < 0.002 * want,
              "u(C)/u(A) = 2(1 + nu_u)/3, the compressibility nu_u < 0.5 leaves behind");
        // And the constrained moduli themselves, as the manual's two statements about one soil.
        const double MuC = constrained(kEu, nu_u);
        check(std::fabs(-kQ * kH / MuC - uC) < 0.02 * std::fabs(uC),
              "the (C) column settles by -q H / M_u with the UNDRAINED pair, within 2%");
    }
    std::printf("  the two descriptions converge as nu_u -> 0.5, which is the only place they\n"
                "  are the same statement -- and 0.5 exactly is singular, so they never meet.\n");

    // ------------------------------------------------------------------ total stress, measured --
    std::printf("\n== a submerged column: total stress, and no buoyancy subtracted ==\n");
    bool ok_tc = false, ok_dr = false;
    const auto TC = run(submerged(m::Drainage::UndrainedC), ok_tc);
    const auto DR = run(submerged(m::Drainage::Drained), ok_dr);
    check(ok_tc && ok_dr, "both submerged columns solve");
    if (ok_tc && ok_dr) {
        const double gw = katai::app::kGammaWater;
        double worst_tc = 0.0, worst_dr = 0.0, worst_h = 0.0;
        for (int n = 0; n < TC.mesh.node_count; ++n) {
            const double d = kH - TC.mesh.y[n];
            if (d < 0.5 || d > kH - 0.5) continue;          // skip the surface/base rows
            const double sv_total = -20.0 * d, sv_eff = -(20.0 - gw) * d;
            worst_tc = std::fmax(worst_tc, std::fabs(TC.stress.stress[n](1) - sv_total) /
                                               std::fabs(sv_total));
            worst_h = std::fmax(worst_h, std::fabs(TC.stress.stress[n](0) - 0.8 * sv_total) /
                                             std::fabs(0.8 * sv_total));
            worst_dr = std::fmax(worst_dr, std::fabs(DR.stress.stress[n](1) - sv_eff) /
                                               std::fabs(sv_eff));
        }
        const int base = 0;
        std::printf("  Undrained (C): sigma_v = -gamma_sat (H - y)      worst error %.2e relative\n"
                    "                 sigma_h = K0 sigma_v (K0 = 0.80)  worst error %.2e\n"
                    "  Drained      : sigma'_v = -gamma' (H - y)        worst error %.2e\n",
                    worst_tc, worst_h, worst_dr);
        check(worst_tc < 1e-6, "the (C) column carries the TOTAL overburden, unbuoyed");
        check(worst_h < 1e-6, "and its K0 multiplies that total stress (MMM 2.7)");
        check(worst_dr < 1e-6, "while the same column declared Drained carries the buoyant one");
        // The water table is still reported as a field -- it exists, the model just does not
        // subtract it from this material. Said in the diagnostics rather than left to be guessed.
        double pore_base = 0.0;
        for (int n = 0; n < TC.mesh.node_count; ++n)
            if (TC.mesh.y[n] < 1e-6) pore_base = std::fmax(pore_base, TC.pore[n]);
        std::printf("  reported water-table field at the base: %.2f kPa (hydrostatic %.2f) -- "
                    "displayed, not subtracted\n", pore_base, gw * kH);
        check(std::fabs(pore_base - gw * kH) < 1e-6, "the water table itself is unchanged by the "
                                                     "material's drainage type");
        (void)base;
        bool told = false;
        for (const auto& d : TC.diagnostics) told |= std::string(d.code) == "K2D-M003";
        check(told, "the run says the stresses of this material are total (K2D-M003)");
    }

    // ---------------------------------------------------------------------------- the refusals --
    std::printf("\n== what Undrained (C) is refused for ==\n");
    {   // Hardening Soil: PLAXIS offers (C) for Linear Elastic and Mohr-Coulomb only.
        m::Project pr = column(m::Drainage::UndrainedC, kEu, 0.495, 0.495);
        pr.materials[0].model = m::SoilModel::HardeningSoil;
        pr.materials[0].c = 60.0; pr.materials[0].phi = 0.0;
        const auto rep = katai::io::validate_project(pr);
        check(!rep.ok(), "Hardening Soil + Undrained (C) is refused by the input contract");
        bool ok_hs = false;
        const auto R = run(pr, ok_hs);
        check(!ok_hs && R.message.find("Undrained (C)") != std::string::npos,
              "and by the constitutive catalogue, by name: " + R.message.substr(0, 60) + "...");
    }
    {   // Consolidation: there are no pore pressures in it to consolidate.
        m::Project pr = submerged(m::Drainage::UndrainedC);
        pr.materials[0].kx = pr.materials[0].ky = 1e-3;
        m::Phase con;
        con.name = "Consolidate";
        con.type = m::PhaseType::Consolidation;
        con.duration = 10.0;
        pr.phases.push_back(con);
        const auto M = katai::app::mesh_from_project(pr);
        check(M.ok, "the consolidation model meshes");
        if (M.ok) {
            const auto res = katai::app::solve_phases(pr, M.mesh,
                katai::app::initial_phase_from(pr.initial_procedure), nullptr, nullptr, {});
            const bool refused = res.size() >= 2 && !res.back().ok &&
                                 res.back().message.find("Undrained (C)") != std::string::npos;
            check(refused, "a consolidation phase on an Undrained (C) material is refused: " +
                               (res.size() >= 2 ? res.back().message.substr(0, 70) : std::string()) +
                               "...");
        }
    }

    std::printf(g_failures ? "\nFAILED (%d)\n" : "\nAll checks passed.\n", g_failures);
    return g_failures ? 1 : 0;
}
