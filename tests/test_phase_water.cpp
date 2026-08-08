// Staged dewatering: a phase may lower the water table, and the soil must respond the way
// Terzaghi says it does.
//
// This is the second capability the DGGT / Schweiger anchored-excavation benchmark asked for
// and did not find. Its second construction stage lowers the water inside the pit from -3 m to
// -17.9 m BEFORE the first cubic metre is dug, and that step is not decoration: it is what
// first loads the wall, and the settlement it causes is part of the published answer. Until
// now the phreatic surface was a property of the PROJECT, so every phase ran at the same water
// level and an excavation could not be dewatered at all.
//
// The oracle is the effective-stress principle, which is exactly what a dewatering step tests.
// Lowering the water table by dh in a laterally confined column raises the effective vertical
// stress by gamma_w*dh at every point that was below the OLD level and stays below the surface,
// so the settlement of a linear-elastic column is the integral of that increment over its
// oedometric stiffness:
//
//     ds = (1/E_oed) * INTEGRAL over the column of d(sigma'_v) dy
//
// With the level dropping from y_hi to y_lo over a column of height H on a rigid base, the
// increment is gamma_w*(y_hi - y_lo) everywhere below y_lo, and it tapers linearly from that
// value at y_lo to zero at y_hi. The integral is therefore
//
//     ds = gamma_w*(y_hi - y_lo)/E_oed * (y_lo + (y_hi - y_lo)/2)
//
// -- an area of a rectangle plus a triangle, with nothing in it that the finite element code
// could have chosen. The unit weight above the water is the same as below in this fixture, so
// the ONLY thing the dewatering changes is the buoyancy, which is the point.
//
// verify: KV-CST-003
//   oracle:   closed_form
//   source:   Terzaghi's effective-stress principle (sigma' = sigma - u) applied to a change of phreatic level in a laterally confined column; the per-phase water condition itself follows the PLAXIS staged-construction contract ("water conditions per phase"), and the case that required it is the DGGT / Schweiger triple-anchored excavation, which dewaters before excavating
//   locator:  ds = gamma_w (y_hi - y_lo) / E_oed * (y_lo + (y_hi - y_lo)/2) with E_oed = E(1-nu)/((1+nu)(1-2nu)), gamma_w = 9.81 kN/m3, measured at the top of the column between the phase that holds the high water level and the phase that lowers it (stated in full)
//   quantity: settlement of the surface of a laterally confined weightless-buoyancy column when a staged phase lowers the phreatic surface from y_hi to y_lo [m]
//   expected: the closed form above
//   band:     2%, as asserted below -- MEASURED +0.01% (0.016580 m against 0.016579 m): the effective-stress field is linear in y, so the tri6 column represents it exactly and the residual is the load-stepping tolerance. The pore pressure at the base is asserted to 1e-6 kPa against the LOWERED level, and a control phase without the override moves 6.0e-07 m -- four orders below the dewatering settlement and inside this tree's nil-phase identity tolerance

#include <katai/analysis/constants.hpp>
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

constexpr double kW = 2.0, kH = 12.0;        // column [m]
constexpr double kE = 2.0e4, kNu = 0.3;      // elastic soil
constexpr double kHi = 10.0, kLo = 3.0;      // phreatic level before / after [m]

m::Project column() {
    m::Project pr;
    pr.name = "staged dewatering";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, kW};
    pr.wy = {kHi, kHi};

    m::Material s;
    s.name = "Sand";
    s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu;
    // The same total unit weight above and below the water: then lowering the table changes the
    // EFFECTIVE stress and nothing else, which is the quantity under test.
    s.gamma_unsat = 20.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

}  // namespace

int main() {
    std::printf("== staged dewatering (KV-CST-003) ==\n");

    m::Project pr = column();
    m::Phase dewater;
    dewater.name = "Lower the water table";
    dewater.water_override = true;
    dewater.wx = {0.0, kW};
    dewater.wy = {kLo, kLo};
    pr.phases.push_back(dewater);

    const katai::io::ValidationReport rep = katai::io::validate_project(pr);
    for (const auto& i : rep.issues)
        std::printf("      %s %s: %s\n", i.severity == katai::io::Severity::Error ? "[error]  "
                                                                                  : "[warning]",
                    i.path.c_str(), i.message.c_str());
    check(rep.ok(), "the project with a per-phase water line validates");

    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + dewatering phases converged");
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return 1;
    }

    // The initial K0 phase establishes the geostatic state and must not displace: the column is
    // in equilibrium with the water at its original level.
    check(res[0].max_disp < 1e-9, "the K0 initial phase does not displace");

    const double Eoed = kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));
    const double dh = kHi - kLo;
    const double want = katai::core::kGammaWater * dh / Eoed * (kLo + 0.5 * dh);

    int top = 0;
    double bd = 1e300;
    for (int n = 0; n < res[1].mesh.node_count; ++n) {
        const double d = std::hypot(res[1].mesh.x[n] - 0.5 * kW, res[1].mesh.y[n] - kH);
        if (d < bd) { bd = d; top = n; }
    }
    const double got = -res[1].disp[top * 2 + 1];
    std::printf("  water %.1f m -> %.1f m: settlement %.6f m vs closed form %.6f m (%+.2f%%)\n",
                kHi, kLo, got, want, 100.0 * (got - want) / want);
    check(std::fabs(got - want) / want < 0.02,
          "the dewatering settlement matches the effective-stress closed form within 2%");

    // The pore pressure the phase reports must be the one it was told to use -- otherwise the
    // settlement could be right for the wrong reason.
    double p_at_base = 0.0;
    for (int n = 0; n < res[1].mesh.node_count; ++n)
        if (res[1].mesh.y[n] < 1e-9) p_at_base = std::fmax(p_at_base, res[1].pore[n]);
    const double want_pore = katai::core::kGammaWater * kLo;
    std::printf("  pore pressure at the base %.4f kPa vs %.4f kPa for the LOWERED level\n",
                p_at_base, want_pore);
    check(std::fabs(p_at_base - want_pore) < 1e-6,
          "the phase reports the pore pressure of its own water level");

    // And the control: without the override the phase runs at the project's level, so nothing
    // moves. A per-phase field that quietly applied everywhere would pass every check above.
    m::Project same = column();
    m::Phase nil;
    nil.name = "No water change";
    same.phases.push_back(nil);
    const auto M2 = katai::app::mesh_from_project(same);
    const auto res2 = katai::app::solve_phases(same, M2.mesh,
                                               katai::app::initial_phase_from(same.initial_procedure));
    check(res2.size() == 2 && res2[1].ok, "the control model solves");
    if (res2.size() == 2 && res2[1].ok) {
        std::printf("  control phase (no override) max|u| = %.3e m (dewatering moved %.3e m)\n",
                    res2[1].max_disp, got);
        // The threshold is this tree's own nil-phase identity tolerance (1e-6 m, the tripwire the
        // staged chain is pinned against elsewhere), not zero: an unchanged staged phase still
        // re-solves and carries a numerical residual. Measured 6.0e-07 m -- four orders of
        // magnitude below the settlement the dewatering causes, which is the comparison that
        // matters here.
        check(res2[1].max_disp < 1e-6,
              "a phase without the override changes nothing beyond the nil-phase residual -- "
              "the water level is per phase, not global");
        check(res2[1].max_disp < 0.01 * got,
              "and it is orders of magnitude below the dewatering it must not be confused with");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
