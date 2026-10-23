// test_axisym_safety -- a factor of safety in an axisymmetric model: a circular footing on clay.
//
// Until 2026-10-22 a Safety phase was refused in an axisymmetric model, so a circular footing, a
// tank or a shaft had no factor of safety. Every trial of the strength-reduction search is a
// nonlinear solve with the phase's options, and the r-z kinematics now reach them all.
//
// Oracle (independent path): on a weightless Tresca clay the collapse pressure of a footing is
// proportional to the undrained strength, so a footing carrying a pressure q has a factor of safety
// q_ult / q -- reached by strength reduction at constant load. The same model is loaded to
// collapse by increasing the pressure at full strength: a second path to the same number that
// shares no step with the strength-reduction search. A sanity bound against the published
// bearing capacity of a smooth circular footing on Tresca clay (Shield 1955; Cox, Eason & Hopkins
// 1961: q_ult = 5.69 s_u) is kept loose, since a flexible uniform pressure on a finite mesh is
// not the rigid punch of the closed solution.
//
// verify: KV-FND-017
//   oracle:   independent_path
//   source:   the proportionality of a Tresca collapse load to the undrained strength (a weightless soil), so that strength reduction at a constant load q finds q_ult/q; the circular smooth-footing bearing capacity q_ult = 5.69 s_u (Cox, A.D., Eason, G. & Hopkins, H.G. (1961). Axially symmetric plastic deformations in soils. Phil. Trans. R. Soc. Lond. A 254, 1-45) as a sanity bound
//   locator:  the factor of safety of a Safety phase after a Plastic phase that placed q, against the load multiplier at collapse of the same model loaded at full strength, both on the same 15-noded mesh
//   quantity: factor of safety [-]; collapse pressure over s_u [-]
//   expected: FoS x q between the collapse pressure found by loading and 3% above it; the collapse pressure within 8% of 5.69 s_u and clear of the plane-strain strip's 5.14
//   band:     one-sided 0..3% between the two paths (the load path stops at its last converged increment, a lower bound; the strength search at its plateau; measured +2.2% on a 0.5 m mesh); 8% on the published value (a flexible pressure, a finite mesh)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

static m::Project footing(bool axisymmetric, double q, bool safety) {
    constexpr double R = 6.0, Hd = 6.0, B = 1.0, su = 20.0;
    m::Project pr; pr.axisymmetric = axisymmetric;
    m::Material clay; clay.name = "Clay"; clay.model = m::SoilModel::MohrCoulomb;
    clay.E = 1.0e4; clay.nu = 0.45; clay.c = su; clay.phi = 0.0; clay.psi = 0.0;
    clay.gamma_unsat = clay.gamma_sat = 0.0; clay.e_init = 0.5; clay.tension_cutoff = false;
    pr.materials.push_back(clay);
    m::SoilPolygon P; P.name = "Clay"; P.material = 0;
    P.x = {0, R, R, 0}; P.y = {0, 0, Hd, Hd};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Footing";
    L.x1 = 0; L.y1 = Hd; L.x2 = B; L.y2 = Hd; L.qx1 = L.qx2 = 0.0; L.qy1 = L.qy2 = -q; L.coarseness = 0.25;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};
    m::Phase load; load.name = "Load"; load.type = m::PhaseType::Plastic; load.load_active = {1}; load.load_steps = 60;
    pr.phases.push_back(load);
    if (safety) {
        m::Phase s; s.name = "Safety"; s.type = m::PhaseType::Safety; s.load_active = {1};
        pr.phases.push_back(s);
    }
    return pr;
}

int main() {
    std::printf("A factor of safety in an axisymmetric model: a circular footing on Tresca clay\n\n");
    constexpr double su = 20.0;
    // (1) collapse by loading at full strength: a pressure well above the capacity, and the
    //     multiplier the phase reaches.
    const double q_big = 6.5 * su;   // just above the capacity: 60 steps of 0.11 s_u
    const m::Project p1 = footing(true, q_big, false);
    const auto M = katai::app::mesh_from_project(p1, 0.6, 15);
    check(M.ok, "the footing meshes");
    if (!M.ok) return 1;
    const auto r1 = katai::app::solve_phases(p1, M.mesh, InitialPhase::K0Procedure);
    const bool collapsed = r1.size() == 2 && !r1[1].ok && r1[1].load_factor > 0.1;
    check(collapsed, "loaded at full strength, the footing stops short of the applied pressure");
    if (r1.size() == 2) std::printf("  load path: %s\n", r1[1].message.substr(0, 160).c_str());
    const double q_ult = r1.size() == 2 ? r1[1].load_factor * q_big : 0.0;
    std::printf("  collapse pressure by loading: %.3f kPa = %.3f s_u (Cox et al. 5.69, Prandtl 5.14)\n", q_ult, q_ult / su);
    check(std::fabs(q_ult / su - 5.69) / 5.69 < 0.08, "the axisymmetric collapse pressure is the circular footing's (5.69 s_u within 8%)");
    check(q_ult / su > 5.14 * 1.04, "...and not the plane-strain strip's (5.14 s_u)");

    // (2) the same footing at half that pressure, then a Safety phase
    const double q = 0.5 * q_ult;
    const m::Project p2 = footing(true, q, true);
    const auto r2 = katai::app::solve_phases(p2, M.mesh, InitialPhase::K0Procedure);
    const bool solved = r2.size() == 3 && r2[1].ok && r2[2].ok && r2[2].fos > 0.0;
    check(solved, "the loaded footing and its Safety phase solve in axisymmetry");
    if (!solved) {
        if (!r2.empty()) std::printf("  (%s)\n", r2.back().message.substr(0, 300).c_str());
    } else {
        const double fos = r2[2].fos;
        std::printf("  factor of safety %.4f%s; FoS x q = %.3f kPa against %.3f kPa by loading (%+.2f%%)\n", fos,
                    r2[2].fos_lower_bound ? " (lower bound)" : "", fos * q, q_ult, 100.0 * (fos * q - q_ult) / q_ult);
        // The load path stops at its last converged increment, so it is a LOWER bound of the
        // collapse pressure; the strength search follows its path to the plateau. One-sided band.
        check(!r2[2].fos_lower_bound && fos * q >= 0.995 * q_ult && fos * q <= 1.03 * q_ult,
              "strength reduction finds the collapse pressure the loading found (the load path a lower bound, within 3%)");
    }
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
