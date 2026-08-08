// Two staged-construction options that decide what a phase IS, not how it is solved.
//
// A phase in PLAXIS carries more than a list of what is switched on. It carries Sum-Mstage -- how
// much of the change is actually applied -- and "Ignore und. behaviour", which lets a phase treat
// undrained soil as drained. Both are ordinary practice: an excavation lift taken half way to see
// where the wall is going, an initial state established without generating excess pore pressures
// that the in-situ soil never had. Neither could be asked for in a .k2d until now, and the
// consequence was not a missing convenience: a partial stage had to be faked by editing geometry,
// and a phase that should have been drained had to have its materials retyped -- which changes the
// material everywhere, including in the phases where it must stay undrained.
//
// Both are verified here by IDENTITIES rather than by tables, because both have an exact answer.
//
// verify: KV-EXC-002
//   oracle:   closed_form
//   source:   linear superposition: for a linear-elastic soil the staged-construction ramp is a linear operator on the configuration imbalance, so applying a fraction of it must give exactly that fraction of the response -- the definition of PLAXIS's Sum-Mstage as "the proportion of the unbalanced force that has been applied" (PLAXIS 2D 2025.1 Reference Manual, staged construction), stated as a testable equality
//   locator:  an excavation on a linear-elastic soil, the digging phase run at Sum-Mstage = 1, 0.5 and 0.25, reading the heave of the excavation floor; and the same phase run in two halves (0.5 followed by the remainder) against the whole stage in one
//   quantity: vertical displacement of the excavation floor [m]
//   expected: u(0.5) = 0.5 u(1) and u(0.25) = 0.25 u(1) to solver precision on a linear material; and a stage taken in two steps arrives where the single step does
//   band:     1e-12 relative on the proportionality (a linear system solved directly twice), 1e-9 m on the two-step path -- both asserted below and measured, not inherited
//
// verify: KV-CST-005
//   oracle:   independent_path
//   source:   PLAXIS 2D 2025.1 Reference Manual, "Ignore und. behaviour (A,B)": the undrained response of Undrained (A)/(B) clusters is temporarily excluded so that no excess pore pressure is generated, while the strength parameters stay as declared. The independent path is the same model with the material DECLARED drained: if the option means what the manual says, the two must be the same calculation
//   locator:  a loaded undrained column solved (a) with drainage = Undrained (A) and the phase ignoring undrained behaviour, and (b) with drainage = Drained, same everything else
//   quantity: settlement of the column top [m], with the steady pore field read alongside to show which part of the water response the option touches
//   expected: bit-identical displacements between the ignoring run and the declared-drained run; the undrained run must differ substantially (otherwise the switch has nothing to switch); the hydrostatic field is untouched
//   band:     exact -- 0.0e+00 between the two routes, asserted below. Measured: undrained (A) 0.001030 m, ignoring 0.029714 m, declared drained 0.029714 m -- the undrained column settles 3.5% of the drained one, which is the excess pore pressure the option removes

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

// A 20 x 10 m block of linear-elastic soil in two layers: the upper one is dug out in phase 1.
// Linear elastic on purpose -- the identity being tested is linearity itself, and a plastic soil
// would answer a different (and much less sharp) question.
m::Project excavation() {
    m::Project pr;
    pr.name = "KV-EXC-002 partial stage";
    pr.x_min = 0; pr.x_max = 20; pr.y_min = 0; pr.y_max = 10;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;

    m::Material s;
    s.name = "Soil";
    s.model = m::SoilModel::LinearElastic;
    s.E = 2.0e4; s.nu = 0.3;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);

    m::SoilPolygon lower;
    lower.name = "Lower";
    lower.material = 0;
    lower.x = {0, 20, 20, 0};
    lower.y = {0, 0, 6, 6};
    lower.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(lower);

    m::SoilPolygon upper;
    upper.name = "Dug out";
    upper.material = 0;
    upper.x = {0, 20, 20, 0};
    upper.y = {6, 6, 10, 10};
    upper.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(upper);

    m::Phase dig;
    dig.name = "Excavate";
    dig.type = m::PhaseType::Plastic;
    dig.poly_active = {1, 0};      // the upper layer goes
    pr.phases.push_back(dig);
    return pr;
}

// Heave of the excavation floor: the node nearest the middle of the new surface.
double floor_heave(const katai::app::SolveResult& R) {
    int best = 0;
    double bd = 1e300;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        const double d = std::hypot(R.mesh.x[n] - 10.0, R.mesh.y[n] - 6.0);
        if (d < bd) { bd = d; best = n; }
    }
    return R.disp[best * 2 + 1];
}

// A 1 x 8 m column of clay under a surcharge -- the standard undrained-versus-drained comparison.
m::Project column(m::Drainage drainage, bool ignore_undrained) {
    m::Project pr;
    pr.name = "KV-CST-005 ignore undrained";
    pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 8;
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, 1.0};
    pr.wy = {8.0, 8.0};
    pr.initial_procedure = m::InitialProcedure::K0Procedure;

    m::Material clay;
    clay.name = "Clay";
    clay.model = m::SoilModel::MohrCoulomb;
    clay.drainage = drainage;
    clay.E = 5.0e3; clay.nu = 0.3;
    clay.c = 20.0; clay.phi = 25.0; clay.psi = 0.0;
    clay.gamma_unsat = 17.0; clay.gamma_sat = 19.0;
    clay.kx = 1e-4; clay.ky = 1e-4;
    pr.materials.push_back(clay);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, 1, 1, 0};
    P.y = {0, 0, 8, 8};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::Load q;
    q.kind = m::LoadKind::Distributed;
    q.name = "Surcharge";
    q.x1 = 0.0; q.y1 = 8.0; q.x2 = 1.0; q.y2 = 8.0;
    q.qy1 = -25.0; q.qy2 = -25.0;
    pr.loads.push_back(q);

    pr.initial.load_active = {0};
    m::Phase load;
    load.name = "Surcharge";
    load.type = m::PhaseType::Plastic;
    load.load_active = {1};
    load.ignore_undrained = ignore_undrained;
    pr.phases.push_back(load);
    return pr;
}

}  // namespace

int main() {
    std::printf("== a partial staged change (KV-EXC-002) ==\n");

    const m::Project full = excavation();
    check(katai::io::validate_project(full).ok(), "the excavation validates");
    const auto M = katai::app::mesh_from_project(full);
    check(M.ok, "excavation meshed");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }

    const auto run = [&M](const m::Project& pr) {
        return katai::app::solve_phases(pr, M.mesh,
                                        katai::app::initial_phase_from(pr.initial_procedure),
                                        nullptr, nullptr, {});
    };

    const auto whole = run(full);
    check(whole.size() == 2 && whole[1].ok, "the whole stage solves");
    if (whole.size() != 2 || !whole[1].ok) return 1;
    const double u_full = floor_heave(whole[1]);
    std::printf("  Sum-Mstage 1.00 -> floor heave %.9f m\n", u_full);

    // For a linear material the ramp is a linear operator, so a fraction of the stage must give
    // exactly that fraction of the response. This is the sharpest statement the option admits,
    // and it is the one that would fail if the fraction scaled the wrong vector -- the loads but
    // not the imbalance, say, or the imbalance but not the structural weight.
    const double fractions[2] = {0.5, 0.25};
    for (double frac : fractions) {
        m::Project part = full;
        part.phases[0].sum_mstage = frac;
        check(katai::io::validate_project(part).ok(),
              "a partial stage validates at " + std::to_string(frac));
        const auto res = run(part);
        check(res.size() == 2 && res[1].ok, "the partial stage solves");
        if (res.size() != 2 || !res[1].ok) return 1;
        const double u = floor_heave(res[1]);
        const double rel = std::fabs(u - frac * u_full) / std::fabs(u_full);
        std::printf("  Sum-Mstage %.2f -> floor heave %.9f m   (%.2f x the whole stage, error "
                    "%.2e)\n", frac, u, u / u_full, rel);
        check(rel < 1e-12, "a fraction of the stage gives exactly that fraction of the response");
        // And the run says so, rather than leaving a half-dug pit looking like a finished one.
        bool told = false;
        for (const auto& d : res[1].diagnostics) told |= std::string(d.code) == "K2D-A007";
        check(told, "the run reports that the stage was only partly applied (K2D-A007)");
    }

    // The other half of the same idea: a stage taken in two steps must arrive where the single
    // step does. This is what makes a partial stage a construction step rather than a scaling
    // trick -- the next phase continues from the partly-changed state.
    m::Project halves = full;
    halves.phases[0].sum_mstage = 0.5;
    m::Phase rest = halves.phases[0];
    rest.name = "Finish the dig";
    rest.sum_mstage = 1.0;
    halves.phases.push_back(rest);
    const auto two_step = run(halves);
    check(two_step.size() == 3 && two_step[2].ok, "the two-step excavation solves");
    if (two_step.size() == 3 && two_step[2].ok) {
        // Phase displacements are reported relative to the phase's own start, so the two steps
        // are added to compare with the single step (the convention pinned by KV-CST-002).
        const double u_two = floor_heave(two_step[1]) + floor_heave(two_step[2]);
        std::printf("  half then the rest -> %.9f m (one step: %.9f m, difference %.2e m)\n",
                    u_two, u_full, std::fabs(u_two - u_full));
        check(std::fabs(u_two - u_full) < 1e-9,
              "a stage taken in two steps arrives where the single step does");
    }

    // Out of range is refused, not clamped: 1.4 means the author expected something else.
    m::Project bad = full;
    bad.phases[0].sum_mstage = 1.4;
    check(!katai::io::validate_project(bad).ok(), "a target above 1 is refused");
    m::Project bad_initial = full;
    bad_initial.initial.sum_mstage = 0.5;
    check(!katai::io::validate_project(bad_initial).ok(),
          "and a partial INITIAL phase is refused: scaled gravity is not a construction step");

    // ------------------------------------------------------------------------------------
    std::printf("\n== ignoring undrained behaviour (KV-CST-005) ==\n");

    const m::Project und = column(m::Drainage::Undrained, false);
    const m::Project ign = column(m::Drainage::Undrained, true);
    const m::Project dry = column(m::Drainage::Drained, false);
    check(katai::io::validate_project(ign).ok(), "the ignore-undrained project validates");

    const auto CM = katai::app::mesh_from_project(und);
    check(CM.ok, "column meshed");
    if (!CM.ok) return 1;
    const auto run_col = [&CM](const m::Project& pr) {
        return katai::app::solve_phases(pr, CM.mesh,
                                        katai::app::initial_phase_from(pr.initial_procedure),
                                        nullptr, nullptr, {});
    };

    const auto R_und = run_col(und);
    const auto R_ign = run_col(ign);
    const auto R_dry = run_col(dry);
    check(R_und.size() == 2 && R_und[1].ok, "the undrained column solves");
    check(R_ign.size() == 2 && R_ign[1].ok, "the ignore-undrained column solves");
    check(R_dry.size() == 2 && R_dry[1].ok, "the drained column solves");
    if (R_und.size() != 2 || R_ign.size() != 2 || R_dry.size() != 2) return 1;

    const auto top_settlement = [](const katai::app::SolveResult& R) {
        int best = 0;
        double bd = 1e300;
        for (int n = 0; n < R.mesh.node_count; ++n) {
            const double d = std::hypot(R.mesh.x[n] - 0.5, R.mesh.y[n] - 8.0);
            if (d < bd) { bd = d; best = n; }
        }
        return -R.disp[best * 2 + 1];
    };
    const double s_und = top_settlement(R_und[1]);
    const double s_ign = top_settlement(R_ign[1]);
    const double s_dry = top_settlement(R_dry[1]);
    std::printf("  undrained (A)                 settlement %.9f m\n", s_und);
    std::printf("  undrained (A), phase ignores  settlement %.9f m\n", s_ign);
    std::printf("  material declared drained     settlement %.9f m\n", s_dry);
    std::printf("  ignore-undrained vs drained   %.3e m\n", std::fabs(s_ign - s_dry));
    check(s_ign == s_dry,
          "ignoring undrained behaviour IS the drained calculation, bit for bit");
    check(std::fabs(s_und - s_dry) / s_dry > 0.05,
          "and it is not the undrained one: the two differ by more than 5%, so the switch has "
          "something to switch");

    // Where the excess pore pressure went, said precisely. `SolveResult::pore` is the STEADY
    // (hydrostatic) field, which both runs share and which this option does not touch -- 8 m of
    // water is 78.5 kPa either way. The excess lives in the effective-stress state, and its
    // absence is exactly what the settlement identity above measures: an undrained column barely
    // moves under load because the water carries it, and this one moved the drained distance.
    double pp_ign = 0.0, pp_und = 0.0;
    for (double p : R_ign[1].pore) pp_ign = std::fmax(pp_ign, std::fabs(p));
    for (double p : R_und[1].pore) pp_und = std::fmax(pp_und, std::fabs(p));
    std::printf("  steady pore field (unchanged by the option): ignoring %.2f kPa, undrained "
                "%.2f kPa\n", pp_ign, pp_und);
    check(std::fabs(pp_ign - pp_und) < 1e-9 && std::fabs(pp_ign - 9.81 * 8.0) < 0.1,
          "the hydrostatic pore field is the same 8 m of water in both, as it must be");
    std::printf("  the undrained run settles %.1f%% of the drained one -- the excess pore "
                "pressure the option removes\n", 100.0 * s_und / s_dry);

    bool told = false;
    for (const auto& d : R_ign[1].diagnostics) told |= std::string(d.code) == "K2D-A008";
    check(told, "the run reports which material was solved drained (K2D-A008)");
    bool quiet = true;
    for (const auto& d : R_und[1].diagnostics) quiet &= std::string(d.code) != "K2D-A008";
    check(quiet, "and says nothing when the option is off");

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
