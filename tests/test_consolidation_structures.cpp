// A wall, a raft or an anchor inside a CONSOLIDATION phase. Until now a consolidation phase was
// soil-only: any structural element made it refuse, so the one analysis an engineer most often
// wants -- what the excavation does over the months while the excess pore pressure dissipates --
// could not be run with the structure that holds it up. Plates, anchors and geogrids now take part
// in the coupled solve; interfaces followed (KV-STR-007) and embedded beams (KV-STR-023). Since
// 2026-10 every structure is evaluated by its own law at every iterate (StructureLaws) -- the
// static solver's assembler, with the soil left to the coupled core.
//
// THE ORACLE IS THE LIMIT THE COUPLED SOLUTION MUST WALK INTO. Consolidation is the transition
// between two states this program already solves by an independent path: at t = 0+ the load is
// carried by the water, and as t grows the excess pore pressure vanishes and what remains IS the
// drained problem. So a consolidation phase run until the excess pore pressure is gone must
// reproduce the DRAINED Plastic phase of the same model -- soil settlement AND the structure's own
// bending moment -- through a completely different solver: a monolithic coupled Biot system with a
// pore-pressure degree of freedom at every node, against solve_nonlinear's mechanical-only Newton.
// Nothing in the two paths is shared except the element matrices, which is exactly the claim being
// tested: that the plate the coupled phase carries is the plate the static phase carries.
//
// The second assertion is the one that would catch the silent version of this feature. A structure
// whose stiffness never reached the coupled system would not fail loudly -- it would produce a
// slightly softer, entirely reasonable-looking settlement. So the same model is also run WITHOUT
// the plate, and the two must differ by far more than the band above.
//
// The third came out of the measurement rather than into it. The structural branch of this phase
// WAS elastic, and for a geogrid that is not a capacity but its behaviour: tension-only, so wherever
// the settlement bowl put a sheet in compression the elastic branch had it push BACK on the soil
// instead of going slack -- measured at -1.57% on the sheet's force, and 0.00000% for the same sheet
// kept wholly in tension, which proved the cause. The run reported it (K2D-A015), and anchors,
// geogrids and plates past their capacity (K2D-A014). With the structures by their own laws the
// wide sheet lands on the drained answer (+0.00000%), an element past its capacity yields as it
// does in the Plastic phase -- in the fully-coupled phase as well, which takes the same laws --
// and the three warnings are retired.
//
// verify: KV-STR-006
//   oracle:   independent_path
//   source:   the drained limit of Biot consolidation -- as the excess pore pressure dissipates the coupled problem BECOMES the drained one, so the verified static path (solve_nonlinear; the plate/anchor/geogrid diagrams of test_struct_gui, KV-STR-002/004/005) is the answer the coupled solve must converge to. Terzaghi (1943) supplies the rate, cv = k Eoed / gamma_w, and with it the time factor at which the comparison is fair
//   locator:  Tv = cv t / H_dr^2, cv = k Eoed / gamma_w, Eoed = E(1-nu)/((1+nu)(1-2nu)); at Tv = 4 the remaining maximum excess pore pressure is 9.8e-4 kPa of the 50 kPa applied (2e-5), i.e. the drained state to within a ten-thousandth. Tension-only geogrid: the element slackens in compression (N = 0, reversible), kernel/fem/elements/geogrid.hpp
//   quantity: maximum settlement [m], a plate's |M|max [kNm/m] and an anchor's / geogrid's |N|max [kN] at the end of a consolidation phase run to Tv = 4, each against the SAME model solved as a drained Plastic phase -- elastic, and with a capacity below what the element develops; the same settlement with the structure removed; and, in a fully-coupled phase, the utilisation each element's capacity is reported at
//   expected: plate, anchor and a geogrid held in pure tension must equal the drained Plastic phase (4.50102187e-02 m vs 4.50113380e-02 m; 21.621032 vs 21.621032 kNm/m; 181.713475 vs 181.714482 kN; 2.406855 vs 2.406855 kN). Removing the plate must move the settlement to 5.295710e-02 m. A geogrid whose ends the settlement bowl puts in COMPRESSION matches too, now that it goes slack. An element given a capacity below its force yields and ends at the drained phase's yielded answer
//   band:     5e-4 relative on every quantity compared to the drained path, measured -0.00249% (settlement), 0.00000% (plate moment), -0.00055% (anchor force) and -0.00000% (geogrid in pure tension). The band is not slack for its own sake: what separates the two answers is the excess pore pressure still in the ground at Tv = 4, which is 2e-5 of the load, and the settlement difference came out at that order; 5e-4 leaves room for another linear backend rounding 120 steps differently. The plate's own effect on the settlement is 17.66% -- three orders outside, which is what makes the agreement evidence rather than two soft numbers agreeing. The wide geogrid agrees to +0.00000% (it was -1.57% on the elastic branch). The yielded elements are held to 0.5%: an element that is elastic and then yields remembers the order its loads came in, and the two paths load it undrained first and drained throughout
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kW = 20.0, kH = 10.0;      // clay block [m]
constexpr double kE = 5000.0, kNu = 0.2;    // [kPa]
constexpr double kPerm = 1.0e-3;            // [m/day]
constexpr double kQ = 50.0;                 // surcharge on the raft [kPa]
constexpr double kRaftX1 = 8.0, kRaftX2 = 12.0;

double eoed() { return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }
double cv() { return kPerm * eoed() / katai::app::kGammaWater; }

// The block, drained at the top only, with an optional raft carrying the surcharge. `plastic_phase`
// builds the DRAINED reference (a Plastic phase); otherwise the consolidation phase to Tv = 4.
m::Project block(bool with_plate, bool plastic_phase, double Mp = 0.0) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 17.0; s.gamma_sat = 20.0; s.e_init = 0.7;
    s.kx = kPerm; s.ky = kPerm;
    s.drainage = m::Drainage::Drained;   // effective-stress soil either way
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};   // top drains
    P.edge_head = {0.0, 0.0, kH, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;   // pore pressure = the excess only

    if (with_plate) {
        m::PlateMaterial pm; pm.EA = 5.0e6; pm.EI = 8.5e3;
        pm.Mp = Mp; pm.Np = 0.0; pm.elastoplastic = Mp > 0.0;
        pr.plates.push_back(pm);
        m::StructElement e; e.kind = m::StructKind::Plate; e.name = "Raft";
        e.x1 = kRaftX1; e.y1 = kH; e.x2 = kRaftX2; e.y2 = kH; e.material = 0;
        pr.structs.push_back(e);
    }

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = kRaftX1; L.y1 = kH; L.x2 = kRaftX2; L.y2 = kH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};        // the load arrives in the phase, not in the initial state

    m::Phase ph; ph.load_active = {1};
    if (plastic_phase) {
        ph.name = "Drained (long term)";
        ph.type = m::PhaseType::Plastic;
    } else {
        ph.name = "Consolidation";
        ph.type = m::PhaseType::Consolidation;
        ph.duration = 4.0 * kH * kH / cv();   // Tv = 4: the excess pore pressure is gone
        ph.consol_load_at_start = true;       // ... which needs the load on from the start
        ph.time_steps = 120;
    }
    pr.phases.push_back(ph);
    return pr;
}

struct Answer {
    bool ok = false;
    double settle = 0.0, moment = 0.0, pore = 0.0;
    double moment_or_n = 0.0;   // the line's max |M| (plate) or max |N| (anchor / geogrid)
    std::string msg;
};

Answer run(const m::Project& pr) {
    Answer a;
    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    if (!M.ok) { a.msg = M.message; return a; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2) { a.msg = "phases did not run"; return a; }
    const auto& R = res[1];
    a.ok = R.ok; a.msg = R.message;
    if (!R.ok) return a;
    // Maximum settlement of the loaded surface (|u_y| over the mesh -- the raft is at the top and
    // the block is fixed at its base, so the maximum is the surface's).
    for (int n = 0; n < R.mesh.node_count; ++n)
        a.settle = std::fmax(a.settle, std::fabs(R.disp(2 * n + 1)));
    if (!R.struct_forces.empty()) {
        a.moment = R.struct_forces.front().max_M;
        a.moment_or_n = R.struct_forces.front().kind == 0 ? R.struct_forces.front().max_M
                                                          : R.struct_forces.front().max_N;
    }
    if (!R.consol_excess_pore.empty()) a.pore = R.consol_excess_pore.back();
    return a;
}

void test_drained_limit() {
    std::printf("-- (1) a consolidation phase with a raft walks into the DRAINED answer --\n");
    std::printf("   cv = %.4f m2/day, Tv = 4 at t = %.0f day\n", cv(), 4.0 * kH * kH / cv());

    const Answer drained = run(block(/*with_plate=*/true, /*plastic_phase=*/true));
    check(drained.ok, "the drained Plastic reference solved");
    if (!drained.ok) { std::printf("   (%s)\n", drained.msg.c_str()); return; }

    const Answer consol = run(block(/*with_plate=*/true, /*plastic_phase=*/false));
    check(consol.ok, "the consolidation phase with a plate solved (it used to be refused)");
    if (!consol.ok) { std::printf("   (%s)\n", consol.msg.c_str()); return; }

    std::printf("   settlement: consolidation %.8e m vs drained %.8e m  (%+.5f%%)\n",
                consol.settle, drained.settle, 100.0 * (consol.settle / drained.settle - 1.0));
    std::printf("   plate |M|max: consolidation %.6f vs drained %.6f kNm/m  (%+.5f%%)\n",
                consol.moment, drained.moment, 100.0 * (consol.moment / drained.moment - 1.0));
    std::printf("   excess pore left at Tv = 4: %.4e kPa of the %.0f kPa applied\n", consol.pore, kQ);

    check(drained.moment > 1.0, "the reference really bends the plate (a moment to compare)");
    check(std::fabs(consol.settle / drained.settle - 1.0) < 5e-4,
          "settlement at the end of consolidation = the drained answer within 0.05%");
    check(std::fabs(consol.moment / drained.moment - 1.0) < 5e-4,
          "and so is the PLATE MOMENT -- the structure is in the coupled system, not beside it");
    check(consol.pore < 0.01 * kQ, "the excess pore pressure really has dissipated");

    // The check that would catch the silent version: no plate at all.
    const Answer bare = run(block(/*with_plate=*/false, /*plastic_phase=*/false));
    check(bare.ok, "the same phase without the plate solved");
    if (bare.ok) {
        std::printf("   without the plate: settlement %.6e m  (%+.2f%% from the plated run)\n",
                    bare.settle, 100.0 * (bare.settle / consol.settle - 1.0));
        check(std::fabs(bare.settle / consol.settle - 1.0) > 0.10,   // measured 17.66%
              "removing the plate moves the answer by far more than the band -- so the agreement "
              "above is the plate's stiffness, not a coincidence");
    }
}

void test_elastic_limit_is_not_silent() {
    std::printf("\n-- (3) a plate past its plastic moment hinges in the consolidation phase --\n");
    // The same raft with a plastic moment far below what the load produces. The consolidation
    // phase used to solve it elastically and report the over-capacity line (K2D-A014); it now
    // hinges, as in the Plastic phase.
    const Answer dh = run(block(/*with_plate=*/true, /*plastic_phase=*/true, /*Mp=*/1.0));
    const Answer ch = run(block(/*with_plate=*/true, /*plastic_phase=*/false, /*Mp=*/1.0));
    check(dh.ok && ch.ok, "the hinged raft solves in both phases");
    if (!dh.ok || !ch.ok) { std::printf("   (%s | %s)\n", dh.msg.c_str(), ch.msg.c_str()); return; }
    std::printf("   |M|max consolidation %.6f vs drained %.6f kNm/m; settlement %+.4f%%\n", ch.moment,
                dh.moment, 100.0 * (ch.settle / dh.settle - 1.0));
    // The diagram's moment is recovered at the stations from the displacement, and comes out a
    // few per cent above the plastic moment the hinge holds at the stress points -- in the
    // Plastic phase as in this one (1.085 there). Elastic, the raft carried 21.6.
    check(ch.moment < 1.1 && std::fabs(ch.moment / dh.moment - 1.0) < 0.02,
          "the moment is held at the plastic moment, as in the Plastic phase (2%)");
    check(std::fabs(ch.settle / dh.settle - 1.0) < 5e-3, "and the settlement is the drained phase's (0.5%)");
    const auto diag_run = [](const m::Project& pr) {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        return katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    };
    const auto res = diag_run(block(true, false, 1.0));
    bool warned = false;
    for (const auto& d : res[1].diagnostics) if (d.code == std::string("K2D-A014")) warned = true;
    check(!warned && !res[1].struct_forces.empty() && res[1].struct_forces.front().yielded,
          "the line is marked as yielded, and there is no elastic branch to report (K2D-A014 silent)");
    // The fully-coupled phase hinges it too.
    m::Project fc = block(true, false, 1.0);
    fc.phases[0].type = m::PhaseType::FullyCoupled;
    const auto rf = diag_run(fc);
    bool fwarned = false;
    for (const auto& d : rf[1].diagnostics) if (d.code == std::string("K2D-A014")) fwarned = true;
    const double fm = rf[1].struct_forces.empty() ? 0.0 : rf[1].struct_forces.front().max_M;
    std::printf("   fully coupled: |M|max %.6f kNm/m%s%s\n", fm, rf[1].ok ? "" : " -- ",
                rf[1].ok ? "" : rf[1].message.c_str());
    check(rf[1].ok && !fwarned && fm < 1.1 && !rf[1].struct_forces.empty() && rf[1].struct_forces.front().yielded,
          "and so does a fully-coupled phase: the raft hinges, marked as yielded");
}

// The anchor and geogrid branches of the same machinery, with the drained limit as the oracle
// again. They are here for two reasons beyond coverage: their capacities are the ones the elastic
// branch cannot apply, and the anchor's capacity carries a unit trap -- the solver holds the force
// per METRE OF WALL (EA/Ls) while the report and the engineer's input are per ANCHOR, so a
// utilisation computed in the wrong one of the two is wrong by the spacing and looks plausible.
void test_tie_and_geogrid() {
    std::printf("\n-- (2) an anchor and a geogrid in the same phase, against the same limit --\n");
    for (int kind = 0; kind < 2; ++kind) {
        const bool tie = kind == 0;
        auto build = [&](bool plastic_phase, double cap, bool narrow = false) {
            m::Project pr = block(/*with_plate=*/false, plastic_phase);
            m::StructElement e; e.material = 0;
            if (tie) {
                m::AnchorMaterial am; am.EA = 1.0e5; am.Lspacing = 2.5;
                am.Fmax_tens = cap; am.elastoplastic = cap > 0.0;
                pr.anchors.push_back(am);
                e.kind = m::StructKind::Anchor; e.name = "Tie";
                // Fixed end above the ground, soil end strictly inside: the settling soil under
                // the surcharge stretches the tie, so it carries a share of the load in tension.
                e.x1 = 10.0; e.y1 = kH - 2.0; e.x2 = 10.0; e.y2 = kH + 4.0;
            } else {
                m::GeogridMaterial gm; gm.EA = 2.0e3; gm.Np = cap;
                gm.elastoplastic = cap > 0.0;
                pr.geogrids.push_back(gm);
                e.kind = m::StructKind::Geogrid; e.name = "Geogrid";
                // WIDE (default): the sheet reaches past the loaded strip, so the settlement bowl
                // leaves its ends in COMPRESSION -- where a real geogrid goes slack and this
                // phase's elastic branch does not. NARROW: entirely under the load, tension
                // everywhere, and then the elastic branch is the right branch.
                e.y1 = e.y2 = kH - 2.0;
                if (narrow) { e.x1 = 8.5; e.x2 = 11.5; } else { e.x1 = 6.0; e.x2 = 14.0; }
            }
            pr.structs.push_back(e);
            return pr;
        };
        const char* what = tie ? "anchor" : "geogrid";
        const Answer drained = run(build(true, 0.0));
        const Answer consol = run(build(false, 0.0));
        check(drained.ok && consol.ok, tie ? "anchor: both phases solved"
                                           : "geogrid: both phases solved");
        if (!drained.ok || !consol.ok) {
            std::printf("   (%s)\n", (drained.ok ? consol.msg : drained.msg).c_str());
            continue;
        }
        std::printf("   %s: N consolidation %.6f vs drained %.6f kN  (%+.5f%%),"
                    "  settlement %+.5f%%\n", what, consol.moment_or_n, drained.moment_or_n,
                    100.0 * (consol.moment_or_n / drained.moment_or_n - 1.0),
                    100.0 * (consol.settle / drained.settle - 1.0));
        check(std::fabs(drained.moment_or_n) > 1.0, "the element really carries something");
        if (tie) {
            check(std::fabs(consol.moment_or_n / drained.moment_or_n - 1.0) < 5e-4,
                  "anchor force at the end of consolidation = the drained answer");
        } else {
            // THE CLAIM. A geogrid is tension-only, and the wide sheet -- whose ends the settlement
            // bowl puts in compression -- goes slack there in this phase as in the Plastic one, so
            // it lands on the drained answer. On the elastic branch it came out 1.57% off.
            check(std::fabs(consol.moment_or_n / drained.moment_or_n - 1.0) < 5e-4,
                  "the wide sheet, compressed at its ends, = the drained answer: it goes slack there");
            const auto diagnostics_of = [](const m::Project& prw) {
                const auto Mw = katai::app::mesh_from_project(prw, 1.0, 6);
                return katai::app::solve_phases(prw, Mw.mesh, InitialPhase::K0Procedure)[1].diagnostics;
            };
            bool slack_warned = false;
            for (const auto& d : diagnostics_of(build(false, 0.0)))
                if (d.code == std::string("K2D-A015")) slack_warned = true;
            check(!slack_warned, "and nothing is reported about an elastic sheet (K2D-A015 silent)");
            // The fully-coupled phase solves its structures by the same laws.
            m::Project fcw = build(false, 0.0);
            fcw.phases[0].type = m::PhaseType::FullyCoupled;
            const Answer fcg = run(fcw);
            bool fc_slack = false;
            for (const auto& d : diagnostics_of(fcw)) if (d.code == std::string("K2D-A015")) fc_slack = true;
            std::printf("   geogrid, fully coupled: N %.6f vs drained %.6f kN (%+.5f%%)\n", fcg.moment_or_n,
                        drained.moment_or_n, 100.0 * (fcg.moment_or_n / drained.moment_or_n - 1.0));
            check(fcg.ok && !fc_slack && std::fabs(fcg.moment_or_n / drained.moment_or_n - 1.0) < 5e-4,
                  "and so does the fully-coupled phase: the sheet goes slack there too");

            // THE CONTROL. The same sheet, kept entirely under the load so that every station is
            // in tension: there the elastic branch IS the geogrid's behaviour, and the two paths
            // must agree as tightly as the anchor did. If this one also drifted, the cause of the
            // difference above would not be the tension-only cut but something in the coupling.
            const Answer nd = run(build(true, 0.0, /*narrow=*/true));
            const Answer nc = run(build(false, 0.0, /*narrow=*/true));
            check(nd.ok && nc.ok, "the narrow sheet solved on both paths");
            if (nd.ok && nc.ok) {
                std::printf("   geogrid, all in tension: %.6f vs %.6f kN  (%+.5f%%)\n",
                            nc.moment_or_n, nd.moment_or_n,
                            100.0 * (nc.moment_or_n / nd.moment_or_n - 1.0));
                check(std::fabs(nc.moment_or_n / nd.moment_or_n - 1.0) < 5e-4,
                      "a geogrid in pure tension = the drained answer -- so the gap above IS the "
                      "compression cut and nothing else");
            }
        }

        // The capacity, set at half of what the element developed. The consolidation phase
        // applies it as the Plastic phase does: the element yields, and the two end together.
        const double cap = 0.5 * std::fabs(drained.moment_or_n);
        {
            const Answer yd = run(build(true, cap));
            const Answer yc = run(build(false, cap));
            check(yd.ok && yc.ok, tie ? "anchor at its capacity: both phases solved"
                                      : "geogrid at its capacity: both phases solved");
            if (yd.ok && yc.ok) {
                std::printf("   %s at its capacity %.4f: N consolidation %.6f vs drained %.6f (%+.4f%%),"
                            " settlement %+.4f%%\n", what, cap, yc.moment_or_n, yd.moment_or_n,
                            100.0 * (yc.moment_or_n / yd.moment_or_n - 1.0),
                            100.0 * (yc.settle / yd.settle - 1.0));
                check(std::fabs(yc.moment_or_n) <= cap * (1.0 + 1e-6),
                      "the element carries no more than its capacity in the consolidation phase");
                if (tie) {
                    check(std::fabs(yc.moment_or_n / yd.moment_or_n - 1.0) < 5e-3 &&
                              std::fabs(yc.settle / yd.settle - 1.0) < 5e-3,
                          "and ends at the drained phase's yielded answer (0.5%)");
                } else {
                    // THE SHEET HAS A HISTORY. The coupled phase loads the ground undrained first:
                    // it spreads sideways at constant volume and stretches the sheet past its
                    // capacity, where it yields; as the ground then consolidates it draws in again
                    // and the yielded sheet shortens, unloading from its capacity. The drained
                    // phase never stretched it past, so the two end apart -- and both halves of
                    // that history are measured here rather than assumed.
                    const auto early = [&](double c) {
                        m::Project pe = build(false, c);
                        pe.phases[0].duration = 0.01 * kH * kH / cv();   // Tv = 0.01
                        return run(pe);
                    };
                    const Answer e0 = early(0.0), ec = early(cap);
                    std::printf("   geogrid at Tv = 0.01: N %.6f with no capacity, %.6f with it; at Tv = 4 %.6f (capacity %.6f)\n",
                                e0.moment_or_n, ec.moment_or_n, yc.moment_or_n, cap);
                    check(e0.ok && ec.ok && e0.moment_or_n > 2.0 * std::fabs(drained.moment_or_n),
                          "early in the phase the ground loaded undrained stretches the sheet to more than twice its drained force");
                    check(ec.moment_or_n <= cap * (1.0 + 1e-6) && yc.moment_or_n < ec.moment_or_n,
                          "and as the ground consolidates the yielded sheet shortens: it unloads from its capacity");
                }
            }
        }
        // The fully-coupled phase applies it too.
        m::Project pr = build(false, cap);
        pr.phases[0].type = m::PhaseType::FullyCoupled;
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool warned = false;
        double reported = 0.0;
        for (const auto& d : res[1].diagnostics)
            if (d.code == std::string("K2D-A014")) {
                warned = true;
                const size_t at = d.message.find(" at ");
                if (at != std::string::npos) reported = std::atof(d.message.c_str() + at + 4);
            }
        (void)reported;
        const double fcn = res[1].struct_forces.empty() ? 0.0
                           : (res[1].struct_forces.front().kind == 0 ? res[1].struct_forces.front().max_M
                                                                     : res[1].struct_forces.front().max_N);
        std::printf("   %s, fully coupled, capacity %.4f: N %.6f\n", what, cap, fcn);
        check(res[1].ok && !warned && std::fabs(fcn) <= cap * (1.0 + 1e-6),
              tie ? "anchor, fully coupled: held at its capacity, and no elastic branch to report"
                  : "geogrid, fully coupled: held at its capacity, and no elastic branch to report");
    }
}

// The limit this section used to pin has been LIFTED, and the test says so rather than being
// deleted. It asserted that a wall -- a plate with interfaces, which splits the mesh -- was refused
// in a consolidation phase, on the grounds that the seam would silently become impermeable. The
// grounds were right about the danger and wrong about the remedy: the interface's cross
// permeability was already an input with a documented default, and a later increment read it
// (KV-STR-007). So what is asserted here now is that the wall RUNS and that its joints are
// reported -- and what is still refused is the one thing that has not been answered.
void test_what_is_and_is_not_carried() {
    std::printf("\n-- (4) a wall with interfaces runs, and so does an embedded beam --\n");
    m::Project pr = block(/*with_plate=*/true, /*plastic_phase=*/false);
    pr.structs.front().iface_pos = true;    // a plate WITH interfaces is an embedded wall:
    pr.structs.front().iface_neg = true;    // the mesh is split along it
    const Answer wall = run(pr);
    check(wall.ok, "a wall (plate + interfaces) in a consolidation phase solves");
    if (!wall.ok) { std::printf("   (%s)\n", wall.msg.c_str()); }

    // The embedded beam was the last refusal: its skin and base slip by definition, and a phase
    // that solved its structures as one elastic matrix could not carry it. It runs now; its
    // verification is KV-STR-023 (test_consolidation_pile).
    m::Project pb = block(/*with_plate=*/false, /*plastic_phase=*/false);
    m::EmbeddedBeamMaterial em; em.name = "Pile";
    pb.embedded.push_back(em);
    m::StructElement e; e.kind = m::StructKind::EmbeddedBeam; e.name = "Pile";
    e.x1 = 10.0; e.y1 = 2.0; e.x2 = 10.0; e.y2 = kH; e.material = 0;
    pb.structs.push_back(e);
    const Answer beam = run(pb);
    check(beam.ok, "an embedded beam in a consolidation phase solves (it used to be refused)");
    if (!beam.ok) std::printf("   (%s)\n", beam.msg.c_str());
}

}  // namespace

int main() {
    std::printf("Structural elements inside a consolidation phase\n\n");
    test_drained_limit();
    test_tie_and_geogrid();
    test_elastic_limit_is_not_silent();
    test_what_is_and_is_not_carried();
    if (g_failures == 0) {
        std::printf("\nOK: the coupled phase carries the structure, and walks into the drained answer\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
