// Measurement instrument for the CONVERGENCE CRITERIA FAMILY (a `study_` target:
// EXCLUDE_FROM_ALL, not wired into CTest).
//
// The question it was written for is N-2's acceptance, stated before the code existed: is there
// a case in this tree that SATISFIES the global force criterion and FAILS a local one? The
// prediction on record was the Hardening Soil oedometer at its 1e-2 default -- the tolerance
// this project inherited from the same source the criteria come from. A prediction that cannot
// fail is not worth making, so the run prints every member of the family for every case and
// leaves the reading to the numbers.
//
// usage: study_convergence_family [case]
//   oedo     KV-CST-002, the Hardening Soil oedometer, over a sweep of tolerances
//   footing  a Mohr-Coulomb strip footing walked towards its limit load (CSP must fall)
//   struct   the corpus interface and embedded-beam cases (Eq. 9-8 and Eq. 9-9)
//   elastic  a linear-elastic block: every local count must be zero and CSP exactly 1
//   all      all of the above (default)
#include <katai/analysis/results.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

// A Mohr-Coulomb strip footing. The point of including it is CSP: a body that plastifies
// through must show the parameter falling, and if it does not, the criterion built on it is
// decoration.
m::Project mc_footing(double q, m::SoilModel model) {
    m::Project pr;
    pr.name = "MC footing";
    pr.x_min = 0.0; pr.x_max = 20.0;
    pr.y_min = 0.0; pr.y_max = 10.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Clay";
    s.model = model;
    s.drainage = m::Drainage::Drained;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.E = 1.0e4; s.nu = 0.3; s.c = 20.0; s.phi = 20.0; s.psi = 0.0;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Ground";
    P.material = 0;
    P.x = {0, 20.0, 20.0, 0};
    P.y = {0, 0, 10.0, 10.0};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Footing";
    L.x1 = 0; L.y1 = 10.0; L.x2 = 2.0; L.y2 = 10.0;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -q;
    pr.loads.push_back(L);

    pr.initial.load_active = {0};
    m::Phase step;
    step.name = "Apply footing load";
    step.load_active = {1};
    pr.phases.push_back(step);
    return pr;
}

void print_header() {
    std::printf("%-26s %7s %10s %10s %5s %5s %12s %7s  %s\n", "case", "CSP", "force err",
                "tolerated", "lf", "iters", "max|u| [m]", "time",
                "local criteria (inaccurate / points, worst error)");
    std::printf("%s\n", std::string(170, '-').c_str());
}

// Run one project and print the whole family for its LAST phase.
void run(const char* label, m::Project pr, double tol) {
    if (tol > 0.0) {
        pr.initial.tolerance = tol;
        for (auto& ph : pr.phases) ph.tolerance = tol;
    }
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("%-26s  MESH FAILED: %s\n", label, M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty()) { std::printf("%-26s  no results\n", label); return; }
    const katai::core::SolveResult& r = res.back();
    const auto& c = r.convergence;
    if (!c.measured) {
        std::printf("%-26s  NOT MEASURED (%s)\n", label,
                    r.ok ? "solved without an accepted iterate?" : r.message.c_str());
        return;
    }
    // The answer itself, so that the cost of a stricter stopping rule can be read next to what
    // it buys: an extra iteration that does not move the number is a cost with no benefit, and
    // one that moves it by a per cent is the reason the rule exists.
    double umax = 0.0;
    for (int i = 0; i < r.disp.size(); ++i) umax = std::max(umax, std::fabs(r.disp[i]));
    std::printf("%-26s %7.5f %10.3e %10.3e %5.3f %5d %12.9f %6.1fs  plastic %4d/%-5d "
                "worst %8.2e %-4s | nl-el %4d/%-5d (of %5d el) worst %8.2e %-4s%s\n",
                label, c.csp, c.force_error, c.tolerated, r.load_factor, r.iterations, umax,
                r.timings.total,
                c.plastic_inaccurate, c.plastic_points, c.worst_plastic_error,
                c.plastic_points_ok() ? "ok" : "FAIL",
                c.nl_elastic_inaccurate, c.nl_elastic_points, c.elastic_points,
                c.worst_nl_elastic_error, c.nl_elastic_ok() ? "ok" : "FAIL",
                r.ok ? "" : "  [phase did not converge]");
    // The structural members of the family print only where the model has them, so a case
    // without an interface or a pile does not carry a row of zeros pretending to be a check.
    if (c.iface_points > 0)
        std::printf("      interfaces + skin springs: %d/%d inaccurate, worst %8.2e  %s\n",
                    c.iface_inaccurate, c.iface_points, c.worst_iface_error,
                    c.iface_points_ok() ? "ok" : "FAIL");
    if (c.feet > 0)
        std::printf("      embedded-beam feet (%d): force error %8.2e of %8.2e tolerated  %s\n",
                    c.feet, c.foot_force_error,
                    c.kFootToleranceFactor * c.tolerated, c.foot_ok() ? "ok" : "FAIL");
    if (c.has_moment)
        std::printf("      moment residual: %8.2e of %8.2e tolerated  %s\n", c.moment_error,
                    c.tolerated, c.moment_ok() ? "ok" : "FAIL");
    if (c.force_ok() && !c.local_ok())
        std::printf("    ^^ GLOBAL MET, LOCAL FAILED -- the case N-2 predicted exists.\n");
}

}  // namespace

int main(int argc, char** argv) {
    const std::string which = argc > 1 ? argv[1] : "all";
    print_header();

    if (which == "oedo" || which == "all") {
        // The corpus file itself, so the case is the one the record publishes numbers for.
        m::Project pr;
        std::string err;
        const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
        if (!m::load_project(path, pr, &err, nullptr)) {
            std::printf("FAIL load %s: %s\n", path.c_str(), err.c_str());
            return 1;
        }
        // The sweep is the point: the tolerance is what the LOCAL criteria are measured against
        // too, so a case that fails locally at 1e-2 and passes at 1e-6 says exactly where the
        // default sits relative to the accuracy it implies.
        const double tols[] = {1e-2, 1e-3, 1e-4, 1e-6};
        for (double t : tols) {
            char lbl[64];
            std::snprintf(lbl, sizeof lbl, "HS oedometer tol=%.0e", t);
            run(lbl, pr, t);
        }
    }

    if (which == "footing" || which == "all") {
        const double qs[] = {100.0, 300.0, 600.0, 900.0};
        for (double q : qs) {
            char lbl[64];
            std::snprintf(lbl, sizeof lbl, "MC footing q=%.0f kPa", q);
            run(lbl, mc_footing(q, m::SoilModel::MohrCoulomb), 0.0);
        }
        // The decisive experiment. On this case the GLOBAL error falls by five orders of
        // magnitude between the first iterate and the last while the LOCAL one falls by a
        // factor of three. Two readings explain that, and they call for opposite actions: the
        // local error is simply slower to converge (tighten the tolerance and it follows), or
        // it has stalled on something iteration cannot fix (tighten the tolerance and it does
        // not move). Only a sweep can say which.
        std::printf("\n-- q=300 kPa, tolerance swept: does the LOCAL error follow the global "
                    "one down? --\n");
        const double tols[] = {1e-2, 1e-4, 1e-6, 1e-8, 1e-10};
        for (double t : tols) {
            char lbl[64];
            std::snprintf(lbl, sizeof lbl, "MC footing tol=%.0e", t);
            run(lbl, mc_footing(300.0, m::SoilModel::MohrCoulomb), t);
        }
    }

    // The two corpus cases that exercise the STRUCTURAL members of the family: an interface
    // (Eq. 9-8) and an embedded beam with a foot (Eq. 9-9). Both are cases the record already
    // publishes numbers for, so what the new criteria say about them can be read against a
    // known answer rather than against nothing.
    if (which == "struct" || which == "all") {
        struct Case { const char* file; const char* label; };
        const Case cases[] = {
            {"kv-str-002-plaxis-sliding-block.k2d", "sliding block (interface)"},
            {"kv-str-004-axial-pile-capacity.k2d", "axial pile (foot + skin)"},
        };
        for (const Case& cs : cases) {
            m::Project pr;
            std::string err;
            if (!m::load_project(std::string(KATAI_CORPUS_DIR) + "/" + cs.file, pr, &err,
                                 nullptr)) {
                std::printf("FAIL load %s: %s\n", cs.file, err.c_str());
                continue;
            }
            run(cs.label, pr, 0.0);
        }
    }

    // The RATE-DEPENDENT case, and it is here because it is the one that refuted the argument
    // for binding these criteria. On the Hardening Soil oedometer, requiring them at the shipped
    // tolerance lands on the converged answer; on this Soft Soil Creep column every duration
    // moved AWAY from the published creep law when they were required. Two readings are
    // possible and they call for opposite actions -- the criteria are wrong for a viscous point,
    // or the looser run was accidentally closer to a law the model does not exactly obey -- and
    // only a tolerance sweep with the criteria OFF can say which. If the global-only answer
    // walks to the same place as the bound one, the bound one is right and the case's band is
    // what has to move.
    if (which == "creep" || which == "all") {
        m::Project pr;
        std::string err;
        if (m::load_project(std::string(KATAI_CORPUS_DIR) + "/kv-cst-010-soft-soil-creep-column.k2d",
                            pr, &err, nullptr)) {
            std::printf("\n-- Soft Soil Creep column: where does the answer go as the global "
                        "tolerance tightens? --\n");
            // mu* ln(1 + t/tau), tau = 1 day, over a 4 m column: the law the case is verified
            // against, evaluated here so the walk can be read against it rather than against
            // the previous run.
            const double mu = 0.008, H = 4.0;
            for (double days : {1.0, 10.0, 100.0}) {
                m::Project at = pr;
                for (auto& ph : at.phases)
                    if (ph.duration > 0.0) ph.duration = days;
                std::printf("   t = %6.0f d   closed form mu* ln(1 + t/tau) = %.9f m\n", days,
                            mu * std::log(1.0 + days / 1.0) * H / 4.0);
                for (double t : {1e-2, 1e-4, 1e-6, 1e-8}) {
                    char lbl[64];
                    std::snprintf(lbl, sizeof lbl, "  SSC %.0f d, tol=%.0e", days, t);
                    run(lbl, at, t);
                }
            }
        } else {
            std::printf("FAIL load creep column: %s\n", err.c_str());
        }
    }

    // The control: linear elasticity. Nothing is plastic, nothing has a stress-dependent
    // modulus, so every local count must be zero and CSP must be exactly 1. A criterion that
    // cannot report "nothing to report" will report noise everywhere else.
    if (which == "elastic" || which == "all")
        run("LE block (control)", mc_footing(50.0, m::SoilModel::LinearElastic), 0.0);
    return 0;
}
