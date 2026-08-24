// The convergence criteria family (Scientific Manual §9.1): what this tree measures at the
// accepted iterate, and what it would cost to be bound by it.
//
// Until this landed the solver checked ONE thing -- a global force residual against a fixed
// scale -- and said "converged". The tests below pin the four statements that make the family
// worth having, and each of them is written so that it can fail: a control that must report
// nothing, a parameter that must fall, a separation that must be large, and a cost that must be
// paid in iterations rather than in the answer.
#include <katai/analysis/results.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::core::SolveResult;

namespace {

int failures = 0;
void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++failures;
}

// A Mohr-Coulomb strip footing on a weightful clay. Loading it harder walks the body from
// almost elastic to a full mechanism, which is the only way to watch CSP do what it is for.
m::Project footing(double q, m::SoilModel model) {
    m::Project pr;
    pr.name = "footing";
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

struct Run {
    bool ok = false;
    katai::core::NewtonResult::Convergence c;
    double umax = 0.0;
    int iterations = 0;
    double load_factor = 0.0;
};

// enforce_local: the tri-state of NumericalControls -- +1 forces the local criteria on, -1
// forces them OFF, which is what the comparisons below need now that ON is the default.
Run solve(m::Project pr, double tol, bool enforce_local) {
    if (tol > 0.0) {
        pr.initial.tolerance = tol;
        for (auto& ph : pr.phases) ph.tolerance = tol;
    }
    const auto M = katai::app::mesh_from_project(pr);
    Run out;
    if (!M.ok) return out;
    katai::app::NumericalControls nc;
    nc.enforce_local_criteria = enforce_local ? 1 : -1;
    const auto res = katai::app::solve_phases(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, nullptr, nc);
    if (res.empty()) return out;
    const SolveResult& r = res.back();
    out.ok = r.ok;
    out.c = r.convergence;
    out.iterations = r.iterations;
    out.load_factor = r.load_factor;
    for (int i = 0; i < r.disp.size(); ++i) out.umax = std::max(out.umax, std::fabs(r.disp[i]));
    return out;
}

// ---------------------------------------------------------------------------------------
// verify: KV-NUM-013
//   oracle:   independent_path
//   source:   PLAXIS 2D Scientific Manual 2025.1 sec 9.1 -- the deformation convergence criteria as a FAMILY: Eq. 9-1 the CSP-normalised global force error, Eq. 9-3/9-4 the moment residual, Eq. 9-5 the local error at a plastic soil stress point and Eq. 9-7 at a non-linear elastic one, with the equilibrium/constitutive stress pair of Fig. 9-1 and Eq. 9-6; and PLAXIS 2D Reference Manual 2025.1 sec 7.13.8, Eq. 7-21/7-22 for the same global check and for CSP itself. THE TWO MANUALS PRINT CSP AS RECIPROCALS OF EACH OTHER: Scientific Eq. 9-2 gives elastic energy over total, Reference Eq. 7-22 gives total over elastic. Only the Reference orientation matches the behaviour BOTH manuals state in words -- "when the solution is fully elastic the Stiffness is equal to unity, whereas at failure the Stiffness approaches zero" -- and the uses built on it there (arc-length engages below CSP 0.5, collapse is reported below 0.015). The inverted form is >= 1 and GROWS without bound as a mechanism forms, which would make Eq. 9-1 loosen towards failure instead of tightening. Implemented from Eq. 7-22, deliberately, and this case is what holds it there
//   locator:  three fixtures on one strip-footing geometry: a LINEAR ELASTIC control that must report nothing at all, the same geometry in Mohr-Coulomb walked from a nearly elastic load to a full mechanism, and the checked-in tests/corpus/kv-cst-002-hs-oedometer.k2d solved three ways -- at the shipped tolerance, at the shipped tolerance with the local criteria binding, and at a tolerance four decades tighter with only the global one
//   quantity: CSP, the global force error, the local error at each soil stress point and the counts of inaccurate points, and the settlement each stopping rule stops at [-; -; m]
//   expected: the elastic control reports CSP equal to 1 with no plastic and no non-linear elastic points and every criterion met; CSP falls monotonically as the footing load rises and passes below the 0.015 the source reports collapse at; a case exists that meets the global criterion by three or more orders of magnitude while failing a local one; and binding the local criteria at the shipped tolerance lands on the settlement that four decades of extra global tolerance produce, for fewer iterations than that tolerance costs
//   band:     CSP within 1e-12 of 1 on the elastic control (the two energy sums reach the same number by different arithmetic) and strictly below 1 wherever a plastic point exists; the local/global separation asserted at 1000x, measured 63000x; the settlements agreed to 0.05%, measured 0.0048%, against a global-only run at the same tolerance that is asserted to differ by more than 0.05% and measures 0.179%
// The separation is the point. On the footing at q = 300 kPa the global force error falls from
// 1.3e-2 to 3.7e-7 between the first iterate and the last -- five orders of magnitude -- while
// the worst local error falls from 6.5e-2 to 2.3e-2, a factor of three. Both vanish at the exact
// solution (the elastic control measures 1.4e-15 and 0.0), so the local error is an error and
// not a measure of nonlinearity; it simply lags the force residual by three to four decades,
// and a solver that stops on the force residual alone stops a step early in local terms.
void test_convergence_family() {
    // (1) The control. Nothing is plastic and no modulus depends on stress, so every local
    // count must be zero and CSP must be one. A criterion that cannot report "nothing to
    // report" will report noise everywhere else.
    const Run le = solve(footing(50.0, m::SoilModel::LinearElastic), 0.0, false);
    check(le.ok, "elastic control solves");
    check(le.c.measured, "elastic control measured the family");
    // Not an exact equality, and the reason is worth stating: CSP's numerator is the work
    // against the stress the material law RETURNED (built through the Lame form inside the
    // elastic predictor) and its denominator is the work against D^e times the same strain
    // (built as a matrix product). On a linear-elastic point those are the same number reached
    // two ways, so agreement to round-off is a small consistency check rather than a tautology.
    std::printf("     elastic control: CSP = %.17g (1 - CSP = %.3e)\n", le.c.csp,
                1.0 - le.c.csp);
    check(std::fabs(le.c.csp - 1.0) < 1e-12,
          "elastic control: CSP is 1 to round-off (the two energy sums agree)");
    check(le.c.plastic_points == 0, "elastic control: no plastic points");
    check(le.c.nl_elastic_points == 0, "elastic control: no stress-dependent elastic points");
    check(le.c.worst_plastic_error == 0.0 && le.c.worst_nl_elastic_error == 0.0,
          "elastic control: no local error of either kind");
    check(le.c.all_ok(), "elastic control: every criterion satisfied");

    // (2) CSP must FALL as the body plastifies, and keep falling to the collapse value the
    // source names. This is the guard on the manual conflict: implemented from the Scientific
    // Manual's printed Eq. 9-2 instead, every one of these would be >= 1 and the sequence would
    // RISE -- so this check does not merely confirm a number, it convicts the other reading.
    const Run f100 = solve(footing(100.0, m::SoilModel::MohrCoulomb), 0.0, false);
    const Run f300 = solve(footing(300.0, m::SoilModel::MohrCoulomb), 0.0, false);
    const Run f900 = solve(footing(900.0, m::SoilModel::MohrCoulomb), 0.0, false);
    std::printf("     CSP: q=100 %.5f  q=300 %.5f  q=900 %.5f (lf %.3f)\n",
                f100.c.csp, f300.c.csp, f900.c.csp, f900.load_factor);
    check(f100.c.csp < 1.0 && f300.c.csp < f100.c.csp && f900.c.csp < f300.c.csp,
          "CSP falls monotonically as the footing load rises (never rises: Eq. 7-22, not 9-2)");
    check(f900.c.csp < 0.015,
          "CSP below 0.015 at a load the model cannot carry (the source's collapse value)");
    check(f900.load_factor < 1.0, "and that load is indeed not carried");

    // (3) The separation. The global criterion is met by orders of magnitude while a local one
    // is not -- the case N-2 was written to find, and the reason one check is not enough.
    check(f300.ok && f300.c.force_ok(), "q=300 converged on the global force criterion");
    const double separation = f300.c.worst_plastic_error / f300.c.force_error;
    std::printf("     q=300: force error %.3e, worst local %.3e, separation %.0fx\n",
                f300.c.force_error, f300.c.worst_plastic_error, separation);
    check(separation > 1000.0,
          "the worst local error exceeds the global force error by more than 1000x");
    check(!f300.c.local_ok(),
          "and the local criterion is FAILED on a run the global one calls converged");
}

// The cost of being bound by the family, measured on the case the record publishes numbers for.
// Two readings were possible before it was run and they call for opposite decisions: the local
// criteria are over-strict (they would buy nothing and cost iterations), or the global-only stop
// is premature (they would buy the converged answer). The measurement decided it.
void test_local_criteria_cost() {
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        check(false, "load kv-cst-002-hs-oedometer.k2d: " + err);
        return;
    }
    // 1e-2 is what this tree ships for the Hardening Soil family, from the same source the
    // criteria come from. 1e-6 is four decades tighter and is what the record uses when it wants
    // the answer rather than a run.
    const Run loose = solve(pr, 1e-2, false);
    const Run bound = solve(pr, 1e-2, true);
    const Run tight = solve(pr, 1e-6, false);
    check(loose.ok && bound.ok && tight.ok, "all three oedometer runs converged");

    const double d_bound = std::fabs(bound.umax - tight.umax) / tight.umax;
    const double d_loose = std::fabs(loose.umax - tight.umax) / tight.umax;
    std::printf("     settlement: loose %.9f (%.3f%%)  bound %.9f (%.4f%%)  tight %.9f\n"
                "     iterations: loose %d  bound %d  tight %d\n",
                loose.umax, 100.0 * d_loose, bound.umax, 100.0 * d_bound, tight.umax,
                loose.iterations, bound.iterations, tight.iterations);

    // Both directions, so neither half can pass vacuously: the global-only stop at the shipped
    // tolerance MISSES the converged answer, and binding the local criteria at the same
    // tolerance FINDS it.
    check(d_loose > 5e-4,
          "the shipped tolerance, global criterion only, stops more than 0.05% short");
    check(d_bound < 5e-4,
          "binding the local criteria at the SAME tolerance lands within 0.05% of the "
          "four-decades-tighter answer");
    check(bound.iterations < tight.iterations,
          "and it gets there in fewer iterations than tightening the global tolerance costs");

    // What it costs where it buys nothing: the criteria must not be free to ignore, but they
    // must also not be paid for twice. The loose run is the one that was wrong, so the extra
    // iterations are the price of the correction, not overhead.
    check(bound.iterations > loose.iterations,
          "the correction is paid for in iterations (it is not free)");
}

// The STRUCTURAL half of the family (Eq. 9-8, Eq. 9-9), on the two corpus cases that have the
// elements to exercise it. Both already publish a number in the verification matrix, so what the
// new criteria say about them can be read against a known answer.
void test_structural_criteria() {
    m::Project block, pile;
    std::string err;
    const std::string dir = std::string(KATAI_CORPUS_DIR) + "/";
    if (!m::load_project(dir + "kv-str-002-plaxis-sliding-block.k2d", block, &err, nullptr) ||
        !m::load_project(dir + "kv-str-004-axial-pile-capacity.k2d", pile, &err, nullptr)) {
        check(false, "load the interface and pile corpus files: " + err);
        return;
    }

    // The sliding block is the interface case: an elastic block pushed until the joint under it
    // slips. Nothing in the soil yields, so the ONLY local criterion with anything to say is the
    // interface one -- which makes it the clean test that Eq. 9-8 is wired at all.
    const Run b = solve(block, 0.0, false);
    check(b.ok, "sliding block solves");
    std::printf("     sliding block: %d/%d interface points inaccurate, worst %.2e (force %.3e)\n",
                b.c.iface_inaccurate, b.c.iface_points, b.c.worst_iface_error, b.c.force_error);
    check(b.c.iface_points > 0, "the interface's slipping points are counted");
    check(b.c.plastic_points == 0,
          "...and the soil contributes none, so this is the interface criterion alone");
    check(b.c.force_ok() && !b.c.iface_points_ok(),
          "the global force criterion is met while the interface criterion is not");
    const Run b2 = solve(block, 0.0, true);
    check(b2.ok && b2.c.iface_points_ok(),
          "binding the local criteria settles the interface points");
    check(b2.iterations > b.iterations, "...and that costs iterations");

    // The pile is the foot case. Eq. 9-9 is a ratio over every foot in the model rather than a
    // count, and it is tolerated at five times the tolerated error -- the source's factor, which
    // this checks is actually applied rather than quietly rounded to the same bar as the rest.
    const Run p = solve(pile, 0.0, false);
    check(p.ok, "axial pile solves");
    std::printf("     axial pile: foot error %.2e of %.2e tolerated; skin+interface %d/%d\n",
                p.c.foot_force_error, p.c.kFootToleranceFactor * p.c.tolerated,
                p.c.iface_inaccurate, p.c.iface_points);
    check(p.c.feet == 1, "the pile's foot is seen");
    check(p.c.iface_points > 0,
          "the skin coupling springs are counted with the interface points, as the source does");
    check(p.c.foot_ok(), "the foot force is in balance on a run that reached full load");
    // Two-sided: the allowance must be the source's five times, not one.
    katai::core::NewtonResult::Convergence probe = p.c;
    probe.foot_force_error = 3.0 * probe.tolerated;
    check(probe.foot_ok(), "a foot error at 3x the tolerated error is allowed (the 5x factor)");
    probe.foot_force_error = 7.0 * probe.tolerated;
    check(!probe.foot_ok(), "...and one at 7x is not");
}

}  // namespace

int main() {
    test_convergence_family();
    test_structural_criteria();
    test_local_criteria_cost();
    std::printf("%s\n", failures == 0 ? "ALL PASS" : "FAILURES");
    return failures == 0 ? 0 : 1;
}
