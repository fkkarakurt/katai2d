// The excess pore pressure as a state of the phase chain.
//
// An undrained phase generates an excess pore pressure; the phases after it must find that pressure
// where it was left, whatever they are. Until 2026-09 they did not. The pressure was re-derived in
// every phase as (Kw/n) times the volume change SINCE THE INITIAL STATE, whatever had produced it:
//
//   - a phase that ignores undrained behaviour loaded the soil drained, and the next undrained phase
//     read that drained compression as a pore pressure -- a phase that changed nothing heaved 28.7 mm;
//   - after a gravity-loading initial phase that ignored undrained behaviour, a downward undrained
//     load heaved the column 41 mm and left 6% of its effective stress;
//   - the ignoring phase itself deleted the pressure generated before it, and the soil consolidated
//     at once in a phase that has no time;
//   - after a consolidation the next undrained phase undid it (linear skeleton: all of it), because
//     the consolidation's own settlement was in the volume change and its computed pressure field
//     was dropped;
//   - and the same flag took the undrained partial factor away from an Undrained (B) soil in a
//     design phase that ignores undrained behaviour: a footing on it carried 11% more.
//
// Each has an exact answer on a laterally confined column, because a confined column is a 1-D
// problem the elements reproduce exactly: M' = E (1 - nu) / ((1 + nu)(1 - 2 nu)) drained,
// M_u = M' + Kw/n undrained, Kw/n = 3 (nu_u - nu) / ((1 - 2 nu_u)(1 + nu)) K'. A phase that changes
// nothing moves nothing; a load q taken undrained settles q H / M_u and moves the effective stress by
// -q M' / M_u; one taken while undrained behaviour is ignored settles q H / M' and leaves the pressure
// already in the ground as it was.
//
// verify: KV-CST-017
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley, New York -- the effective stress principle and the one-dimensional confined column; Biot, M.A. (1941). General theory of three-dimensional consolidation. J. Appl. Phys. 12, 155-164, for the undrained limit of the coupled column (the pore fluid adds its bulk stiffness Kw/n to the skeleton); the undrained relations of isotropic linear elasticity as recorded in docs/references/effective-stress-formulation.md; KATAI 2D input contract (docs/k2d-format.md, `ignoreund`: no water stiffness and no new excess pore pressure, the pressure generated earlier kept)
//   locator:  a 1 x 8 m column, sides on rollers, base fixed, water table at the surface, E = 5000 kN/m2, nu = 0.3, nu_u = 0.495 (so Kw/n = 45 K' = 187500 kN/m2, M' = 6730.77 kN/m2, M_u = 194230.77 kN/m2), surcharges q1 = 25 and q2 = 15 kN/m2 on the top. Stated in full, per chain after a K0 initial phase unless stated: (A) q1 ignoring -> nil undrained: 0; (B) q1 ignoring -> q2 undrained: q2 H / M_u; (C) q1 undrained -> nil ignoring: 0; (D) q1 undrained -> q2 ignoring: q2 H / M'; (E) gravity-loading initial phase ignoring -> q1 undrained: q1 H / M_u and d(sigma'_yy) = -q1 M' / M_u; (F) q1 undrained -> consolidation to dissipation -> nil undrained: 0, then q2 undrained: q2 H / M_u, for a linear-elastic and a Mohr-Coulomb (elastic range) skeleton; (G) q1 undrained -> consolidation for 5 days -> nil undrained: 0, then consolidation to dissipation: the drained total q1 H / M' over the four phases; (H) q1 undrained -> fully coupled to dissipation -> nil undrained: 0; (I) the checked-in KV-FND-010 footing on an Undrained (B) Tresca clay in an EC7 DA3 phase: the collapse load factor with undrained behaviour ignored equals the one without (phi = 0, so the water cannot change the capacity; only the partial factor could)
//   quantity: settlement of the column top in each phase [m] (phase increments); vertical effective stress at mid-depth [kN/m2]; collapse load factor [-]
//   expected: undrained q1 1.030e-3 m, q2 6.178e-4 m; drained q1 2.9714e-2 m, q2 1.7829e-2 m; a nil phase 0; (I) equal factors
//   band:     nil phases 1e-9 m absolute (the defects this record exists for measured 2.87e-2 m in (A), (C) and the linear-elastic (F), and 9.9e-4 m in the Mohr-Coulomb (F); (E) read 4.11e-2 m of heave for 1.03e-3 m of settlement; all now 0 exactly, the phase starting in equilibrium); closed-form increments 1e-6 relative (measured 3e-14 at worst); the drained total after consolidation 1e-3 relative (measured -7.3e-5: the phase ends with 3e-3 kN/m2 of 25 undissipated); (I) the two factors within one step of the collapse search, 1/320 (measured equal, 0.471875; before, 0.52500 against 0.47187). All as asserted below

#include <katai/io/project_io.hpp>
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
    std::fflush(stdout);
}

constexpr double kE = 5000.0, kNu = 0.3, kNuU = 0.495, kH = 8.0, kQ1 = 25.0, kQ2 = 15.0;

double m_drained() { return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }
double kw_over_n() {
    const double k_eff = kE / (3.0 * (1.0 - 2.0 * kNu));
    return 3.0 * (kNuU - kNu) / ((1.0 - 2.0 * kNuU) * (1.0 + kNu)) * k_eff;
}
double m_undrained() { return m_drained() + kw_over_n(); }

// The column, its two surcharges switched off; `model` 0 = linear elastic, 1 = Mohr-Coulomb.
m::Project column(int model, m::Drainage drainage = m::Drainage::Undrained) {
    m::Project pr;
    pr.name = "KV-CST-017 excess pore chain";
    pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = kH;
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, 1.0};
    pr.wy = {kH, kH};
    pr.initial_procedure = m::InitialProcedure::K0Procedure;

    m::Material clay;
    clay.name = "Clay";
    clay.model = model == 0 ? m::SoilModel::LinearElastic : m::SoilModel::MohrCoulomb;
    clay.drainage = drainage;
    clay.E = kE; clay.nu = kNu; clay.nu_u = kNuU;
    clay.c = 20.0; clay.phi = 25.0; clay.psi = 0.0;
    clay.gamma_unsat = 17.0; clay.gamma_sat = 19.0;
    clay.kx = 1e-4; clay.ky = 1e-4;
    pr.materials.push_back(clay);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, 1, 1, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    for (double q : {kQ1, kQ2}) {
        m::Load L;
        L.kind = m::LoadKind::Distributed;
        L.name = q == kQ1 ? "q1" : "q2";
        L.x1 = 0.0; L.y1 = kH; L.x2 = 1.0; L.y2 = kH;
        L.qy1 = L.qy2 = -q;
        pr.loads.push_back(L);
    }
    pr.initial.load_active = {0, 0};
    return pr;
}

enum class Step { Q1, Q1Ignoring, Q2, Q2Ignoring, Nil, NilIgnoring, Consolidate, ConsolidatePart,
                  FullyCoupled };

void add(m::Project& pr, Step s) {
    m::Phase ph;
    ph.type = m::PhaseType::Plastic;
    // Loads carry over: the last phase's activity, then this step's change.
    ph.load_active = pr.phases.empty() ? std::vector<char>{0, 0} : pr.phases.back().load_active;
    switch (s) {
        case Step::Q1: ph.name = "q1"; ph.load_active[0] = 1; break;
        case Step::Q1Ignoring: ph.name = "q1 ignoring"; ph.load_active[0] = 1; ph.ignore_undrained = true; break;
        case Step::Q2: ph.name = "q2"; ph.load_active[1] = 1; break;
        case Step::Q2Ignoring: ph.name = "q2 ignoring"; ph.load_active[1] = 1; ph.ignore_undrained = true; break;
        case Step::Nil: ph.name = "nil"; break;
        case Step::NilIgnoring: ph.name = "nil ignoring"; ph.ignore_undrained = true; break;
        case Step::Consolidate:
            ph.name = "consolidate"; ph.type = m::PhaseType::Consolidation;
            ph.duration = 5000.0; ph.time_steps = 200; break;
        case Step::ConsolidatePart:
            ph.name = "consolidate 5 days"; ph.type = m::PhaseType::Consolidation;
            ph.duration = 5.0; ph.time_steps = 50; break;
        case Step::FullyCoupled:
            ph.name = "fully coupled"; ph.type = m::PhaseType::FullyCoupled;
            ph.duration = 5000.0; ph.time_steps = 200; break;
    }
    pr.phases.push_back(ph);
}

std::vector<katai::app::SolveResult> solve(const m::Project& pr) {
    std::vector<katai::app::SolveResult> none;
    const auto rep = katai::io::validate_project(pr);
    if (!rep.ok()) {
        for (const auto& i : rep.issues) std::printf("      invalid: %s %s\n", i.path.c_str(), i.message.c_str());
        return none;
    }
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return none;
    return katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure),
                                    nullptr, nullptr, {});
}

bool all_ok(const std::vector<katai::app::SolveResult>& R, size_t n, const std::string& what) {
    bool ok = R.size() == n;
    for (const auto& r : R) ok = ok && r.ok;
    check(ok, what + " solves");
    if (!ok)
        for (const auto& r : R)
            if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    return ok;
}

// Settlement of the column top in the phase (positive down).
double settle(const katai::app::SolveResult& R) {
    int best = 0;
    double bd = 1e300;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        const double d = std::hypot(R.mesh.x[n] - 0.5, R.mesh.y[n] - kH);
        if (d < bd) { bd = d; best = n; }
    }
    return -R.disp[best * 2 + 1];
}

// Vertical effective stress at mid-depth, averaged over the nodes on that line.
double syy_mid(const katai::app::SolveResult& R) {
    double s = 0.0;
    int k = 0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.y[n] - 0.5 * kH) < 1e-9) { s += R.stress.stress[n](1); ++k; }
    return k ? s / k : 0.0;
}

void expect_nil(const katai::app::SolveResult& R, const std::string& what) {
    const double u = settle(R);
    std::printf("  %-58s %+.3e m\n", what.c_str(), u);
    check(std::fabs(u) < 1e-9, what + ": a phase that changes nothing moves nothing");
}

void expect_close(double got, double want, double rel, const std::string& what) {
    const double e = (got - want) / want;
    std::printf("  %-58s %.9e (closed form %.9e, %+.2e)\n", what.c_str(), got, want, e);
    check(std::fabs(e) < rel, what);
}

}  // namespace

int main() {
    const double Mp = m_drained(), Mu = m_undrained();
    std::printf("M' = %.4f  Kw/n = %.4f  M_u = %.4f kN/m2\n", Mp, kw_over_n(), Mu);
    const double und1 = kQ1 * kH / Mu, und2 = kQ2 * kH / Mu, dr1 = kQ1 * kH / Mp, dr2 = kQ2 * kH / Mp;

    std::printf("\n== (A) a drained stage leaves no pore pressure behind ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1Ignoring); add(pr, Step::Nil);
        const auto R = solve(pr);
        if (all_ok(R, 3, "(A)")) {
            expect_close(settle(R[1]), dr1, 1e-6, "(A) q1 ignoring undrained behaviour settles q1 H / M'");
            expect_nil(R[2], "(A) the undrained phase after it");
            const double ds = syy_mid(R[2]) - syy_mid(R[1]);
            std::printf("  (A) change of sigma'_yy at mid-depth %+.3e kN/m2\n", ds);
            check(std::fabs(ds) < 1e-6, "(A) and the effective stress stays where the drained stage put it");
        }
    }

    std::printf("\n== (B) ...and the next undrained load is answered from there ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1Ignoring); add(pr, Step::Q2);
        const auto R = solve(pr);
        if (all_ok(R, 3, "(B)"))
            expect_close(settle(R[2]), und2, 1e-6, "(B) q2 undrained settles q2 H / M_u");
    }

    std::printf("\n== (C) ignoring undrained behaviour keeps the pressure already there ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1); add(pr, Step::NilIgnoring);
        const auto R = solve(pr);
        if (all_ok(R, 3, "(C)")) {
            expect_close(settle(R[1]), und1, 1e-6, "(C) q1 undrained settles q1 H / M_u");
            expect_nil(R[2], "(C) the ignoring phase after it");
        }
    }

    std::printf("\n== (D) ...and loads drained on top of it ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1); add(pr, Step::Q2Ignoring);
        const auto R = solve(pr);
        if (all_ok(R, 3, "(D)"))
            expect_close(settle(R[2]), dr2, 1e-6, "(D) q2 while ignoring settles q2 H / M'");
    }

    std::printf("\n== (E) gravity loading that ignored undrained behaviour ==\n");
    {
        m::Project pr = column(1);
        pr.initial_procedure = m::InitialProcedure::GravityLoading;
        pr.initial.ignore_undrained = true;
        add(pr, Step::Q1);
        const auto R = solve(pr);
        if (all_ok(R, 2, "(E)")) {
            expect_close(settle(R[1]), und1, 1e-6, "(E) q1 undrained settles q1 H / M_u");
            expect_close(syy_mid(R[1]) - syy_mid(R[0]), -kQ1 * Mp / Mu, 1e-6,
                         "(E) and moves sigma'_yy by -q1 M' / M_u");
        }
    }

    for (int model : {0, 1}) {
        const std::string tag = model == 0 ? "linear elastic" : "Mohr-Coulomb";
        std::printf("\n== (F) after a consolidation, %s skeleton ==\n", tag.c_str());
        m::Project pr = column(model);
        add(pr, Step::Q1); add(pr, Step::Consolidate); add(pr, Step::Nil); add(pr, Step::Q2);
        const auto R = solve(pr);
        if (all_ok(R, 5, "(F) " + tag)) {
            expect_close(settle(R[1]) + settle(R[2]), dr1, 1e-3,
                         "(F) " + tag + ": undrained + consolidation reach the drained q1 H / M'");
            expect_nil(R[3], "(F) " + tag + ": the undrained phase after it");
            expect_close(settle(R[4]), und2, 1e-6, "(F) " + tag + ": then q2 undrained settles q2 H / M_u");
        }
    }

    std::printf("\n== (G) a consolidation stopped part-way ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1); add(pr, Step::ConsolidatePart); add(pr, Step::Nil); add(pr, Step::Consolidate);
        const auto R = solve(pr);
        if (all_ok(R, 5, "(G)")) {
            expect_nil(R[3], "(G) the undrained phase between the two consolidations");
            expect_close(settle(R[1]) + settle(R[2]) + settle(R[3]) + settle(R[4]), dr1, 1e-3,
                         "(G) and the remainder still dissipates to the drained q1 H / M'");
        }
    }

    std::printf("\n== (H) after a fully coupled phase ==\n");
    {
        m::Project pr = column(1);
        add(pr, Step::Q1); add(pr, Step::FullyCoupled); add(pr, Step::Nil);
        const auto R = solve(pr);
        if (all_ok(R, 4, "(H)")) expect_nil(R[3], "(H) the undrained phase after it");
    }

    std::printf("\n== (I) an Undrained (B) soil keeps gamma_cu when undrained behaviour is ignored ==\n");
    {
        m::Project base;
        std::string err;
        const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-fnd-010-prandtl-strip-footing.k2d";
        check(m::load_project(path, base, &err, nullptr), "the KV-FND-010 footing loads");
        base.materials[0].drainage = m::Drainage::UndrainedB;
        base.phases[0].design_approach = m::DesignApproach::EC7_DA3;
        m::Project ign = base;
        ign.phases[0].ignore_undrained = true;
        const auto Ra = solve(base);
        const auto Rb = solve(ign);
        if (Ra.size() == 2 && Rb.size() == 2) {
            const double a = Ra[1].load_factor, b = Rb[1].load_factor;
            std::printf("  collapse load factor: undrained %.6f, ignoring %.6f (ratio %.5f; the "
                        "partial factors differ by %.5f)\n", a, b, b / a, 1.4 / 1.25);
            check(std::fabs(b - a) <= 1.0 / 320.0 + 1e-12,
                  "(I) the capacity does not depend on the water, so neither may the factor");
        } else {
            check(false, "(I) both footing runs return their phases");
        }
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
