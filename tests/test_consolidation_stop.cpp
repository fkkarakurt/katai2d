// A consolidation phase that ends when the ground has consolidated, not when a duration guessed in
// advance runs out (the `cstop` key of a Consolidation phase in docs/k2d-format.md: minimum excess
// pore pressure / degree of consolidation). The design question has the shape "how
// long until it has settled out?", and until this the only thing that could be said to the phase
// was a number of days -- so the answer had to be obtained by guessing a span, reading the curve,
// and guessing again.
//
// Three things are verified here. The second is the reason this test exists at all; the third the
// test found on its own, and it was a crash.
//
//  (1) THE MACHINERY. The march changes its time step as it goes, and it does that by RESTARTING
//      the fixed-dt core, which is the verified one (KV-CON-001/002). That is only legitimate if a
//      restart is exact: one run of N steps must equal k runs of N/k chained through the state.
//      Asserted below for both cores, to 1e-12 relative -- if it ever stops holding, every stop
//      criterion is measuring a different problem from the one Terzaghi validated.
//
//  (2) THE DEFINITION. "Degree of consolidation" names two different numbers. Classically it is
//      the settlement ratio, settlement over final settlement; KATAI's input contract (`cdeg`)
//      defines the stop target instead as the PRESSURE ratio, the maximum excess pore pressure
//      over the maximum initial excess pore pressure p_max/p_max,initial. The two
//      are 21.6% apart in TIME on the one problem where both are known in closed form, so a build
//      that implements the name instead of the definition would be wrong by a fifth of the answer
//      and would look entirely reasonable doing it. The test pins the pressure ratio AGAINST the
//      settlement ratio: it asserts the stop time matches Tv = 1.0311 and that it is nowhere near
//      Tv = 0.8481, which is where the settlement ratio reaches the same 90%.
//
//  (3) A REFUSED SOLVE IS AN ANSWER. Driving the same column past its Mohr-Coulomb strength makes
//      the coupled tangent singular; the verifying backend refuses to return a vector that does not
//      satisfy the system, and both coupled cores used to let that refusal escape and ABORT THE
//      PROCESS (0xC0000409). They now end the time step with it, the way the static solver has
//      since it began checking its answers. A solver may report that it could not solve; it may
//      never answer by terminating.
//
// verify: KV-CON-003
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley -- one-dimensional consolidation; the stop criterion and the definition of its target are KATAI 2D input contract (docs/k2d-format.md, cstop / cdeg / cminp), which defines the degree of consolidation as the maximum excess pore pressure over the maximum initial excess pore pressure, not the settlement ratio the name suggests. Series derivation recorded in docs/references/consolidation-formulation.md
//   locator:  for a column with uniform initial excess pore pressure u0 drained at one end, u(Z,T)/u0 = sum_j (2/M) sin(M Z) exp(-M^2 T), M = (2j+1) pi/2, Z measured from the drained end, T = cv t / H_dr^2, cv = k Eoed / gamma_w; the MAXIMUM over the column is the value at the impermeable end, u(1,T)/u0 = sum_j (2/M) (-1)^j exp(-M^2 T). The settlement ratio is the other series, U(T) = 1 - sum_j (2/M^2) exp(-M^2 T). The reference the ratio is taken against has its own closed form: the undrained pressure a surcharge q generates in the confined column, p_u = q / (1 + n Eoed / Kw)
//   quantity: the time [day] at which a consolidation phase asked to stop at 90% degree of consolidation stops; the same time for a phase asked to stop at 1 kPa maximum excess pore pressure under the same 10 kPa surcharge (the same 10% of what it generates, so the same instant); the reference pressure of the ratio [kPa]; and the exactness of a restart of the fixed-dt core, which is what lets the march change its step at all
//   expected: u(1,T)/u0 = 0.1 at T = 1.031105, i.e. t = 14.566 day for this column (H_dr = 12 m, cv = k Eoed / gamma_w = 0.1 * 1000 / 9.81 = 10.194 m2/day) -- and NOT t = 11.980 day, which is where the SETTLEMENT ratio reaches the same 90% (T = 0.848085), 21.6% earlier. The settlement ratio at the stop is 93.9%, not 90%. The reference is p_u = 10 / (1 + 0.3333 * 1000 / 2e6) = 9.99833 kPa. A restart is exact to round-off
//   band:     2.5% on the stop time; measured +1.45% (14.776 vs 14.566), and the error is the log-time march's, not the mesh's: a fine equal-step grid gives 14.60 day on 99, 169 AND 453 nodes alike, while the march's error halves with each doubling of its steps-per-level (8 -> +5.5%, 16 -> +3.1%, 32 -> +1.5%, 64 -> +0.8%, 128 -> +0.4%; the shipped 32 is that measurement's choice). The other three are tight because they are identities rather than discretisations: the two criteria must agree to 0.5% (measured 0.00%), the elastoplastic path with the same criterion to 0.5% (measured 0.00%), the reference to 0.01 kPa (measured 0.00000%), and a restarted march to 1e-12 relative on displacement and pore pressure (measured 6e-16)
#include <katai/analysis/consolidation.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;
using katai::core::DofMap;
using katai::core::GaussState;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Permeability;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;
namespace linsolve = katai::linsolve;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// Maximum excess pore pressure in the column / u0: the value at the impermeable end (the profile
// increases monotonically away from the drainage boundary). This is the series the pressure-ratio
// degree of consolidation (`cdeg`) is a ratio of.
double terzaghi_pmax(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 200; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += (2.0 / M) * ((j % 2 == 0) ? 1.0 : -1.0) * std::exp(-M * M * Tv);
    }
    return s;
}
// The OTHER definition: the average degree of settlement.
double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 200; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += (2.0 / (M * M)) * std::exp(-M * M * Tv);
    }
    return 1.0 - s;
}
// Time factor at which a monotone decreasing f reaches `target` (bisection; the series are smooth).
double time_factor_at(double (*f)(double), double target) {
    double lo = 1e-4, hi = 8.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (f(mid) > target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

katai::core::ConsolidationSolveFactory nonsym_factory() {
    return [](const katai::math::CsrMatrix& A) {
        std::shared_ptr<linsolve::DirectSolver> s =
            linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(A);
        return std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(
            [s](const Eigen::VectorXd& b) { return s->solve(b); });
    };
}

// ---------------------------------------------------------------------------------------------
// (1) The restart identity. The whole adaptive march rests on it.
// ---------------------------------------------------------------------------------------------
struct Column {
    Mesh mesh; DofMap dofs; std::vector<char> drained;
    Column(Mesh m, DofMap d, std::vector<char> dr)
        : mesh(std::move(m)), dofs(std::move(d)), drained(std::move(dr)) {}
};
Column make_column(double W, double H, int nx, int ny) {
    RectangularDomain domain{0.0, 0.0, W, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri15(domain, nx, ny);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;
    return Column(std::move(mesh), std::move(dofs), std::move(drained));
}

void test_restart_identity() {
    std::printf("-- (1) restarting the core is exact: one run of 40 steps = four runs of 10 --\n");
    constexpr double W = 0.2, H = 1.0, E = 1000.0, nu = 0.0, u0 = 10.0;
    constexpr double k = 1.0e-4, gamma_w = 10.0, kw_over_n = 1.0e9;
    constexpr double dt = 2.0;
    constexpr int kSteps = 40, kChunks = 4;

    Column c = make_column(W, H, 1, 8);
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> p0(c.mesh.node_count, u0);

    for (int pass = 0; pass < 2; ++pass) {
        const bool plastic = pass == 1;
        // Mohr-Coulomb on the second pass: the elastoplastic core carries a COMMITTED Gauss state
        // between steps, and it is that state -- not just the pore field -- a restart must hand on.
        // c' huge: Mohr-Coulomb that never yields here. The identity being asserted is that the
        // COMMITTED Gauss state is carried across a restart, and it is carried whether or not it
        // is plastic; a column that yields under this u0 is laterally confined and cannot deform
        // into the yield, which drives the tangent singular -- that case is its own test below.
        MaterialModel mc; mc.type = MaterialType::MohrCoulomb;
        mc.youngs_modulus = E; mc.poisson_ratio = nu;
        mc.cohesion = 1.0e6; mc.friction_angle = 0.44; mc.dilatancy_angle = 0.0;   // [rad] ~ 25 deg
        const std::vector<MaterialModel> mm =
            plastic ? std::vector<MaterialModel>{mc}
                    : std::vector<MaterialModel>{{MaterialType::LinearElastic, E, nu}};

        Eigen::VectorXd whole_disp;
        Eigen::VectorXd whole_pore;
        std::vector<GaussState> whole_state;
        if (plastic) {
            const auto r = katai::core::solve_consolidation_plastic(
                c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n, c.drained, {}, p0, dt, kSteps, {},
                nullptr, nonsym_factory());
            check(r.converged, "MC: the whole run converged");
            whole_disp = r.series.displacement.back();
            whole_pore = r.series.pore.back();
            whole_state = r.committed;
        } else {
            const auto r = katai::core::solve_consolidation(c.mesh, c.dofs, mm, perm, gamma_w,
                                                            kw_over_n, c.drained, p0, dt, kSteps);
            whole_disp = r.displacement.back();
            whole_pore = r.pore.back();
        }

        // The same march, restarted three times. Displacement accumulates (each run starts from
        // zero); the pore field and the committed state ARE the state, and are handed on.
        Eigen::VectorXd part_disp = Eigen::VectorXd::Zero(whole_disp.size());
        std::vector<double> p(p0);
        std::vector<GaussState> state;
        Eigen::VectorXd part_pore;
        bool ok = true;
        for (int chunk = 0; chunk < kChunks; ++chunk) {
            if (plastic) {
                const auto r = katai::core::solve_consolidation_plastic(
                    c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n, c.drained, state, p, dt,
                    kSteps / kChunks, {}, nullptr, nonsym_factory());
                ok = ok && r.converged;
                part_disp += r.series.displacement.back();
                part_pore = r.series.pore.back();
                state = r.committed;
            } else {
                const auto r = katai::core::solve_consolidation(c.mesh, c.dofs, mm, perm, gamma_w,
                                                                kw_over_n, c.drained, p, dt,
                                                                kSteps / kChunks);
                part_disp += r.displacement.back();
                part_pore = r.pore.back();
            }
            for (int n = 0; n < c.mesh.node_count; ++n) p[n] = part_pore(n);
        }
        check(ok, plastic ? "MC: every chunk converged" : "LE: every chunk ran");

        const double dref = std::fmax(whole_disp.cwiseAbs().maxCoeff(), 1e-30);
        const double pref = std::fmax(whole_pore.cwiseAbs().maxCoeff(), 1e-30);
        const double derr = (whole_disp - part_disp).cwiseAbs().maxCoeff() / dref;
        const double perr = (whole_pore - part_pore).cwiseAbs().maxCoeff() / pref;
        double serr = 0.0;
        for (size_t g = 0; g < whole_state.size(); ++g)
            for (int i = 0; i < 3; ++i)
                serr = std::fmax(serr, std::fabs(whole_state[g].stress(i) - state[g].stress(i)));
        std::printf("   %s: displacement %.3e, pore %.3e (relative)%s\n", plastic ? "MC" : "LE",
                    derr, perr, plastic ? "" : "");
        if (plastic) std::printf("   MC: committed effective stress %.3e kPa (absolute)\n", serr);
        check(derr < 1e-12, plastic ? "MC: restarted displacement = one run"
                                    : "LE: restarted displacement = one run");
        check(perr < 1e-12, plastic ? "MC: restarted pore field = one run"
                                    : "LE: restarted pore field = one run");
        if (plastic) check(serr < 1e-9, "MC: restarted committed effective stress = one run");
    }
}

// A refused linear solve is an ANSWER. Found by this test on its first honest run: a confined
// column dissipating an excess pore pressure larger than its Mohr-Coulomb strength drives the
// coupled tangent singular, the verifying backend refuses to return a vector that does not satisfy
// the system -- and the elastoplastic Biot core let that refusal escape the call stack and ABORT
// THE PROCESS (0xC0000409). The static solver has caught the same refusal since it began checking
// its answers; both coupled cores now do too. A solver may report that it could not solve; it may
// never answer by terminating, because a crash is the one result no engineering judgement can be
// applied to.
void test_refused_solve_is_an_answer() {
    std::printf("\n-- (2) a singular tangent is REFUSED, not aborted --\n");
    constexpr double W = 0.2, H = 1.0, E = 1000.0, nu = 0.0, u0 = 10.0;
    constexpr double k = 1.0e-4, gamma_w = 10.0, kw_over_n = 1.0e9;
    Column c = make_column(W, H, 1, 8);
    // c' = 2 kPa, phi = 25 deg: with nu = 0 the horizontal effective stress stays zero while the
    // vertical one climbs to u0 = 10 kPa as the pore pressure dissipates, so the stress path runs
    // straight through the Mohr-Coulomb envelope (failure would need sigma_1 <= 3.14 c' = 6.3 kPa).
    MaterialModel mc; mc.type = MaterialType::MohrCoulomb;
    mc.youngs_modulus = E; mc.poisson_ratio = nu;
    mc.cohesion = 2.0; mc.friction_angle = 0.44; mc.dilatancy_angle = 0.0;
    const std::vector<MaterialModel> mm = {mc};
    const std::vector<Permeability> perm = {{k, k}};
    const std::vector<double> p0(c.mesh.node_count, u0);

    const auto r = katai::core::solve_consolidation_plastic(
        c.mesh, c.dofs, mm, perm, gamma_w, kw_over_n, c.drained, {}, p0, 2.0, 40, {},
        nullptr, nonsym_factory());
    std::printf("   converged = %d, steps recorded = %zu (the run RETURNED)\n",
                (int)r.converged, r.series.times.size());
    check(!r.converged, "the core reports that it could not solve, instead of terminating");
    check(!r.series.times.empty(), "and it hands back the steps it did complete");
}

// ---------------------------------------------------------------------------------------------
// (3) The stop criteria through the phase path, against the closed form.
// ---------------------------------------------------------------------------------------------
constexpr double kW = 1.0, kH = 12.0;      // column: single (top) drainage -> H_dr = H
constexpr double kE = 1000.0, kNu = 0.0;   // nu = 0 -> Eoed = E
constexpr double kPerm = 0.1;              // [m/day]
constexpr double kQ = 10.0;                // surcharge [kPa]

// The project of KV-CON-002 (a laterally confined column, drained at the top, surcharged at t = 0+)
// with the consolidation phase configured by the caller.
m::Project column_project(m::SoilModel model, bool load_in_phase = true) {
    m::Project pr;
    m::Material s; s.model = model;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 16.0; s.gamma_sat = 18.0; s.e_init = 0.5;
    s.c = 50.0; s.phi = 25.0; s.psi = 0.0;   // confined compression stays elastic -> = Terzaghi
    s.kx = kPerm; s.ky = kPerm;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};   // top drains (p = 0)
    P.edge_head = {0.0, 0.0, kH, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;   // pore = excess only (no hydrostatic background)

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 0; L.y1 = kH; L.x2 = kW; L.y2 = kH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};                       // surcharge OFF in the initial K0 phase

    m::Phase consol; consol.name = "Consolidation";
    consol.type = m::PhaseType::Consolidation;
    consol.load_active = {load_in_phase ? 1 : 0};        // ON -> applied at t = 0+
    pr.phases.push_back(consol);
    return pr;
}

std::vector<katai::core::SolveResult> run(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr, 0.4, 6);   // the KV-CON-002 mesh
    if (!M.ok) { check(false, "column meshed"); return {}; }
    return katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
}

void test_stop_criteria() {
    const double cv = kPerm * kE / katai::app::kGammaWater;   // nu = 0 -> Eoed = E
    const double s_inf = kQ * kH / kE;
    const double Tv_pressure = time_factor_at(terzaghi_pmax, 0.1);        // the cdeg definition
    const double Tv_settle = time_factor_at([](double T) { return 1.0 - terzaghi_U(T); }, 0.1);
    const double t_pressure = Tv_pressure * kH * kH / cv;
    const double t_settle = Tv_settle * kH * kH / cv;
    std::printf("\n-- (3) 90%% degree of consolidation: WHICH 90%%? --\n");
    std::printf("   pressure ratio (cdeg)     : Tv = %.6f -> t = %.4f day\n", Tv_pressure, t_pressure);
    std::printf("   settlement ratio (classic): Tv = %.6f -> t = %.4f day\n", Tv_settle, t_settle);
    std::printf("   the two are %.1f%% apart in time; the settlement ratio at the pressure-90%% "
                "instant is %.4f\n", 100.0 * (t_pressure / t_settle - 1.0), terzaghi_U(Tv_pressure));
    check(t_pressure > 1.15 * t_settle,
          "the closed forms themselves disagree by more than 15% -- the definition is not a detail");

    // --- degree of consolidation, 90% --------------------------------------------------------
    m::Project pr = column_project(m::SoilModel::LinearElastic);
    pr.phases[0].consol_stop = m::ConsolStop::DegreeOfConsolidation;
    pr.phases[0].consol_degree = 90.0;
    pr.phases[0].duration = 1.0;    // deliberately absurd: it must NOT be used
    pr.phases[0].time_steps = 25;
    const auto res = run(pr);
    check(res.size() == 2, "initial + consolidation phases ran");
    if (res.size() != 2 || !res[1].ok) {
        if (res.size() == 2) std::printf("   (%s)\n", res[1].message.c_str());
        check(false, "the phase stopped on its target");
        return;
    }
    const auto& C = res[1];
    const double t_stop = C.consol_time.back();
    std::printf("   stopped at t = %.4f day after %d steps (oracle %.4f, %+.2f%%)\n",
                t_stop, C.iterations, t_pressure, 100.0 * (t_stop / t_pressure - 1.0));
    std::printf("   excess pore %.4f kPa of the %.4f kPa generated -> degree %.4f\n",
                C.consol_excess_pore.back(), C.consol_pore_reference, C.consol_degree_reached);
    check(C.consol_stop == katai::core::ConsolidationStop::DegreeOfConsolidation,
          "the result carries WHICH ending it was asked for");
    check(C.consol_stop_met, "the result says the target was reached");
    check(std::fabs(t_stop / t_pressure - 1.0) < 0.025,
          "the stop time matches the closed-form PRESSURE ratio within 2.5%");
    check(t_stop > 1.15 * t_settle,
          "and it is NOT the settlement-ratio time: the definition is the pressure one");
    check(C.consol_degree_reached >= 0.90 && C.consol_degree_reached < 0.92,
          "the degree reached is just past the target, not far past it (the crossing is localised)");
    // The denominator of the ratio is a physical quantity with its own closed form: the undrained
    // pressure a surcharge generates in a confined column with a compressible pore fluid,
    // p_u = q / (1 + n Eoed / Kw). It is verified here because everything the criterion means
    // depends on it -- and because reading it off the transient is exactly where a too-small first
    // time step used to put a 12.8% drainage-boundary overshoot into the answer.
    const double n_por = 0.5 / 1.5;                       // e / (1 + e), the file's e_init = 0.5
    const double p_undrained = kQ / (1.0 + n_por * kE / 2.0e6);   // Kw = 2e6 kPa (the solver's water)
    std::printf("   reference %.5f kPa against the closed form %.5f (%+.4f%%)\n",
                C.consol_pore_reference, p_undrained,
                100.0 * (C.consol_pore_reference / p_undrained - 1.0));
    check(std::fabs(C.consol_pore_reference - p_undrained) < 0.01,
          "the reference of the ratio IS the undrained pressure the surcharge generated");
    // The settlement at that instant is NOT 90% of the final settlement -- it is 93.6%. A build
    // that had implemented the name would report this number as the target it met.
    const double U_settle_at_stop = C.consol_settlement.back() / s_inf;
    const double Tv_stop = cv * t_stop / (kH * kH);
    std::printf("   settlement ratio at the stop: %.4f (closed form at its own Tv = %.4f: %.4f)"
                " -- not 0.90\n", U_settle_at_stop, Tv_stop, terzaghi_U(Tv_stop));
    check(std::fabs(U_settle_at_stop - terzaghi_U(Tv_stop)) < 0.02,
          "the settlement ratio at the stop matches Terzaghi's OTHER series");
    check(U_settle_at_stop > 0.92, "and it is well above the 90% the phase was asked for");
    check(C.message.find("PRESSURE ratio") != std::string::npos,
          "the message carries the definition next to the number");

    // --- minimum excess pore pressure, 1 kPa -------------------------------------------------
    // 1 kPa is 10% of the 10 kPa this stage generates, so on THIS problem the two criteria ask
    // for the same instant. They are different code paths and different questions; agreeing here
    // is a fact about the problem, and it is worth pinning because it can only hold if both are
    // measuring the same pressure.
    m::Project pr2 = column_project(m::SoilModel::LinearElastic);
    pr2.phases[0].consol_stop = m::ConsolStop::MinExcessPore;
    pr2.phases[0].consol_min_pore = 0.1 * kQ;
    const auto res2 = run(pr2);
    if (res2.size() == 2 && res2[1].ok) {
        const double t2 = res2[1].consol_time.back();
        std::printf("\n   minimum excess pore 1 kPa: t = %.4f day (%.2f%% from the 90%% run)\n",
                    t2, 100.0 * (t2 / t_stop - 1.0));
        check(std::fabs(t2 / t_stop - 1.0) < 0.005,
              "1 kPa left of 10 kPa generated is the same instant as 90% dissipated");
        check(res2[1].consol_stop == katai::core::ConsolidationStop::MinExcessPore,
              "and the result says which of the two was asked for");
    } else {
        check(false, "the minimum-excess-pore phase solved");
        if (res2.size() == 2) std::printf("   (%s)\n", res2[1].message.c_str());
    }

    // --- the same target through the elastoplastic core ---------------------------------------
    // Confined compression under this surcharge stays elastic (c' = 50 kPa), so the coupled Newton
    // path must land on the same instant. It is the path a real consolidation phase takes.
    m::Project pr3 = column_project(m::SoilModel::MohrCoulomb);
    pr3.phases[0].consol_stop = m::ConsolStop::DegreeOfConsolidation;
    pr3.phases[0].consol_degree = 90.0;
    const auto res3 = run(pr3);
    if (res3.size() == 2 && res3[1].ok) {
        const double t3 = res3[1].consol_time.back();
        std::printf("\n   Mohr-Coulomb (coupled Newton): t = %.4f day (%.2f%% from the LE run)\n",
                    t3, 100.0 * (t3 / t_stop - 1.0));
        check(std::fabs(t3 / t_stop - 1.0) < 0.005,
              "the elastoplastic path stops at the same instant on an elastic stress path");
    } else {
        check(false, "the Mohr-Coulomb phase solved");
        if (res3.size() == 2) std::printf("   (%s)\n", res3[1].message.c_str());
    }
}

// ---------------------------------------------------------------------------------------------
// (4) What happens when the target is NOT reached, and when it was never real.
// ---------------------------------------------------------------------------------------------
void test_honest_endings() {
    std::printf("\n-- (4) a target that is not reached is refused, not reported --\n");
    m::Project pr = column_project(m::SoilModel::LinearElastic);
    pr.phases[0].consol_stop = m::ConsolStop::DegreeOfConsolidation;
    pr.phases[0].consol_degree = 90.0;
    pr.phases[0].consol_max_steps = 8;      // one level of the march: nowhere near 90%
    const auto res = run(pr);
    check(res.size() == 2, "the phase ran");
    if (res.size() == 2) {
        std::printf("   (%s)\n", res[1].message.c_str());
        check(!res[1].ok, "the phase REFUSED rather than reporting the time it happened to stop at");
        check(res[1].message.find("did not reach its target") != std::string::npos,
              "and the message says the target was not reached");
        check(res[1].message.find("degree of consolidation of") != std::string::npos,
              "and states how far it got, in the definition it was asked in");
    }

    std::printf("\n-- (5) a target that was already met at t = 0 is a warning, not a settlement --\n");
    // The surcharge is left OFF in the consolidation phase, so the stage changes nothing and
    // generates no excess pore pressure. The criterion is then true at the first time step, and the
    // time that comes out is the size of that step -- which is exactly the number a report would
    // otherwise present as "consolidated in 0.0004 days".
    m::Project pr2 = column_project(m::SoilModel::LinearElastic, /*load_in_phase=*/false);
    pr2.phases[0].consol_stop = m::ConsolStop::DegreeOfConsolidation;
    pr2.phases[0].consol_degree = 90.0;
    const auto res2 = run(pr2);
    check(res2.size() == 2, "the no-change phase ran");
    if (res2.size() == 2) {
        bool warned = false;
        for (const auto& d : res2[1].diagnostics)
            if (d.code == std::string("K2D-A013")) warned = true;
        std::printf("   solved=%d, peak excess pore = %.3e kPa, diagnostics = %zu\n",
                    (int)res2[1].ok, res2[1].consol_pore_reference, res2[1].diagnostics.size());
        for (const auto& d : res2[1].diagnostics)
            if (d.code == std::string("K2D-A013")) std::printf("   (%s)\n", d.message.c_str());
        check(res2[1].ok, "the phase still solves (the state is legitimate, if degenerate)");
        check(warned, "K2D-A013: it says the target was met because nothing generated pressure");
    }
}

}  // namespace

int main() {
    std::printf("Consolidation ended on a state: the cstop stop criteria\n\n");
    test_restart_identity();
    test_refused_solve_is_an_answer();
    test_stop_criteria();
    test_honest_endings();
    if (g_failures == 0) {
        std::printf("\nOK: the march restarts exactly, and it stops where the closed form says\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
