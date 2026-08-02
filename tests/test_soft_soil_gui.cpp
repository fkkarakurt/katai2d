// SOFT SOIL through the FULL GUI compute path (Stage 3: build_problem mapping + K0 seeding with
// OCR/POP + honest guards). Oracles are the model's own 1-D ln-law closed forms on a laterally
// confined column with self-weight sigma_v0(z) = gamma z (depth z from the surface):
//   NC surcharge q:      u_top = lambda* * I(q),   I(q) = (H+a) ln(H+a) - H ln H - a ln a,  a = q/gamma
//                        (eps_v(z) = lambda* ln((sigma_v0+q)/sigma_v0) integrated over depth; on the
//                        K0NC line p is proportional to sigma_v, so the ratio is the same in p);
//   OC (OCR / POP) q:    the same integral on the kappa* line while sigma_v stays below the
//                        preconsolidation -- the app path must SEED p_p from OCR/POP for this.
// Bands are honest: the kernel floors p' at 1 kPa near the (stress-free) surface and a thin
// near-surface strip re-yields under OCR (q > (OCR-1) sigma_v0 there), so OC checks use class
// bands and a DISCRIMINATION ratio (NC/OC settlement), not exact equality.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kGammaSoil = 17.0, kH = 8.0, kW = 4.0;
constexpr double kLam = 0.10, kKap = 0.02;

m::Project ss_column() {
    m::Project pr;
    m::Material s;
    s.name = "Soft clay";
    s.model = m::SoilModel::SoftSoil;
    s.lam_star = kLam; s.kap_star = kKap; s.nu_ur = 0.15;
    s.c = 1.0; s.phi = 25.0; s.psi = 0.0;
    s.k0nc_auto = true;
    s.gamma_unsat = kGammaSoil; s.gamma_sat = 19.0;
    pr.materials.push_back(s);
    m::SoilPolygon L; L.material = 0;
    L.x = {0, kW, kW, 0}; L.y = {0, 0, kH, kH};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(L);
    pr.has_water = false;
    return pr;
}

void add_surcharge_phase(m::Project& pr, double q) {
    m::Load ld; ld.kind = m::LoadKind::Distributed; ld.name = "q";
    ld.x1 = 0; ld.y1 = kH; ld.x2 = kW; ld.y2 = kH; ld.qy1 = -q; ld.qy2 = -q; ld.qx1 = ld.qx2 = 0;
    pr.loads.push_back(ld);
    pr.initial.load_active = {0};
    m::Phase ph; ph.name = "Surcharge"; ph.poly_active = {1}; ph.load_active = {1};
    pr.phases.push_back(ph);
}

double uy_top(const katai::app::SolveResult& R) {
    // Surface settlement: mean uy of the top-edge nodes (uniform in 1D, the mean averages noise).
    double s = 0.0; int n = 0;
    for (int i = 0; i < R.mesh.node_count; ++i)
        if (R.mesh.y[i] > kH - 1e-6) { s += R.disp[i * 2 + 1]; ++n; }
    return n > 0 ? s / n : 0.0;
}

// 1-D ln-law depth integral: I(a) = (H+a) ln(H+a) - H ln H - a ln a.
double depth_integral(double a) {
    return (kH + a) * std::log(kH + a) - kH * std::log(kH) - a * std::log(a);
}
// NC (virgin) surcharge: on the K0NC line p is proportional to sigma_v and the elastoplastic
// oedometric relation is eps_1 = lambda* ln(sigma_v/sigma_v0) -> u = lambda* I(q/gamma).
double lnlaw_integral(double q) { return depth_integral(q / kGammaSoil); }
// ELASTIC (kappa*) surcharge -- the OC/POP oracle. On the elastic path the lateral ratio of the
// INCREMENT is nu/(1-nu) (not K0NC) and Eoed_el = K + 4G/3 = 3K(1-nu)/(1+nu) with K = p'/kappa*.
// Writing p along the path as p = C sigma_v0 + B dsigma_v, B = (1+nu)/(3(1-nu)), C = (1+2K0NC)/3,
// the identity (3(1-nu)/(1+nu)) * B = 1 collapses the integral to
//   eps_1(z) = kappa* ln(1 + (B/C) q / sigma_v0(z))   ->   u = kappa* I(r q / gamma), r = B/C.
double elastic_integral(double q, double nu, double k0nc) {
    const double B = (1.0 + nu) / (3.0 * (1.0 - nu)), C = (1.0 + 2.0 * k0nc) / 3.0;
    return depth_integral((B / C) * q / kGammaSoil);
}
const double kK0nc = 1.0 - std::sin(25.0 * 3.14159265358979323846 / 180.0);

// Run K0 initial + surcharge phase; return (initial max|u|, phase surface settlement).
bool run_column(const m::Project& pr, double& u0, double& utop, std::string& msg) {
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    if (!M.ok) { msg = "mesh failed"; return false; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2 || !res[0].ok || !res[1].ok) {
        msg = res.empty() ? "no result" : res.back().message;
        return false;
    }
    u0 = res[0].max_disp;
    utop = uy_top(res[1]);
    return true;
}

void test_nc_surcharge() {
    std::printf("-- NC column: K0 seed (pp = f_bar(sigma0)) + surcharge on the lambda* line --\n");
    m::Project pr = ss_column();
    add_surcharge_phase(pr, 50.0);
    double u0 = 0, ut = 0; std::string msg;
    const bool ok = run_column(pr, u0, ut, msg);
    check(ok, "K0 + surcharge phases converged");
    if (!ok) { std::printf("   (%s)\n", msg.c_str()); return; }
    std::printf("   initial max|u| = %.3e m (admissible NC seed)\n", u0);
    check(u0 < 5e-3, "K0 phase is (near-)undisturbed: the ss_initial_pp seed is admissible");
    const double u_ex = -kLam * lnlaw_integral(50.0);
    std::printf("   surcharge: u_top = %.4f m (1-D ln-law %.4f, err %+.1f%%)\n", ut, u_ex,
                100.0 * (ut - u_ex) / std::fabs(u_ex));
    check(std::fabs(ut - u_ex) < 0.05 * std::fabs(u_ex),
          "NC settlement matches the lambda* depth integral (5%)");
}

void test_ocr_memory() {
    std::printf("-- OCR = 4: the app path must seed p_p from OCR (kappa*-class response) --\n");
    const double q = 10.0;
    double u0 = 0, unc = 0, uoc = 0; std::string msg;
    m::Project nc = ss_column(); add_surcharge_phase(nc, q);
    bool ok = run_column(nc, u0, unc, msg);
    check(ok, "NC reference column converged");
    if (!ok) { std::printf("   (%s)\n", msg.c_str()); return; }

    m::Project oc = ss_column();
    oc.materials[0].oc_mode = 1; oc.materials[0].OCR = 4.0;
    add_surcharge_phase(oc, q);
    ok = run_column(oc, u0, uoc, msg);
    check(ok, "OCR = 4 column converged");
    if (!ok) { std::printf("   (%s)\n", msg.c_str()); return; }

    // OCR now raises the automatic K0 too (an audit fix; the PLAXIS Ref formula):
    //   K0_oc = K0nc*OCR - nu_ur/(1-nu_ur)*(OCR-1)  -> the initial lateral coefficient C changes too.
    const double k0oc = kK0nc * 4.0 - 0.15 / 0.85 * 3.0;
    const double uk = -kKap * elastic_integral(q, 0.15, k0oc);
    std::printf("   u_NC = %.4f  u_OCR4 = %.4f  (elastic kappa* form %.4f, K0_oc = %.3f)  "
                "ratio = %.2f\n", unc, uoc, uk, k0oc, unc / uoc);
    check(unc / uoc > 3.0, "preconsolidation memory: NC settles > 3x more than OCR = 4");
    check(std::fabs(uoc - uk) < 0.15 * std::fabs(uk),
          "OCR = 4 response matches the elastic kappa* integral (15%; thin surface strip re-yields)");
}

void test_pop_memory() {
    std::printf("-- POP = 100 kPa: previously ignored silently, now an equivalent-ratio seed --\n");
    const double q = 10.0;
    double u0 = 0, up = 0; std::string msg;
    m::Project pp = ss_column();
    pp.materials[0].oc_mode = 2; pp.materials[0].POP = 100.0;
    add_surcharge_phase(pp, q);
    const bool ok = run_column(pp, u0, up, msg);
    check(ok, "POP = 100 column converged");
    if (!ok) { std::printf("   (%s)\n", msg.c_str()); return; }
    // q = 10 stays inside the preconsolidation EVERYWHERE (sigma_v0 + q < sigma_v0 + POP), so the
    // whole depth responds elastically -- a tighter band than the OCR case is justified (the
    // remaining deviation is the kernel's p' >= 1 kPa floor near the stress-free surface).
    const double uk = -kKap * elastic_integral(q, 0.15, kK0nc);
    std::printf("   u_POP = %.4f (elastic kappa* form %.4f, err %+.1f%%)\n", up, uk,
                100.0 * (up - uk) / std::fabs(uk));
    check(std::fabs(up - uk) < 0.10 * std::fabs(uk),
          "POP response matches the elastic kappa* integral (10%)");
}

void test_undrained_refusal() {
    std::printf("-- honest guard: Soft Soil + Undrained (A) is refused, not silently drained --\n");
    m::Project pr = ss_column();
    pr.materials[0].drainage = m::Drainage::Undrained;
    add_surcharge_phase(pr, 10.0);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    bool refused = false;
    for (const auto& r : res)
        if (!r.ok && r.message.find("Soft Soil") != std::string::npos) refused = true;
    check(refused, "refused with an explicit Soft Soil + Undrained message");
}

void test_dynamic_stiffness_smoke() {
    std::printf("-- dynamic (linear) after K0: per-Gauss E from K_ur = p'/kappa* --\n");
    m::Project pr = ss_column();
    m::Phase dyn; dyn.name = "Shake"; dyn.type = m::PhaseType::Dynamic; dyn.poly_active = {1};
    dyn.duration = 0.25; dyn.time_steps = 25;
    dyn.seismic_amp = 0.5; dyn.seismic_freq = 2.0;
    pr.phases.push_back(dyn);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    const bool ok = res.size() == 2 && res[0].ok && res[1].ok;
    check(ok, "K0 + dynamic phases ran");
    if (!ok) { for (const auto& r : res) if (!r.ok) std::printf("   (%s)\n", r.message.c_str()); return; }
    std::printf("   model f1 estimate = %.3f Hz (from the committed-stress SS stiffness)\n",
                res[1].dyn_model_f1);
    check(res[1].dyn_model_f1 > 0.1 && res[1].dyn_model_f1 < 100.0,
          "f1 estimated from the stress-dependent SS stiffness is finite and physical");
}

// ===================== SOFT SOIL CREEP (Stage 3 GUI wiring) =====================
constexpr double kMu = 0.005;

m::Project ssc_column() {
    m::Project pr = ss_column();
    pr.materials[0].name = "Creeping clay";
    pr.materials[0].model = m::SoilModel::SoftSoilCreep;
    pr.materials[0].mu_star = kMu;
    return pr;
}

// SSC is an ISOTACH model: "primary" compression is creep at a high rate. So the
// GUI wiring is pinned by TIME CLASSES on the same column:
//   T = 0 (timeless stage)  -> elastic (kappa*) + MC only: far stiffer than lambda*;
//   T = tau = 1 day         -> the 1-day NC isotach IS the lambda* line (Eq 11-14's
//                              definition), so u matches Soft Soil's NC integral;
//   T = 100 vs 1 day        -> depth-UNIFORM secondary creep, du = mu* H ln(100).
// The phase's 'Time interval [day]' drives all of this through the solver's
// SumMstage-proportional time distribution (Stage 2 contract) -- exactly the
// plumbing this test exists to pin.
void test_ssc_time_classes() {
    std::printf("-- SSC time classes: T = 0 (kappa*), T = tau (~lambda* isotach), creep tail --\n");
    const double q = 50.0;
    double u_seed = -1.0;
    // Load WITH time: the phase's 'Time interval' is distributed over the load
    // increments in Delta-lambda proportion (the Stage 2 SumMstage contract) --
    // the PLAXIS-realistic staged use. (The pathological alternative -- a TIMELESS
    // full load, then a pure hold -- leaves the state overstressed by p_eq/pp up
    // to ~50 with a creep rate ~ ratio^((lam*-kap*)/mu*) = ratio^16, a shock the
    // one-increment plastic phase honestly fails on; holding periods belong in a
    // CONSOLIDATION phase, where pore pressure raises p' gradually and time is
    // sub-stepped. Stated in the mu* help text.)
    auto run_T = [&](double days, double& ut) {
        m::Project pr = ssc_column();
        add_surcharge_phase(pr, q);
        pr.phases.back().duration = days;
        double u0 = 0;
        std::string msg;
        const bool ok = run_column(pr, u0, ut, msg);
        if (!ok) std::printf("   (T=%g day: %s)\n", days, msg.c_str());
        if (ok) u_seed = std::max(u_seed, u0);
        return ok;
    };
    double u_0d = 0, u_1d = 0, u_100d = 0;
    check(run_T(0.0, u_0d) && run_T(1.0, u_1d) && run_T(100.0, u_100d),
          "SSC columns converged for T = 0, 1 and 100 days");
    if (u_seed < 0.0) return;
    std::printf("   initial max|u| = %.3e m (SSC K0 seed)\n", u_seed);
    check(u_seed < 5e-3, "K0 seed admissible for SSC (initial phase near-undisturbed)");
    const double u_lam = -kLam * lnlaw_integral(q);   // Soft Soil's NC (lambda*) class
    std::printf("   u(T=0) = %.4f  u(T=1d) = %.4f  u(T=100d) = %.4f  (lambda* form %.4f)\n",
                u_0d, u_1d, u_100d, u_lam);
    // Loading over tau = 1 day approaches the 1-day NC isotach = the lambda* line;
    // the SumMstage distribution leaves the end state slightly YOUNGER than 1 day
    // (later increments arrive late with little time to relax), so the settlement
    // sits somewhat above the closed form -- measured -18%; a 25% class band is the
    // honest statement. The 100-vs-1-day DIFFERENCE cancels that youth deficit and
    // pins the mu* machinery exactly.
    check(std::fabs(u_1d - u_lam) < 0.25 * std::fabs(u_lam),
          "loading over tau = 1 day lands near the lambda* NC isotach (25% class)");
    check(std::fabs(u_0d) < 0.5 * std::fabs(u_lam),
          "a TIMELESS stage is far stiffer (isotach: no time -> no 'primary' creep)");
    const double d_sec = u_100d - u_1d;
    const double d_ex = -kMu * kH * std::log(100.0);
    std::printf("   creep tail u(100d) - u(1d) = %.4f m (mu* H ln100 = %.4f, err %+.1f%%)\n",
                d_sec, d_ex, 100.0 * (d_sec - d_ex) / std::fabs(d_ex));
    check(std::fabs(d_sec - d_ex) < 0.10 * std::fabs(d_ex),
          "the 100 vs 1 day difference is the depth-uniform mu* ln(t) creep (10%)");
}

// HOLDING periods for SSC. The consolidation route was tried first and did NOT
// verify: on a fast-draining column the settlement DECREASED with time-step
// refinement (T=100 day: -0.519 at 50 steps, -0.506 at 100, diverged at 200;
// T=1 day: -0.442 at 10 steps, -0.403 at 25, diverged at 50) and the creep tail
// came out 37% short of the exact mu* H ln(100) oracle the verified plastic-phase
// path meets at +0.1%. A recommendation that cannot be verified must not ship:
// SSC x Consolidation is now REFUSED (pinned below) and the working alternative is
// pinned instead -- load WITH time, then hold in a timed NIL Plastic phase (no
// creep shock: the state sits on its isotach when the hold starts, rates are
// moderate). Isotach uniqueness says the chain must land where the single ramped
// run lands.
void test_ssc_consolidation_hold() {
    std::printf("-- SSC holding: timed-NIL chain verified; consolidation refused --\n");
    // Two-phase hold: timed loading (T=1, verified) + a timed NIL hold (T=99).
    // No shock: the state sits on the isochrone, the creep rate is moderate — it must land
    // where the single-phase T=100 ramped run does (isochrone uniqueness).
    auto run_hold_chain = [&](double t_load, double t_hold, double& ut) {
        m::Project pr = ssc_column();
        add_surcharge_phase(pr, 50.0);
        pr.phases.back().duration = t_load;
        m::Phase hold;
        hold.name = "hold";
        hold.type = m::PhaseType::Plastic;
        hold.poly_active = {1};
        hold.load_active = {1};
        hold.duration = t_hold;
        pr.phases.push_back(hold);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.size() != 3 || !res[0].ok || !res[1].ok || !res[2].ok) {
            std::printf("   (chain %g+%g day: %s)\n", t_load, t_hold,
                        res.empty() ? "no result" : res.back().message.c_str());
            return false;
        }
        ut = uy_top(res[1]) + uy_top(res[2]);
        return true;
    };
    double u_chain = 0, u_ramp = 0;
    {
        m::Project pr = ssc_column();
        add_surcharge_phase(pr, 50.0);
        pr.phases.back().duration = 100.0;
        double u0 = 0;
        std::string msg;
        check(run_column(pr, u0, u_ramp, msg), "single-phase T = 100 day reference converged");
    }
    check(run_hold_chain(1.0, 99.0, u_chain),
          "load-with-time + timed NIL hold chain converged (no creep shock: the "
          "state sits on its isotach when the hold starts)");
    // The chain ends at a TRUE age of ~100 days, so it should sit on the ideal
    // 100-day isotach: u = lambda* I(q/gamma) + mu* H ln(100/tau). The single ramped
    // run is genuinely YOUNGER (the SumMstage distribution hands later increments
    // little time), so it settles LESS -- that ordering is the physics, not an error.
    const double u_ideal = -kLam * lnlaw_integral(50.0) - kMu * kH * std::log(100.0);
    std::printf("   hold chain u(1d + 99d) = %.4f m (ideal 100-day isotach %.4f, err %+.1f%%); "
                "ramped u(100d) = %.4f m\n",
                u_chain, u_ideal, 100.0 * (u_chain - u_ideal) / std::fabs(u_ideal), u_ramp);
    check(std::fabs(u_chain - u_ideal) < 0.05 * std::fabs(u_ideal),
          "the aged chain sits on the ideal 100-day isotach closed form (5%)");
    check(std::fabs(u_chain) > std::fabs(u_ramp),
          "the ramped run is younger and settles less (SumMstage youth, physical)");
    // SSC x Consolidation: honest refusal until the Biot x creep interaction verifies.
    {
        m::Project pr = ssc_column();
        m::Phase ph;
        ph.name = "consolidate";
        ph.type = m::PhaseType::Consolidation;
        ph.poly_active = {1};
        ph.duration = 10.0;
        ph.time_steps = 10;
        pr.phases.push_back(ph);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool refused = false;
        for (const auto& r : res)
            if (!r.ok && r.message.find("Soft Soil Creep") != std::string::npos) refused = true;
        check(refused, "SSC + consolidation refused with an explicit message (unverified path)");
    }
}

void test_ssc_undrained_refusal() {
    std::printf("-- honest guard: Soft Soil Creep + Undrained (A) refused too --\n");
    m::Project pr = ssc_column();
    pr.materials[0].drainage = m::Drainage::Undrained;
    add_surcharge_phase(pr, 10.0);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    bool refused = false;
    for (const auto& r : res)
        if (!r.ok && r.message.find("Soft Soil") != std::string::npos) refused = true;
    check(refused, "refused with an explicit Soft Soil (Creep) + Undrained message");
}

void test_ssc_dynamic_smoke() {
    std::printf("-- SSC dynamic (linear) after K0: per-Gauss E from K_ur = p'/kappa* --\n");
    m::Project pr = ssc_column();
    m::Phase dyn; dyn.name = "Shake"; dyn.type = m::PhaseType::Dynamic; dyn.poly_active = {1};
    dyn.duration = 0.25; dyn.time_steps = 25;
    dyn.seismic_amp = 0.5; dyn.seismic_freq = 2.0;
    pr.phases.push_back(dyn);
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    const bool ok = res.size() == 2 && res[0].ok && res[1].ok;
    check(ok, "K0 + SSC dynamic phases ran");
    if (!ok) { for (const auto& r : res) if (!r.ok) std::printf("   (%s)\n", r.message.c_str()); return; }
    std::printf("   model f1 estimate = %.3f Hz (from the committed-stress SSC stiffness)\n",
                res[1].dyn_model_f1);
    check(res[1].dyn_model_f1 > 0.1 && res[1].dyn_model_f1 < 100.0,
          "f1 estimated from the stress-dependent SSC stiffness is finite and physical");
}

}  // namespace

int main() {
    std::printf("SOFT SOIL through the GUI compute path -- Stage 3 V&V\n\n");
    test_nc_surcharge();
    std::printf("\n");
    test_ocr_memory();
    std::printf("\n");
    test_pop_memory();
    std::printf("\n");
    test_undrained_refusal();
    std::printf("\n");
    test_dynamic_stiffness_smoke();
    std::printf("\n");
    test_ssc_time_classes();
    std::printf("\n");
    test_ssc_consolidation_hold();
    std::printf("\n");
    test_ssc_undrained_refusal();
    std::printf("\n");
    test_ssc_dynamic_smoke();
    if (g_failures == 0) {
        std::printf("\nOK: Soft Soil (+ Creep) is wired end-to-end (mapping, K0 + OCR/POP p_p "
                    "seeding, phase time -> creep, honest guards, dynamic stiffness)\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
