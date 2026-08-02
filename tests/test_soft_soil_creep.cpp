// SOFT SOIL CREEP -- Stage 1 material-point V&V (PLAXIS MMM sec 11; locked formulation in
// docs/references/soft-soil-creep-formulation.md). Oracles are the model's DEFINING closed forms,
// derived independently in the formulation doc (sec 5) and none share code with the kernel:
//   (a) CONSTANT isotropic stress, NC start: eps_v^c(t) = mu* ln(1 + t/tau) EXACTLY
//       (Buisman/Garlanger secondary compression, slope mu* per ln-cycle);
//   (b) initial creep-rate ratio between OCR = 2 and NC equals OCR^-beta, beta = (lam*-kap*)/mu*
//       (the overconsolidation kill-switch that MMM 11.11 warns about);
//   (c) isotropic stress RELAXATION (zero total strain): p(t) = p0 (1 + (lam*/(kap* tau)) t)^(-mu*/lam*)
//       and p_p grows exactly by (p0/p)^(kap*/(lam*-kap*)) -- both pinned;
//   (d) 24-hour staged oedometer: end-of-day states advance on the lam* NC line (the tau = 1 day
//       DEFINITION of normal consolidation, Eq 11-13/14) and develop K0NC laterally;
//   (e) FAST drained triaxial: failure ON the Mohr-Coulomb line (creep has no failure of its own;
//       MC is checked AFTER the creep update, manual sec 11.7).
// STAGE 2 -- FE constitutive wiring with TIME: dt enters integrate_point as a trailing
// parameter (0 = no creep, the bit-identical legacy path); the static solver distributes the
// phase time_interval over increments by d(lambda). The FE section pins the WIRING (the kernel
// is pinned by (a)-(e) above):
//   (g1) plane-strain integrate_point chain (strain + time) == the principal kernel chain;
//   (g2) axisymmetric branch == the same kernel chain;
//   (g3) tri15 column BVP: compress by prescribed settlement (with time), then HOLD the
//        displacement over further timed solves -- every Gauss point relaxes exactly like the
//        kernel chain (assembler + Newton + FD tangent + time bookkeeping end-to-end).
#include <katai/materials/soft_soil_creep.hpp>

#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace katai::core::softsoilcreep;
namespace ss = katai::core::softsoil;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

Params clay() {
    Params P;
    P.lam_star = 0.10;
    P.kap_star = 0.02;
    P.mu_star = 0.005;   // lam*/mu* = 20, beta = 16
    P.nu_ur = 0.15;
    P.c = 1.0;
    P.phi = 25.0 * kPi / 180.0;
    P.psi = 0.0;
    P.K0nc = 1.0 - std::sin(P.phi);
    return P;
}

// Hold the ISOTROPIC stress constant at p0 while time runs: per step, solve the isotropic strain
// increment (secant) so that the returned mean stress stays at p0. Total eps_v then equals the
// accumulated creep strain (the elastic state does not change).
struct CreepRun { double eps_v; Eigen::Vector3d sig; double pp; };
CreepRun hold_isotropic(const Params& P, double p0, double pp0, double t_end) {
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = pp0, eps_v = 0.0, t = 0.0;
    double de_guess = 0.0;
    while (t < t_end - 1e-12) {
        const double dt = std::min(t < 2.0 ? 0.02 : 0.49, t_end - t);
        double de = de_guess;   // warm start
        Eigen::Vector3d s_new = sig; double pp_new = pp, dev_new = 0.0;
        for (int it = 0; it < 60; ++it) {
            const auto r = ssc_step(P, sig, pp, Eigen::Vector3d::Constant(de / 3.0), dt);
            const double res = r.sig.mean() - p0;
            s_new = r.sig; pp_new = r.pp; dev_new = de;
            if (std::fabs(res) < 1e-8 * p0) break;
            const double h = 1e-9 + 1e-6 * std::fabs(de);
            const auto r2 = ssc_step(P, sig, pp, Eigen::Vector3d::Constant((de + h) / 3.0), dt);
            const double drde = (r2.sig.mean() - r.sig.mean()) / h;
            if (!std::isfinite(drde) || std::fabs(drde) < 1e-30) break;
            de -= res / drde;
        }
        sig = s_new; pp = pp_new; eps_v += dev_new; t += dt;
        de_guess = dev_new;
    }
    return {eps_v, sig, pp};
}

void test_constant_stress_creep() {
    std::printf("-- (a) constant isotropic stress: eps_v^c = mu* ln(1 + t/tau) (exact form) --\n");
    const Params P = clay();
    const double p0 = 100.0;
    for (double t_end : {2.0, 100.0}) {
        const auto r = hold_isotropic(P, p0, p0, t_end);
        const double ex = P.mu_star * std::log(1.0 + t_end / P.tau_day);
        std::printf("   t = %6.1f d: eps_v = %.6f (exact %.6f, err %+.2f%%)   pp = %.2f\n", t_end,
                    r.eps_v, ex, 100.0 * (r.eps_v - ex) / ex, r.pp);
        check(std::fabs(r.eps_v - ex) < 0.01 * ex,
              "secondary compression follows mu* ln(1+t/tau) within 1%");
    }
}

void test_ocr_rate_ratio() {
    std::printf("-- (b) initial creep-rate ratio = OCR^-beta (the overconsolidation kill-switch) --\n");
    const Params P = clay();
    const double p0 = 100.0, dt = 1e-4;
    const Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    const auto r_nc = ssc_step(P, sig, p0, Eigen::Vector3d::Zero(), dt);
    const auto r_oc = ssc_step(P, sig, 2.0 * p0, Eigen::Vector3d::Zero(), dt);
    const double beta = (P.lam_star - P.kap_star) / P.mu_star;
    const double ratio = r_nc.devc / r_oc.devc, ex = std::pow(2.0, beta);
    std::printf("   L_NC = %.3e  L_OCR2 = %.3e  ratio = %.0f (2^beta = %.0f, err %+.2f%%)\n",
                r_nc.devc, r_oc.devc, ratio, ex, 100.0 * (ratio - ex) / ex);
    check(std::fabs(ratio - ex) < 0.03 * ex, "rate ratio equals OCR^beta within 3%");
}

void test_relaxation() {
    std::printf("-- (c) isotropic relaxation (zero total strain): closed-form p(t) and pp(t) --\n");
    const Params P = clay();
    const double p0 = 100.0, t_end = 100.0;
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0, t = 0.0;
    while (t < t_end - 1e-12) {
        const double dt = std::min(t < 1.0 ? 0.01 : 0.1, t_end - t);
        const auto r = ssc_step(P, sig, pp, Eigen::Vector3d::Zero(), dt);
        sig = r.sig; pp = r.pp; t += dt;
    }
    const double p_ex =
        p0 * std::pow(1.0 + (P.lam_star / (P.kap_star * P.tau_day)) * t_end, -P.mu_star / P.lam_star);
    const double pp_ex = p0 * std::pow(p0 / p_ex, P.kap_star / (P.lam_star - P.kap_star));
    std::printf("   p(100d) = %.3f (exact %.3f, err %+.2f%%)   pp = %.3f (exact %.3f)\n",
                sig.mean(), p_ex, 100.0 * (sig.mean() - p_ex) / p_ex, pp, pp_ex);
    check(std::fabs(sig.mean() - p_ex) < 0.005 * p_ex,
          "relaxed stress follows the closed form within 0.5%");
    check(std::fabs(pp - pp_ex) < 0.005 * pp_ex,
          "aged pre-consolidation follows (p0/p)^(kap/(lam-kap)) within 0.5%");
}

// Oedometer driver: solve the axial strain increment d1 (lateral = 0) so sigma_1 tracks a target.
struct OedoState { Eigen::Vector3d sig; double pp; double eps1; };
void oedo_to(const Params& P, OedoState& st, double s1_target, double time, int nstep) {
    const double dt = time / nstep;
    double d_guess = 0.0;
    for (int i = 0; i < nstep; ++i) {
        const double s_tgt = st.sig(0) + (s1_target - st.sig(0)) * 1.0;   // hold/track target
        double d = d_guess;
        Eigen::Vector3d s_new = st.sig; double pp_new = st.pp, d_new = 0.0;
        for (int it = 0; it < 60; ++it) {
            const auto r = ssc_step(P, st.sig, st.pp, Eigen::Vector3d(d, 0.0, 0.0), dt);
            const double res = r.sig(0) - s_tgt;
            s_new = r.sig; pp_new = r.pp; d_new = d;
            if (std::fabs(res) < 1e-7 * std::max(1.0, s_tgt)) break;
            const double h = 1e-9 + 1e-6 * std::fabs(d);
            const auto r2 = ssc_step(P, st.sig, st.pp, Eigen::Vector3d(d + h, 0.0, 0.0), dt);
            const double drdd = (r2.sig(0) - r.sig(0)) / h;
            if (!std::isfinite(drdd) || std::fabs(drdd) < 1e-30) break;
            double step = -res / drdd;
            step = std::clamp(step, -0.02, 0.02);
            d += step;
        }
        st.sig = s_new; st.pp = pp_new; st.eps1 += d_new;
        d_guess = d_new;
    }
}

void test_staged_oedometer() {
    std::printf("-- (d) 24-hour staged oedometer: end-of-day states on the lam* NC line + K0NC --\n");
    const Params P = clay();
    const ss::Params S = P.ss();
    // NC oedometric start at sigma_v = 100 (K0NC laterals, pp = p_eq: exactly on the surface).
    OedoState st;
    st.sig = Eigen::Vector3d(100.0, P.K0nc * 100.0, P.K0nc * 100.0);
    st.pp = ss::ss_initial_pp(S, st.sig, 1.0);
    st.eps1 = 0.0;
    // Let the start state age one day so it is a true end-of-day NC point (stage-0 reference).
    // Stage 1 still carries a seating transient (the aged state is slightly overconsolidated),
    // so the lam* pin reads the STEADY stages (2-3); stage 1 is printed for the record.
    oedo_to(P, st, 100.0, 1.0, 40);
    double s_prev = st.sig(0), e_prev = st.eps1;
    double slope_min = 1e9, slope_max = -1e9;
    int stage = 0;
    for (double target : {150.0, 225.0, 337.5}) {
        oedo_to(P, st, target, 0.001, 25);   // fast load (~86 s)
        oedo_to(P, st, target, 1.0, 40);     // hold 24 h
        const double slope = (st.eps1 - e_prev) / std::log(st.sig(0) / s_prev);
        std::printf("   stage end: sigma_v = %.1f  eps_1 = %.5f  slope = %.5f  K0 = %.4f\n",
                    st.sig(0), st.eps1, slope, st.sig(1) / st.sig(0));
        if (++stage >= 2) {
            slope_min = std::min(slope_min, slope);
            slope_max = std::max(slope_max, slope);
        }
        s_prev = st.sig(0); e_prev = st.eps1;
    }
    check(slope_min > 0.98 * P.lam_star && slope_max < 1.02 * P.lam_star,
          "steady end-of-day compression slope = lam* within 2% (tau = 1 day defines the NC line)");
    // K0NC pin: the STAGED path's end-of-day K0 sits ABOVE K0NC by design (the fast load drops
    // K0 elastically, the hold recovers it -- the end-of-day value is the upper turning point).
    // The M(K0NC) calibration is pinned on the CONTINUOUS (CRS) path instead: in the isotache
    // framework a constant-rate oedometer holds the same stress RATIO at any rate.
    for (int i = 0; i < 600; ++i) {
        double d = 2e-4;
        Eigen::Vector3d s_new = st.sig; double pp_new = st.pp;
        const auto r = ssc_step(P, st.sig, st.pp, Eigen::Vector3d(d, 0.0, 0.0), 0.02);
        s_new = r.sig; pp_new = r.pp;
        st.sig = s_new; st.pp = pp_new; st.eps1 += d;
    }
    const double K0crs = st.sig(1) / st.sig(0);
    std::printf("   CRS segment: sigma_v = %.1f  K0 = %.4f (K0NC = %.4f, err %+.2f%%)\n",
                st.sig(0), K0crs, P.K0nc, 100.0 * (K0crs - P.K0nc) / P.K0nc);
    check(std::fabs(K0crs - P.K0nc) < 0.03 * P.K0nc,
          "constant-rate oedometer develops K0NC within 3% (M(K0NC) behavioural pin, creep flow)");
    check(std::fabs(st.sig(1) - st.sig(2)) < 1e-6 * st.sig(0),
          "the two lateral stresses stay tied on the creep path (corner handling)");
}

void test_fast_triaxial_mc() {
    std::printf("-- (e) fast drained triaxial: failure ON the Mohr-Coulomb line --\n");
    const Params P = clay();
    const double p0 = 100.0;
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = 4.0 * p0;   // overconsolidated so cap creep stays negligible during the fast test
    const int n = 600;
    const double d1 = 0.20 / n, dt = 1e-8;   // ~5 ms total: creep contribution negligible
    double dlat = -d1 * 0.4;
    for (int i = 0; i < n; ++i) {
        Eigen::Vector3d s_new = sig;
        double pp_new = pp;
        for (int it = 0; it < 60; ++it) {
            const auto r = ssc_step(P, sig, pp, Eigen::Vector3d(d1, dlat, dlat), dt);
            const double res = r.sig(2) - p0;
            s_new = r.sig; pp_new = r.pp;
            if (std::fabs(res) < 1e-7 * p0) break;
            const double h = 1e-7;
            const auto r2 = ssc_step(P, sig, pp, Eigen::Vector3d(d1, dlat + h, dlat + h), dt);
            const double drds = (r2.sig(2) - r.sig(2)) / h;
            if (!std::isfinite(drds) || std::fabs(drds) < 1e-30) { dlat = -0.5 * d1; continue; }
            double step = -res / drds;
            step = std::clamp(step, -5.0 * d1, 5.0 * d1);
            dlat += step;
            if (!std::isfinite(dlat)) dlat = -0.5 * d1;
        }
        sig = s_new; pp = pp_new;
    }
    const double q = sig(0) - sig(2);
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi);
    const double q_exact = (p0 * (1 + sphi) + 2 * P.c * cphi) / (1 - sphi) - p0;
    std::printf("   q = %.3f (MC closed form %.3f, err %+.2f%%)   sigma3 = %.3f\n", q, q_exact,
                100.0 * (q - q_exact) / q_exact, sig(2));
    check(std::fabs(sig(2) - p0) < 0.01 * p0, "the mixed control held the cell pressure");
    check(std::fabs(q - q_exact) < 0.005 * q_exact,
          "fast triaxial failure sits ON the Mohr-Coulomb line (< 0.5%)");
}

// ================================ STAGE 2: FE wiring ================================
namespace fe {
using namespace katai::core;

MaterialModel ssc_model() {
    MaterialModel m;
    m.type = MaterialType::SoftSoilCreep;
    m.ssc = clay();
    return m;
}

// Kernel oracle: n_c oedometric compression steps (d, dt_c) from an NC start, then n_h timed
// holds (0, dt_h) -- the FE paths below must reproduce this chain exactly.
struct KChain { Eigen::Vector3d sig; double pp; };
KChain kernel_chain(const Params& P, double p0, double d, double dt_c, int n_c, double dt_h,
                    int n_h) {
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0;
    for (int i = 0; i < n_c; ++i) {
        const auto r = ssc_step(P, sig, pp, Eigen::Vector3d(d, 0.0, 0.0), dt_c);
        sig = r.sig; pp = r.pp;
    }
    for (int i = 0; i < n_h; ++i) {
        const auto r = ssc_step(P, sig, pp, Eigen::Vector3d::Zero(), dt_h);
        sig = r.sig; pp = r.pp;
    }
    return {sig, pp};
}

double rel(double a, double b) { return std::fabs(a - b) / std::max(1.0, std::fabs(b)); }

void test_g1_integrate_point_identity() {
    std::printf("-- (g1) plane-strain integrate_point (strain + time) == kernel chain --\n");
    const MaterialModel m = ssc_model();
    const double p0 = 100.0, d = 1e-3, dt_c = 0.01, dt_h = 0.5;
    const int n_c = 40, n_h = 40;
    GaussState gs;
    gs.stress = Eigen::Vector3d(-p0, -p0, 0.0); gs.stress_zz = -p0; gs.pp = p0;
    Eigen::Matrix3d Dt;
    for (int i = 0; i < n_c; ++i) {
        GaussState tr;
        integrate_point(m, gs, Eigen::Vector3d(0.0, -d, 0.0), tr, Dt, TangentMode::kConsistent,
                        dt_c);
        gs = tr;
    }
    for (int i = 0; i < n_h; ++i) {
        GaussState tr;
        integrate_point(m, gs, Eigen::Vector3d::Zero(), tr, Dt, TangentMode::kConsistent, dt_h);
        gs = tr;
    }
    const KChain k = kernel_chain(clay(), p0, d, dt_c, n_c, dt_h, n_h);
    std::printf("   FE: syy = %.6f  sxx = %.6f  szz = %.6f  pp = %.4f\n", -gs.stress(1),
                -gs.stress(0), -gs.stress_zz, gs.pp);
    std::printf("   kernel: s1 = %.6f  slat = %.6f  pp = %.4f\n", k.sig(0), k.sig(1), k.pp);
    check(rel(-gs.stress(1), k.sig(0)) < 1e-8, "FE vertical stress == kernel major (1e-8)");
    check(rel(-gs.stress(0), k.sig(1)) < 1e-8, "FE lateral stress == kernel lateral (1e-8)");
    check(rel(-gs.stress_zz, k.sig(2)) < 1e-8, "FE out-of-plane == kernel lateral (1e-8)");
    check(rel(gs.pp, k.pp) < 1e-8, "FE pp == kernel pp (1e-8)");
}

void test_g2_axisym_identity() {
    std::printf("-- (g2) axisymmetric integrate_point (strain + time) == kernel chain --\n");
    const MaterialModel m = ssc_model();
    const double p0 = 100.0, d = 1e-3, dt_c = 0.01, dt_h = 0.5;
    const int n_c = 40, n_h = 40;
    GaussState gs;
    gs.stress = Eigen::Vector3d(-p0, -p0, 0.0); gs.stress_zz = -p0; gs.pp = p0;
    Eigen::Matrix4d Dt;
    for (int i = 0; i < n_c; ++i) {
        GaussState tr;
        integrate_point_axisym(m, gs, Eigen::Vector4d(0.0, -d, 0.0, 0.0), tr, Dt,
                               TangentMode::kConsistent, dt_c);
        gs = tr;
    }
    for (int i = 0; i < n_h; ++i) {
        GaussState tr;
        integrate_point_axisym(m, gs, Eigen::Vector4d::Zero(), tr, Dt, TangentMode::kConsistent,
                               dt_h);
        gs = tr;
    }
    const KChain k = kernel_chain(clay(), p0, d, dt_c, n_c, dt_h, n_h);
    check(rel(-gs.stress(1), k.sig(0)) < 1e-8, "axisym axial == kernel major (1e-8)");
    check(rel(-gs.stress(0), k.sig(1)) < 1e-8, "axisym radial == kernel lateral (1e-8)");
    check(rel(-gs.stress_zz, k.sig(2)) < 1e-8, "axisym HOOP == kernel lateral (1e-8)");
    check(rel(gs.pp, k.pp) < 1e-8, "axisym pp == kernel pp (1e-8)");
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& kmat, const Eigen::VectorXd& r) {
    auto s = katai::linsolve::make_direct_solver(katai::linsolve::MatrixType::RealNonsymmetric);
    s->factorize(kmat);
    return s->solve(r);
}

void test_g3_relaxation_bvp() {
    std::printf("-- (g3) tri15 column BVP: timed compression, then held-displacement "
                "relaxation --\n");
    const double p0 = 100.0, H = 1.0, W = 0.5;
    const double u_step = 0.004, dt_c = 0.01;   // 10 x 0.004 m over 0.1 day
    const double dt_h = 0.5;                    // then 10 x 0.5 day holds
    const int n_c = 10, n_h = 10;
    using katai::geometry::RectangularDomain;
    const RectangularDomain dom{0.0, 0.0, W, H, 0};
    katai::mesh::Mesh mesh = katai::mesh::generate_structured_tri15(dom, 1, 2);

    DofMap dofs(mesh.node_count, 2);
    std::vector<double> presc(2 * mesh.node_count, 0.0);
    for (int nn : mesh.bottom_nodes) dofs.fix_node_component(nn, 1);
    for (int nn : mesh.left_nodes) dofs.fix_node_component(nn, 0);
    for (int nn : mesh.right_nodes) dofs.fix_node_component(nn, 0);
    for (int nn : mesh.top_nodes) dofs.fix_node_component(nn, 1);
    dofs.finalize();
    for (int nn : mesh.top_nodes) presc[dofs.global_dof(nn, 1)] = -u_step;

    const std::vector<MaterialModel> mm = {ssc_model()};
    const size_t ng = (size_t)mesh.element_count * Tri15Element::kGaussCount;
    std::vector<GaussState> state(ng);
    for (auto& g : state) { g.stress = Eigen::Vector3d(-p0, -p0, 0.0); g.stress_zz = -p0; g.pp = p0; }

    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), (Eigen::Index)presc.size());
    const Eigen::VectorXd pv_zero = Eigen::VectorXd::Zero((Eigen::Index)presc.size());
    katai::core::NewtonOptions opt;
    opt.load_steps = 1; opt.max_iterations = 40; opt.tolerance = 1e-9;
    bool all_conv = true;
    for (int call = 0; call < n_c; ++call) {
        opt.time_interval = dt_c;
        const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, opt, state,
                                                    {}, {}, pv);
        all_conv = all_conv && r.converged;
        state = r.gauss_states;
    }
    for (int call = 0; call < n_h; ++call) {
        opt.time_interval = dt_h;
        const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f0, solve_unsym, opt, state,
                                                    {}, {}, pv_zero);
        all_conv = all_conv && r.converged;
        state = r.gauss_states;
    }
    check(all_conv, "all timed compression + hold solves converged (Newton on the FD tangent)");

    const KChain k = kernel_chain(clay(), p0, u_step / H, dt_c, n_c, dt_h, n_h);
    double dev_syy = 0.0, dev_sxx = 0.0, dev_pp = 0.0;
    for (const auto& g : state) {
        dev_syy = std::max(dev_syy, rel(-g.stress(1), k.sig(0)));
        dev_sxx = std::max(dev_sxx, rel(-g.stress(0), k.sig(1)));
        dev_pp = std::max(dev_pp, rel(g.pp, k.pp));
    }
    std::printf("   Gauss syy = %.4f (kernel %.4f)  max rel dev: syy %.2e  sxx %.2e  pp %.2e\n",
                -state[0].stress(1), k.sig(0), dev_syy, dev_sxx, dev_pp);
    check(dev_syy < 1e-6, "every Gauss vertical stress == kernel chain (1e-6)");
    check(dev_sxx < 1e-6, "every Gauss lateral stress == kernel chain (1e-6)");
    check(dev_pp < 1e-6, "every Gauss pp == kernel chain (1e-6)");
    // And the holds genuinely RELAXED the column (the time bookkeeping did real work): compare
    // against the same chain WITHOUT the hold phases.
    const KChain k_nohold = kernel_chain(clay(), p0, u_step / H, dt_c, n_c, 0.0, 0);
    std::printf("   relaxation over %g d: sigma_v %.4f -> %.4f\n", n_h * dt_h, k_nohold.sig(0),
                k.sig(0));
    check(-state[0].stress(1) < 0.99 * k_nohold.sig(0),
          "held-displacement column relaxed by > 1% (time actually flowed through the solver)");
}

}  // namespace fe

}  // namespace

int main() {
    std::printf("SOFT SOIL CREEP (PLAXIS MMM sec 11) -- Stage 1 material-point closed-form V&V\n\n");
    test_constant_stress_creep();
    std::printf("\n");
    test_ocr_rate_ratio();
    std::printf("\n");
    test_relaxation();
    std::printf("\n");
    test_staged_oedometer();
    std::printf("\n");
    test_fast_triaxial_mc();
    std::printf("\nSTAGE 2 -- FE constitutive wiring (time through integrate_point)\n\n");
    fe::test_g1_integrate_point_identity();
    std::printf("\n");
    fe::test_g2_axisym_identity();
    std::printf("\n");
    fe::test_g3_relaxation_bvp();
    if (g_failures == 0) {
        std::printf("\nOK: Buisman/Garlanger secondary slope, OCR^-beta rate, closed-form "
                    "relaxation, the 1-day NC line, MC failure, and the FULL FE wiring "
                    "(plane strain + axisym + timed BVP) all reproduce their defining forms\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
