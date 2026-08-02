// SOFT SOIL model -- Stage 1 material-point V&V (PLAXIS MMM sec 10; locked formulation in
// docs/references/soft-soil-formulation.md). Oracles are the model's own DEFINING closed forms,
// none of which share code with the kernel's return mapping:
//   (a) isotropic VIRGIN loading follows eps_v = lam* ln(p/p0) (the ln-law that defines lam*);
//   (b) unloading follows the kap* line and does NOT move p_p (memory preserved);
//   (c) reloading past p_p REJOINS the virgin line (pre-consolidation memory);
//   (d) on the virgin isotropic line p == p_p identically (the cap top rides the stress);
//   (e) OEDOMETER primary loading develops sigma_h/sigma_v == K0NC -- the BEHAVIOURAL pin of the
//       M(K0NC) formula (MMM: "M can be chosen such that a known value of K0nc is matched");
//   (f) drained TRIAXIAL compression fails exactly on the Mohr-Coulomb line (M is NOT the failure
//       criterion in Soft Soil -- the manual is explicit about that).
// STAGE 2 -- FE constitutive wiring (integrate_point / ss_return_core). The kernel is pinned by
// (a)-(f); these pin the WIRING around it (Voigt trial, spectral decomposition, rank pairing,
// coaxial reconstruction, FD tangent, assembler/solver):
//   (g1) plane-strain integrate_point oedometer chain == the principal kernel chain (identity)
//        and develops K0NC through the FULL FE path;
//   (g2) frame indifference: the same path rotated 30 deg in-plane, stresses rotated back == (g1);
//   (g3) axisymmetric integrate_point (hoop = real strain) == the same principal kernel chain;
//   (g4) FE BVP: tri15 confined column (solve_nonlinear, prescribed top settlement, 10 committed
//        steps) reproduces the SAME kernel chain at every Gauss point -- assembler + Newton +
//        FD-tangent wiring, with the NC start immediately on the cap (plastic from step 1).
#include <katai/materials/soft_soil.hpp>

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

using namespace katai::core::softsoil;

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
    P.nu_ur = 0.15;
    P.c = 1.0;
    P.phi = 25.0 * kPi / 180.0;
    P.psi = 0.0;
    P.K0nc = 1.0 - std::sin(P.phi);   // Jaky
    return P;
}

// Drive an isotropic strain path in n steps of volumetric increment dv (positive = compression).
struct IsoRun { Eigen::Vector3d sig; double pp, eps_v; };
IsoRun iso_drive(const Params& P, Eigen::Vector3d sig, double pp, double dv_total, int n) {
    double eps_v = 0.0;
    for (int i = 0; i < n; ++i) {
        const Eigen::Vector3d de = Eigen::Vector3d::Constant(dv_total / n / 3.0);
        const auto r = ss_step(P, sig, pp, de);
        sig = r.sig; pp = r.pp; eps_v += dv_total / n;
    }
    return {sig, pp, eps_v};
}

void test_isotropic_lnlaw() {
    std::printf("-- (a,b,c,d) isotropic ln-law: virgin lam*, unload kap*, memory --\n");
    const Params P = clay();
    const double p0 = 50.0;
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0;   // normally consolidated start (on the cap top)

    // (a) virgin: compress to eps_v = lam* ln(4) -> p should reach 4 p0.
    const double ev1 = P.lam_star * std::log(4.0);
    auto r1 = iso_drive(P, sig, pp, ev1, 400);
    const double p1 = r1.sig.mean();
    std::printf("   virgin: p = %.4f (target %.4f)   rel err %.2e\n", p1, 4 * p0,
                std::fabs(p1 - 4 * p0) / (4 * p0));
    check(std::fabs(p1 - 4 * p0) < 1e-6 * 4 * p0,
          "virgin isotropic loading follows eps_v = lam* ln(p/p0) (exponential-mean return: exact)");
    // (d) on the virgin line the cap top rides the stress: p == p_p.
    check(std::fabs(r1.pp - p1) < 0.005 * p1, "on the virgin line p == p_p (the cap top rides)");

    // (b) unload to half the strain of one decade: slope kap*, p_p frozen.
    const double pp_before = r1.pp;
    const double ev2 = -P.kap_star * std::log(2.0);   // elastic: p should halve
    auto r2 = iso_drive(P, r1.sig, r1.pp, ev2, 200);
    const double p2 = r2.sig.mean();
    std::printf("   unload: p = %.4f (target %.4f)   pp = %.4f (was %.4f)\n", p2, p1 / 2, r2.pp,
                pp_before);
    check(std::fabs(p2 - p1 / 2) < 0.002 * p1, "unloading follows the kap* line (exact exponential)");
    check(r2.pp == pp_before, "unloading does NOT move p_p (memory preserved, bitwise)");

    // (c) reload past p_p: total eps_v at p = 8 p0 equals the VIRGIN value (path independence of
    // the ln-law memory): eps_v = lam* ln(8) + (elastic loop closes exactly).
    const double ev3 = P.kap_star * std::log(2.0) + P.lam_star * std::log(2.0);
    auto r3 = iso_drive(P, r2.sig, r2.pp, ev3, 400);
    const double p3 = r3.sig.mean();
    std::printf("   reload: p = %.4f (virgin target %.4f)\n", p3, 8 * p0);
    check(std::fabs(p3 - 8 * p0) < 0.01 * 8 * p0,
          "reloading rejoins the virgin line (pre-consolidation memory, < 1%)");
}

void test_oedometer_k0nc() {
    std::printf("-- (e) oedometer primary loading -> sigma_h/sigma_v == K0NC (M-formula pin) --\n");
    const Params P = clay();
    const double p0 = 50.0;
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0;
    // 1D compression: eps = (d, 0, 0), well past the initial state so the stress path converges
    // to the K0NC line.
    const int n = 800;
    const double d = 0.12 / n;
    for (int i = 0; i < n; ++i) {
        const auto r = ss_step(P, sig, pp, Eigen::Vector3d(d, 0.0, 0.0));
        sig = r.sig; pp = r.pp;
    }
    const double K0 = sig(1) / sig(0);   // lateral / axial
    std::printf("   sigma_v = %.2f, sigma_h = %.2f -> K0 = %.4f (target K0NC = %.4f, err %+.2f%%)\n",
                sig(0), sig(1), K0, P.K0nc, 100 * (K0 - P.K0nc) / P.K0nc);
    check(sig(0) > 4 * p0, "the oedometer really loaded far into primary compression (teeth)");
    check(std::fabs(K0 - P.K0nc) < 0.02 * P.K0nc,
          "primary 1D compression develops K0NC within 2% (M(K0NC) behavioural pin)");
}

void test_triaxial_mc_failure() {
    std::printf("-- (f) drained triaxial: failure ON the Mohr-Coulomb line (M is not failure) --\n");
    const Params P = clay();
    const double p0 = 100.0;   // cell pressure
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0;
    // Mixed control: drive eps1, iterate eps2 = eps3 so sigma3 stays at the cell pressure.
    const int n = 600;
    const double d1 = 0.20 / n;
    double dlat = -d1 * 0.4;   // WARM-STARTED across steps: near failure the flow is nearly
                               // isochoric (psi = 0 -> eps_lat ~ -eps1/2) and a cold guess makes
                               // the secant lose the cell pressure exactly when MC activates.
    for (int i = 0; i < n; ++i) {
        Eigen::Vector3d s_new = sig;
        double pp_new = pp;
        for (int it = 0; it < 60; ++it) {
            const auto r = ss_step(P, sig, pp, Eigen::Vector3d(d1, dlat, dlat));
            const double res = r.sig(2) - p0;
            s_new = r.sig; pp_new = r.pp;
            if (std::fabs(res) < 1e-7 * p0) break;
            // FD step must sit ABOVE the kernel's own numerical-gradient noise floor (its internal
            // stress perturbations are ~1e-6*scale): 1e-7 strain moves sigma3 by ~1e-3 kPa, well
            // resolved; 1e-9 was inside the noise and produced NaN.
            const double h = 1e-7;
            const auto r2 = ss_step(P, sig, pp, Eigen::Vector3d(d1, dlat + h, dlat + h));
            const double drds = (r2.sig(2) - r.sig(2)) / h;
            if (!std::isfinite(drds) || std::fabs(drds) < 1e-30) { dlat = -0.5 * d1; continue; }
            double step = -res / drds;
            step = std::clamp(step, -5.0 * d1, 5.0 * d1);   // keep the secant on the rails
            dlat += step;
            if (!std::isfinite(dlat)) dlat = -0.5 * d1;
        }
        sig = s_new; pp = pp_new;
    }
    const double q = sig(0) - sig(2);
    const double sphi = std::sin(P.phi), cphi = std::cos(P.phi);
    const double s1_f = (p0 * (1 + sphi) + 2 * P.c * cphi) / (1 - sphi);
    const double q_exact = s1_f - p0;
    std::printf("   q = %.3f (MC closed form %.3f, err %+.2f%%)   sigma3 = %.3f (cell %.1f)\n",
                q, q_exact, 100 * (q - q_exact) / q_exact, sig(2), p0);
    check(std::fabs(sig(2) - p0) < 0.01 * p0, "the mixed control held the cell pressure");
    check(std::fabs(q - q_exact) < 0.002 * q_exact,
          "triaxial failure sits ON the Mohr-Coulomb line (< 0.2%)");
}

// ================================ STAGE 2: FE wiring ================================
namespace fe {
using namespace katai::core;

MaterialModel ss_model() {
    MaterialModel m;
    m.type = MaterialType::SoftSoil;
    m.ssoil = clay();
    return m;
}

// Principal kernel oedometer chain (compression-positive), n steps of d: the ORACLE for every
// FE path below (itself pinned by (a)-(f)).
struct KChain { Eigen::Vector3d sig; double pp; };
KChain kernel_oedometer(const Params& P, double p0, double d, int n) {
    Eigen::Vector3d sig = Eigen::Vector3d::Constant(p0);
    double pp = p0;
    for (int i = 0; i < n; ++i) {
        const auto r = ss_step(P, sig, pp, Eigen::Vector3d(d, 0.0, 0.0));
        sig = r.sig; pp = r.pp;
    }
    return {sig, pp};
}

double rel(double a, double b) { return std::fabs(a - b) / std::max(1.0, std::fabs(b)); }

void test_g1_integrate_point_identity() {
    std::printf("-- (g1) plane-strain integrate_point oedometer == principal kernel chain --\n");
    const MaterialModel m = ss_model();
    const double p0 = 50.0;
    const int n = 200;
    const double d = 5e-4;   // total eps_yy = 0.10, deep into primary compression
    GaussState gs;
    gs.stress = Eigen::Vector3d(-p0, -p0, 0.0); gs.stress_zz = -p0; gs.pp = p0;
    Eigen::Matrix3d Dt;
    for (int i = 0; i < n; ++i) {
        GaussState tr;
        integrate_point(m, gs, Eigen::Vector3d(0.0, -d, 0.0), tr, Dt);
        gs = tr;
    }
    const KChain k = kernel_oedometer(m.ssoil, p0, d, n);
    std::printf("   FE: syy = %.6f  sxx = %.6f  szz = %.6f  pp = %.4f\n", -gs.stress(1),
                -gs.stress(0), -gs.stress_zz, gs.pp);
    std::printf("   kernel: s1 = %.6f  slat = %.6f  pp = %.4f\n", k.sig(0), k.sig(1), k.pp);
    check(rel(-gs.stress(1), k.sig(0)) < 1e-8, "FE vertical stress == kernel major (1e-8)");
    check(rel(-gs.stress(0), k.sig(1)) < 1e-8, "FE lateral stress == kernel lateral (1e-8)");
    check(rel(-gs.stress_zz, k.sig(2)) < 1e-8, "FE out-of-plane == kernel lateral (1e-8)");
    check(rel(gs.pp, k.pp) < 1e-8, "FE pp == kernel pp (1e-8)");
    const double K0 = gs.stress(0) / gs.stress(1);
    std::printf("   K0 through FE path = %.4f (K0NC = %.4f)\n", K0, m.ssoil.K0nc);
    check(std::fabs(K0 - m.ssoil.K0nc) < 0.02 * m.ssoil.K0nc,
          "K0NC develops through the FULL FE path (2%)");
}

void test_g2_frame_indifference() {
    std::printf("-- (g2) frame indifference: 30-deg rotated path, stresses rotated back --\n");
    const MaterialModel m = ss_model();
    const double p0 = 50.0, th = 30.0 * kPi / 180.0;
    const double c = std::cos(th), s = std::sin(th);
    const int n = 150;
    const double d = 6e-4;
    GaussState ga, gb;
    ga.stress = Eigen::Vector3d(-p0, -p0, 0.0); ga.stress_zz = -p0; ga.pp = p0;
    gb = ga;   // isotropic start is rotation-invariant
    // eps' = Q eps Q^T for eps = diag(0, -d): exx' = -d s^2, eyy' = -d c^2, gxy' = 2 d s c.
    const Eigen::Vector3d dea(0.0, -d, 0.0);
    const Eigen::Vector3d deb(-d * s * s, -d * c * c, 2.0 * d * s * c);
    Eigen::Matrix3d Dt;
    for (int i = 0; i < n; ++i) {
        GaussState tr;
        integrate_point(m, ga, dea, tr, Dt); ga = tr;
        integrate_point(m, gb, deb, tr, Dt); gb = tr;
    }
    // Rotate B's stress back: S = Q^T S' Q.
    Eigen::Matrix2d Q, Sp;
    Q << c, -s, s, c;
    Sp << gb.stress(0), gb.stress(2), gb.stress(2), gb.stress(1);
    const Eigen::Matrix2d Sb = Q.transpose() * Sp * Q;
    std::printf("   back-rotated: sxx %.6f (ref %.6f)  syy %.6f (ref %.6f)  sxy %.2e\n",
                Sb(0, 0), ga.stress(0), Sb(1, 1), ga.stress(1), Sb(0, 1));
    check(rel(Sb(0, 0), ga.stress(0)) < 1e-6, "rotated-back sxx matches (1e-6)");
    check(rel(Sb(1, 1), ga.stress(1)) < 1e-6, "rotated-back syy matches (1e-6)");
    check(std::fabs(Sb(0, 1) - ga.stress(2)) < 1e-6 * p0, "rotated-back sxy matches (1e-6 p0)");
    check(rel(gb.stress_zz, ga.stress_zz) < 1e-6, "sigma_zz frame-invariant (1e-6)");
    check(rel(gb.pp, ga.pp) < 1e-6, "pp frame-invariant (1e-6)");
}

void test_g3_axisym_identity() {
    std::printf("-- (g3) axisymmetric integrate_point oedometer == principal kernel chain --\n");
    const MaterialModel m = ss_model();
    const double p0 = 50.0;
    const int n = 200;
    const double d = 5e-4;
    GaussState gs;
    gs.stress = Eigen::Vector3d(-p0, -p0, 0.0); gs.stress_zz = -p0; gs.pp = p0;   // [r,z,rz]+theta
    Eigen::Matrix4d Dt;
    for (int i = 0; i < n; ++i) {
        GaussState tr;
        integrate_point_axisym(m, gs, Eigen::Vector4d(0.0, -d, 0.0, 0.0), tr, Dt);
        gs = tr;
    }
    const KChain k = kernel_oedometer(m.ssoil, p0, d, n);
    check(rel(-gs.stress(1), k.sig(0)) < 1e-8, "axisym axial stress == kernel major (1e-8)");
    check(rel(-gs.stress(0), k.sig(1)) < 1e-8, "axisym radial stress == kernel lateral (1e-8)");
    check(rel(-gs.stress_zz, k.sig(2)) < 1e-8, "axisym HOOP stress == kernel lateral (1e-8)");
    check(rel(gs.pp, k.pp) < 1e-8, "axisym pp == kernel pp (1e-8)");
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& kmat, const Eigen::VectorXd& r) {
    auto s = katai::linsolve::make_direct_solver(katai::linsolve::MatrixType::RealNonsymmetric);
    s->factorize(kmat);
    return s->solve(r);
}

void test_g4_fe_oedometer_bvp() {
    std::printf("-- (g4) FE BVP: tri15 confined column vs the same kernel chain --\n");
    const double p0 = 50.0, H = 1.0, W = 0.5, u_step = 0.01;   // 10 x 0.01 -> eps_yy = 0.10
    const int ncalls = 10;
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

    const std::vector<MaterialModel> mm = {ss_model()};
    const size_t ng = (size_t)mesh.element_count * Tri15Element::kGaussCount;
    std::vector<GaussState> state(ng);
    for (auto& g : state) { g.stress = Eigen::Vector3d(-p0, -p0, 0.0); g.stress_zz = -p0; g.pp = p0; }

    const Eigen::VectorXd f0 = Eigen::VectorXd::Zero(dofs.equation_count());
    Eigen::Map<const Eigen::VectorXd> pv(presc.data(), (Eigen::Index)presc.size());
    bool all_conv = true;
    // 10 sequential single-step solves with the state carried forward: the committed step sizes
    // are then EXACTLY uniform (no adaptive sub-stepping ambiguity) and the kernel chain with
    // d = u_step/H per step is the exact oracle. NC start -> plastic (cap) from the FIRST step,
    // so the FD tangent path is what Newton actually runs on.
    for (int call = 0; call < ncalls; ++call) {
        const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f0, solve_unsym,
                                                    {1, 40, 1e-9}, state, {}, {}, pv);
        all_conv = all_conv && r.converged;
        state = r.gauss_states;
    }
    check(all_conv, "all 10 prescribed-settlement steps converged (Newton on the FD tangent)");

    const KChain k = kernel_oedometer(clay(), p0, u_step / H, ncalls);
    double dev_syy = 0.0, dev_sxx = 0.0, dev_pp = 0.0, uni = 0.0;
    for (const auto& g : state) {
        dev_syy = std::max(dev_syy, rel(-g.stress(1), k.sig(0)));
        dev_sxx = std::max(dev_sxx, rel(-g.stress(0), k.sig(1)));
        dev_pp = std::max(dev_pp, rel(g.pp, k.pp));
        uni = std::max(uni, std::fabs(g.stress(1) - state[0].stress(1)));
    }
    std::printf("   Gauss syy = %.4f (kernel %.4f)  max rel dev: syy %.2e  sxx %.2e  pp %.2e\n",
                -state[0].stress(1), k.sig(0), dev_syy, dev_sxx, dev_pp);
    std::printf("   uniformity max|syy - syy0| = %.2e\n", uni);
    check(dev_syy < 1e-6, "every Gauss vertical stress == kernel chain (1e-6)");
    check(dev_sxx < 1e-6, "every Gauss lateral stress == kernel chain (1e-6)");
    check(dev_pp < 1e-6, "every Gauss pp == kernel chain (1e-6)");
    check(uni < 1e-8 * p0, "stress field uniform across the column (wiring sanity)");
    const double K0 = state[0].stress(0) / state[0].stress(1);
    std::printf("   K0 through the FULL BVP = %.4f (K0NC = %.4f)\n", K0, clay().K0nc);
    check(std::fabs(K0 - clay().K0nc) < 0.02 * clay().K0nc, "K0NC develops in the BVP (2%)");
}

}  // namespace fe

}  // namespace

int main() {
    std::printf("SOFT SOIL (PLAXIS MMM sec 10) -- Stage 1 material-point closed-form V&V\n\n");
    test_isotropic_lnlaw();
    std::printf("\n");
    test_oedometer_k0nc();
    std::printf("\n");
    test_triaxial_mc_failure();
    std::printf("\nSTAGE 2 -- FE constitutive wiring\n\n");
    fe::test_g1_integrate_point_identity();
    std::printf("\n");
    fe::test_g2_frame_indifference();
    std::printf("\n");
    fe::test_g3_axisym_identity();
    std::printf("\n");
    fe::test_g4_fe_oedometer_bvp();
    if (g_failures == 0) {
        std::printf("\nOK: ln-law (lam*/kap*), pre-consolidation memory, K0NC via M, MC failure, "
                    "and the FULL FE wiring (plane strain + axisym + rotated frame + BVP) all "
                    "reproduce the kernel/closed forms\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
