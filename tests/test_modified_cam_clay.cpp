// MODIFIED CAM-CLAY: the material point against the model's own closed forms, the finite-element
// wrapper against the material point, and a confined column through the whole program.
//
// The model (Roscoe & Burland 1968) is the Soft Soil core read with the von Mises q, the
// critical-state M = 6 sin(phi_cs)/(3 - sin(phi_cs)), no cohesion and no Mohr-Coulomb surface
// (soft_soil.hpp, Params::cam_clay). None of the oracles below shares code with the return
// mapping: they are the model's equations, solved in closed form or integrated along the stress
// variable p with a step a thousand times finer than the test's.
//
// verify: KV-CST-023
//   oracle:   closed_form
//   source:   Roscoe, K.H. & Burland, J.B. (1968). On the generalised stress-strain behaviour of wet clay. In Engineering Plasticity, Cambridge University Press, 535-609 -- the modified Cam-Clay ellipse, its associated flow and its hardening; Wood, D.M. (1990). Soil Behaviour and Critical State Soil Mechanics. Cambridge University Press -- the undrained strength su = (M/2) p0 (OCR/2)^Lambda, Lambda = (lambda - kappa)/lambda, and the triaxial paths
//   locator:  lambda* = 0.10, kappa* = 0.02, nu = 0.2, phi_cs = 25 deg (M = 0.984), p0 = 100 kPa, isotropically consolidated: (a) undrained triaxial compression at OCR 1, 1.5 and 8 -- on the plastic part of the path p_p = p_c0 (p0/p)^(kappa*/(lambda* - kappa*)) and q = M sqrt(p (p_p - p)), and at 40% axial strain the state is at the critical state p_f = p0 (OCR/2)^Lambda, q_f = M p_f; (b) the same paths' axial strain against the integral of the elastic and plastic shear strain over p; (c) drained triaxial compression at OCR 1, q = 3 (p - p0), its axial and volumetric strain against the same integral; (d) the axisymmetric finite-element integration of (a) against the material point; (e) a confined column of the material through the whole program, K0 procedure on the model's own K0 and a surcharge: the ratio of the horizontal to the vertical effective stress increment against the closed-form K0 of one-dimensional compression
//   quantity: (a) p_p and q on the path, p and q at the end [kPa]; (b), (c) axial and volumetric strain [-]; (d) the stresses and p_p [kPa]; (e) the material point's K0 at two step sizes, and d sigma'_h / d sigma'_v in the column [-]
//   expected: (a) the closed forms; (b), (c) the integrals; (d) the material point; (e) K0 = (3 - eta)/(3 + 2 eta), eta the root of eta kappa*/(4.5 r) + 2 eta (lambda* - kappa*)/(M^2 - eta^2) = (2/3) lambda*, r = (1 - 2 nu)/(1 + nu) -- 0.72331
//   band:     (a) 1e-9 relative for the volume constraint (measured 3e-14), 1e-6 for q on the ellipse (2e-9), 1e-4 for the end state against p_f and M p_f (it is reached: 57.435 = 57.435 kPa at OCR 1); (b), (c) 0.5% (measured 0.23%, 0.13% and 0.14%, at the start of the plastic path, where the strain is smallest; the return is backward Euler, first order in the step); (d) 1e-8; (e) 5e-4 for the material point at both step sizes (measured -4e-5 at 1e-4 and -1.9e-4 at 1e-3: what is left is the elastic shear modulus held at the trial's p over the step) and 0.2% for the column (measured 0.08%). The return along the TRIAL's normal, which Soft Soil uses, put the column at 0.74252 (+2.66%) and stopped the undrained path short of the critical state; the normal of the returned state removes both
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>
#include <katai/materials/soft_soil.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using namespace katai::core;
namespace m = katai::model;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

constexpr double kLam = 0.10, kKap = 0.02, kNu = 0.2, kPhiDeg = 25.0, kP0 = 100.0;

softsoil::Params mcc() {
    softsoil::Params P;
    P.cam_clay = true;
    P.lam_star = kLam;
    P.kap_star = kKap;
    P.nu_ur = kNu;
    P.phi = kPhiDeg * kPi / 180.0;
    return P;
}
double M_cs() { const double s = std::sin(kPhiDeg * kPi / 180.0); return 6.0 * s / (3.0 - s); }
double G_at(double p) { return 1.5 * (p / kKap) * (1.0 - 2.0 * kNu) / (1.0 + kNu); }
double q_of(const Eigen::Vector3d& s) {
    const Eigen::Vector3d d = s - Eigen::Vector3d::Constant(s.mean());
    return std::sqrt(1.5 * d.squaredNorm());
}

// ---- the material point, driven in triaxial compression (compression-positive principals) ----
struct Pt { double ea, ev, p, q, pp; };
std::vector<Pt> undrained(double ocr, double ea_end, int n) {
    const softsoil::Params P = mcc();
    Eigen::Vector3d s = Eigen::Vector3d::Constant(kP0);
    double pp = ocr * kP0;
    std::vector<Pt> out{{0.0, 0.0, kP0, 0.0, pp}};
    const double d = ea_end / n;
    for (int i = 0; i < n; ++i) {
        const auto r = softsoil::ss_step(P, s, pp, Eigen::Vector3d(d, -0.5 * d, -0.5 * d));
        if (r.unstable) { std::printf("      (unstable at step %d)\n", i); return {}; }
        s = r.sig; pp = r.pp;
        out.push_back({d * (i + 1), 0.0, s.mean(), q_of(s), pp});
    }
    return out;
}
// Drained: the lateral strain of each step is found so the cell pressure stays at p0.
std::vector<Pt> drained(double ea_end, int n) {
    const softsoil::Params P = mcc();
    Eigen::Vector3d s = Eigen::Vector3d::Constant(kP0);
    double pp = kP0, ev = 0.0;
    std::vector<Pt> out{{0.0, 0.0, kP0, 0.0, pp}};
    const double d = ea_end / n;
    double e = -0.5 * d;
    for (int i = 0; i < n; ++i) {
        auto lateral = [&](double el) {
            return softsoil::ss_step(P, s, pp, Eigen::Vector3d(d, el, el)).sig(1) - kP0;
        };
        double g = lateral(e);
        for (int it = 0; it < 50 && std::fabs(g) > 1e-10 * kP0; ++it) {
            const double h = 1e-7 * std::max(std::fabs(e), d);
            const double dg = (lateral(e + h) - g) / h;
            e -= g / dg;
            g = lateral(e);
        }
        const auto r = softsoil::ss_step(P, s, pp, Eigen::Vector3d(d, e, e));
        s = r.sig; pp = r.pp; ev += d + 2.0 * e;
        out.push_back({d * (i + 1), ev, s.mean(), q_of(s), pp});
    }
    return out;
}

// ---- the oracles: the model's equations integrated over p (Simpson, fine) -------------------
// Undrained: eps_v = 0, so p_p(p) and q(p) are closed forms and eps_a = eps_q. The elastic shear
// strain int dq/(3G) is taken by parts -- G is proportional to p, so it is q/(3G) + int q/(3 G p) dp
// -- because at the normally consolidated start, the tip of the ellipse, dq/dp is infinite; the
// elastic part before first yield, q_y/(3 G(p0)), is the boundary term at p0. The plastic shear
// strain is int -kappa*/p * 2 eta/(M^2 - eta^2) dp, from p0 to p.
double undrained_ea_at(double ocr, double p_to) {
    const double M = M_cs(), pc0 = ocr * kP0, ex = kKap / (kLam - kKap);
    auto pc = [&](double p) { return pc0 * std::pow(kP0 / p, ex); };
    auto q = [&](double p) { return M * std::sqrt(std::max(0.0, p * (pc(p) - p))); };
    auto integrand = [&](double p) {
        const double eta = q(p) / p;
        return q(p) / (3.0 * G_at(p) * p) - kKap / p * 2.0 * eta / (M * M - eta * eta);
    };
    const int n = 20000;
    const double a = kP0, b = p_to, h = (b - a) / n;
    double sum = integrand(a) + integrand(b);
    for (int i = 1; i < n; ++i) sum += integrand(a + i * h) * (i % 2 ? 4.0 : 2.0);
    return q(p_to) / (3.0 * G_at(p_to)) + sum * h / 3.0;
}
// Drained (OCR 1): q = 3 (p - p0), on the ellipse p_p = p + q^2/(M^2 p);
// eps_v = kappa* ln(p/p0) + (lambda* - kappa*) ln(p_p/p_p0), eps_q = int dq/(3G) + deps_v^p 2 eta/(M^2 - eta^2),
// eps_a = eps_q + eps_v/3.
void drained_strains_at(double p_to, double& ea, double& ev) {
    const double M = M_cs();
    auto q = [&](double p) { return 3.0 * (p - kP0); };
    auto pc = [&](double p) { return p + q(p) * q(p) / (M * M * p); };
    auto integrand = [&](double p) {
        const double h = 1e-6 * p;
        const double dlnpc = (std::log(pc(p + h)) - std::log(pc(p - h))) / (2.0 * h);
        const double eta = q(p) / p;
        return 3.0 / (3.0 * G_at(p)) + (kLam - kKap) * dlnpc * 2.0 * eta / (M * M - eta * eta);
    };
    const int n = 20000;
    const double a = kP0, b = p_to, h = (b - a) / n;
    double sum = integrand(a) + integrand(b);
    for (int i = 1; i < n; ++i) sum += integrand(a + i * h) * (i % 2 ? 4.0 : 2.0);
    const double eq = sum * h / 3.0;
    ev = kKap * std::log(p_to / kP0) + (kLam - kKap) * std::log(pc(p_to) / kP0);
    ea = eq + ev / 3.0;
}
// One-dimensional compression on the normally consolidated line: the stress ratio eta at which
// the total strain is uniaxial, eps_q / eps_v = 2/3.
double k0_closed_form() {
    const double M = M_cs(), r = (1.0 - 2.0 * kNu) / (1.0 + kNu);
    auto g = [&](double eta) {
        return eta * kKap / (4.5 * r) + 2.0 * eta * (kLam - kKap) / (M * M - eta * eta) - 2.0 / 3.0 * kLam;
    };
    double lo = 1e-9, hi = M * (1.0 - 1e-12);
    for (int i = 0; i < 200; ++i) { const double mid = 0.5 * (lo + hi); (g(mid) > 0.0 ? hi : lo) = mid; }
    const double eta = 0.5 * (lo + hi);
    return (3.0 - eta) / (3.0 + 2.0 * eta);
}

// One-dimensional compression at the material point, from an isotropic state, until the stress
// ratio settles: the model's own K0 of normal consolidation.
double point_k0(double step) {
    const softsoil::Params P = mcc();
    Eigen::Vector3d s = Eigen::Vector3d::Constant(10.0);
    double pp = 10.0;
    const int n = (int)std::lround(2.0 / step);   // to 2.0 of axial strain: some ten decades of p
    for (int i = 0; i < n; ++i) {
        const auto r = softsoil::ss_step(P, s, pp, Eigen::Vector3d(step, 0.0, 0.0));
        s = r.sig; pp = r.pp;
    }
    return s(1) / s(0);
}

void test_undrained() {
    std::printf("-- (a), (b) undrained triaxial compression at OCR 1, 1.5 and 8 --\n");
    const double M = M_cs(), Lam = (kLam - kKap) / kLam, ex = kKap / (kLam - kKap);
    for (double ocr : {1.0, 1.5, 8.0}) {
        const auto path = undrained(ocr, 0.40, 8000);
        check(!path.empty(), "OCR " + std::to_string(ocr).substr(0, 4) + ": the path is integrated without an instability");
        if (path.empty()) continue;
        double worst_pp = 0.0, worst_q = 0.0, worst_ea = 0.0, worst_at = 0.0;
        for (const Pt& s : path) {
            if (s.q * s.q / (M * M * s.p) + s.p - s.pp < -1e-6 * s.pp) continue;   // inside the ellipse: elastic
            worst_pp = std::max(worst_pp, rel(s.pp, ocr * kP0 * std::pow(kP0 / s.p, ex)));
            worst_q = std::max(worst_q, rel(s.q, M * std::sqrt(s.p * (s.pp - s.p))));
            const double eta = s.q / s.p;
            if (eta < 0.9 * M && s.ea > 0.002) {
                const double r = rel(s.ea, undrained_ea_at(ocr, s.p));
                if (r > worst_ea) worst_at = s.ea;
                worst_ea = std::max(worst_ea, r);
            }
        }
        const Pt& f = path.back();
        const double pf = kP0 * std::pow(ocr / 2.0, Lam);
        std::printf("  OCR %.1f: p_p on the volume constraint %.1e, q on the ellipse %.1e, eps_a against the integral %.2e (at %.4f);"
                    " end p %.3f (p_f %.3f), q %.3f (M p_f %.3f), su %.3f kPa\n",
                    ocr, worst_pp, worst_q, worst_ea, worst_at, f.p, pf, f.q, M * pf, 0.5 * f.q);
        const std::string tag = "OCR " + std::to_string(ocr).substr(0, 3) + ": ";
        check(worst_pp < 1e-9, tag + "p_p follows the undrained volume constraint (1e-9)");
        check(worst_q < 1e-6, tag + "q stays on the ellipse (1e-6)");
        // The critical state is an asymptote: the end state is compared with the closed-form
        // path at the p it reached, and its distance from the critical state is printed.
        check(rel(f.p, pf) < 1e-4 && rel(f.q, M * pf) < 1e-4,
              tag + "the path ends at the critical state, su = (M/2) p0 (OCR/2)^Lambda (1e-4)");
        check(worst_ea < 5e-3, tag + "the axial strain is the integral of the model's shear strain (0.5%)");
    }
}

void test_drained() {
    std::printf("-- (c) drained triaxial compression, normally consolidated --\n");
    const double M = M_cs();
    const auto path = drained(0.30, 6000);
    double worst_ea = 0.0, worst_ev = 0.0, worst_path = 0.0, worst_at = 0.0;
    for (const Pt& s : path) {
        if (s.ea < 0.002) continue;
        worst_path = std::max(worst_path, std::fabs(s.q - 3.0 * (s.p - kP0)) / kP0);
        if (s.q / s.p > 0.9 * M) continue;
        double ea, ev;
        drained_strains_at(s.p, ea, ev);
        if (rel(s.ea, ea) > worst_ea) worst_at = s.ea;
        worst_ea = std::max(worst_ea, rel(s.ea, ea));
        worst_ev = std::max(worst_ev, rel(s.ev, ev));
    }
    const Pt& f = path.back();
    const double pf = 3.0 * kP0 / (3.0 - M);
    std::printf("  q = 3 (p - p0) to %.1e; eps_a %.2e (worst at eps_a %.4f), eps_v %.2e against the integral; end q/(M p) = %.5f, p %.2f (p_f %.2f)\n",
                worst_path, worst_ea, worst_at, worst_ev, f.q / (M * f.p), f.p, pf);
    check(worst_path < 1e-8, "the stress path is the drained one, q = 3 (p - p0)");
    check(worst_ea < 5e-3 && worst_ev < 5e-3, "axial and volumetric strain are the integrals of the model (0.5%)");
    bool rising = true;
    for (size_t k = 1; k < path.size(); ++k)
        rising = rising && path[k].q / path[k].p >= path[k - 1].q / path[k - 1].p - 1e-12;
    check(rising && f.q / (M * f.p) < 1.0, "the stress ratio rises towards the critical state and stays below it");
}

void test_wrapper_identity() {
    std::printf("-- (d) the axisymmetric finite-element integration == the material point --\n");
    MaterialModel mm;
    mm.type = MaterialType::SoftSoil;
    mm.ssoil = mcc();
    GaussState gs;
    gs.stress = Eigen::Vector3d(-kP0, -kP0, 0.0);   // [r, z, rz] + theta, tension-positive
    gs.stress_zz = -kP0;
    gs.pp = 1.5 * kP0;
    const int n = 400;
    const double d = 0.10 / n;
    Eigen::Matrix4d Dt;
    for (int i = 0; i < n; ++i) {
        GaussState tr;
        integrate_point_axisym(mm, gs, Eigen::Vector4d(0.5 * d, -d, 0.0, 0.5 * d), tr, Dt);
        gs = tr;
    }
    const auto path = undrained(1.5, 0.10, n);
    const Pt& f = path.back();
    const Eigen::Vector3d s(-gs.stress(1), -gs.stress(0), -gs.stress_zz);
    std::printf("  p %.9f / %.9f, q %.9f / %.9f, p_p %.9f / %.9f\n", s.mean(), f.p, q_of(s), f.q, gs.pp, f.pp);
    check(rel(s.mean(), f.p) < 1e-8 && rel(q_of(s), f.q) < 1e-8 && rel(gs.pp, f.pp) < 1e-8,
          "p, q and p_p agree with the material point (1e-8)");
    check(std::fabs(gs.stress(0) - gs.stress_zz) < 1e-8 * kP0, "the radial and hoop stresses stay equal");
}

// The model has no tension cut-off, whatever the material's switch says (it is on by default): at an
// overconsolidation where the elastic path reaches a tensile principal stress, an undrained test still ends on
// the critical state line, at p_f = p0 (R/2)^Lambda.
void test_no_cutoff() {
    std::printf("-- the critical state of a heavily overconsolidated clay, with the default switch --\n");
    MaterialParams p;
    p.lam_star = 0.10; p.kap_star = 0.02; p.nu_ur = 0.2; p.phi_rad = 25.0 * 3.14159265358979323846 / 180.0;
    p.tension_cutoff = true;
    p.tensile_strength = 0.0;
    const MaterialModel mm = find_model("ModifiedCamClay")->build(p);
    check(!mm.tension_cutoff, "the model takes no cut-off");
    const double sphi = std::sin(p.phi_rad), M = 6.0 * sphi / (3.0 - sphi), Lam = 1.0 - p.kap_star / p.lam_star;
    constexpr double kR = 20.0;
    GaussState gs;
    gs.stress = Eigen::Vector3d(-kP0, -kP0, 0.0);
    gs.stress_zz = -kP0;
    gs.pp = kR * kP0;
    const int n = 6000;
    const double d = 0.30 / n;
    Eigen::Matrix4d Dt;
    for (int i = 0; i < n; ++i) {
        GaussState tr;
        integrate_point_axisym(mm, gs, Eigen::Vector4d(0.5 * d, -d, 0.0, 0.5 * d), tr, Dt);
        gs = tr;
    }
    const Eigen::Vector3d s(-gs.stress(1), -gs.stress(0), -gs.stress_zz);
    const double pf = kP0 * std::pow(0.5 * kR, Lam);
    std::printf("  OCR %.0f: p %.3f, q %.3f against p_f %.3f, M p_f %.3f\n", kR, s.mean(), q_of(s), pf, M * pf);
    check(rel(s.mean(), pf) < 2e-3 && rel(q_of(s), M * pf) < 2e-3, "it ends on the critical state line (0.2%)");
}

void test_column() {
    std::printf("-- (e) a confined column through the whole program --\n");
    const double K0 = k0_closed_form();
    m::Project pr;
    pr.name = "Cam-Clay column";
    pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 4;
    pr.has_water = false;
    pr.mesh.elem_size = 0.5;
    m::Material s; s.name = "Clay"; s.model = m::SoilModel::ModifiedCamClay;
    s.lam_star = kLam; s.kap_star = kKap; s.nu_ur = kNu; s.phi = kPhiDeg; s.c = 0.0; s.psi = 0.0;
    s.gamma_unsat = s.gamma_sat = 10.0;
    s.k0_auto = false; s.k0 = K0;
    s.tension_cutoff = false;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Clay"; P.material = 0;
    P.x = {0, 1, 1, 0}; P.y = {0, 0, 4, 4};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 0; L.y1 = 4; L.x2 = 1; L.y2 = 4; L.qy1 = L.qy2 = -100.0;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};
    m::Phase ph; ph.name = "Surcharge"; ph.type = m::PhaseType::Plastic; ph.load_active = {1};
    pr.phases.push_back(ph);
    const auto M = katai::app::mesh_from_project(pr);
    const auto res = M.ok ? katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure)
                          : std::vector<katai::app::SolveResult>{};
    const bool ok = res.size() == 2 && res[0].ok && res[1].ok;
    for (const auto& r : res) if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    check(ok, "the column solves, K0 procedure and a surcharge");
    if (!ok) return;
    double worst = 0.0, ratio_mid = 0.0;
    for (int k = 0; k < M.mesh.node_count; ++k) {
        const double y = M.mesh.y[k];
        if (y < 0.5 || y > 3.5) continue;   // away from the fixed base and the loaded top
        const auto& a = res[0].stress.stress[k];
        const auto& b = res[1].stress.stress[k];
        const double ratio = (b(0) - a(0)) / (b(1) - a(1));
        worst = std::max(worst, rel(ratio, K0));
        if (std::fabs(y - 2.0) < 1e-9) ratio_mid = ratio;
    }
    const double k_fine = point_k0(1e-4), k_coarse = point_k0(1e-3);
    std::printf("  closed-form K0 %.5f; the material point's %.5f (steps of 1e-4) and %.5f (1e-3); d sigma'_h / d sigma'_v at mid-height %.5f, worst %.2e\n",
                K0, k_fine, k_coarse, ratio_mid, worst);
    check(rel(k_fine, K0) < 5e-4 && rel(k_coarse, K0) < 5e-4,
          "the material point compresses at the closed-form K0, at either step (5e-4)");
    check(worst < 2e-3, "the column compresses at the model's one-dimensional K0 (0.2%)");
}
}  // namespace

int main() {
    std::printf("== Modified Cam-Clay (KV-CST-023) ==\n");
    test_undrained();
    test_drained();
    test_wrapper_identity();
    test_no_cutoff();
    test_column();
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
