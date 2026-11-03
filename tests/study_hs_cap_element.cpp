// Study (not in the suite): the undrained element tests of test_hs_multiaxial with the cap's M and beta chosen on the
// command line, against the published deviators, mean stresses and first-yield points.
//   study_hs_cap_element [M=<value>|Mfac=<factor>] [beta=<value>|betafac=<factor>|betadoc=1] [paths=1]
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/registry.hpp>

#include <chrono>
#include <cmath>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace katai::core;

namespace {
constexpr double kPi = 3.14159265358979323846;
struct Clay {
    const char* name;
    double E50, Eoed, Eur, c, phi, g07, G0, k0nc, ocr, pop;
    double k0, sv;
    double pc_ref, qc_ref, pe_ref, qe_ref, qc3_ref, qe3_ref;
    double kink_c, kink_e;
};

bool g_nold = false;
int g_seed = 0;   // 0 the present stress, 1 the pre-consolidation state, 2 none
MaterialModel build(const Clay& k) {
    MaterialParams p;
    p.E50_ref = k.E50; p.Eoed_ref = k.Eoed; p.Eur_ref = k.Eur; p.m = 1.0; p.p_ref = 100;
    p.Rf = 0.9; p.nu_ur = 0.15; p.c = k.c; p.phi_rad = k.phi * kPi / 180; p.psi_rad = 0;
    p.G0_ref = k.G0; p.gamma07 = k.g07; p.k0nc_auto = false; p.k0nc = k.k0nc;
    p.tension_cutoff = true; p.tensile_strength = 0;
    p.E = k.Eur; p.nu = 0.15;
    return find_model(k.G0 > 0 ? "HSsmall" : "HardeningSoil")->build(p);
}

GaussState seed(const MaterialModel& m, double sv, double sh, const Clay& k) {
    GaussState gs;
    gs.stress = Eigen::Vector3d(-sh, -sv, 0.0);
    gs.stress_zz = -sh;
    const auto pe = hs_small_strain_params(m.hs, 0.0);
    const HsSeed s = k.pop > 0 ? hs_seed_from_history(pe, sh, sv, 0.0, sh, 2, 1.0, k.pop)
                               : hs_seed_from_history(pe, sh, sv, 0.0, sh, 1, k.ocr, 0.0);
    gs.pp = s.pp;
    gs.gamma_p = s.gamma_p;
    if (g_seed == 1) {
        const double svc = k.pop > 0 ? sv + k.pop : k.ocr * sv;
        gs.gamma_p = hs_seed_from_history(pe, k.k0nc * svc, svc, 0.0, k.k0nc * svc, 1, 1.0, 0.0).gamma_p;
    }
    if (g_seed == 2) gs.gamma_p = 0.0;
    return gs;
}

struct Path { double q3, q10, p10, kink; std::vector<std::array<double, 3>> pts; };

// E50 as a drained triaxial test reads it: from an isotropic state at sigma3 = p_ref (normally consolidated, or far
// inside the cap), sigma3 held, the secant at half the failure deviator.
double measured_e50(HardeningSoilParams p, bool nc, double e50_step = 0.0) {
    p.G0_ref = 0.0;
    const double s3 = p.p_ref, qf = p.q_failure(s3);
    Eigen::Vector3d sig(s3, s3, s3);
    double gp = 0.0, pp = nc ? hs_initial_pp(p, sig) : 100.0 * s3;
    double e1 = 0.0, q_prev = 0.0, e_prev = 0.0;
    const double de1 = 2.0e-5 * (qf / (e50_step > 0.0 ? e50_step : p.E50(s3))) / 0.01 * (e50_step > 0.0 ? 0.1 : 1.0);
    for (int i = 0; i < 20000; ++i) {
        // de3 so that the lateral stress stays
        double de3 = -0.2 * de1;
        HsIntegrated r{};
        for (int it = 0; it < 8; ++it) {
            r = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, de3, de3));
            const double err = r.stress(1) - s3;
            if (std::fabs(err) < 1e-8 * s3) break;
            const HsIntegrated h = hs_integrate(p, sig, gp, pp, Eigen::Vector3d(de1, de3 + 1e-3 * de1, de3 + 1e-3 * de1));
            const double slope = (h.stress(1) - r.stress(1)) / (1e-3 * de1);
            de3 -= err / slope;
        }
        sig = r.stress; gp = r.gamma_p; pp = r.pp;
        e1 += de1;
        const double q = sig(0) - sig(1);
        if (q >= 0.5 * qf) {
            const double e50 = e_prev + (0.5 * qf - q_prev) / (q - q_prev) * (e1 - e_prev);
            return 0.5 * qf / e50;
        }
        q_prev = q; e_prev = e1;
    }
    return 0.0;
}

void report_e50(const char* name, double E50, double Eoed, double Eur, double m, double c, double phi, double psi, double nu, double k0nc) {
    MaterialParams q;
    q.E50_ref = E50; q.Eoed_ref = Eoed; q.Eur_ref = Eur; q.m = m; q.p_ref = 100;
    q.Rf = 0.9; q.nu_ur = nu; q.c = c; q.phi_rad = phi * kPi / 180; q.psi_rad = psi * kPi / 180;
    q.k0nc_auto = false; q.k0nc = k0nc; q.tension_cutoff = true; q.tensile_strength = 0;
    q.E = Eur; q.nu = nu;
    MaterialModel mm = find_model("HardeningSoil")->build(q);
    const double nc = measured_e50(mm.hs, true), oc = measured_e50(mm.hs, false);
    double lo = 1.0, hi = 3.0;
    for (int i = 0; i < 25; ++i) {
        const double mid = 0.5 * (lo + hi);
        HardeningSoilParams t = mm.hs; t.E50_ref = E50 * mid;
        if (measured_e50(t, true) < E50) lo = mid; else hi = mid;
    }
    std::printf("E50 %s: asked %.0f, read with the cap %.0f (%.3f), far from it %.0f (%.3f); factor on E_i for the asked value %.3f; M %.3f\n",
                name, E50, nc, nc / E50, oc, oc / E50, 0.5 * (lo + hi), mm.hs.cap_alpha);
}

Path undrained(const MaterialModel& m, const Clay& k, int dir, bool keep) {
    GaussState gs = seed(m, k.sv, k.k0 * k.sv, k);
    const int n = 5000;
    const double d = 0.10 / n;
    const double p0 = k.sv * (1.0 + 2.0 * k.k0) / 3.0;
    Path out{0, 0, 0, 0, {}};
    double pv = 0.0;
    for (int i = 1; i <= n; ++i) {
        GaussState tr = gs;
        Eigen::Matrix4d T;
        const Eigen::Vector4d de(dir * d / 2, -dir * d, 0.0, dir * d / 2);
        integrate_point_axisym(m, gs, de, tr, T, TangentMode::kContinuum);
        gs = tr;
        const double sr = -gs.stress(0), sz = -gs.stress(1), sth = -gs.stress_zz;
        const double p = (sr + sz + sth) / 3.0;
        out.q10 = std::fabs(sz - sr);
        out.p10 = p;
        if (i == n * 3 / 10) out.q3 = out.q10;
        if (pv == 0.0 && i == 3) pv = p;
        if (out.kink == 0.0 && pv > 0.0 && std::fabs(p - pv) > 0.012 * pv) out.kink = out.q10;
        if (keep && i % 10 == 0) out.pts.push_back({100.0 * i * d, out.q10, p});
    }
    (void)p0;
    return out;
}
}  // namespace

int main(int argc, char** argv) {
    std::map<std::string, double> a;
    for (int i = 1; i < argc; ++i) {
        const char* eq = std::strchr(argv[i], '=');
        if (eq) a[std::string(argv[i], eq - argv[i])] = std::atof(eq + 1);
    }
    const bool paths = a.count("paths") > 0;
    g_nold = a.count("nold") > 0;
    if (a.count("seed")) g_seed = (int)a["seed"];
    const Clay clays[] = {
        {"Clay4", 15400, 8300, 30800, 1, 32.0, 4e-4, 39230, 0.53, 1.00, 32, 0.61, 126.5, 80.4, 107.1, 81.0, 75.6, 106.3, 74.5, 83.3, 61.0},
        {"Clay5", 13100, 6600, 26200, 0, 30.5, 6e-4, 40960, 0.53, 1.00, 32, 0.58, 186.3, 117.7, 144.6, 118.6, 103.5, 142.1, 101.7, 111.8, 84.0},
        {"Clay6", 13100, 6600, 26200, 0, 30.5, 8e-4, 36060, 0.53, 1.14, 0, 0.57, 248.8, 154.7, 190.2, 156.6, 136.4, 186.7, 133.9, 148.1, 109.0},
        {"Clay7", 10500, 5285, 21000, 0, 30.5, 9e-4, 37460, 0.53, 1.14, 0, 0.56, 408.0, 252.4, 309.5, 255.3, 221.7, 301.3, 216.1, 250.1, 180.0},
    };
    if (a.count("berlin")) {
        // Berlin sand III, drained and undrained triaxial from an isotropic 200 kPa, as published with the reference curves.
        MaterialParams q;
        q.E50_ref = 105000; q.Eoed_ref = 105000; q.Eur_ref = 315000; q.m = 0.55; q.p_ref = 100;
        q.Rf = 0.9; q.nu_ur = 0.2; q.c = 1.0; q.phi_rad = 38.0 * kPi / 180; q.psi_rad = 6.0 * kPi / 180;
        q.k0nc_auto = false; q.k0nc = 0.38; q.tension_cutoff = true; q.tensile_strength = 0;
        q.E = 315000; q.nu = 0.2;
        MaterialModel mm = find_model("HardeningSoil")->build(q);
        HardeningSoilParams pp = mm.hs;
        std::printf("M %.3f beta %.3e f %.4f\n", pp.cap_alpha, pp.cap_beta, pp.ei_scale);
        if (a.count("e50fac")) pp.E50_ref *= a["e50fac"];
        if (a.count("M")) pp.cap_alpha = a["M"];
        if (a.count("beta")) pp.cap_beta = a["beta"];
        if (a.count("betafac")) pp.cap_beta *= a["betafac"];
        if (a.count("f")) pp.ei_scale = a["f"];
        const double s3 = 200.0;
        for (int drained = 1; drained >= 0; --drained) {
            Eigen::Vector3d sig(s3, s3, s3);
            double gp = 0.0, ppv = hs_initial_pp(pp, sig);
            const double de1 = 2.0e-5;
            double e1 = 0.0, e3 = 0.0;
            std::printf("%s: eps1 eps_q s1 s3 q p eps_v\n", drained ? "drained" : "undrained");
            int next = 0;
            const double marks[] = {0.0005, 0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05};
            for (int i = 0; i < 2500; ++i) {
                double de3 = drained ? -0.2 * de1 : -0.5 * de1;
                HsIntegrated r{};
                if (drained) {
                    for (int it = 0; it < 8; ++it) {
                        r = hs_integrate(pp, sig, gp, ppv, Eigen::Vector3d(de1, de3, de3));
                        const double err = r.stress(1) - s3;
                        if (std::fabs(err) < 1e-8 * s3) break;
                        const HsIntegrated h = hs_integrate(pp, sig, gp, ppv, Eigen::Vector3d(de1, de3 + 1e-3 * de1, de3 + 1e-3 * de1));
                        de3 -= err / ((h.stress(1) - r.stress(1)) / (1e-3 * de1));
                    }
                } else {
                    r = hs_integrate(pp, sig, gp, ppv, Eigen::Vector3d(de1, de3, de3));
                }
                sig = r.stress; gp = r.gamma_p; ppv = r.pp;
                e1 += de1; e3 += de3;
                if (next < 13 && e1 >= marks[next] - 1e-12) {
                    std::printf("  %.4f %.5f %.1f %.1f %.1f %.1f %+.5f\n", e1, (2.0 / 3.0) * (e1 - e3), sig(0), sig(1), sig(0) - sig(1),
                                (sig(0) + 2 * sig(1)) / 3, e1 + 2 * e3);
                    ++next;
                }
            }
        }
        return 0;
    }
    if (a.count("lim")) {
        // How much of the asked E50 the cap leaves room for, as the shear surface is made stiffer.
        for (double k0nc : {0.53, 0.45}) {
            for (double ratio : {0.40, 0.45, 0.50, 0.53, 0.60, 0.70, 0.80, 1.00}) {
                MaterialParams q;
                q.E50_ref = 13100; q.Eoed_ref = 13100 * ratio; q.Eur_ref = 26200; q.m = 1.0; q.p_ref = 100;
                q.Rf = 0.9; q.nu_ur = 0.15; q.c = 0; q.phi_rad = 30.5 * kPi / 180; q.psi_rad = 0;
                q.k0nc_auto = false; q.k0nc = k0nc; q.tension_cutoff = true; q.tensile_strength = 0;
                q.E = 26200; q.nu = 0.15;
                MaterialModel mm = find_model("HardeningSoil")->build(q);
                std::printf("K0nc %.2f Eoed/E50 %.2f: M %.3f; E50 read / asked at f =", k0nc, ratio, mm.hs.cap_alpha);
                for (double f : {1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 100.0}) {
                    HardeningSoilParams t = mm.hs;
                    t.E50_ref = 13100 * f;
                    std::printf(" %.3f", measured_e50(t, true, 13100.0) / 13100.0);
                }
                std::printf("\n");
            }
        }
        return 0;
    }
    if (a.count("e50")) {
        for (const Clay& k : clays) report_e50(k.name, k.E50, k.Eoed, k.Eur, 1.0, k.c, k.phi, 0.0, 0.15, k.k0nc);
        report_e50("L3 clay", 4000, 3300, 12000, 1.0, 1.0, 25.0, 0.0, 0.15, 0.5774);
        report_e50("L3 sand", 40000, 40000, 120000, 0.5, 0.0, 32.0, 2.0, 0.2, 0.4701);
        report_e50("L2 clay", 5600, 5000, 20000, 1.0, 10.0, 25.0, 0.0, 0.2, 0.5774);
        report_e50("L7 silt", 20000, 20000, 60000, 0.5, 1.0, 30.0, 0.0, 0.2, 0.5);
        report_e50("L7 sand", 30000, 30000, 90000, 0.5, 0.0, 34.0, 4.0, 0.2, 0.4408);
        report_e50("L7 loam", 12000, 8000, 36000, 0.8, 5.0, 29.0, 0.0, 0.2, 0.5152);
        return 0;
    }
    if (a.count("stress")) {
        // The calibration over a spread of parameters: every one must end finite, in range, and in reasonable time.
        std::mt19937 rng(20261102);
        const auto u = [&](double lo, double hi) { return lo + (hi - lo) * std::uniform_real_distribution<double>(0.0, 1.0)(rng); };
        int bad = 0, ceiling = 0, floor_ = 0;
        double slowest = 0.0, total = 0.0;
        const int n = (int)a["stress"];
        for (int i = 0; i < n; ++i) {
            MaterialParams q;
            q.E50_ref = std::exp(u(std::log(1e3), std::log(2e5)));
            q.Eoed_ref = q.E50_ref * u(0.2, 2.0);
            q.Eur_ref = q.E50_ref * u(1.5, 10.0);
            q.m = u(0.0, 1.0); q.p_ref = 100; q.Rf = u(0.6, 0.99); q.nu_ur = u(0.1, 0.35);
            const double phi = u(5.0, 45.0);
            q.phi_rad = phi * kPi / 180; q.c = u(0.0, 1.0) < 0.3 ? 0.0 : u(0.0, 50.0);
            q.psi_rad = u(0.0, 1.0) < 0.5 ? 0.0 : u(0.0, std::max(0.0, std::min(phi - 5.0, 15.0))) * kPi / 180;
            q.k0nc_auto = u(0.0, 1.0) < 0.5; q.k0nc = u(0.3, 0.9); q.tension_cutoff = true; q.tensile_strength = 0;
            q.E = q.Eur_ref; q.nu = q.nu_ur;
            const auto t0 = std::chrono::steady_clock::now();
            const MaterialModel mm = find_model("HardeningSoil")->build(q);
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            slowest = std::max(slowest, sec); total += sec;
            const double f = mm.hs.ei_scale, al = mm.hs.cap_alpha, be = mm.hs.cap_beta;
            const bool ok = std::isfinite(f) && f >= 1.0 && f <= HardeningSoilParams::kMaxEiScale + 1e-9 && std::isfinite(al) && al > 0.0 && std::isfinite(be) && be > 0.0;
            if (f >= HardeningSoilParams::kMaxEiScale) ++ceiling;
            if (f <= 1.0) ++floor_;
            if (!ok || sec > 20.0) {
                ++bad;
                std::printf("BAD %d E50 %.0f Eoed/E50 %.2f Eur/E50 %.2f m %.2f phi %.1f c %.1f psi %.1f nu %.2f Rf %.2f K0nc %.2f: f %.3f M %.3f beta %.3e, %.1f s\n", i, q.E50_ref, q.Eoed_ref / q.E50_ref, q.Eur_ref / q.E50_ref, q.m, phi, q.c, q.psi_rad * 180 / kPi, q.nu_ur, q.Rf, mm.hs.K0nc, f, al, be, sec);
            }
        }
        std::printf("%d sets: %d bad, %d at the ceiling, %d at f = 1; mean %.2f s, slowest %.2f s\n", n, bad, ceiling, floor_, total / n, slowest);
        return bad ? 1 : 0;
    }
    if (a.count("shallow")) {
        // The three shallow, more overconsolidated clays of the same source; their start stresses from the curves at zero strain.
        const Clay shallow[] = {
            {"Clay1", 7000, 3002, 27000, 1, 32.0, 4e-4, 77480, 0.53, 2.00, 0, 0.70, 49.0, 36.7, 49.4, 37.2, 34.5, 48.7, 34.2, 0.0, 0.0},
            {"Clay2", 7000, 3400, 33700, 1, 32.0, 4e-4, 39700, 0.47, 1.65, 0, 0.70, 64.3, 45.1, 60.3, 46.6, 43.0, 58.2, 41.9, 0.0, 0.0},
            {"Clay3", 7000, 3068, 24250, 1, 32.0, 4e-4, 40210, 0.53, 1.00, 32, 0.63, 84.6, 57.2, 75.9, 57.6, 52.9, 74.8, 52.2, 0.0, 0.0},
        };
        for (const Clay& k : shallow) {
            MaterialModel m = build(k);
            const Path c = undrained(m, k, +1, false), e = undrained(m, k, -1, false);
            std::printf("%s f %.3f M %.3f | C q3 %.1f/%.1f (%+.1f%%) q10 %.1f/%.1f (%+.1f%%) p10 %.1f/%.1f | E q3 %.1f/%.1f (%+.1f%%) q10 %.1f/%.1f (%+.1f%%) p10 %.1f/%.1f\n",
                        k.name, m.hs.ei_scale, m.hs.cap_alpha, c.q3, k.qc3_ref, 100 * (c.q3 / k.qc3_ref - 1), c.q10, k.qc_ref, 100 * (c.q10 / k.qc_ref - 1), c.p10, k.pc_ref,
                        e.q3, k.qe3_ref, 100 * (e.q3 / k.qe3_ref - 1), e.q10, k.qe_ref, 100 * (e.q10 / k.qe_ref - 1), e.p10, k.pe_ref);
        }
        return 0;
    }
    double s2 = 0.0;
    int n2 = 0;
    for (const Clay& k : clays) {
        MaterialModel m = build(k);
        const double M0 = m.hs.cap_alpha, b0 = m.hs.cap_beta;
        if (g_nold) m.hs.dilatancy_cut = true;
        if (a.count("e50fac")) m.hs.E50_ref *= a["e50fac"];
        const double nu = k.E50 > 0 ? 0.15 : 0.15;
        const double r = (k.Eur / k.Eoed) * k.k0nc / ((1.0 + 2.0 * k.k0nc) * (1.0 - 2.0 * nu));
        const double Ks = k.Eur / (3.0 * (1.0 - 2.0 * nu));
        const double bdoc = 100.0 * (r - 1.0) / Ks;
        if (a.count("M")) m.hs.cap_alpha = a["M"];
        if (a.count("Mfac")) m.hs.cap_alpha = M0 * a["Mfac"];
        if (a.count("beta")) m.hs.cap_beta = a["beta"];
        if (a.count("betafac")) m.hs.cap_beta = b0 * a["betafac"];
        if (a.count("betadoc")) m.hs.cap_beta = bdoc;
        const Path c = undrained(m, k, +1, paths), e = undrained(m, k, -1, paths);
        const double err[] = {c.q3 / k.qc3_ref - 1, c.q10 / k.qc_ref - 1, c.p10 / k.pc_ref - 1,
                              e.q3 / k.qe3_ref - 1, e.q10 / k.qe_ref - 1, e.p10 / k.pe_ref - 1};
        for (double x : err) { s2 += x * x; ++n2; }
        std::printf("%s f %.3f M %.3f (cal %.3f) beta %.3e (cal %.3e doc %.3e) | C q3 %+.1f%% q10 %+.1f%% p10 %+.1f%% kink %.1f/%.1f | "
                    "E q3 %+.1f%% q10 %+.1f%% p10 %+.1f%% kink %.1f/%.1f\n",
                    k.name, m.hs.ei_scale, m.hs.cap_alpha, M0, m.hs.cap_beta, b0, bdoc, 100 * err[0], 100 * err[1], 100 * err[2], c.kink, k.kink_c,
                    100 * err[3], 100 * err[4], 100 * err[5], e.kink, k.kink_e);
        if (paths) {
            for (const auto& q : c.pts) std::printf("PATH %s C %.3f %.3f %.3f\n", k.name, q[0], q[1], q[2]);
            for (const auto& q : e.pts) std::printf("PATH %s E %.3f %.3f %.3f\n", k.name, q[0], q[1], q[2]);
        }
    }
    std::printf("rms of the six quantities over the four clays: %.2f%%\n", 100.0 * std::sqrt(s2 / n2));
    return 0;
}
