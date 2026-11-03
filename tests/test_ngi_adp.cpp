// NGI-ADP: the material point against the model's own closed forms.
//
// The model (ngi_adp.hpp; docs/references/ngi-adp-formulation.md) is an anisotropic undrained
// shear strength model in total stress: three strengths -- plane strain active, plane strain
// passive, direct simple shear -- and the shear strain at failure in three tests. None of the
// oracles below shares code with the return mapping: they are the strengths the model is built
// to reach, and the strain at which a test run from the in-situ state reaches them.
//
// verify: KV-CST-024
//   oracle:   closed_form
//   source:   Grimstad, G., Andresen, L. & Jostad, H.P. (2012). NGI-ADP: Anisotropic shear strength model for clay. International Journal for Numerical and Analytical Methods in Geomechanics 36(4), 483-497 -- the yield surface translated from the in-situ state, hardening to the three strengths at the three failure strains; Davis, E.H. & Christian, J.T. (1971). Bearing capacity of anisotropic cohesive soil. J. Soil Mech. Found. Div. ASCE 97(SM5), 753-769 -- the plane-strain ellipse the surface reduces to at failure
//   locator:  (a)-(d) at the point: suA = 50 kPa, suP/suA = 0.4, suDSS/suA = 0.7, G = 500 suA, nu_u = 0.495, gf = 1.5 / 4 / 3 % (compression, extension, simple shear), an in-situ state sigma_v0 = 100, sigma_h0 = 60 kPa (compression positive); strain-controlled paths at isochoric strain: plane strain active and passive, triaxial compression and extension, direct simple shear; (f) through the program: a 1 x 1 m plane strain block (tri6, 0.25 m), smooth base and side, its top pushed down and pulled up 0.05 m; (e) the corpus' Prandtl strip footing (KV-FND-010, B = 1.2 m, su = 10 kPa, q = 60 kPa) with its material made isotropic NGI-ADP (suP = suDSS = suA, G = 500 su); (g) the same footing on anisotropic NGI-ADP with suP = 0.4 suA and suDSS = (suA + suP)/2, where the plane strain failure ellipse is a circle about a shifted centre
//   quantity: the shear stress at 20% strain [kPa]; the total shear strain at which the hardening reaches 1 [-]; (f) the mean (sigma_yy - sigma_xx)/2 over the block [kPa]; (e) the footing's collapse pressure [kPa]
//   expected: plane strain (sigma_yy - sigma_xx)/2 = suA and -suP; triaxial suC,TX = 0.99 suA and -suE,TX = -(S x - Da)/2, x = sqrt(3/(4 H(1))); simple shear t_xy = suDSS; each test, run from the in-situ state, reaches its strength at the total shear strain entered for it; (f) suA and -suP; (e) Prandtl (2 + pi) su = 51.42 kPa, and what the same mesh carries on a Tresca clay (Mohr-Coulomb, phi = 0); (g) (2 + pi)(suA + suP)/2 = 35.99 kPa -- the constant shift of the centre is no part of equilibrium, a surcharge at the free surface that comes back out under the footing (the circular case of Davis & Christian 1971) -- and a Tresca clay of that mean strength on the same mesh
//   band:     1e-6 relative on the strengths (measured to the printed 6 digits); 1e-3 relative on the failure strains (the step of the path, 1e-5, interpolated; measured 0.000% in triaxial compression and extension, -0.017% in simple shear, where the hardening runs at a Lode angle that moves slightly while the normal stresses settle); (f) 1e-4 (measured to the printed 6 digits); (e) 3% on Prandtl (measured -0.81%) and 1% on the Tresca run (-0.73%); (g) 3% on the closed form (measured +0.03%) and 1% on the Tresca run (0.000%). The collapse pressures resolve to the load step the run ended on (51.000 = 0.85 x 60, 36.000 = 0.60 x 60), so the measured differences are as fine as the stepping, not finer
#include <katai/jobs/driver.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/materials/ngi_adp.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace adp = katai::core::adp;
namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }

constexpr double kSuA = 50.0;

adp::Params params(double r_tx) {
    adp::Params P;
    P.suA = kSuA;
    P.suP = 0.4 * kSuA;
    P.suDSS = 0.7 * kSuA;
    P.G = 500.0 * kSuA;
    P.nu = 0.495;
    P.gfC = 0.015;
    P.gfE = 0.04;
    P.gfDSS = 0.03;
    P.r_tx = r_tx;
    return P;
}

adp::State in_situ() {
    adp::State s;
    s.sig << 60.0, 100.0, 60.0, 0.0;
    s.s0 = s.sig;
    return s;
}

struct Path {
    double tau_end = 0.0;      // the shear stress the path ends at
    double gamma_fail = -1.0;  // the total shear strain (its size) at which the hardening reached 1; -1 not reached
    bool ok = true;
};

// A strain-controlled path along `dir` (per unit of the shear strain measure), `steps` steps of `h`.
// `tau` reads the path's shear stress, `gamma` its shear strain from the strain increment.
template <class Tau, class Gam>
Path drive(const adp::Params& P, const adp::Vec4& dir, double h, int steps, Tau tau, Gam gamma) {
    Path out;
    adp::State s = in_situ();
    const adp::Derived d = adp::derive(P);
    double gC, gE, gD;
    adp::plastic_failure_strains(P, d, s.s0, gC, gE, gD);
    double g_total = 0.0;
    for (int i = 0; i < steps; ++i) {
        adp::State n;
        const adp::Vec4 de = h * dir;
        if (!adp::integrate(P, s, de, n)) { out.ok = false; return out; }
        const double g_next = g_total + std::fabs(gamma(de));
        // kappa = 1 is reached where gp passes the failure strain of the path's orientation
        const adp::Mod m = adp::modified(n.sig, n.s0, 1.0, d);
        const double gpf = adp::failure_strain(m, gC, gE, gD);
        if (out.gamma_fail < 0.0 && n.gp >= gpf) {
            const double t = (gpf - s.gp) / (n.gp - s.gp);
            out.gamma_fail = g_total + t * (g_next - g_total);
        }
        s = n;
        g_total = g_next;
    }
    out.tau_end = tau(s.sig);
    return out;
}

}  // namespace

// Direct simple shear before the peak. The failure strain is read at the orientation of the modified deviator, which moves
// with the hardening itself; read at a k that is not the answer, the curve was 5.7% low at 5% of the failure strain, where
// triaxial compression and extension agree to 0.1%. Symbols of Ukritchon & Boonyatee (2015), Geotech. Eng. J. 46(1), Fig. 13.
void test_dss_curve() {
    std::printf("-- the direct simple shear curve before its peak --\n");
    adp::Params P;
    P.suA = 60.0; P.suP = 24.0; P.suDSS = 42.0; P.G = 200.0 * 60.0; P.nu = 0.495;
    P.gfC = 0.04; P.gfE = 0.15; P.gfDSS = 0.07; P.r_tx = 0.99;
    adp::State s;
    s.sig << -0.0, -0.0, -0.0, 0.0;
    s.s0 = s.sig;
    const double g_end = 0.0105, h = g_end / 3000.0;
    double g = 0.0;
    double at35 = 0.0, at70 = 0.0, at105 = 0.0;
    for (int i = 0; i < 3000; ++i) {
        adp::State n;
        adp::Vec4 de; de << 0.0, 0.0, 0.0, h;
        if (!adp::integrate(P, s, de, n)) { check(false, "the path integrates"); return; }
        s = n;
        g += h;
        if (i == 999) at35 = std::fabs(s.sig(3));
        if (i == 1999) at70 = std::fabs(s.sig(3));
        if (i == 2999) at105 = std::fabs(s.sig(3));
    }
    std::printf("   tau at gamma 0.35, 0.70, 1.05%%: %.3f %.3f %.3f against 14.405 21.314 26.024\n", at35, at70, at105);
    check(rel(at35, 14.4048) < 0.012 && rel(at70, 21.3135) < 0.01 && rel(at105, 26.0244) < 0.01,
          "on the published curve (1%, the symbols' own scatter)");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // what it got through is in the log if it runs out of time
    std::printf("== NGI-ADP at the material point (KV-CST-024) ==\n");
    test_dss_curve();
    const auto tau_v = [](const adp::Vec4& s) { return 0.5 * (s(1) - s(0)); };
    const auto tau_xy = [](const adp::Vec4& s) { return s(3); };
    const auto g_v = [](const adp::Vec4& de) { return de(1) - de(0); };   // eps1 - eps3
    const auto g_xy = [](const adp::Vec4& de) { return de(3); };
    adp::Vec4 ps_a, ps_p, tx_c, tx_e, dss;
    ps_a << -0.5, 0.5, 0.0, 0.0;        // plane strain, vertical compression: eps_yy - eps_xx = 1
    ps_p << 0.5, -0.5, 0.0, 0.0;
    tx_c << -1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, 0.0;   // eps1 - eps3 = 1, isochoric
    tx_e << 1.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0, 0.0;
    dss << 0.0, 0.0, 0.0, 1.0;
    {
        const adp::Params P = params(0.99);
        const adp::Derived d = adp::derive(P);
        std::printf("      derived: a1 = %.6f, H(1) = %.6f, suC,TX = %.4f, suE,TX = %.4f kPa\n", d.a1, d.H1, d.suC_tx, d.suE_tx);
        const Path a = drive(P, ps_a, 1e-4, 2000, tau_v, g_v);
        const Path p = drive(P, ps_p, 1e-4, 2000, tau_v, g_v);
        const Path c = drive(P, tx_c, 1e-4, 2000, tau_v, g_v);
        const Path e = drive(P, tx_e, 1e-4, 2000, tau_v, g_v);
        const Path s = drive(P, dss, 1e-4, 2000, tau_xy, g_xy);
        std::printf("      plane strain active %.6f (suA %.4f), passive %.6f (-suP %.4f)\n", a.tau_end, P.suA, p.tau_end, -P.suP);
        std::printf("      triaxial compression %.6f (suC,TX %.4f), extension %.6f (-suE,TX %.4f)\n", c.tau_end, d.suC_tx, e.tau_end, -d.suE_tx);
        std::printf("      direct simple shear %.6f (suDSS %.4f)\n", s.tau_end, P.suDSS);
        check(a.ok && p.ok && c.ok && e.ok && s.ok, "every path integrates");
        check(rel(a.tau_end, P.suA) < 1e-6, "plane strain active: (sigma_yy - sigma_xx)/2 = suA");
        check(rel(p.tau_end, -P.suP) < 1e-6, "plane strain passive: (sigma_yy - sigma_xx)/2 = -suP");
        check(std::fabs(d.suC_tx / P.suA - 0.99) < 1e-12 && rel(c.tau_end, d.suC_tx) < 1e-6,
              "triaxial compression: suC,TX = 0.99 suA");
        check(rel(e.tau_end, -d.suE_tx) < 1e-6, "triaxial extension: -suE,TX");
        check(rel(s.tau_end, P.suDSS) < 1e-6, "direct simple shear: t_xy = suDSS");
    }
    {
        // The strain at failure: a test run from the in-situ state reaches its strength at the
        // total shear strain entered for it.
        const adp::Params P = params(0.99);
        const Path c = drive(P, tx_c, 1e-5, 6000, tau_v, g_v);
        const Path e = drive(P, tx_e, 1e-5, 6000, tau_v, g_v);
        const Path s = drive(P, dss, 1e-5, 6000, tau_xy, g_xy);
        std::printf("      failure strains: triaxial compression %.6f (gf_C %.4f), extension %.6f (gf_E %.4f), simple shear %.6f (gf_DSS %.4f)\n",
                    c.gamma_fail, P.gfC, e.gamma_fail, P.gfE, s.gamma_fail, P.gfDSS);
        check(c.ok && c.gamma_fail > 0.0 && rel(c.gamma_fail, P.gfC) < 1e-3,
              "triaxial compression fails at the total shear strain entered");
        check(e.ok && e.gamma_fail > 0.0 && rel(e.gamma_fail, P.gfE) < 1e-3,
              "triaxial extension fails at the total shear strain entered");
        check(s.ok && s.gamma_fail > 0.0 && rel(s.gamma_fail, P.gfDSS) < 1e-3,
              "direct simple shear fails at the total shear strain entered");
    }
    // ---- through the whole program ----
    const auto adp_material = [](double su, double ratioP, double ratioDSS) {
        m::Material s;
        s.name = "NGI-ADP clay"; s.model = m::SoilModel::NgiAdp; s.drainage = m::Drainage::UndrainedC;
        s.gamma_unsat = s.gamma_sat = 0.0;
        s.c = su; s.phi = 0.0; s.psi = 0.0;
        s.nu = 0.495;
        s.E = 2.0 * (1.0 + s.nu) * 500.0 * su;   // G = 500 suA
        s.adp_suP = ratioP; s.adp_suDSS = ratioDSS;
        s.adp_gfC = 1.5; s.adp_gfE = 4.0; s.adp_gfDSS = 3.0;
        s.k0_auto = false; s.k0 = 1.0;
        return s;
    };
    {
        // (f) A plane strain block pushed down and pulled up by its top, smooth on the base and
        // the left side: (sigma_yy - sigma_xx)/2 reaches suA and -suP.
        const auto block = [&](double uy) {
            m::Project pr;
            pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 1;
            pr.has_water = false;
            pr.initial_procedure = m::InitialProcedure::GravityLoading;
            pr.mesh.elem_size = 0.25;
            pr.mesh.auto_refine = false;
            pr.materials.push_back(adp_material(kSuA, 0.4, 0.7));
            m::SoilPolygon g;
            g.name = "Block"; g.material = 0; g.x = {0, 1, 1, 0}; g.y = {0, 0, 1, 1};
            g.edge_bc = {(int)m::BCType::VerticallyFixed, 0, 0, (int)m::BCType::HorizontallyFixed};
            pr.polygons.push_back(g);
            m::PrescribedDisp D;
            D.name = "Top"; D.x1 = 0; D.y1 = 1; D.x2 = 1; D.y2 = 1; D.set_uy = true; D.uy = uy;
            pr.disps.push_back(D);
            pr.initial.disp_active = {0};
            m::Phase ph = pr.initial;
            ph.name = "Shear"; ph.disp_active = {1};
            pr.phases.push_back(ph);
            const auto M = katai::app::mesh_from_project(pr);
            double tau = 0.0;
            bool ok = M.ok;
            if (ok) {
                const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
                ok = res.size() == 2 && res[1].ok;
                if (ok) {
                    for (const auto& s : res[1].stress.stress) tau += 0.5 * (s(0) - s(1));   // compression positive
                    tau /= (double)res[1].stress.stress.size();
                }
            }
            return std::make_pair(ok, tau);
        };
        const auto a = block(-0.05), p = block(0.05);
        std::printf("      FE block: pushed %.6f (suA %.4f), pulled %.6f (-suP %.4f) kPa\n", a.second, kSuA, p.second, -0.4 * kSuA);
        check(a.first && rel(a.second, kSuA) < 1e-4, "through the program: a plane strain block pushed down reaches suA");
        check(p.first && rel(p.second, -0.4 * kSuA) < 1e-4, "...and pulled up, -suP");
    }
    {
        // (e) A strip footing on isotropic NGI-ADP clay (suP = suDSS = suA) is a footing on a
        // Tresca clay: the plane strain failure is at w = 0, where the rounding is exact. The
        // corpus' Prandtl case (KV-FND-010) with its material replaced.
        m::Project pr;
        std::string err;
        const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-fnd-010-prandtl-strip-footing.k2d";
        const bool loaded = m::load_project(path, pr, &err, nullptr);
        check(loaded, "the Prandtl case loads");
        if (loaded) {
            const auto run = [](const m::Project& p) {
                const auto M = katai::app::mesh_from_project(p);
                if (!M.ok) return -1.0;
                const auto res = katai::app::solve_phases(p, M.mesh, katai::app::initial_phase_from(p.initial_procedure));
                return res.size() == 2 ? res[1].load_factor : -1.0;
            };
            const double q = 60.0, su = pr.materials[0].c;
            const double lf_mc = run(pr);
            m::Project pa = pr;
            pa.materials[0] = adp_material(su, 1.0, 1.0);
            const double lf_adp = run(pa);
            const double prandtl = (2.0 + 3.14159265358979323846) * su;
            std::printf("      strip footing: NGI-ADP isotropic %.3f kPa, Tresca (Mohr-Coulomb, phi = 0) %.3f kPa, Prandtl (2 + pi) su %.3f kPa\n",
                        lf_adp * q, lf_mc * q, prandtl);
            check(lf_adp > 0.0 && rel(lf_adp * q, prandtl) < 0.03, "a footing on isotropic NGI-ADP clay carries (2 + pi) su, within 3%");
            check(lf_adp > 0.0 && lf_mc > 0.0 && rel(lf_adp, lf_mc) < 0.01,
                  "...and what the same mesh carries on a Tresca clay, within 1%");

            // (g) ANISOTROPIC, AND STILL CLOSED FORM. With suDSS = (suA + suP)/2 the plane strain
            // failure ellipse is a circle of radius (suA + suP)/2 about a centre shifted by
            // (suA - suP)/2 in (sigma_yy - sigma_xx)/2. A constant deviatoric shift is no part of
            // equilibrium: at the free surface it acts as a surcharge, under the footing it comes
            // back out, and the Prandtl field gives q_u = (2 + pi) (suA + suP)/2 -- the bearing
            // capacity on the mean of the active and passive strengths (the circular case of Davis &
            // Christian 1971). Held against the closed form and against a Tresca clay of that mean
            // strength on the same mesh.
            const double ratioP = 0.4;
            m::Project pn = pr;
            pn.materials[0] = adp_material(su, ratioP, 0.5 * (1.0 + ratioP));
            const double lf_an = run(pn);
            m::Project pt = pr;
            pt.materials[0].c = 0.5 * (1.0 + ratioP) * su;
            const double lf_tm = run(pt);
            const double dc = (2.0 + 3.14159265358979323846) * 0.5 * (1.0 + ratioP) * su;
            std::printf("      anisotropic (suP = 0.4 suA, suDSS = 0.7 suA): %.3f kPa, Tresca at the mean strength %.3f kPa, (2 + pi)(suA + suP)/2 %.3f kPa\n",
                        lf_an * q, lf_tm * q, dc);
            check(lf_an > 0.0 && rel(lf_an * q, dc) < 0.03,
                  "anisotropic clay with suDSS = (suA + suP)/2 carries (2 + pi)(suA + suP)/2, within 3%");
            check(lf_an > 0.0 && lf_tm > 0.0 && rel(lf_an, lf_tm) < 0.01,
                  "...and what the same mesh carries on a Tresca clay of that mean strength, within 1%");
        }
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
