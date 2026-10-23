// test_hoek_brown_safety -- strength reduction of the Hoek-Brown criterion, and a factor of safety
// with a Hoek-Brown rock.
//
// Until 2026-10-22 a Safety phase on a model with a Hoek-Brown material was refused: the criterion
// has no c' and phi' to divide. Strength reduction is defined for it inside the criterion (Benz,
// Schwab, Kauther & Vermeer 2008): at each sigma'_3 the envelope's tangent is a Mohr-Coulomb line
// (c_t, phi_t), reduced as Mohr-Coulomb is reduced, c_t / Msf and tan(phi_t) / Msf, and the reduced
// envelope passes through the reduced line at that sigma'_3. The code implements the published
// closed form f_red = f / eta(f', Msf) in a rationalised version.
//
// Oracle (closed form, constructed independently): the tangent of the envelope at a compressive
// confinement sigma_3 (compression positive) has the slope K_p = 1 + g, g = -d fbar / d sigma'_3,
// hence sin(phi_t) = g / (2 + g), and the intercept 2 c_t sqrt(K_p) = fbar - g sigma_3; the reduced
// line evaluated at the same sigma_3 is the reduced envelope. (1) The code's reduced envelope and
// its slope against that construction, for two rock masses, across the confining range, at
// several factors. (2) A homogeneous block of rock under a confining stress p3 and a major stress
// p1 has its factor of safety where the reduced envelope passes through (p3, p1): found here by
// bisection on the construction, and by the Safety phase on the finite element model.
//
// verify: KV-CST-021
//   oracle:   closed_form
//   source:   Benz, T., Schwab, R., Kauther, R.A. & Vermeer, P.A. (2008). A Hoek-Brown criterion with intrinsic material strength factorization. Int. J. Rock Mech. Min. Sci. 45(2), 210-222 -- strength reduction of the Hoek-Brown envelope through its local Mohr-Coulomb tangent; the tangent construction evaluated in this file
//   locator:  the reduced envelope fbar_red(sigma'_3) and its slope against the reduced tangent line at the same sigma'_3; the factor of safety of a homogeneous block under (p3, p1) against the factor at which the reduced envelope passes through that state
//   quantity: reduced deviator at failure [kPa]; its slope [-]; factor of safety [-]
//   expected: the reduced envelope equal to the construction within 1e-10 relative, Msf = 1 bit for bit the criterion; the slope within 1e-6 of a central difference; the factor of safety within 1% of the construction's, and a Safety phase that is no longer refused
//   band:     1e-10 / 1e-6 (arithmetic); 1% on the factor of safety (the search's resolution on a block that fails all at once)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/materials/hoek_brown.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
namespace hb = katai::core::hoekbrown;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

// The reduced deviator at compressive confinement s3c > 0 (compression positive), by the tangent
// construction alone -- no eta, no code from hoek_brown.hpp beyond the unreduced envelope.
static double reduced_by_tangent(double s3c, double M, const hb::Params& P, const hb::Constants& C) {
    const double sci = std::fabs(P.sigci);
    const double B = C.mb * (s3c / sci) + C.s;
    const double f = sci * std::pow(B, C.a);
    const double g = C.a * C.mb * std::pow(B, C.a - 1.0);   // d f / d s3c
    const double Kp = 1.0 + g;
    const double phi = std::asin(g / (2.0 + g));
    const double c = (f - g * s3c) / (2.0 * std::sqrt(Kp));
    const double phr = std::atan(std::tan(phi) / M), cr = c / M;
    const double Kpr = (1.0 + std::sin(phr)) / (1.0 - std::sin(phr));
    return (Kpr - 1.0) * s3c + 2.0 * cr * std::sqrt(Kpr);
}

int main() {
    std::printf("Hoek-Brown strength reduction, and a factor of safety in rock\n\n");
    // (1) the reduced envelope against the construction
    {
        double worst = 0.0, worst_slope = 0.0;
        bool identity = true;
        const double rocks[2][4] = {{30000.0, 12.0, 45.0, 0.0}, {80000.0, 25.0, 75.0, 0.5}};   // sigci, mi, GSI, D
        for (const auto& r : rocks) {
            hb::Params P; P.sigci = r[0]; P.mi = r[1]; P.gsi = r[2]; P.D = r[3];
            const hb::Constants C = hb::constants_of(P);
            for (double s3c : {0.0, 10.0, 200.0, 1500.0, 8000.0}) {
                const double x = -s3c;
                identity = identity && hb::fbar(x, P, C) == hb::fbar_full(x, P, C);
                for (double M : {1.1, 1.5, 2.5, 4.0}) {
                    hb::Params Q = P; Q.msf = M;
                    const double got = hb::fbar(x, Q, C), want = reduced_by_tangent(s3c, M, P, C);
                    worst = std::max(worst, std::fabs(got - want) / want);
                    const double h = 1e-4 * (1.0 + s3c);
                    const double fd = (hb::fbar(x + h, Q, C) - hb::fbar(x - h, Q, C)) / (2.0 * h);
                    if (s3c > 0.0)   // at sigma_3 = 0 the central difference straddles the uni-axial point
                        worst_slope = std::max(worst_slope, std::fabs(hb::dfbar_dx(x, Q, C) - fd) / std::fabs(fd));
                }
            }
        }
        std::printf("  reduced envelope against the tangent construction: worst %.2e; slope against a central difference: worst %.2e\n",
                    worst, worst_slope);
        check(identity, "Msf = 1 is the criterion bit for bit");
        check(worst < 1e-10, "the reduced envelope is the reduced tangent line at every confinement and factor");
        check(worst_slope < 1e-6, "its slope is the derivative of the reduced envelope");
    }

    // (2) a homogeneous block under (p3, p1), and its factor of safety
    {
        constexpr double W = 2.0, H = 2.0, p3 = 400.0, Msf_target = 1.6;
        m::Material rock; rock.name = "Rock"; rock.model = m::SoilModel::HoekBrown;
        rock.E = 2.0e6; rock.nu = 0.3; rock.sig_ci = 20000.0; rock.mi = 10.0; rock.gsi = 50.0; rock.hb_D = 0.0;
        rock.psi = 0.0; rock.gamma_unsat = rock.gamma_sat = 0.0; rock.e_init = 0.5; rock.tension_cutoff = false;
        hb::Params P; P.sigci = rock.sig_ci; P.mi = rock.mi; P.gsi = rock.gsi; P.D = rock.hb_D;
        const hb::Constants C = hb::constants_of(P);
        // the major stress at which the reduced envelope (at Msf_target) passes through the state
        const double p1 = p3 + reduced_by_tangent(p3, Msf_target, P, C);
        // and, as a check of the bisection the test relies on, the factor recovered from it
        double lo = 1.0, hi = 5.0;
        for (int it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (p3 + reduced_by_tangent(p3, mid, P, C) > p1) lo = mid; else hi = mid;
        }
        const double want = 0.5 * (lo + hi);
        m::Project pr;
        pr.materials.push_back(rock);
        m::SoilPolygon B; B.name = "Rock"; B.material = 0;
        B.x = {0, W, W, 0}; B.y = {0, 0, H, H};
        B.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::Free, (int)m::BCType::Free,
                     (int)m::BCType::HorizontallyFixed};
        pr.polygons.push_back(B);
        pr.has_water = false;
        m::Load top; top.kind = m::LoadKind::Distributed; top.name = "p1";
        top.x1 = 0; top.y1 = H; top.x2 = W; top.y2 = H; top.qx1 = top.qx2 = 0.0; top.qy1 = top.qy2 = -p1;
        m::Load side; side.kind = m::LoadKind::Distributed; side.name = "p3";
        side.x1 = W; side.y1 = 0; side.x2 = W; side.y2 = H; side.qx1 = side.qx2 = -p3; side.qy1 = side.qy2 = 0.0;
        pr.loads = {top, side};
        pr.initial.load_active = {0, 0};
        m::Phase load; load.name = "Load"; load.type = m::PhaseType::Plastic; load.load_active = {1, 1}; load.load_steps = 10;
        m::Phase safety; safety.name = "Safety"; safety.type = m::PhaseType::Safety; safety.load_active = {1, 1};
        pr.phases = {load, safety};
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 3 && res[1].ok && res[2].ok && res[2].fos > 0.0;
        check(ok, "a Safety phase on a Hoek-Brown rock runs (it was refused)");
        if (!ok) {
            if (!res.empty()) std::printf("  (%s)\n", res.back().message.substr(0, 300).c_str());
        } else {
            bool noted = false;
            for (const auto& d : res[2].diagnostics) noted = noted || d.code == "K2D-A037";
            std::printf("  block under p3 = %.0f, p1 = %.1f kPa: factor of safety %.4f%s, construction %.4f (%+.2f%%)\n",
                        p3, p1, res[2].fos, res[2].fos_lower_bound ? " (lower bound)" : "", want,
                        100.0 * (res[2].fos - want) / want);
            check(!res[2].fos_lower_bound && std::fabs(res[2].fos - want) / want < 0.01,
                  "the factor of safety is the one at which the reduced envelope passes through the state");
            check(noted, "and the phase says how the Hoek-Brown strength was reduced (K2D-A037)");
        }
    }
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
