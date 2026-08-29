// A circular tunnel unloaded into a Hoek-Brown rock mass, against the closed-form solution.
//
// This case exists because it is the one that FOUND three defects the material-point tests could
// not see, and it is here so that they cannot come back. All three had the same shape -- a model
// that had been made reachable, in a place where its name was not written:
//
//   1. `integrate_point_axisym` had no Hoek-Brown branch and no `default:`, so a rock material in
//      axisymmetry left the trial stress and the tangent exactly as the caller passed them. Now
//      also guarded generally by test_material_axisym_coverage.
//   2. The tangent handed to the solver was the ELASTIC operator. That integrates the material
//      correctly and every material-point test passes with it; it does not iterate an
//      ill-conditioned boundary value problem to equilibrium. The phase stalled once the plastic
//      annulus passed about 4% of the tunnel radius, and load steps did not buy it back.
//   3. The solver's hybrid-tangent retry -- fail on the cheap tangent, re-enter THAT increment
//      with the stronger one -- was armed by a flag that tested for one model BY NAME, so the
//      consistent tangent added in (2) was never reached. The tell was that the run stalled at
//      the identical load factor with the new tangent as without it, and at 40 steps as at 150.
//
// So what this test asserts is deliberately two things: that the unloading CONVERGES at all (the
// regression guard for 2 and 3), and that where it converges it lands on the closed form (the
// verification). Either one alone would have passed at some point during the repair.
//
// THE MODEL IS A RADIAL STRIP IN AXISYMMETRY. In axisymmetry the three normal stresses are
// (sigma_r, sigma_z, sigma_theta); map the tunnel onto it with r the radius, theta the hoop and z
// the TUNNEL AXIS, and holding u_z = 0 on both horizontal edges is exactly the tunnel's
// plane-strain condition. The cross-section is then one-dimensional and needs no polygonal
// approximation of a circle.
//
// verify: KV-CST-016
//   oracle:   closed_form
//   source:   the elasto-plastic solution for a circular opening in a Hoek-Brown medium under hydrostatic in-situ stress (Carranza-Torres & Fairhurst 1999, IJRMMS 36(6) 777-809; generalised a >= 1/2 in Carranza-Torres 2004, IJRMMS 41(S1) 629-639), derived from radial equilibrium in this file rather than quoted, and checked by reducing to their scaled forms at a = 1/2. Rock mass and stress state: case D1 of Table A1.1 in Hoek, Carranza-Torres, Diederichs & Corkum (2008), Integration of geotechnical and structural design in tunnelling, 56th Annual Geotechnical Engineering Conference, University of Minnesota
//   locator:  ln(r/R) = 1/(m_b(1-a)) [ (m_b sigma_r/sigma_ci + s)^(1-a) - (m_b p_i/sigma_ci + s)^(1-a) ] in the plastic annulus, with the elastic-plastic boundary stress from 2(sigma_0 - sigma_r_ep) = sigma_ci (m_b sigma_r_ep/sigma_ci + s)^a and the elastic branch sigma_r = sigma_0 - (sigma_0 - sigma_r_ep)(R_pl/r)^2
//   quantity: the radial effective stress profile sigma_r(r) around a tunnel unloaded to a support pressure below the critical one [kPa], and the critical support pressure itself
//   expected: the closed-form profile through both the plastic annulus and the elastic field, and yielding that starts at the critical pressure and not before
//   band:     1.5% on the radial stress on both sides of the elastic-plastic boundary -- MEASURED at 0.59% inside the plastic annulus and 0.72% outside it, on 0.5 m tri6 elements with the outer boundary at 40 R. The wall node is reported but NOT asserted: there the closed-form stress is the applied support pressure itself and a recovered nodal stress at a free boundary is the weakest number a displacement formulation produces (measured -0.10%, which is 0.003 MPa on a 28 MPa problem). Refining to 0.4 m in the published run of the same problem (python/examples/rock_tunnel_ground_reaction.py) moves the wall closure by 0.04% and improves the plastic-annulus stress, which is what a discretisation error does and a model error does not
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// Case D1 of Table A1.1: a graphitic phyllite at ~1100 m depth (Yacambu-Quibor).
constexpr double kSigCi = 50.0e3;    // kPa, INTACT rock
constexpr double kMi = 7.0, kGsi = 48.0, kD = 0.0;
constexpr double kE = 7.5e6, kNu = 0.25;
constexpr double kSig0 = 28.0e3;     // kPa, hydrostatic in situ
constexpr double kR = 2.5;           // m, tunnel radius
// The outer radius has to be far enough that the elastic field outside the plastic annulus is
// effectively the infinite one the closed form assumes. It is also a CONVERGENCE parameter, which
// is less obvious: at 20 R this same problem, on a FINER mesh, could not be taken past a support
// pressure of 4 MPa, while at 40 R a coarser one walks down to 0.06 MPa. A pressure boundary that
// close leaves the annulus breathing against a spring the solution does not have.
constexpr double kB = 100.0;         // m, outer radius (40 R)
constexpr double kH = 2.0;           // m, strip height

double mb_of() { return kMi * std::exp((kGsi - 100.0) / (28.0 - 14.0 * kD)); }
double s_of() { return std::exp((kGsi - 100.0) / (9.0 - 3.0 * kD)); }
double a_of() { return 0.5 + (std::exp(-kGsi / 15.0) - std::exp(-20.0 / 3.0)) / 6.0; }

// Radial stress at the elastic-plastic boundary: 2(sigma_0 - x) = sigma_ci (m_b x/sigma_ci + s)^a.
// It is also the critical support pressure, because there the plastic radius equals R.
double sigma_r_ep() {
    const double mb = mb_of(), s = s_of(), a = a_of();
    const auto g = [&](double x) {
        return 2.0 * (kSig0 - x) - kSigCi * std::pow(mb * x / kSigCi + s, a);
    };
    double lo = 0.0, hi = kSig0;
    if (g(lo) < 0.0) return -1.0;      // never yields, even unsupported
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (g(mid) > 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

double plastic_radius(double p_i) {
    const double mb = mb_of(), s = s_of(), a = a_of(), xep = sigma_r_ep();
    if (xep < 0.0 || p_i >= xep) return kR;
    const double k = 1.0 / (mb * (1.0 - a));
    return kR * std::exp(k * (std::pow(mb * xep / kSigCi + s, 1.0 - a) -
                              std::pow(mb * p_i / kSigCi + s, 1.0 - a)));
}

// Compression positive.
double sigma_r_closed(double r, double p_i) {
    const double mb = mb_of(), s = s_of(), a = a_of();
    const double rpl = plastic_radius(p_i), xep = sigma_r_ep();
    if (r >= rpl) return kSig0 - (kSig0 - xep) * (rpl / r) * (rpl / r);
    const double k = 1.0 / (mb * (1.0 - a));
    const double v = std::pow(mb * p_i / kSigCi + s, 1.0 - a) + std::log(r / kR) / k;
    return (std::pow(v, 1.0 / (1.0 - a)) - s) * kSigCi / mb;
}

// The support pressures the excavation is let down through [kPa]. Coarse while the rock is still
// elastic, closer together once the plastic annulus is open -- seven stages is enough to take the
// rock well past the critical pressure, which is the point of the regression guard, without paying
// for the twenty a published ground reaction curve wants.
//
// The MESH matters as much as the staging here, and the first version of this test got it wrong:
// at 1 m elements the plastic annulus at p_i = 4 MPa is 1.53 m thick, so it is spanned by about
// one and a half elements. The run did not stall -- it reported that the body had SOFTENED INTO A
// MECHANISM, with the stiffness parameter down to 0.025, which is what an under-resolved
// perfectly-plastic annulus genuinely becomes. The engine was right and the fixture was wrong;
// 0.5 m elements put four across it.
// Six stages stop at 6 MPa, where the plastic annulus is already 43% of the tunnel radius --
// far past the critical pressure, and far past where all three defects showed themselves (the
// worst of them stalled at R_pl/R = 1.04). Going further is what the published ground reaction
// curve in python/examples is for; a suite test buys the guard, not the whole curve.
const std::vector<double> kLevels = {20.0e3, 14.0e3, 12.0e3, 10.0e3, 8.0e3, 6.0e3, 4.5e3, 3.0e3};
constexpr double kElemSize = 0.5;

m::Project tunnel() {
    m::Project pr;
    pr.axisymmetric = true;
    m::Material rock; rock.name = "Graphitic phyllite"; rock.model = m::SoilModel::HoekBrown;
    rock.E = kE; rock.nu = kNu;
    rock.sig_ci = kSigCi; rock.mi = kMi; rock.gsi = kGsi; rock.hb_D = kD;
    rock.psi = 0.0; rock.sig_psi = 0.0;
    // Weightless in all but name: this problem is a uniform far field, not a gravity gradient,
    // and the strip is 2 m tall so any real unit weight would add 0.2 per mille of sigma_0 as an
    // axial disturbance. It is not exactly zero because the K0 initial phase wants a weight.
    rock.gamma_unsat = rock.gamma_sat = 0.001;
    // The schema's default tension cut-off is left ON at sigma_t = 0. In THIS problem it is
    // inactive either way -- every stress around a tunnel under 28 MPa of ground is compressive,
    // so the closed form derived above applies unchanged -- but it is what a user gets without
    // touching the box, so it is what the verification runs.
    rock.tension_cutoff = true;
    rock.tensile_strength = 0.0;
    pr.materials.push_back(rock);

    m::SoilPolygon P; P.material = 0;
    P.x = {kR, kB, kB, kR};
    P.y = {0.0, 0.0, kH, kH};
    // bottom and top vertically fixed => u_z = 0 => plane strain along the tunnel axis;
    // both radial faces free, the pressures act there.
    P.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::VerticallyFixed, (int)m::BCType::Free};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed};
    P.edge_head = {0, 0, 0, 0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::Load far; far.kind = m::LoadKind::Distributed; far.name = "In-situ stress";
    far.x1 = kB; far.y1 = 0.0; far.x2 = kB; far.y2 = kH;
    far.qx1 = far.qx2 = -kSig0; far.qy1 = far.qy2 = 0.0;   // inward, -r
    pr.loads.push_back(far);

    // THE SUPPORT PRESSURE THAT NEVER COMES OFF. The relief pieces below sum to
    // sigma_0 - kLevels.back(), so without this one the tunnel would start at that reduced
    // pressure and end at ZERO rather than at the last level -- the stages would be right and
    // every pressure they are compared against would be wrong. (That is exactly what the first
    // version of this fixture did, and it looked like a convergence problem: the run failed on
    // the last stage because the last stage was taking the wall to nothing.)
    m::Load residual; residual.kind = m::LoadKind::Distributed;
    residual.name = "Support that stays";
    residual.x1 = kR; residual.y1 = 0.0; residual.x2 = kR; residual.y2 = kH;
    residual.qx1 = residual.qx2 = +kLevels.back();
    residual.qy1 = residual.qy2 = 0.0;
    pr.loads.push_back(residual);

    double previous = kSig0;
    for (double target : kLevels) {
        m::Load piece; piece.kind = m::LoadKind::Distributed;
        piece.name = "Relief";
        piece.x1 = kR; piece.y1 = 0.0; piece.x2 = kR; piece.y2 = kH;
        piece.qx1 = piece.qx2 = +(previous - target);      // outward, +r
        piece.qy1 = piece.qy2 = 0.0;
        pr.loads.push_back(piece);
        previous = target;
    }

    const int nload = (int)pr.loads.size();
    pr.initial.load_active.assign(nload, 1);               // the tunnel starts supported
    for (size_t k = 0; k < kLevels.size(); ++k) {
        m::Phase ph; ph.type = m::PhaseType::Plastic;
        ph.name = "Relief";
        ph.load_active.assign(nload, 1);
        for (size_t j = 0; j <= k; ++j) ph.load_active[2 + j] = 0;   // pieces let down so far
        ph.load_steps = 50;
        ph.max_iterations = 1500;
        pr.phases.push_back(ph);
    }
    return pr;
}

}  // namespace

int main() {
    std::printf("A tunnel unloaded into a Hoek-Brown rock mass\n\n");
    const double mb = mb_of(), s = s_of(), a = a_of();
    std::printf("   m_b = %.4f   s = %.6f   a = %.4f   (the source's Table A1.1 prints "
                "1.093, 0.0031, 0.507)\n", mb, s, a);
    check(std::fabs(mb - 1.093) < 5e-4 && std::fabs(s - 0.0031) < 5e-5 &&
              std::fabs(a - 0.507) < 5e-4,
          "the rock-mass constants are the published ones to the digits printed");

    const double xep = sigma_r_ep();
    std::printf("   critical support pressure = %.3f MPa   (%.1f%% of the in-situ stress)\n",
                xep / 1000.0, 100.0 * xep / kSig0);
    check(xep > 0.0 && xep < kSig0, "the rock yields before the tunnel is fully unloaded");
    std::printf("   plastic radius at the last stage: %.3f R\n\n",
                plastic_radius(kLevels.back()) / kR);

    const auto pr = tunnel();
    const auto M = katai::app::mesh_from_project(pr, kElemSize, 6);
    if (!M.ok) { std::printf("FAIL: mesh: %s\n", M.message.c_str()); return 1; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);

    // (1) THE REGRESSION GUARD. Every stage must converge. With the elastic tangent, or with the
    // consistent tangent the hybrid retry could not reach, the run stopped partway and reported a
    // stalled Newton search -- which is also what a genuinely hard problem looks like, so the
    // assertion is on the load factor and not on the message.
    bool all_converged = res.size() == kLevels.size() + 1;
    for (size_t i = 1; i < res.size(); ++i) {
        std::printf("   stage %zu (p_i = %5.1f MPa): ok=%d  lambda=%.3f  %s\n", i,
                    kLevels[i - 1] / 1000.0, (int)res[i].ok, res[i].load_factor,
                    res[i].ok ? "" : res[i].message.c_str());
        all_converged = all_converged && res[i].ok && res[i].load_factor > 0.999;
    }
    check(all_converged, "every unloading stage reaches equilibrium (it did not, three times)");
    if (!all_converged) {
        std::fprintf(stderr, "\n%d check(s) failed\n", g_failures + 1);
        return 1;
    }

    // (2) THE VERIFICATION. The radial stress profile, through the plastic annulus and out into
    // the elastic field, against the closed form at the same support pressure.
    const double p_i = kLevels.back();
    const auto& last = res.back();
    const double rpl = plastic_radius(p_i);
    std::printf("\n   radial stress at p_i = %.1f MPa   (closed-form R_pl = %.3f R)\n",
                p_i / 1000.0, rpl / kR);
    std::printf("      r/R      FE [MPa]   closed form [MPa]     difference\n");
    double worst_plastic = 0.0, worst_elastic = 0.0;
    int n_plastic = 0, n_elastic = 0;
    for (int i = 0; i < last.mesh.node_count; ++i) {
        if (std::fabs(last.mesh.y[i] - 0.5 * kH) > 1e-6) continue;
        const double r = last.mesh.x[i];
        if (r > 6.0 * kR) continue;
        const double fe = -last.stress.stress[i](0);      // compression positive
        const double cf = sigma_r_closed(r, p_i);
        const double rel = (fe - cf) / std::max(std::fabs(cf), 1.0);
        // The wall node is excluded from the band: there the closed-form stress is the applied
        // support pressure itself and a recovered NODAL stress at a free boundary is the weakest
        // number a displacement formulation produces. It is reported, not asserted.
        if (r > kR + 1e-9) {
            if (r < rpl) { worst_plastic = std::max(worst_plastic, std::fabs(rel)); ++n_plastic; }
            else { worst_elastic = std::max(worst_elastic, std::fabs(rel)); ++n_elastic; }
        }
        std::printf("   %6.2f   %11.3f   %17.3f   %+12.2f%%%s\n", r / kR, fe / 1000.0,
                    cf / 1000.0, 100.0 * rel, r <= kR + 1e-9 ? "   (wall node)" : "");
    }
    std::printf("\n   worst deviation: %.2f%% inside the plastic annulus (%d node(s)), "
                "%.2f%% outside it (%d)\n", 100.0 * worst_plastic, n_plastic,
                100.0 * worst_elastic, n_elastic);
    // AN ASSERTION NEEDS SOMETHING TO ASSERT ON. An earlier version of this test stopped at a
    // support pressure whose plastic annulus reached only 1.427 R, while the nearest interior
    // node sat at 1.50 R -- so not one sample fell inside it, "the worst deviation in the
    // plastic annulus" was zero by vacuum, and the band below could not have failed however
    // wrong the model was. The last stage now goes far enough that the annulus contains nodes,
    // and this check is what says so.
    check(n_plastic >= 2, "the profile actually samples the plastic annulus");
    check(worst_elastic < 0.015,
          "the elastic field is the closed form's, within 1.5%");
    check(worst_plastic < 0.015,
          "and so is the plastic annulus, within 1.5%");

    if (g_failures == 0) {
        std::printf("\nOK: the tunnel converges and lands on the closed-form solution\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
