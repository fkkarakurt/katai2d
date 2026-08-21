// Measurement instrument for the Hardening Soil stress-point integrator (a `study_` target:
// EXCLUDE_FROM_ALL, not wired into CTest). One run per process, so the solver's own trace
// (KATAI_NL_DEBUG) belongs to exactly one case.
//
// It is what KV-NUM-009's 2026-08-20 re-measurement was made with, and it is kept so the axes it
// sweeps stay sweepable: `point` walks the integrator alone over a fixed strain at a chosen outer
// step count and tolerance, `law` walks the corpus material's oedometer against the closed form
// with no FE at all, `triax` walks the Berlin Sand drained triaxial to its failure plateau,
// `calib` times and reports the cap calibration, and the numeric form runs the KV-CST-002
// boundary-value case over one of three stress ranges.
//
// The question it was written for: KV-NUM-009's stress-range signature printed 100->200 kPa as
// +0.6419% on the MKL composition and +0.4583% on portable/Eigen, while the other two ranges
// agreed bit-for-bit. Accumulated round-off does not behave like that; a discrete event in one
// run does. It was the substep subdivision, which had no tolerance of its own -- the two
// compositions now agree to fifteen significant figures on all three ranges.
//
// usage: study_hs_integration <range 0|1|2> <tolerance> <staged load steps> [seating steps]
#include <katai/materials/hardening_soil_plastic.hpp>
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace m = katai::model;

// Material-point mode: the integrator alone, no FE, which is the only place the error control's
// own behaviour can be watched without a run taking minutes.
// usage: study_hs_integration point [outer steps] [stol]
// The SAME total oedometric strain is walked in N outer calls at a fixed integration tolerance.
// If the integrator's substepping controls the error, the answer must not depend on N -- that is
// the whole claim an error-controlled scheme makes. Any dependence left is error the substeps
// cannot see, and the outer step is the only other axis in this walk.
static int point_mode(int nouter, double stol) {
    using katai::core::HardeningSoilParams;
    using katai::core::HsIntegrated;
    HardeningSoilParams p;
    p.p_ref = 100.0; p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.nu_ur = 0.2; p.friction = 30.0 * 3.14159265358979323846 / 180.0;
    p.cap_alpha = 0.6; p.cap_beta = 2.0e-3; p.Rf = 0.9;

    const double eps_total = 0.02;          // fixed path: 2% vertical strain, lateral restrained
    const double de1 = eps_total / nouter;
    Eigen::Vector3d sig(5.0, 5.0, 5.0);
    double gp = 0.0, pp = 5.0;
    long long acc = 0, rej = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < nouter; ++i) {
        const HsIntegrated r = katai::core::hs_integrate(p, sig, gp, pp,
                                                        Eigen::Vector3d(de1, 0.0, 0.0), stol);
        acc += r.nsub; rej += r.saturated;
        if (!r.stress.allFinite()) { std::printf("  NON-FINITE at step %d\n", i); return 2; }
        sig = r.stress; gp = r.gamma_p; pp = r.pp;
    }
    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    std::printf("point: outer=%-5d stol=%-8g sigma1=%.9f  K0=%.9f  gamma_p=%.9e  "
                "substeps=%lld rejected=%lld (%.1f per call)  [%.0f ms]\n",
                nouter, stol, sig(0), sig(2) / sig(0), gp, acc, rej,
                (double)acc / nouter, ms);
    return 0;
}

// usage: study_hs_integration law [stol] [strain step]
// KV-CST-002's material, walked as an oedometer at the STRESS POINT, with no FE anywhere. The
// closed form the case is judged against assumes the model reproduces Eoed_ref (sigma/p_ref)^m;
// this asks whether it does, so that the deviation the FE run reports can be attributed to the
// model or to the FE.
static int law_mode(double stol, double de1) {
    using katai::core::HardeningSoilParams;
    using katai::core::HsIntegrated;
    HardeningSoilParams p;                       // tests/corpus/kv-cst-002-hs-oedometer.k2d
    p.p_ref = 100.0; p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.nu_ur = 0.2; p.Rf = 0.9;
    p.friction = 35.0 * 3.14159265358979323846 / 180.0;
    p.dilatancy = 5.0 * 3.14159265358979323846 / 180.0;
    p.cohesion = 0.0;
    katai::core::hs_calibrate_cap(p, 1.0 - std::sin(p.friction));   // k0nc_auto in the file

    const double s0 = 0.5;
    Eigen::Vector3d sig(s0, s0, s0);
    double gp = 0.0, pp = s0, eps = 0.0;
    const double marks[5] = {50.0, 100.0, 200.0, 400.0, 800.0};
    double eps_at[5] = {0, 0, 0, 0, 0};
    int next = 0;
    for (int i = 0; i < 4000000 && next < 5; ++i) {
        const HsIntegrated r = katai::core::hs_integrate(p, sig, gp, pp,
                                                        Eigen::Vector3d(de1, 0.0, 0.0), stol);
        if (!r.stress.allFinite()) { std::printf("  NON-FINITE\n"); return 2; }
        while (next < 5 && r.stress(0) >= marks[next]) {
            // linear in strain between the two samples, so the mark is not a step-size artefact
            const double f = (marks[next] - sig(0)) / (r.stress(0) - sig(0));
            eps_at[next++] = eps + f * de1;
        }
        sig = r.stress; gp = r.gamma_p; pp = r.pp; eps += de1;
    }
    std::printf("law: stol=%-8g de=%-8g alpha=%.4f beta=%.4e  K0=%.4f\n",
                stol, de1, p.cap_alpha, p.cap_beta, sig(2) / sig(0));
    const double h = 4.0;   // the case's column height, so the numbers read as settlements
    for (int i = 0; i + 1 < 5; ++i) {
        const double a = marks[i], b = marks[i + 1];
        const double got = h * (eps_at[i + 1] - eps_at[i]);
        const double want = h * std::pow(p.p_ref, p.m) / p.Eoed_ref *
                            (std::pow(b, 1.0 - p.m) - std::pow(a, 1.0 - p.m)) / (1.0 - p.m);
        std::printf("   %3.0f -> %3.0f kPa   model %.9f m   closed form %.9f m   %+.4f%%\n",
                    a, b, got, want, 100.0 * (got - want) / want);
    }
    return 0;
}

// usage: study_hs_integration triax [max axial strain] [stol]
// test_hs_berlin's drained triaxial, with the strain limit opened up: the shipped check asks
// whether q has reached the qf plateau by 5% axial strain, and the question that decides whether
// the check or the model is wrong is whether it reaches it AT ALL.
static int triax_mode(double eps_max, double stol) {
    using katai::core::HardeningSoilParams;
    using katai::core::HsIntegrated;
    HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 105e3; p.Eur_ref = 315e3; p.Eoed_ref = 105e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38 * 3.14159265358979323846 / 180;
    p.dilatancy = 6 * 3.14159265358979323846 / 180; p.cohesion = 1.0; p.Rf = 0.9;
    katai::core::hs_calibrate_cap(p, 0.38);
    const double sigma3 = 200.0;
    const double qf = p.q_failure(sigma3);
    Eigen::Vector3d sig(sigma3, sigma3, sigma3);
    double gp = 0, pp = sigma3, eps1 = 0, epsv = 0;
    double q_at_5 = 0.0, epsv_at_5 = 0.0;
    for (int s = 0; s < 200000; ++s) {
        double dlat = 0.0; HsIntegrated r;
        for (int it = 0; it < 30; ++it) {
            r = katai::core::hs_integrate(p, sig, gp, pp,
                                          Eigen::Vector3d(5e-6, dlat, dlat), stol);
            const double g = r.stress(2) - sigma3;
            if (std::fabs(g) <= 1e-8 * (1 + sigma3)) break;
            dlat -= g / (r.tangent(2, 1) + r.tangent(2, 2));
        }
        epsv += 5e-6 + 2.0 * dlat;
        sig = r.stress; gp = r.gamma_p; pp = r.pp; eps1 += 5e-6;
        if (q_at_5 == 0.0 && eps1 >= 0.05) { q_at_5 = sig(0) - sig(2); epsv_at_5 = epsv; }
        if (eps1 > eps_max) break;
    }
    const double q = sig(0) - sig(2);
    std::printf("triax: stol=%-8g qf=%.2f | at 5%%: q=%.2f (%.3f%% of qf) eps_v=%+.5f"
                " | at %.1f%%: q=%.2f (%.3f%% of qf) eps_v=%+.5f\n",
                stol, qf, q_at_5, 100.0 * q_at_5 / qf, epsv_at_5,
                100.0 * eps_max, q, 100.0 * q / qf, epsv);
    return 0;
}

// Calibration mode: what the cap calibration alone costs, and what it produces. It is a pure
// function of the material, so if it is slow it is slow once per material -- but "once" is per
// solve today, which is why it is timed here.
static int calib_mode() {
    using katai::core::HardeningSoilParams;
    HardeningSoilParams p;
    p.p_ref = 100.0; p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.nu_ur = 0.2; p.friction = 30.0 * 3.14159265358979323846 / 180.0;
    p.Rf = 0.9;
    const auto t0 = std::chrono::steady_clock::now();
    katai::core::hs_calibrate_cap(p, 1.0 - std::sin(p.friction));
    const double ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t0).count();
    double e = 0.0, k0 = 0.0;
    katai::core::hs_oedometer_probe(p, e, k0);
    std::printf("calib: alpha=%.9f beta=%.9e  ->  Eoed(p_ref)=%.6f (target %.1f)  "
                "K0=%.6f (target %.6f)  [%.0f ms]\n",
                p.cap_alpha, p.cap_beta, e, p.Eoed_ref, k0, 1.0 - std::sin(p.friction), ms);
    return 0;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::string(argv[1]) == "point")
        return point_mode(argc > 2 ? std::atoi(argv[2]) : 400,
                          argc > 3 ? std::atof(argv[3]) : 0.0);
    if (argc > 1 && std::string(argv[1]) == "law")
        return law_mode(argc > 2 ? std::atof(argv[2]) : 0.0,
                        argc > 3 ? std::atof(argv[3]) : 2.0e-5);
    if (argc > 1 && std::string(argv[1]) == "triax")
        return triax_mode(argc > 2 ? std::atof(argv[2]) : 0.20,
                          argc > 3 ? std::atof(argv[3]) : 0.0);
    if (argc > 1 && std::string(argv[1]) == "calib") return calib_mode();
    const int idx = argc > 1 ? std::atoi(argv[1]) : 1;
    const double tol = argc > 2 ? std::atof(argv[2]) : 1e-6;
    const int steps = argc > 3 ? std::atoi(argv[3]) : 160;
    const int seat = argc > 4 ? std::atoi(argv[4]) : 40;
    const int window = argc > 5 ? std::atoi(argv[5]) : 0;

    const double ranges[3][2] = {{50.0, 100.0}, {100.0, 200.0}, {200.0, 400.0}};

    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-cst-002-hs-oedometer.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        std::printf("FAIL load: %s\n", err.c_str());
        return 1;
    }
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("FAIL mesh: %s\n", M.message.c_str()); return 1; }

    pr.loads[0].qy1 = pr.loads[0].qy2 = -ranges[idx][0];
    pr.loads[1].qy1 = pr.loads[1].qy2 = -(ranges[idx][1] - ranges[idx][0]);
    pr.initial.load_steps = seat;  pr.initial.tolerance = tol;  pr.initial.max_iterations = 500;
    pr.phases[0].load_steps = steps; pr.phases[0].tolerance = tol; pr.phases[0].max_iterations = 500;

    katai::app::NumericalControls nc;
    nc.line_search_window = window;
    const auto res = katai::app::solve_phases(
        pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure), nullptr, nullptr, nc);
    if (res.size() != 2 || !res[1].ok) {
        std::printf("range=%d tol=%g steps=%d  DID NOT SOLVE\n", idx, tol, steps);
        return 2;
    }
    const auto top_of = [](const katai::core::SolveResult& r) {
        int top = 0;
        double best = 1e300;
        for (int n = 0; n < r.mesh.node_count; ++n) {
            const double d = std::hypot(r.mesh.x[n] - 0.5, r.mesh.y[n] - 4.0);
            if (d < best) { best = d; top = n; }
        }
        return -r.disp[top * 2 + 1];
    };
    const double u_seat = top_of(res[0]);   // where the divergence is born, if it is born early
    const double u = top_of(res[1]);

    // The closed form of the HS oedometric law over this range (c = 0, so sin(phi) cancels).
    const double kH = 4.0, kEoed = 30000.0, kPref = 100.0, kM = 0.5;
    const double want = std::pow(kPref, kM) / kEoed *
                        (std::pow(ranges[idx][1], 1.0 - kM) - std::pow(ranges[idx][0], 1.0 - kM)) /
                        (1.0 - kM) * kH;

    std::printf("range=%d %3.0f->%3.0f  tol=%g steps=%d seat=%d  u=%.15e  dev=%+.6f%%  "
                "iters(seat/staged)=%d/%d  lf=%.6f  window=%d  u_seat=%.15e\n",
                idx, ranges[idx][0], ranges[idx][1], tol, steps, seat, u,
                100.0 * (u - want) / want, res[0].iterations, res[1].iterations,
                res[1].load_factor, window, u_seat);
    return 0;
}
