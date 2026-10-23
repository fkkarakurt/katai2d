// test_sliding_block_dynamic -- the permanent slip of a block on a frictional joint under a ground
// acceleration pulse, against Newmark's sliding block.
//
// The nonlinear dynamic phase had been checked against its own limits only (the linear path for an
// elastic soil, the quasi-static path for a slow load) until 2026-10-22: nothing tied the
// displacement it accumulates during shaking to an answer known in closed form. That displacement
// is what a seismic slope or retaining-wall assessment reads.
//
// Oracle (closed form; Newmark, N.M. (1965). Effects of earthquakes on dams and embankments.
// Geotechnique 15(2), 139-160): a rigid block resting on a horizontal joint of friction angle
// phi_i, with no cohesion, follows the ground until the ground's acceleration exceeds the joint's
// yield acceleration N = g tan(phi_i), and then slides. For a rectangular pulse of acceleration A
// lasting t0 the block's slip relative to the ground is
//     d = V^2 / (2 N) (1 - N/A),   V = A t0,
// reached at t0 + (A - N) t0 / N, after which the block is at rest on the ground again. The finite
// element model is a stiff base layer fixed at its bottom, a stiff block on it, and an interface
// between the two with the joint's friction. Rayleigh damping of 5% is aimed at the block's own
// frequency, so that its elastic sway dies out; its mass-proportional part is negligible, because
// that part acts on the velocity relative to the base -- for a sliding block, the slip rate -- as a
// viscous drag the closed form does not have. Undamped, the stiff block's ringing ratchets the joint
// on after the block should have stopped (measured +10..+25% of the slip).
//
// THE CASE FOUND TWO DEFECTS (2026-10-22), both now fixed:
//   * the stiffness-proportional damping was built on the full initial stiffness, joint springs
//     included, and a slipping joint met a viscous resistance far above its strength: 5% damping
//     left 0.8% of the closed-form slip. Contact springs are now left out of that matrix;
//   * the soil column above a node on the model's right boundary weighed nothing (a point on a
//     right or top edge is outside every region by the half-open point-in-polygon rule), so the
//     joint's end station started with no normal stress and no friction: the slip came out 8.6%
//     high on 6-noded elements and 3.9% on 15-noded ones, the share of the two end stations.
//
// verify: KV-DYN-004
//   oracle:   closed_form
//   source:   Newmark, N.M. (1965). Effects of earthquakes on dams and embankments. Geotechnique 15(2), 139-160 -- the rigid block on a frictional plane under a rectangular acceleration pulse: d = V^2/(2N) (1 - N/A), N = g tan(phi) for a horizontal plane
//   locator:  the horizontal displacement of the block's top relative to the fixed base at the end of the phase (the block at rest again), against d evaluated in this file
//   quantity: permanent slip of the block [m], for two pulse amplitudes on 6- and 15-noded elements; a pulse below the yield acceleration
//   expected: the slip within 1% of the closed form for A = 2N and A = 3N on both element orders; a pulse below N leaves no slip (under 1% of the 2N slip)
//   band:     1% -- a stiff, not rigid, block and base, and a time step of t0/200; measured -0.06% / -0.46% on both orders

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

namespace m = katai::model;
using katai::app::InitialPhase;

static int g_failures = 0;
static void check(bool ok, const std::string& what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what.c_str());
    if (!ok) ++g_failures;
}

struct Slip { bool ok = false; double d = 0.0; std::string msg; };

static Slip run(double A, double t0, double phi_i_deg, int order) {
    constexpr double W = 4.0, hb = 1.0, hk = 1.0, E = 2.0e7, nu = 0.2, gamma = 20.0;
    constexpr double dt = 1e-3, T = 0.8;
    m::Project pr;
    m::Material base; base.name = "Base"; base.model = m::SoilModel::LinearElastic;
    base.E = E; base.nu = nu; base.gamma_unsat = base.gamma_sat = gamma; base.e_init = 0.5;
    // The joint takes the strength of the soil named for it: c = 0, phi = phi_i, Rinter = 1.
    m::Material joint = base; joint.name = "Joint"; joint.c = 0.0; joint.phi = phi_i_deg; joint.psi = 0.0;
    joint.rinter_rigid = false; joint.Rinter = 1.0;
    m::Material block = base; block.name = "Block";
    pr.materials = {base, block, joint};
    m::SoilPolygon Pb; Pb.name = "Base"; Pb.material = 0;
    Pb.x = {0, W, W, 0}; Pb.y = {0, 0, hb, hb};
    Pb.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free};
    m::SoilPolygon Pk; Pk.name = "Block"; Pk.material = 1;
    Pk.x = {0, W, W, 0}; Pk.y = {hb, hb, hb + hk, hb + hk};
    Pk.edge_bc = {(int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free};
    pr.polygons = {Pb, Pk};
    pr.has_water = false;
    m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Joint";
    j.x1 = 0; j.y1 = hb; j.x2 = W; j.y2 = hb; j.iface_material = 2;
    pr.structs.push_back(j);
    m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Pulse";
    p.seismic_wave = m::SeismicWave::Record; p.seismic_amp = 1.0; p.record_dt = dt;
    const int n_pulse = (int)std::lround(t0 / dt);
    p.accel_record.assign((size_t)std::lround(T / dt) + 1, 0.0);
    for (int i = 0; i < n_pulse; ++i) p.accel_record[(size_t)i] = A;
    // 5% damping aimed at the block's own frequency (~500 Hz) with a negligible mass-proportional
    // part (f1 = 0.01 Hz): the mass term acts on the velocity relative to the base, which for a
    // sliding block is the slip rate itself.
    p.damping_ratio = 0.05; p.rayleigh_f1 = 0.01; p.rayleigh_f2 = 500.0;
    p.duration = T; p.time_steps = (int)std::lround(T / dt);
    p.dynamic_nonlinear = true;
    pr.phases.push_back(p);
    Slip s;
    const auto M = katai::app::mesh_from_project(pr, 0.5, order);
    if (!M.ok) { s.msg = M.message; return s; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2 || !res[1].ok) { s.msg = res.empty() ? "no result" : res.back().message; return s; }
    const auto& R = res[1];
    // mean ux over the block's top minus mean ux over the base's bottom (it moves with the ground)
    double top = 0.0, bot = 0.0; int nt = 0, nb = 0;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        const double y = R.mesh.y[n], ux = R.disp[2 * n];
        if (std::fabs(y - (hb + hk)) < 1e-9) { top += ux; ++nt; }
        if (std::fabs(y) < 1e-9) { bot += ux; ++nb; }
    }
    if (nt == 0 || nb == 0) { s.msg = "no nodes found"; return s; }
    s.ok = true;
    s.d = std::fabs(top / nt - bot / nb);
    return s;
}

int main() {
    std::printf("Newmark's sliding block through the nonlinear dynamic phase\n\n");
    const double g = 9.81, phi = 15.0, t0 = 0.2;
    const double N = g * std::tan(phi * 3.14159265358979323846 / 180.0);
    for (int order : {6, 15}) {
        double d2 = 0.0;
        for (double ratio : {2.0, 3.0}) {
            const double A = ratio * N;
            const double V = A * t0;
            const double want = V * V / (2.0 * N) * (1.0 - N / A);
            const Slip s = run(A, t0, phi, order);
            char label[160];
            std::snprintf(label, sizeof(label), "%d-noded, A = %.0fN = %.3f m/s2 for %.2f s", order, ratio, A, t0);
            if (!s.ok) { check(false, std::string(label) + ": solved (" + s.msg + ")"); continue; }
            std::printf("  %s: slip %.5f m, Newmark %.5f m (%+.2f%%)\n", label, s.d, want, 100.0 * (s.d - want) / want);
            check(std::fabs(s.d - want) / want < 0.01, std::string(label) + ": the permanent slip is Newmark's within 1%");
            if (ratio == 2.0) d2 = s.d;
        }
        const Slip s = run(0.7 * N, t0, phi, order);
        std::printf("  %d-noded, A = 0.7N: slip %.3e m\n", order, s.d);
        check(s.ok && d2 > 0.0 && s.d < 0.01 * d2, std::to_string(order) + "-noded: a pulse below the yield acceleration leaves no slip");
    }
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
