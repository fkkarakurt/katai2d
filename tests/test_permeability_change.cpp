// test_permeability_change -- consolidation with a permeability that falls as the clay compresses,
// against an independent one-dimensional column.
//
// The change of permeability (`ck`: log10(k/k0) = de/ck) is what made a published embankment
// consolidate at its published rate (the end of consolidation moved by a third when it arrived),
// and until 2026-10-22 it was checked only by an identity -- a time scaling -- that a wrong law of
// the same shape would also pass.
//
// The oracle is a finite-difference column written in this file, sharing no code with the solver:
// a linear elastic skeleton (a Mohr-Coulomb soil that stays elastic), so the strain is
// eps = (q - u) / E_oed, and the permeability
//     k = k0 10^(de/ck),    de = (1 + e0) (exp(-eps) - 1)    (the void ratio of the compressed soil),
// read, as the model reads it, from the state at the start of each time step. One-dimensional
// continuity with compressible pore water, in conservative form:
//     (1/E_oed + n/K_w) du/dt = d/dz( (k/gamma_w) du/dz ),
// backward Euler on the finite element phase's own time steps.
//
// verify: KV-CON-006
//   oracle:   independent_path
//   source:   the change of permeability log10(k/k0) = de/ck with de from the volumetric strain (KATAI 2D input contract, docs/k2d-format.md, materials[i].ck); one-dimensional consolidation of a linear skeleton (Terzaghi's equation written with a variable permeability in conservative form)
//   locator:  the finite-difference column in this file: 801 nodes, the permeability of each half-cell from the strain at the start of the step, backward Euler with the finite element phase's time steps
//   quantity: surface settlement at the end of each time step [m], with ck and (as a control) without it
//   expected: the settlement history within 1% of the column with and without ck; the ck run measurably slower than the constant one (the comparison can tell the two apart)
//   band:     1% -- a 0.25 m 6-noded mesh against an 801-node grid, the same time steps

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <katai/analysis/constants.hpp>
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

// Settlement at the end of each step of a column of thickness H, drained at the top (z = 0),
// closed at the base, loaded by q at t = 0.
static std::vector<double> fd_column(double H, double Eoed, double n_over_kw, double k0, double e0,
                                     double ck, double q, double duration, int steps, int N = 800) {
    const double gw = katai::app::kGammaWater, dz = H / N, dt = duration / steps;
    std::vector<double> u(N + 1, 0.0), eps(N + 1, 0.0), out;
    for (int st = 0; st < steps; ++st) {
        // half-cell permeabilities from the strain at the start of the step
        std::vector<double> kh(N, k0);
        if (ck > 0.0)
            for (int i = 0; i < N; ++i) {
                const double e_mid = 0.5 * (eps[i] + eps[i + 1]);
                const double de = (1.0 + e0) * std::expm1(-e_mid);
                kh[i] = k0 * std::pow(10.0, de / ck);
            }
        // Per node, with eps = (q - u)/Eoed (the load is on from the first step, which starts from
        // u = 0 and eps = 0):
        //   (q - u)/Eoed - eps_n - n/Kw (u - u_n) + dt [kr (u_{i+1} - u) - kl (u - u_{i-1})] = 0.
        std::vector<double> a(N + 1, 0.0), b(N + 1, 0.0), c(N + 1, 0.0), r(N + 1, 0.0);
        b[0] = 1.0; r[0] = 0.0;   // drained top
        for (int i = 1; i <= N; ++i) {
            const double kl = kh[i - 1] / (gw * dz * dz);
            const double kr = i < N ? kh[i] / (gw * dz * dz) : 0.0;   // closed base
            const double eps_n = eps[i], u_n = u[i];
            b[i] = -(1.0 / Eoed) - n_over_kw - dt * (kr + kl);
            a[i] = dt * kl;
            c[i] = dt * kr;
            r[i] = -(q / Eoed - eps_n + n_over_kw * u_n);
        }
        std::vector<double> cp(N + 1), rp(N + 1);
        cp[0] = c[0] / b[0]; rp[0] = r[0] / b[0];
        for (int i = 1; i <= N; ++i) {
            const double den = b[i] - a[i] * cp[i - 1];
            cp[i] = c[i] / den;
            rp[i] = (r[i] - a[i] * rp[i - 1]) / den;
        }
        u[N] = rp[N];
        for (int i = N - 1; i >= 0; --i) u[i] = rp[i] - cp[i] * u[i + 1];
        double s = 0.0;
        for (int i = 0; i <= N; ++i) {
            eps[i] = (q - u[i]) / Eoed;
            s += (i == 0 || i == N ? 0.5 : 1.0) * eps[i] * dz;
        }
        out.push_back(s);
    }
    return out;
}

int main() {
    std::printf("Consolidation with a change of permeability against an independent column\n\n");
    constexpr double W = 1.0, H = 5.0, E = 2000.0, nu = 0.3, q = 50.0, k0 = 1.0e-4, e0 = 1.0, ck = 0.1;
    constexpr double T = 2000.0;
    constexpr int steps = 80;
    const double Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
    const double n_over_kw = (e0 / (1.0 + e0)) / 2.0e6;

    auto build = [&](double ck_value) {
        m::Project pr;
        m::Material s; s.name = "Clay"; s.model = m::SoilModel::MohrCoulomb;   // elastic range: c is large
        s.E = E; s.nu = nu; s.c = 1.0e4; s.phi = 0.0; s.psi = 0.0; s.drainage = m::Drainage::Undrained;
        s.gamma_unsat = s.gamma_sat = 18.0; s.e_init = e0; s.kx = s.ky = k0; s.ck = ck_value;
        pr.materials.push_back(s);
        const int C = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                     (int)m::BCType::HorizontallyFixed};
        P.edge_flow = {C, C, Hd, C}; P.edge_head = {0, 0, H, 0};
        pr.polygons.push_back(P);
        pr.has_water = true; pr.wx = {0, W}; pr.wy = {H, H};
        m::Load L; L.kind = m::LoadKind::Distributed; L.name = "q";
        L.x1 = 0; L.y1 = H; L.x2 = W; L.y2 = H; L.qy1 = L.qy2 = -q;
        pr.loads.push_back(L);
        pr.initial.load_active = {0};
        m::Phase ph; ph.name = "Consolidation"; ph.type = m::PhaseType::Consolidation;
        ph.duration = T; ph.time_steps = steps; ph.consol_load_at_start = true; ph.load_active = {1};
        pr.phases.push_back(ph);
        return pr;
    };
    double end_const = 0.0, mid_const = 0.0, mid_ck = 0.0;
    for (double ckv : {0.0, ck}) {
        const m::Project pr = build(ckv);
        const auto M = katai::app::mesh_from_project(pr, 0.25, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const std::string tag = ckv > 0 ? "ck = 0.1" : "constant k";
        check(M.ok && res.size() == 2 && res[1].ok, tag + ": the column solves");
        if (res.size() != 2 || !res[1].ok) { if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str()); continue; }
        const auto& R = res[1];
        const std::vector<double> fd = fd_column(H, Eoed, n_over_kw, k0, e0, ckv, q, T, steps);
        double worst = 0.0;
        for (int s = 0; s < steps && s + 1 < (int)R.consol_settlement.size(); ++s) {
            const double fe = R.consol_settlement[(size_t)s + 1], ref = fd[(size_t)s];
            worst = std::max(worst, std::fabs(fe - ref) / ref);
            if (s == 0 || s == 9 || s == 19 || s == 39 || s == steps - 1)
                std::printf("  %s  t = %6.1f d   FE %.5f m   column %.5f m   (%+.3f%%)\n", tag.c_str(),
                            R.consol_time[(size_t)s + 1], fe, ref, 100.0 * (fe - ref) / ref);
        }
        bool finite = true;
        for (double v : fd) finite = finite && std::isfinite(v);
        std::printf("  %s: worst %.3f%%\n", tag.c_str(), 100.0 * worst);
        check(finite && worst < 0.01, tag + ": the settlement history follows the independent column within 1%");
        const double mid = R.consol_settlement[(size_t)steps / 4 + 1];
        if (ckv > 0) mid_ck = mid; else { mid_const = mid; end_const = R.consol_settlement.back(); }
    }
    std::printf("  at t = %.0f d: %.5f m with constant k, %.5f m with ck (final %.5f m)\n", T / 4, mid_const, mid_ck, end_const);
    check(mid_ck < 0.9 * mid_const, "the falling permeability slows the consolidation measurably (the comparison can tell)");
    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
