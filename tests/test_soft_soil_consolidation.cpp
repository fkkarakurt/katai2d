// test_soft_soil_consolidation -- a Soft Soil clay consolidating under a fill, against an
// independent one-dimensional solution of the same problem.
//
// Until 2026-10-22 consolidation with the Soft Soil model had no independent check: the corpus
// column (KV-CON-002) is linear elastic, and the Soft Soil cases are drained element paths. The
// pair is the everyday case of a soft-ground embankment, and both halves are nonlinear -- the
// clay's compressibility falls as it consolidates, so its coefficient of consolidation rises.
//
// The oracle is a finite-difference solution written in this file, sharing no code with the
// solver. A normally consolidated Soft Soil clay loaded along its K0nc line compresses as
//     eps_v = lambda* ln(sigma'_v / sigma'_v0)        (the model's primary loading line),
// and in one dimension the excess pore pressure u of a column with a drained top and a closed
// base obeys (small strains, constant permeability, compressible pore water: the volume the clay
// loses is the water it expels, less what the water's own compression takes up)
//     d eps_v / dt - (n / K_w) du/dt = -(k / gamma_w) d2u/dz2,     sigma'_v = sigma'_v0 + q - u,
// solved here with backward Euler in time, Newton per step and a fine grid. The finite element
// model is the same column -- 2 m of very permeable, stiff sand on 4 m of clay, water at the
// surface, sides and base closed, a surcharge placed at the start of a consolidation phase --
// and the sand's own elastic compression (a tenth of a millimetre) is added to the oracle.
//
// verify: KV-CON-005
//   oracle:   independent_path
//   source:   Soft Soil primary loading line eps_v = lambda* ln(p'/p'0) on a proportional (K0nc) stress path (the model's definition, KATAI 2D docs/k2d-format.md, Soft Soil lam_star); one-dimensional consolidation with a stress-dependent compressibility m_v = lambda*/sigma'_v (Terzaghi's equation with the strain written in full)
//   locator:  the finite-difference column in this file: 801 nodes over the clay, backward Euler with the finite element phase's own time steps, Newton on the strain at every node
//   quantity: surface settlement of the column at the end of each time step [m]; final settlement against the closed form sum of lambda* ln((sigma'_v0 + q) / sigma'_v0) over the clay
//   expected: the drained end state within 1% of the closed-form final settlement; the settlement history within 1.5% of the finite-difference column at every step, and its degree of consolidation (each against its own end state) within 0.005; the consolidation's end within 0.5% of the drained end state
//   band:     as asserted -- measured -0.48% on the drained end state (the initial K0 of 1 - sin(phi) is not exactly the model's own one-dimensional K0nc, so the path is not quite proportional), -0.38..-0.59% on the history and 0.06% between the two end states; the time steps are shared, so the time error is common to both

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

// ---- the oracle: the clay column by finite differences ---------------------------------------
struct Column {
    double H;            // clay thickness [m]
    double s0_top;       // effective vertical stress at the clay's top before loading [kPa]
    double gamma_eff;    // buoyant unit weight of the clay [kN/m3]
    double lambda;       // lambda*
    double k;            // permeability [m/day]
    double n_over_kw;    // n / K_w [1/kPa]
    double q;            // surcharge [kPa]
};
// Settlement of the clay at the end of each of `steps` equal steps over `duration` days.
static std::vector<double> fd_settlement(const Column& c, double duration, int steps, int N = 800) {
    const double gw = katai::app::kGammaWater;
    const double dz = c.H / N, dt = duration / steps, D = c.k / gw;
    std::vector<double> s0(N + 1), u(N + 1, 0.0), e(N + 1, 0.0);
    for (int i = 0; i <= N; ++i) s0[i] = c.s0_top + c.gamma_eff * dz * i;
    std::vector<double> out;
    for (int st = 0; st < steps; ++st) {
        const std::vector<double> un = u, en = e;
        // Newton on u at time n+1; node 0 drained (u = 0), node N closed (mirror).
        for (int it = 0; it < 50; ++it) {
            std::vector<double> a(N + 1, 0.0), b(N + 1, 0.0), cc(N + 1, 0.0), r(N + 1, 0.0);
            b[0] = 1.0; r[0] = -u[0];
            double rmax = 0.0;
            for (int i = 1; i <= N; ++i) {
                const double sv = s0[i] + c.q - u[i];
                const double eps = c.lambda * std::log(sv / s0[i]);
                const double deps_du = -c.lambda / sv;
                const double up = i < N ? u[i + 1] : u[i - 1];   // closed base: du/dz = 0
                const double lap = (up - 2.0 * u[i] + u[i - 1]) / (dz * dz);
                // residual F = (eps - en)/dt - n/Kw (u - un)/dt + D lap = 0
                const double F = (eps - en[i]) / dt - c.n_over_kw * (u[i] - un[i]) / dt + D * lap;
                rmax = std::max(rmax, std::fabs(F));
                b[i] = deps_du / dt - c.n_over_kw / dt - 2.0 * D / (dz * dz);
                a[i] = D / (dz * dz);
                cc[i] = i < N ? D / (dz * dz) : 0.0;
                if (i == N) a[i] = 2.0 * D / (dz * dz);
                r[i] = -F;
            }
            // Thomas algorithm
            std::vector<double> cp(N + 1), rp(N + 1), du(N + 1);
            cp[0] = cc[0] / b[0]; rp[0] = r[0] / b[0];
            for (int i = 1; i <= N; ++i) {
                const double den = b[i] - a[i] * cp[i - 1];
                cp[i] = cc[i] / den;
                rp[i] = (r[i] - a[i] * rp[i - 1]) / den;
            }
            du[N] = rp[N];
            for (int i = N - 1; i >= 0; --i) du[i] = rp[i] - cp[i] * du[i + 1];
            // Projected step: the effective stress stays positive (a full Newton step from the
            // unloaded state would overshoot the load the pressure can carry).
            double dmax = 0.0;
            for (int i = 0; i <= N; ++i) {
                const double cap = s0[i] + c.q - 1e-6 * (s0[i] + c.q);
                const double un1 = std::min(u[i] + du[i], cap);
                dmax = std::max(dmax, std::fabs(un1 - u[i]));
                u[i] = un1;
            }
            if (dmax < 1e-11 && rmax < 1e-9) break;
        }
        double s = 0.0;
        for (int i = 0; i <= N; ++i) {
            e[i] = i == 0 ? c.lambda * std::log((s0[0] + c.q) / s0[0])
                          : c.lambda * std::log((s0[i] + c.q - u[i]) / s0[i]);
            const double w = (i == 0 || i == N) ? 0.5 : 1.0;
            s += w * e[i] * dz;
        }
        out.push_back(s);
    }
    return out;
}

int main() {
    std::printf("Soft Soil consolidation against an independent one-dimensional column\n\n");
    const double gw = katai::app::kGammaWater;
    constexpr double W = 1.0, Hc = 4.0, Hs = 2.0, q = 40.0;
    constexpr double lam = 0.10, kap = 0.02, k_clay = 1.0e-3, e0 = 1.5;
    constexpr double g_sand = 20.0, g_clay = 16.0, E_sand = 1.0e6, nu_sand = 0.3;
    constexpr double T = 1000.0;
    constexpr int steps = 100;

    m::Project pr;
    m::Material clay; clay.name = "Clay"; clay.model = m::SoilModel::SoftSoil;
    clay.drainage = m::Drainage::Undrained;
    clay.lam_star = lam; clay.kap_star = kap; clay.nu_ur = 0.15; clay.c = 1.0; clay.phi = 25.0; clay.psi = 0.0;
    clay.gamma_unsat = g_clay; clay.gamma_sat = g_clay; clay.e_init = e0; clay.kx = clay.ky = k_clay;
    clay.oc_mode = 1; clay.OCR = 1.0;
    m::Material sand; sand.name = "Sand"; sand.model = m::SoilModel::LinearElastic;
    sand.E = E_sand; sand.nu = nu_sand; sand.gamma_unsat = g_sand; sand.gamma_sat = g_sand; sand.e_init = 0.5;
    sand.kx = sand.ky = 10.0;
    pr.materials = {clay, sand};
    const int C = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
    m::SoilPolygon Pc; Pc.name = "Clay"; Pc.material = 0;
    Pc.x = {0, W, W, 0}; Pc.y = {0, 0, Hc, Hc};
    Pc.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                  (int)m::BCType::HorizontallyFixed};
    Pc.edge_flow = {C, C, C, C}; Pc.edge_head = {0, 0, 0, 0};
    m::SoilPolygon Ps; Ps.name = "Sand"; Ps.material = 1;
    Ps.x = {0, W, W, 0}; Ps.y = {Hc, Hc, Hc + Hs, Hc + Hs};
    Ps.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                  (int)m::BCType::HorizontallyFixed};
    Ps.edge_flow = {C, C, Hd, C}; Ps.edge_head = {0, 0, Hc + Hs, 0};
    pr.polygons = {Pc, Ps};
    pr.has_water = true; pr.wx = {0, W}; pr.wy = {Hc + Hs, Hc + Hs};
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Fill";
    L.x1 = 0; L.y1 = Hc + Hs; L.x2 = W; L.y2 = Hc + Hs; L.qy1 = L.qy2 = -q;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};
    m::Phase ph; ph.name = "Consolidation"; ph.type = m::PhaseType::Consolidation;
    ph.duration = T; ph.time_steps = steps; ph.consol_load_at_start = true;
    ph.load_active = {1};
    pr.phases.push_back(ph);

    const auto M = katai::app::mesh_from_project(pr, 0.25, 6);
    check(M.ok, "the column meshes");
    if (!M.ok) return 1;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 2 && res[0].ok && res[1].ok, "the initial and the consolidation phase solve");
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str());
        return 1;
    }
    const auto& R = res[1];

    // The oracle's column.
    Column col;
    col.H = Hc;
    col.s0_top = Hs * (g_sand - gw);
    col.gamma_eff = g_clay - gw;
    col.lambda = lam;
    col.k = k_clay;
    col.n_over_kw = (e0 / (1.0 + e0)) / 2.0e6;
    col.q = q;
    const std::vector<double> fd = fd_settlement(col, T, steps);
    const double Eoed_sand = E_sand * (1 - nu_sand) / ((1 + nu_sand) * (1 - 2 * nu_sand));
    const double s_sand = q * Hs / Eoed_sand;
    // closed-form final settlement of the clay
    double s_final = 0.0;
    {
        const int N = 4000;
        const double dz = Hc / N;
        for (int i = 0; i <= N; ++i) {
            const double s0 = col.s0_top + col.gamma_eff * dz * i;
            const double w = (i == 0 || i == N) ? 0.5 : 1.0;
            s_final += w * lam * std::log((s0 + q) / s0) * dz;
        }
    }
    // The same column loaded DRAINED in a Plastic phase: the end state of the consolidation without
    // the time in between. It separates the two halves of the comparison below -- the compression
    // law, which this end state tests alone, and the coupling, which only the history can test.
    double s_drained = -1.0;
    {
        m::Project d = pr;
        d.materials[0].drainage = m::Drainage::Drained;
        d.phases[0].type = m::PhaseType::Plastic;
        d.phases[0].load_steps = 40;
        const auto rd = katai::app::solve_phases(d, M.mesh, InitialPhase::K0Procedure);
        check(rd.size() == 2 && rd[1].ok, "the drained Plastic phase solves");
        if (rd.size() == 2 && rd[1].ok)
            for (int n = 0; n < rd[1].mesh.node_count; ++n)
                if (std::fabs(rd[1].mesh.y[n] - (Hc + Hs)) < 1e-9) s_drained = std::max(s_drained, -rd[1].disp[2 * n + 1]);
    }
    const double s_closed = s_final + s_sand;
    std::printf("  drained end state: %.5f m against the closed form %.5f m (%+.3f%%)\n", s_drained, s_closed,
                100.0 * (s_drained - s_closed) / s_closed);
    check(s_drained > 0 && std::fabs(s_drained - s_closed) / s_closed < 0.01,
          "the drained end state is the compression law's, within 1%");

    // The FE settlement history: consol_settlement[k] is the max |u_y| at time consol_time[k]
    // (k = 0 is the phase's start). Compared twice: as settlements, and as degrees of
    // consolidation U = s(t) / s(end state), each against its own end state -- the second takes
    // the compression law's small offset out and leaves the coupling alone.
    const size_t nfe = R.consol_settlement.size();
    check(nfe == (size_t)steps + 1, "one settlement per time step, plus the start (" + std::to_string(nfe) + ")");
    double worst_rel = 0.0, worst_U = 0.0;
    for (int s = 0; s < steps && s + 1 < (int)nfe; ++s) {
        const double fe = R.consol_settlement[(size_t)s + 1];
        const double ref = fd[(size_t)s] + s_sand;
        worst_rel = std::max(worst_rel, std::fabs(fe - ref) / ref);
        const double U_fe = fe / s_drained, U_fd = ref / s_closed;
        worst_U = std::max(worst_U, std::fabs(U_fe - U_fd));
        if (s == 0 || s == 4 || s == 9 || s == 24 || s == 49 || s == steps - 1)
            std::printf("  t = %7.1f d   FE %.5f m   column %.5f m   (%+.3f%%)   U %.4f / %.4f\n",
                        R.consol_time[(size_t)s + 1], fe, ref, 100.0 * (fe - ref) / ref, U_fe, U_fd);
    }
    const double fe_end = R.consol_settlement.back();
    std::printf("  worst: %.3f%% in settlement, %.4f in the degree of consolidation\n", 100.0 * worst_rel, worst_U);
    bool finite = std::isfinite(s_final) && std::isfinite(fe_end) && std::isfinite(s_drained);
    for (double v : fd) finite = finite && std::isfinite(v);
    check(finite, "the column and the finite element history are finite numbers");
    check(finite && worst_rel < 0.015, "the settlement history follows the independent column within 1.5%");
    check(finite && worst_U < 0.005, "...and its degree of consolidation within 0.005 at every step (the coupling)");
    const double U_end = fd.back() / s_final;
    check(U_end > 0.99, "the phase runs long enough to reach the end state (column U = " + std::to_string(U_end) + ")");
    check(finite && std::fabs(fe_end - s_drained) / s_drained < 0.005,
          "the consolidation ends where the drained phase does, within 0.5%");

    std::printf("\n%s (%d failure%s)\n", g_failures ? "FAILED" : "PASSED", g_failures, g_failures == 1 ? "" : "s");
    return g_failures ? 1 : 0;
}
