// Transient groundwater flow in an axisymmetric model: a well pumped at a constant rate from a
// confined aquifer, against Theis.
//
// A transient-flow phase in an axisymmetric model was refused (K2D-S002). Its conductivity and its
// storage are now integrated per radian, like the steady flow's (KV-FLW-012): the drawdown spreads
// from the well as the radial diffusion it is, which in plane strain would be a trench draining.
//
// verify: KV-FLW-014
//   oracle:   closed_form
//   source:   Theis, C.V. (1935). The relation between the lowering of the piezometric surface and the rate and duration of discharge of a well using ground-water storage. Transactions, American Geophysical Union 16(2), 519-524 -- the drawdown around a well pumped at a constant rate Q from a confined aquifer of transmissivity T and storativity S
//   locator:  s(r, t) = Q / (4 pi T) W(u), u = r^2 S / (4 T t), W(u) = E1(u) = -gamma - ln u + sum_k (-1)^(k+1) u^k / (k k!); T = k b, S = S_s b with the specific storage of the phase, S_s = n gamma_w / K_w (n = e / (1 + e), K_w = 2e6 kPa); the well is the inner boundary at r_w = 0.1 m given the outflow q_n = Q / (2 pi r_w b), the outer boundary at R = 300 m holds the initial head
//   quantity: the drawdown at r = 10 m and 30 m after 0.005, 0.01 and 0.02 day [m], in three transient-flow phases each continuing from the head the one before ended with
//   expected: the Theis drawdown; the outer boundary is not yet felt there (its image well adds less than 0.3% at 30 m and 0.02 day)
//   band:     1% of the Theis drawdown at every node within 2 m of either radius (backward Euler in 50-200 uniform steps, tri6); measured 0.42 / 0.17 / 0.10% at 30 m and 0.25 / 0.11 / 0.06% at 10 m after 0.005 / 0.01 / 0.02 day -- the error falls as the steps of the march become small against the time
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kRw = 0.1, kR = 300.0, kB = 10.0, kK = 1.0;   // [m], [m/day]
constexpr double kH0 = 20.0;                                   // initial head, above the top
constexpr double kQ = 50.0;                                    // pumping rate [m3/day]
constexpr double kE = 0.5;                                     // void ratio

double well_function(double u) {   // W(u) = E1(u), the series (u < ~5 here)
    double sum = -0.5772156649015329 - std::log(u), term = 1.0;
    for (int k = 1; k < 200; ++k) {
        term *= -u / k;                    // (-u)^k / k!
        sum -= term / k;
        if (std::fabs(term / k) < 1e-17) break;
    }
    return sum;
}

m::Project aquifer() {
    m::Project pr;
    pr.name = "Theis";
    pr.axisymmetric = true;
    pr.x_min = 0; pr.x_max = kR; pr.y_min = 0; pr.y_max = kB;
    pr.has_water = true; pr.wx = {0.0, kR}; pr.wy = {kH0, kH0};
    pr.mesh.elem_size = 2.0;
    m::Material s; s.name = "Aquifer"; s.model = m::SoilModel::LinearElastic; s.E = 5e4; s.nu = 0.3;
    s.gamma_unsat = 20.0; s.gamma_sat = 20.0; s.kx = s.ky = kK; s.e_init = kE;
    pr.materials.push_back(s);
    const int Cl = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head, Fx = (int)m::FlowBCType::Flux;
    const int Fr = (int)m::BCType::Free, Hz = (int)m::BCType::HorizontallyFixed, Fu = (int)m::BCType::FullyFixed;
    const double qn = -kQ / (2.0 * kPi * kRw * kB);   // out of the soil, into the well
    // Near the well the elements are small (the head varies as ln r there), further out larger.
    m::SoilPolygon near; near.name = "Near"; near.material = 0; near.coarseness = 0.15;
    near.x = {kRw, 6.0, 6.0, kRw}; near.y = {0, 0, kB, kB};
    near.edge_bc = {Fu, Fr, Fr, Hz};
    near.edge_flow = {Cl, Cl, Cl, Fx};
    near.edge_head = {0, 0, 0, 0};
    near.edge_flux = {0, 0, 0, qn};
    m::SoilPolygon far; far.name = "Far"; far.material = 0;
    far.x = {6.0, kR, kR, 6.0}; far.y = {0, 0, kB, kB};
    far.edge_bc = {Fu, Hz, Fr, Fr};
    far.edge_flow = {Cl, Hd, Cl, Cl};
    far.edge_head = {0, kH0, 0, 0};
    pr.polygons = {near, far};
    const double t[] = {0.005, 0.005, 0.01};   // three phases: to 0.005, 0.01 and 0.02 day
    const int steps[] = {50, 50, 100};
    for (int k = 0; k < 3; ++k) {
        m::Phase p; p.name = "Pumping " + std::to_string(k + 1); p.type = m::PhaseType::TransientFlow;
        p.duration = t[k]; p.time_steps = steps[k];
        pr.phases.push_back(p);
    }
    return pr;
}
}  // namespace

int main() {
    std::printf("== a well pumped from a confined aquifer, axisymmetric transient flow (KV-FLW-014) ==\n");
    const m::Project pr = aquifer();
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "the aquifer meshes");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    const bool ok = res.size() == 4 && res[1].ok && res[2].ok && res[3].ok;
    for (const auto& r : res)
        if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    check(ok, "three transient-flow phases solve in axisymmetry");
    if (!ok) return 1;

    const double n = kE / (1.0 + kE), Ss = n * 9.81 / 2.0e6, T = kK * kB, S = Ss * kB;
    const double times[] = {0.005, 0.01, 0.02};
    double worst = 0.0;
    for (int k = 0; k < 3; ++k) {
        const auto& h = res[(size_t)k + 1].flow_head;
        if ((int)h.size() < M.mesh.node_count) { check(false, "the phase hands on its head"); return 1; }
        for (double r : {10.0, 30.0}) {
            // Every node within 2 m of that radius, each against Theis at its own radius.
            double worst_here = 0.0, s_fe = 0.0, s_th = 0.0, u = 0.0;
            int cnt = 0;
            for (int nd = 0; nd < M.mesh.node_count; ++nd) {
                const double rn = M.mesh.x[nd];
                if (std::fabs(rn - r) > 2.0) continue;
                const double un = rn * rn * S / (4.0 * T * times[k]);
                const double sth = kQ / (4.0 * kPi * T) * well_function(un);
                const double sfe = kH0 - h[(size_t)nd];
                const double e = std::fabs(sfe / sth - 1.0);
                if (e >= worst_here) { worst_here = e; s_fe = sfe; s_th = sth; u = un; }
                ++cnt;
            }
            if (cnt == 0) { check(false, "nodes near r = " + std::to_string(r)); continue; }
            worst = std::max(worst, worst_here);
            std::printf("  t = %.3f d, r ~ %2.0f m (%d nodes): worst drawdown %.5f m, Theis %.5f m (u = %.4f)  %.3f%%\n",
                        times[k], r, cnt, s_fe, s_th, u, 100.0 * worst_here);
        }
    }
    check(worst < 0.01, "the drawdown is Theis's at both radii and three times (1%)");
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
