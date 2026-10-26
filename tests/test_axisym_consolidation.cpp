// Consolidation in an axisymmetric model: a column against Terzaghi, and the unit cell of a
// vertical drain against Barron.
//
// A consolidation phase in an axisymmetric model was refused (K2D-S002). Its core is now the
// elastoplastic one under the axisymmetric kinematics: the hoop strain and stress in the skeleton,
// the pore pressure acting on the hoop too, and every integral -- equilibrium, coupling, flow and
// storage -- per radian. Radial consolidation towards a vertical drain is the design question that
// needs it: the unit cell of a drain is a cylinder.
//
// verify: KV-CON-008
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley -- one-dimensional consolidation; Barron, R.A. (1948). Consolidation of fine-grained soils by drain wells. Transactions of the ASCE 113, 718-742 -- radial consolidation of the cylinder of soil around a drain, under equal vertical strain
//   locator:  (A) a column of radius 1 m and height H = 4 m about the axis, nu = 0 and both radial sides on rollers, so the problem is one-dimensional and U(Tv) = 1 - sum_j (2/M^2) exp(-M^2 Tv), M = (2j+1) pi / 2, Tv = cv t / H^2, cv = k Eoed / gamma_w, drained at the top; (B) the unit cell of a drain: an annulus r_w = 0.05 .. r_e = 0.75 m (n = 15), 1 m high, drained at the drain only, U_h = 1 - exp(-8 T_h / mu), mu = n^2/(n^2-1) ln n - (3n^2-1)/(4n^2), T_h = c_h t / d_e^2, d_e = 2 r_e; (C) the column with no flow boundary declared against the same column with the axis declared closed, the top and the outer side open
//   quantity: the degree of consolidation, settlement over its final value, at Tv = 0.5 (A) and when Barron's U_h is 0.5 and 0.9 (B) [-]; the settlement history with and without the declaration (C) [m]
//   expected: (A) Terzaghi's U = 0.7639, and the plane-strain column's; (B) Barron's U_h; (C) the same history -- the axis is a line of symmetry, not one of the sides that drain by default
//   band:     (A) 0.5% (measured -0.197%, the error of the backward-Euler march of 120 steps that KV-NUM-010 measures on the plane column: 0.1965%) and 0.2% against plane strain (measured 0.001%); (B) 3% -- the finite elements carry a flexible load, free strain, and Barron's closed form is the equal-strain one, which he showed lies within a few per cent of it (measured +1.3% at t50, -0.6% at t90); (C) 1e-12 m (measured 0; with the axis draining by default the histories part by 1.6 mm)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <Eigen/Core>

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
constexpr double kE = 1000.0, kK = 0.1, kQ = 10.0;   // nu = 0: Eoed = E [kPa]; [m/day]; [kPa]

double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 200; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += 2.0 / (M * M) * std::exp(-M * M * Tv);
    }
    return 1.0 - s;
}

const int Cl = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
const int Fr = (int)m::BCType::Free, Hz = (int)m::BCType::HorizontallyFixed, Fu = (int)m::BCType::FullyFixed;

// x0..x1 radially, height h; flow on edges bottom, outer, top, inner; `declare` false = none drawn.
m::Project cell(bool axi, double x0, double x1, double h, std::vector<int> flow, bool declare,
                double duration, int steps, double esize) {
    m::Project pr;
    pr.name = "axisymmetric consolidation";
    pr.axisymmetric = axi;
    pr.x_min = 0; pr.x_max = x1; pr.y_min = 0; pr.y_max = h;
    pr.has_water = false;   // the pore pressure is the excess alone
    pr.mesh.elem_size = esize;
    m::Material s; s.name = "Clay"; s.model = m::SoilModel::LinearElastic; s.E = kE; s.nu = 0.0;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = 1.0; s.kx = s.ky = kK;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Clay"; P.material = 0;
    P.x = {x0, x1, x1, x0}; P.y = {0, 0, h, h};
    P.edge_bc = {Fu, Hz, Fr, Hz};
    if (declare) {
        P.edge_flow = flow;
        P.edge_head = {0.0, h, h, h};
    }
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = x0; L.y1 = h; L.x2 = x1; L.y2 = h; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};
    m::Phase c; c.name = "Consolidation"; c.type = m::PhaseType::Consolidation;
    c.load_active = {1}; c.duration = duration; c.time_steps = steps; c.consol_load_at_start = true;
    pr.phases.push_back(c);
    return pr;
}

struct History { bool ok = false; std::vector<double> t, s; };
History solve(const m::Project& pr) {
    History H;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return H; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return H;
    }
    H.t = res[1].consol_time;
    H.s = res[1].consol_settlement;
    H.ok = !H.t.empty();
    return H;
}
double at(const History& h, double t) {   // linear between the recorded steps
    for (size_t k = 1; k < h.t.size(); ++k)
        if (h.t[k] >= t) {
            const double f = (t - h.t[k - 1]) / (h.t[k] - h.t[k - 1]);
            return h.s[k - 1] + f * (h.s[k] - h.s[k - 1]);
        }
    return h.s.back();
}
}  // namespace

int main() {
    std::printf("== consolidation in axisymmetry: Terzaghi and Barron (KV-CON-008) ==\n");
    const double gw = katai::app::kGammaWater, cv = kK * kE / gw;

    // (A) the one-dimensional column, drained at the top
    const double H = 4.0, Tv = 0.5, tA = Tv * H * H / cv, sinf = kQ * H / kE;
    const std::vector<int> topA = {Cl, Cl, Hd, Cl};
    const History A = solve(cell(true, 0.0, 1.0, H, topA, true, tA, 120, 0.4));
    const History Ap = solve(cell(false, 0.0, 1.0, H, topA, true, tA, 120, 0.4));
    check(A.ok && Ap.ok, "(A) the column consolidates, axisymmetric and plane strain");
    if (A.ok && Ap.ok) {
        const double U = A.s.back() / sinf, Up = Ap.s.back() / sinf, Ut = terzaghi_U(Tv);
        std::printf("  (A) U(Tv = 0.5): axisymmetric %.5f, plane strain %.5f, Terzaghi %.5f (%+.3f%%)\n", U, Up,
                    Ut, 100.0 * (U / Ut - 1.0));
        check(std::fabs(U / Ut - 1.0) < 5e-3, "(A) the axisymmetric column is Terzaghi's (0.5%)");
        check(std::fabs(U / Up - 1.0) < 2e-3, "(A) and the plane-strain column's (0.2%)");
    }

    // (B) the unit cell of a drain: radial flow to r_w only. Two phases, to Barron's t50 and t90; the
    // degree of consolidation from the settlement of the top AVERAGED over the cell's area (int u r dr
    // / int r dr) -- the free-strain settlement is not uniform, and its largest value, at the drain,
    // is ahead of the average.
    const double rw = 0.05, re = 0.75, n = re / rw, de = 2.0 * re;
    const double mu = n * n / (n * n - 1.0) * std::log(n) - (3.0 * n * n - 1.0) / (4.0 * n * n);
    const double t50 = std::log(2.0) * mu * de * de / (8.0 * cv), t90 = std::log(10.0) * mu * de * de / (8.0 * cv);
    {
        m::Project pb = cell(true, rw, re, 1.0, {Cl, Cl, Cl, Hd}, true, t50, 120, 0.05);
        m::Phase p2 = pb.phases[0];
        p2.name = "Consolidation to t90"; p2.duration = t90 - t50; p2.time_steps = 160;
        p2.consol_load_at_start = false;
        pb.phases.push_back(p2);
        const auto MB = katai::app::mesh_from_project(pb);
        const auto rb = MB.ok ? katai::app::solve_phases(pb, MB.mesh, katai::app::InitialPhase::K0Procedure)
                              : std::vector<katai::app::SolveResult>{};
        const bool okB = rb.size() == 3 && rb[1].ok && rb[2].ok;
        for (const auto& r : rb) if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
        check(okB, "(B) the unit cell of the drain consolidates in two phases");
        if (okB) {
            // the area average of the top's settlement, by the trapezoid rule over its nodes in r
            auto average_top = [&](const Eigen::VectorXd& u) {
                std::vector<std::pair<double, double>> top;
                for (int k = 0; k < MB.mesh.node_count; ++k)
                    if (std::fabs(MB.mesh.y[k] - 1.0) < 1e-9) top.push_back({MB.mesh.x[k], -u[2 * k + 1]});
                std::sort(top.begin(), top.end());
                double num = 0.0, den = 0.0;
                for (size_t k = 1; k < top.size(); ++k) {
                    const double r0 = top[k - 1].first, r1 = top[k].first, dr = r1 - r0;
                    num += 0.5 * dr * (top[k - 1].second * r0 + top[k].second * r1);
                    den += 0.5 * dr * (r0 + r1);
                }
                return den > 0.0 ? num / den : 0.0;
            };
            const double sB = kQ * 1.0 / kE;
            const double U50 = average_top(rb[1].disp_total) / sB, U90 = average_top(rb[2].disp_total) / sB;
            std::printf("  (B) n = %.0f, mu = %.4f: U(t50 = %.4f d) = %.4f (Barron 0.5), U(t90 = %.4f d) = %.4f (Barron 0.9)\n",
                        n, mu, t50, U50, t90, U90);
            check(std::fabs(U50 - 0.5) < 0.03 * 0.5 && std::fabs(U90 - 0.9) < 0.03 * 0.9,
                  "(B) the drain consolidates the cell at Barron's rate (3%)");
        }
    }

    // (C) no flow boundary declared: the top and the outer side drain, the axis does not
    const std::vector<int> openC = {Cl, Hd, Hd, Cl};   // bottom closed, outer and top open, axis closed
    const History Cd = solve(cell(true, 0.0, 1.0, H, openC, true, tA, 40, 0.4));
    const History Cn = solve(cell(true, 0.0, 1.0, H, openC, false, tA, 40, 0.4));
    check(Cd.ok && Cn.ok, "(C) both columns consolidate");
    if (Cd.ok && Cn.ok) {
        double diff = Cd.s.size() == Cn.s.size() ? 0.0 : 1e300;
        for (size_t k = 0; diff < 1e300 && k < Cd.s.size(); ++k) diff = std::max(diff, std::fabs(Cd.s[k] - Cn.s[k]));
        std::printf("  (C) settlement histories, none declared against the axis declared closed: %.3e m (final %.6f m)\n",
                    diff, Cd.s.back());
        check(diff < 1e-12, "(C) undeclared, the axis is closed and the outer side and top drain");
    }
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
