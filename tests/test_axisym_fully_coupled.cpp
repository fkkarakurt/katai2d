// A fully coupled flow-deformation phase in an axisymmetric model: a column against Terzaghi, the
// unit cell of a vertical drain against Barron, and both against the consolidation phase.
//
// The fully coupled phase was refused in axisymmetry (K2D-S002). Its elastoplastic core now takes
// the consolidation twin's kinematics: the hoop strain and stress in the skeleton, the pore
// pressure on the hoop, and every integral per radian. In a saturated soil the phase IS a
// consolidation (the degree of saturation is one, the Bishop parameter one, the storage the water's
// own), so the consolidation phase of the same model is an independent path to the same history,
// and the closed forms of KV-CON-008 hold for it as they do for the consolidation.
//
// verify: KV-CON-009
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley -- one-dimensional consolidation; Barron, R.A. (1948). Consolidation of fine-grained soils by drain wells. Transactions of the ASCE 113, 718-742 -- radial consolidation of the cylinder of soil around a drain, under equal vertical strain
//   locator:  (A) a column of radius 1 m and height H = 4 m about the axis, nu = 0 and both radial sides on rollers, drained at the top, a fully coupled phase to Tv = 0.5 in 120 steps; (B) the unit cell of a drain, r_w = 0.05 .. r_e = 0.75 m (n = 15), 1 m high, drained at the drain only, fully coupled phases to Barron's t50 and t90; each also as a consolidation phase of the same model
//   quantity: the degree of consolidation, settlement over its final value [-]; the settlement history of the fully coupled phase against the consolidation phase's [m]
//   expected: (A) Terzaghi's U = 0.7639; (B) Barron's U_h = 0.5 and 0.9; and the consolidation phase's settlement at every recorded time
//   band:     (A) 0.5% against Terzaghi (the backward-Euler march, as KV-CON-008); (B) 3% (Barron's equal-strain closed form against free strain, as KV-CON-008); 1e-6 of the final settlement between the two phases (two cores, one saturated problem). Measured: -0.197% against Terzaghi (KV-CON-008's own); Barron +1.3% at t50 and -0.6% at t90; the two phases equal to the last digit, 0.0 m over all 121 times
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

// x0..x1 radially, height h; flow on edges bottom, outer, top, inner.
m::Project cell(m::PhaseType type, double x0, double x1, double h, std::vector<int> flow,
                double duration, int steps, double esize) {
    m::Project pr;
    pr.name = "axisymmetric fully coupled";
    pr.axisymmetric = true;
    pr.x_min = 0; pr.x_max = x1; pr.y_min = 0; pr.y_max = h;
    pr.has_water = false;   // the pore pressure is the excess alone
    pr.mesh.elem_size = esize;
    m::Material s; s.name = "Clay"; s.model = m::SoilModel::LinearElastic; s.E = kE; s.nu = 0.0;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = 1.0; s.kx = s.ky = kK;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Clay"; P.material = 0;
    P.x = {x0, x1, x1, x0}; P.y = {0, 0, h, h};
    P.edge_bc = {Fu, Hz, Fr, Hz};
    P.edge_flow = flow;
    P.edge_head = {0.0, h, h, h};
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = x0; L.y1 = h; L.x2 = x1; L.y2 = h; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};
    m::Phase c; c.name = "Coupled"; c.type = type;
    c.load_active = {1}; c.duration = duration; c.time_steps = steps; c.consol_load_at_start = true;
    pr.phases.push_back(c);
    return pr;
}

struct History { bool ok = false; std::vector<double> t, s; std::vector<Eigen::VectorXd> disp; };
History solve(const m::Project& pr) {
    History H;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return H; }
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    if (res.size() != pr.phases.size() + 1) return H;
    for (size_t k = 1; k < res.size(); ++k) {
        if (!res[k].ok) { std::printf("      (%s)\n", res[k].message.substr(0, 300).c_str()); return H; }
        H.disp.push_back(res[k].disp_total);
    }
    H.t = res[1].consol_time;
    H.s = res[1].consol_settlement;
    H.ok = !H.t.empty();
    return H;
}
}  // namespace

int main() {
    std::printf("== a fully coupled phase in axisymmetry: Terzaghi, Barron, and the consolidation phase (KV-CON-009) ==\n");
    const double gw = katai::app::kGammaWater, cv = kK * kE / gw;

    // (A) the one-dimensional column, drained at the top
    const double H = 4.0, Tv = 0.5, tA = Tv * H * H / cv, sinf = kQ * H / kE;
    const std::vector<int> topA = {Cl, Cl, Hd, Cl};
    const History A = solve(cell(m::PhaseType::FullyCoupled, 0.0, 1.0, H, topA, tA, 120, 0.4));
    const History Ac = solve(cell(m::PhaseType::Consolidation, 0.0, 1.0, H, topA, tA, 120, 0.4));
    check(A.ok && Ac.ok, "(A) the column consolidates, fully coupled and as a consolidation phase");
    if (A.ok && Ac.ok) {
        const double U = A.s.back() / sinf, Ut = terzaghi_U(Tv);
        double diff = A.s.size() == Ac.s.size() ? 0.0 : 1e300;
        for (size_t k = 0; diff < 1e300 && k < A.s.size(); ++k) diff = std::max(diff, std::fabs(A.s[k] - Ac.s[k]));
        std::printf("  (A) U(Tv = 0.5): fully coupled %.5f, Terzaghi %.5f (%+.3f%%); against the consolidation phase, %zu times: max |ds| %.3e m (%.2e of the final)\n",
                    U, Ut, 100.0 * (U / Ut - 1.0), A.s.size(), diff, diff / A.s.back());
        check(std::fabs(U / Ut - 1.0) < 5e-3, "(A) the fully coupled column is Terzaghi's (0.5%)");
        check(diff < 1e-6 * A.s.back(), "(A) and its settlement history is the consolidation phase's (1e-6)");
    }

    // (B) the unit cell of a drain: radial flow to r_w only, two phases to Barron's t50 and t90
    const double rw = 0.05, re = 0.75, n = re / rw, de = 2.0 * re;
    const double mu = n * n / (n * n - 1.0) * std::log(n) - (3.0 * n * n - 1.0) / (4.0 * n * n);
    const double t50 = std::log(2.0) * mu * de * de / (8.0 * cv), t90 = std::log(10.0) * mu * de * de / (8.0 * cv);
    const auto two_phases = [&](m::PhaseType type) {
        m::Project pb = cell(type, rw, re, 1.0, {Cl, Cl, Cl, Hd}, t50, 120, 0.05);
        m::Phase p2 = pb.phases[0];
        p2.name = "To t90"; p2.duration = t90 - t50; p2.time_steps = 160;
        p2.consol_load_at_start = false;
        pb.phases.push_back(p2);
        return pb;
    };
    const m::Project pb = two_phases(m::PhaseType::FullyCoupled);
    const History B = solve(pb);
    const History Bc = solve(two_phases(m::PhaseType::Consolidation));
    check(B.ok && Bc.ok && B.disp.size() == 2 && Bc.disp.size() == 2,
          "(B) the unit cell of the drain consolidates in two phases, fully coupled and as consolidation");
    if (B.ok && Bc.ok && B.disp.size() == 2 && Bc.disp.size() == 2) {
        const auto MB = katai::app::mesh_from_project(pb);
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
        const double U50 = average_top(B.disp[0]) / sB, U90 = average_top(B.disp[1]) / sB;
        const double C50 = average_top(Bc.disp[0]) / sB, C90 = average_top(Bc.disp[1]) / sB;
        std::printf("  (B) n = %.0f, mu = %.4f: U(t50) = %.5f (Barron 0.5, consolidation %.5f), U(t90) = %.5f (Barron 0.9, consolidation %.5f)\n",
                    n, mu, U50, C50, U90, C90);
        check(std::fabs(U50 - 0.5) < 0.03 * 0.5 && std::fabs(U90 - 0.9) < 0.03 * 0.9,
              "(B) the drain consolidates the cell at Barron's rate (3%)");
        check(std::fabs(U50 - C50) < 1e-6 && std::fabs(U90 - C90) < 1e-6,
              "(B) and as the consolidation phase does (1e-6)");
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
