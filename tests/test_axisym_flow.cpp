// Steady groundwater flow in an axisymmetric model: radial flow to a well, against Thiem.
//
// In an axisymmetric model x is the radius and every integral of the flow is taken per radian --
// the conductivity, the nodal flux, a prescribed boundary flux, the seams. Until this the flow
// solve did not know the model was axisymmetric at all: a radial flow came out as the plane flow
// through a slab, a straight head profile where the ground has a logarithmic one, and a discharge
// that is not the well's. This test is that difference made visible and closed.
//
// verify: KV-FLW-012
//   oracle:   closed_form
//   source:   Thiem, G. (1906). Hydrologische Methoden. Leipzig: J.M. Gebhardt -- the steady radial flow to a well fully penetrating a confined aquifer of thickness b and permeability k, between the well radius r_w at head h_w and a radius R at head h_R
//   locator:  h(r) = h_w + (h_R - h_w) ln(r / r_w) / ln(R / r_w); the discharge of the whole ring Q = 2 pi k b (h_R - h_w) / ln(R / r_w); with the inner boundary given the inflow q_n = -Q / (2 pi r_w b) instead of a head, the head there is h_R - Q ln(R / r_w) / (2 pi k b)
//   quantity: the head at every node [m]; the discharge of the whole ring [m3/day]; the head at the well with the discharge prescribed [m]
//   expected: the logarithmic profile, the Thiem discharge and the Thiem drawdown; the same geometry in plane strain gives the straight profile of a slab, so the two kinematics are told apart; a well drawn as a line is refused in axisymmetry (its discharge per metre out of the plane has no meaning there, and on the axis the ring integral would weigh it by a radius of zero)
//   band:     heads within 0.1% of the head difference (measured 0.043%), discharge within 0.5% (measured +0.005%), the well head within 0.5% of the drawdown (measured 0.009%) -- tri6, 0.4 m elements on a 0.5 .. 50 m annulus; in plane strain the same geometry is 2.06 m off the radial profile at r = 5 m
#include <katai/jobs/flow_driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kRw = 0.5, kR = 50.0, kB = 10.0;   // well radius, outer radius, aquifer thickness [m]
constexpr double kK = 1.0;                          // permeability [m/day]
constexpr double kHw = 20.0, kHr = 25.0;            // heads, above the top: confined throughout

// The aquifer as a ring section: edges bottom (closed), outer (head), top (closed), well (`inner`).
m::Project aquifer(bool axisymmetric, int inner, double inner_value) {
    m::Project pr;
    pr.name = "radial flow";
    pr.axisymmetric = axisymmetric;
    pr.x_min = 0; pr.x_max = kR; pr.y_min = 0; pr.y_max = kB;
    pr.has_water = false;
    pr.mesh.elem_size = 0.4;
    m::Material s; s.name = "Aquifer"; s.model = m::SoilModel::LinearElastic; s.E = 1e4; s.nu = 0.3;
    s.kx = s.ky = kK;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Aquifer"; P.material = 0;
    P.x = {kRw, kR, kR, kRw}; P.y = {0, 0, kB, kB};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head,
                   (int)m::FlowBCType::Closed, inner};
    P.edge_head = {0.0, kHr, 0.0, inner == (int)m::FlowBCType::Head ? inner_value : 0.0};
    P.edge_flux = {0.0, 0.0, 0.0, inner == (int)m::FlowBCType::Flux ? inner_value : 0.0};
    pr.polygons.push_back(P);
    return pr;
}
double thiem(double r) { return kHw + (kHr - kHw) * std::log(r / kRw) / std::log(kR / kRw); }
}  // namespace

int main() {
    std::printf("== steady radial flow to a well, axisymmetric (KV-FLW-012) ==\n");
    const double Q = 2.0 * kPi * kK * kB * (kHr - kHw) / std::log(kR / kRw);

    // (1) Heads prescribed at the well and outside: the Thiem profile and discharge.
    const m::Project pr = aquifer(true, (int)m::FlowBCType::Head, kHw);
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "the annulus meshes");
    if (!M.ok) return 1;
    const auto R = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(R.ok, "the axisymmetric flow solves");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return 1; }
    double worst = 0.0;
    for (int n = 0; n < M.mesh.node_count; ++n)
        worst = std::max(worst, std::fabs(R.head[n] - thiem(M.mesh.x[n])));
    std::printf("  heads: worst |h - Thiem| = %.3e m of a %.1f m difference; Q = %.4f against %.4f m3/day (%+.3f%%)\n",
                worst, kHr - kHw, R.discharge, Q, 100.0 * (R.discharge / Q - 1.0));
    std::printf("  %s\n", R.message.c_str());
    check(worst < 1e-3 * (kHr - kHw), "the heads are the logarithmic Thiem profile (0.1% of the difference)");
    check(std::fabs(R.discharge / Q - 1.0) < 5e-3, "the discharge of the whole ring is Thiem's (0.5%)");

    // (2) The same geometry in plane strain is a slab: a straight profile. What the flow gave for
    // an axisymmetric model until now, and far from the radial one.
    const m::Project ps = aquifer(false, (int)m::FlowBCType::Head, kHw);
    const auto Rp = katai::app::solve_groundwater_flow(ps, M.mesh);
    check(Rp.ok, "the plane-strain control solves");
    if (Rp.ok) {
        double lin = 0.0, at5 = 0.0;
        for (int n = 0; n < M.mesh.node_count; ++n) {
            const double h_lin = kHw + (kHr - kHw) * (M.mesh.x[n] - kRw) / (kR - kRw);
            lin = std::max(lin, std::fabs(Rp.head[n] - h_lin));
            if (std::fabs(M.mesh.x[n] - 5.0) < 0.21) at5 = Rp.head[n];
        }
        std::printf("  plane strain: worst |h - straight| = %.3e m; at r = 5 m %.3f (radial %.3f)\n",
                    lin, at5, thiem(5.0));
        check(lin < 1e-6 * (kHr - kHw) && std::fabs(at5 - thiem(5.0)) > 1.0,
              "in plane strain the profile is straight, metres away from the radial one");
    }

    // (3) The discharge prescribed at the well instead of its head: the Thiem drawdown.
    const double qn = -Q / (2.0 * kPi * kRw * kB);   // out of the soil
    const m::Project pf = aquifer(true, (int)m::FlowBCType::Flux, qn);
    const auto Rf = katai::app::solve_groundwater_flow(pf, M.mesh);
    check(Rf.ok, "the flow with the discharge prescribed at the well solves");
    if (Rf.ok) {
        double hw = 0.0; int cnt = 0;
        for (int n = 0; n < M.mesh.node_count; ++n)
            if (std::fabs(M.mesh.x[n] - kRw) < 1e-9) { hw += Rf.head[n]; ++cnt; }
        hw /= std::max(1, cnt);
        std::printf("  well head with Q prescribed: %.5f m against %.5f (drawdown %.3f m)\n", hw, kHw,
                    kHr - kHw);
        check(std::fabs(hw - kHw) < 5e-3 * (kHr - kHw), "the well is drawn down by Thiem's drawdown (0.5%)");
    }

    // (4) A well drawn as a line in an axisymmetric model is refused, with the way to model it.
    m::Project pw = aquifer(true, (int)m::FlowBCType::Head, kHw);
    m::HydroLine w; w.name = "Well"; w.kind = m::HydroKind::Well; w.q = 10.0;
    w.x1 = 10.0; w.y1 = 1.0; w.x2 = 10.0; w.y2 = 9.0;
    pw.hydros.push_back(w);
    const auto Mw = katai::app::mesh_from_project(pw);
    const auto Rw = Mw.ok ? katai::app::solve_groundwater_flow(pw, Mw.mesh) : katai::app::FlowResult{};
    std::printf("  (%s)\n", Rw.message.c_str());
    check(Mw.ok && !Rw.ok && Rw.message.find("axisymmetric") != std::string::npos,
          "a well line in axisymmetry is refused, saying how to model it");

    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
