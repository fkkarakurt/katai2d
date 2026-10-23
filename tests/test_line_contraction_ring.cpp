// A tunnel lining that contracts in the ground: the interaction case of a line contraction.
//
// The two bounding cases -- a plate free to shrink and a plate held at both ends -- are checked by
// python/tests/line_contraction.py. Between them the ground decides how much of the contraction
// is realised, and for a circular ring in elastic ground that has a closed form. The contraction C
// (an area loss) is a free hoop strain eps0 = -C/2 of the ring. The ground resists a uniform
// radial displacement u of a circular cavity of radius R with the stiffness of a thick cylinder
// (Lame) whose outer boundary b is held,
//     k = [2 (lambda + G) / b^2 + 2 G / R^2] / (1/R - R / b^2),        p = k |u|,
// and the ring is in equilibrium with that traction, N = p R, while its hoop strain is
// u / R = eps0 + N / EA. Hence
//     u = eps0 R / (1 + k R^2 / EA),     N = -k R u.
// The model is a quarter of a weightless square of half-width b with the ring on 8 plate segments;
// the square's held sides stand in for the cylinder's held circle, which the correction term
// (about 1% of k here) makes nearly irrelevant.
//
// verify: KV-STR-019
//   oracle:   closed_form
//   source:   the thick-walled cylinder under a uniform radial displacement of its inner face, plane strain (Lame; Timoshenko, S.P. & Goodier, J.N. (1970). Theory of Elasticity, 3rd ed. McGraw-Hill, sec. 30), with the ring's hoop equilibrium N = p R and its hoop strain u / R = eps0 + N / EA
//   locator:  a quarter model, 60 x 60 m of weightless linear elastic ground (E = 1e4 kPa, nu = 0.3), symmetry on the axes and the far sides fully fixed; a 3 m lining of 8 plate segments (EA = 2.5e4 kN/m) installed in an excavated tunnel and then given a contraction of 1%
//   quantity: mean radial displacement of the lining's nodes [m]; mean axial force of the lining [kN/m]
//   expected: within 2% of the closed form for both (u = -7.8 mm, N = +60 kN/m tension)
//   band:     2% (the lining is an 8-segment polygon, not a circle; the far boundary a square, not the held circle)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
}  // namespace

int main() {
    constexpr double R = 3.0, b = 60.0, E = 1.0e4, nu = 0.3, EA = 2.5e4, C = 0.01;
    constexpr int nseg = 8;
    const double G = E / (2.0 * (1.0 + nu)), lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double k = (2.0 * (lam + G) / (b * b) + 2.0 * G / (R * R)) / (1.0 / R - R / (b * b));
    const double eps0 = -0.5 * C;
    const double u_want = eps0 * R / (1.0 + k * R * R / EA);
    const double N_want = -k * R * u_want;

    m::Project pr;
    m::Material s; s.name = "Ground"; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = nu; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    // The quarter square (edges: base y = 0 symmetry, right side held, top held, left x = 0 symmetry)
    m::SoilPolygon G0; G0.name = "Ground"; G0.material = 0;
    G0.x = {0, b, b, 0}; G0.y = {0, 0, b, b};
    G0.edge_bc = {(int)m::BCType::VerticallyFixed, (int)m::BCType::FullyFixed,
                  (int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(G0);
    // The tunnel: the quarter disk inside the ring's polygon (drawn over the ground).
    std::vector<double> rx, ry;
    for (int i = 0; i <= nseg; ++i) {
        const double a = 0.5 * kPi * i / nseg;
        rx.push_back(R * std::cos(a)); ry.push_back(R * std::sin(a));
    }
    m::SoilPolygon T; T.name = "Tunnel"; T.material = 0;
    T.x = {0.0}; T.y = {0.0};
    for (int i = 0; i <= nseg; ++i) { T.x.push_back(rx[i]); T.y.push_back(ry[i]); }
    T.edge_bc.assign(T.x.size(), 0);
    T.edge_bc.front() = (int)m::BCType::VerticallyFixed;      // along y = 0
    T.edge_bc.back() = (int)m::BCType::HorizontallyFixed;     // along x = 0
    pr.polygons.push_back(T);
    pr.x_min = 0; pr.x_max = b; pr.y_min = 0; pr.y_max = b;
    pr.has_water = false;
    m::PlateMaterial pm; pm.name = "Lining"; pm.EA = EA; pm.EI = 1.0e2; pm.nu = 0.0; pm.w = 0.0;
    pr.plates.push_back(pm);
    for (int i = 0; i < nseg; ++i) {
        m::StructElement p; p.kind = m::StructKind::Plate; p.name = "Lining " + std::to_string(i + 1);
        p.material = 0; p.x1 = rx[i]; p.y1 = ry[i]; p.x2 = rx[i + 1]; p.y2 = ry[i + 1];
        pr.structs.push_back(p);
    }
    pr.initial.struct_active.assign(nseg, 0);
    m::Phase dig; dig.name = "Excavate and line"; dig.poly_active = {1, 0};
    dig.struct_active.assign(nseg, 1);
    pr.phases.push_back(dig);
    m::Phase shrink = dig; shrink.name = "Contract"; shrink.struct_contraction.assign(nseg, C);
    pr.phases.push_back(shrink);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "quarter model meshed");
    if (!M.ok) { std::printf("  (%s)\n", M.message.c_str()); return 1; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[1].ok && res[2].ok, "excavation, lining and contraction ran");
    if (res.size() != 3 || !res[2].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str());
        return 1;
    }
    const auto& Rz = res[2];
    double ur = 0.0; int nr = 0;
    for (int n = 0; n < Rz.mesh.node_count; ++n) {
        const double r = std::hypot(Rz.mesh.x[n], Rz.mesh.y[n]);
        bool corner = false;
        for (int i = 0; i <= nseg; ++i) corner = corner || std::hypot(Rz.mesh.x[n] - rx[i], Rz.mesh.y[n] - ry[i]) < 1e-9;
        if (!corner || r < 1e-9) continue;   // the ring's own vertices: on the circle exactly
        ur += (Rz.disp[2 * n] * Rz.mesh.x[n] + Rz.disp[2 * n + 1] * Rz.mesh.y[n]) / r;
        ++nr;
    }
    ur = nr ? ur / nr : 0.0;
    double N = 0.0; int nn = 0;
    for (const auto& f : Rz.struct_forces)
        for (const auto& q : f.stations) { N += q.N; ++nn; }
    N = nn ? N / nn : 0.0;
    std::printf("  ground stiffness k = %.2f kPa/m (2G/R = %.2f), k R^2 / EA = %.3f\n", k, 2.0 * G / R, k * R * R / EA);
    std::printf("  radial displacement %.6f m (closed form %.6f m, %+.2f%%) over %d vertices\n",
                ur, u_want, 100.0 * (ur - u_want) / u_want, nr);
    std::printf("  axial force %.3f kN/m (closed form %.3f kN/m, %+.2f%%) over %d stations\n",
                N, N_want, 100.0 * (N - N_want) / N_want, nn);
    check(std::fabs(ur - u_want) < 0.02 * std::fabs(u_want), "the lining moves inward as far as the ground lets it (2%)");
    check(std::fabs(N - N_want) < 0.02 * std::fabs(N_want), "and carries the hoop force the ground holds it with (2%)");
    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all lining contraction checks passed\n");
    return 0;
}
