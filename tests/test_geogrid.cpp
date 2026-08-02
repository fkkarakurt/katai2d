// Geogrid (P2.4): the 3-node tension-only axial membrane (elements/geogrid.hpp) wired into the
// solver. N = EA*eps, tension only (compression -> slack, N=0, reversible), optional N_p cap
// (plastic). Translational DOFs shared with the soil; no bending / rotation (PLAXIS MMM sec 18.2).
//
// (A) Isolated axial bar (closed form): a single geogrid element fixed at one end, axial end
//     load F -> the quadratic bar reproduces the exact linear field, u_B = F Lg/EA, u_mid = F Lg/2EA.
// (B) Tension-only cutoff (defining behaviour): a geogrid chain anchored at its left end. Pulling
//     the right end (tension) stiffens the response (less deflection than no grid); pushing it
//     (compression) leaves the grid slack -> deflection equals the no-grid case.
// (C) N_p yield: stretching past N_p caps the axial force at N_p and leaves a permanent eps_p.
// (See docs/references/structural-plate-formulation.md sec 8.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/geogrid.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::GeogridElement;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::Structures;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(k);
    return s->solve(r);
}
std::vector<int> sorted_top(const Mesh& mesh) {
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    return top;
}

// (A) Single geogrid bar, one end fixed, axial end load -> exact linear field.
void test_geogrid_axial_tension() {
    constexpr double W = 8.0, H = 2.0;
    const int nx = 8;
    const double EA = 1.0e5, F = 100.0;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, 2);
    const std::vector<int> top = sorted_top(mesh);
    const int A = top[0], mid = top[1], B = top[2];
    const double Lg = mesh.x[B] - mesh.x[A];  // = W/nx = 1

    // Fix everything, then free only ux of mid and B (A stays fixed -> the bar's anchored end).
    DofMap dofs(mesh.node_count, 2);
    for (int n = 0; n < mesh.node_count; ++n) {
        dofs.fix_node_component(n, 1);
        if (n != mid && n != B) dofs.fix_node_component(n, 0);
    }
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    f(dofs.equation(dofs.global_dof(B, 0))) = F;  // pull B in +x (tension)

    Structures st;
    st.geogrids = {GeogridElement{{A, B, mid}, {EA, -1.0}}};
    const MaterialModel m{MaterialType::LinearElastic, 1.0e-3, 0.3};  // negligible soil
    const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                {1, 50, 1e-12}, {}, {}, st);
    check(r.converged, "geogrid axial-tension solve converged");
    const double uB = r.displacement[dofs.global_dof(B, 0)];
    const double um = r.displacement[dofs.global_dof(mid, 0)];
    std::printf("  axial bar: u_B=%.6e (F Lg/EA=%.6e)  u_mid=%.6e (half)\n",
                uB, F * Lg / EA, um);
    check(close(uB, F * Lg / EA, 2e-3), "u_B = F Lg/EA (quadratic bar, exact linear field)");
    check(close(um, 0.5 * F * Lg / EA, 2e-3), "u_mid = F Lg/(2EA)");
}

// Build a geogrid chain along the (x-sorted) top edge. nx soil columns -> nx geogrid elements.
std::vector<GeogridElement> build_top_geogrid(const Mesh& mesh, double EA, double Np) {
    const std::vector<int> top = sorted_top(mesh);
    const int nn = static_cast<int>(top.size());  // 2*nx + 1
    std::vector<GeogridElement> g;
    for (int e = 0; e * 2 + 2 < nn; ++e)
        g.push_back(GeogridElement{{top[2 * e], top[2 * e + 2], top[2 * e + 1]}, {EA, Np}});
    return g;
}

// (B) Tension-only cutoff: anchored-left geogrid chain; tension stiffens, compression goes slack.
double run_chain(double F_signed, bool with_grid) {
    constexpr double W = 8.0, H = 2.0;
    const int nx = 8;
    const double E_soil = 1.0e4, EA = 2.0e5;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, 2);
    const std::vector<int> top = sorted_top(mesh);
    const int right = top.back();

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);  // geogrid left end anchored in x
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    f(dofs.equation(dofs.global_dof(right, 0))) = F_signed;

    Structures st;
    if (with_grid) st.geogrids = build_top_geogrid(mesh, EA, -1.0);
    const MaterialModel m{MaterialType::LinearElastic, E_soil, 0.3};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                {2, 50, 1e-10}, {}, {}, st);
    check(r.converged, "geogrid chain solve converged");
    return r.displacement[dofs.global_dof(right, 0)];
}

void test_geogrid_tension_only() {
    constexpr double F = 400.0;
    const double pull_grid = run_chain(+F, true), pull_none = run_chain(+F, false);   // tension
    const double push_grid = run_chain(-F, true), push_none = run_chain(-F, false);   // compression
    std::printf("  pull(+x): grid=%.6e none=%.6e   push(-x): grid=%.6e none=%.6e\n",
                pull_grid, pull_none, push_grid, push_none);
    check(pull_grid < 0.7 * pull_none, "geogrid in TENSION reduces deflection (carries load)");
    check(close(push_grid, push_none, 0.02), "geogrid in COMPRESSION is slack (~ no grid)");
}

// (C) N_p yield: stretch the isolated bar past N_p -> force capped, permanent eps_p.
void test_geogrid_Np_yield() {
    constexpr double W = 8.0, H = 2.0;
    const int nx = 8;
    const double E_soil = 1.0e4, EA = 1.0e5;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, 2);
    const std::vector<int> top = sorted_top(mesh);
    const int right = top.back();

    auto run = [&](double Np) {
        DofMap dofs(mesh.node_count, 2);
        for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
        for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
        dofs.finalize();
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        f(dofs.equation(dofs.global_dof(right, 0))) = 250.0;  // pull right (+x, tension)
        Structures st; st.geogrids = build_top_geogrid(mesh, EA, Np);
        const MaterialModel m{MaterialType::LinearElastic, E_soil, 0.3};
        const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                    {4, 60, 1e-10}, {}, {}, st);
        check(r.converged, "geogrid N_p solve converged");
        // Recover the axial force at the first geogrid element's first Gauss point.
        const auto& ge = st.geogrids[0];
        katai::core::geogrid::NodeCoords Xe;
        katai::core::geogrid::Dof ug = katai::core::geogrid::Dof::Zero();
        for (int k = 0; k < 3; ++k) {
            Xe(k, 0) = mesh.x[ge.nodes[k]]; Xe(k, 1) = mesh.y[ge.nodes[k]];
            ug(2 * k + 0) = r.displacement[dofs.global_dof(ge.nodes[k], 0)];
            ug(2 * k + 1) = r.displacement[dofs.global_dof(ge.nodes[k], 1)];
        }
        const auto kin = katai::core::geogrid::axial_kin(Xe, katai::core::geogrid::gauss_xi()[0]);
        const double eps = (kin.Be * ug)(0);
        const double ep = r.geogrid_plastic.empty() ? 0.0 : r.geogrid_plastic[0];
        const double N = EA * (eps - ep);
        const double ux = r.displacement[dofs.global_dof(right, 0)];
        return std::make_tuple(N, ep, ux);
    };

    const auto el = run(1.0e9);   // effectively unlimited -> elastic
    const auto pl = run(40.0);    // cap below the elastic force
    const double N_el = std::get<0>(el), N_pl = std::get<0>(pl);
    const double ep_pl = std::get<1>(pl), ux_el = std::get<2>(el), ux_pl = std::get<2>(pl);
    std::printf("  N_p: elastic N=%.3f (eps_p=%.2e)  capped N=%.3f (N_p=40, eps_p=%.3e)  "
                "ux_el=%.4e ux_pl=%.4e\n", N_el, std::get<1>(el), N_pl, ep_pl, ux_el, ux_pl);
    check(N_el > 40.0, "elastic geogrid force exceeds the cap (cap is active in the capped run)");
    check(close(N_pl, 40.0, 0.03), "axial force capped at N_p");
    check(ep_pl > 1e-9, "N_p yield leaves a permanent plastic strain eps_p");
    check(ux_pl > ux_el, "yielding geogrid -> larger end displacement than elastic");
}

} // namespace

int main() {
    test_geogrid_axial_tension();
    test_geogrid_tension_only();
    test_geogrid_Np_yield();
    if (g_failures == 0) {
        std::printf("OK: geogrid verified (axial bar closed form, tension-only cutoff, N_p yield)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
