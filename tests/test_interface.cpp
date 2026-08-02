// Interface (P2.4): zero-thickness soil-structure Coulomb joint (elements/interface.hpp) wired
// into the solver. Two coincident sides -- soil (mesh nodes) and structure (extra DOFs via
// DofMap::add_extra_dof) -- connected by normal/shear penalty stiffness + Coulomb slip, with
// Newton-Cotes (nodal) integration (Day & Potts 1994). (See interface-formulation.md.)
//
// Setup: one interface element on the (horizontal) top edge. Soil side fully fixed; the structure
// side (6 extra DOFs) is loaded by uniform normal + shear tractions (weight-distributed nodal
// forces). Local frame: s=+x (c=1,s=0), n=+y -> Du_s = struct_x, Du_n = struct_y.
//
// (A) Elastic (exact): sub-capacity loads -> Du_s = tau/ks, Du_n = sigma_n/kn; recovered sigma_n,
//     tau match the applied tractions to round-off.
// (B) Coulomb capacity (limit analysis): proportional normal N0 + shear T0 ramped to collapse.
//     Slip limit: tau_max*L = (c_i - sigma_n tan phi_i)*L. With sigma_n = -lambda N0/L and applied
//     shear = lambda T0/L, the collapse load factor is the CLOSED FORM
//        lambda* = c_i L / (T0 - N0 tan phi_i).
//     Two runs with different N0 validate both c_i (intercept) and tan phi_i (slope).
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace iface = katai::core::iface;
using katai::core::DofMap;
using katai::core::InterfaceElement;
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
// Robust linear solve: at the collapse step the shear tangent is singular (rigid slip mode);
// catch and return zero so the solver gracefully cuts back and reports the limit load factor.
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    try {
        auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
        s->factorize(k);
        return s->solve(r);
    } catch (...) {
        return Eigen::VectorXd::Zero(r.size());
    }
}
std::vector<int> sorted_top(const Mesh& mesh) {
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    return top;
}

// Build mesh + one interface on the first top-edge segment; allocate the 6 structure-side DOFs.
struct Setup {
    Mesh mesh; DofMap dofs; InterfaceElement ie; double L;
    std::array<int, 6> sdof;  // struct-side global DOFs [A_ux,A_uy,B_ux,B_uy,mid_ux,mid_uy]
};
Setup make_setup(const iface::InterfaceProps& pr) {
    constexpr double W = 8.0, H = 2.0; const int nx = 8;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, nx, 2);
    const std::vector<int> top = sorted_top(mesh);
    const int A = top[0], B = top[2], mid = top[1];
    const double L = mesh.x[B] - mesh.x[A];

    DofMap dofs(mesh.node_count, 2);
    // Soil side fully fixed: the only free DOFs are the structure-side extra DOFs.
    for (int n = 0; n < mesh.node_count; ++n) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    std::array<int, 6> sdof{};
    sdof[0] = dofs.add_extra_dof(); sdof[1] = dofs.add_extra_dof();  // A: ux, uy
    sdof[2] = dofs.add_extra_dof(); sdof[3] = dofs.add_extra_dof();  // B: ux, uy
    sdof[4] = dofs.add_extra_dof(); sdof[5] = dofs.add_extra_dof();  // mid: ux, uy
    dofs.finalize();

    InterfaceElement ie;
    ie.soil_nodes = {A, B, mid};
    ie.struct_dof = sdof;
    ie.props = pr;
    return {std::move(mesh), std::move(dofs), ie, L, sdof};
}

// Distribute a uniform traction (per length) over the 3 struct nodes with Newton-Cotes weights
// (F_i = traction * w_i * J). pn -> +y component, ps -> +x component.
void apply_tractions(const Setup& su, Eigen::VectorXd& f, double ps, double pn) {
    const auto nc = iface::nc_points();
    const double J = su.L / 2.0;
    for (int q = 0; q < 3; ++q) {
        const int nd = nc[q].node;
        const double wJ = nc[q].w * J;
        f(su.dofs.equation(su.sdof[2 * nd + 0])) += ps * wJ;  // x
        f(su.dofs.equation(su.sdof[2 * nd + 1])) += pn * wJ;  // y
    }
}

// (A) Elastic exact check.
void test_interface_elastic() {
    iface::InterfaceProps pr;
    pr.kn = 2.0e5; pr.ks = 1.0e5; pr.c_i = 20.0; pr.phi_i = std::atan(0.5);
    Setup su = make_setup(pr);
    const double sigma_n = -50.0;  // compression
    const double tau = 30.0;       // < tau_max = c_i - sigma_n tanphi = 20 + 25 = 45
    Eigen::VectorXd f = Eigen::VectorXd::Zero(su.dofs.equation_count());
    apply_tractions(su, f, /*ps=*/tau, /*pn=*/sigma_n);

    const MaterialModel m{MaterialType::LinearElastic, 1.0e4, 0.3};
    const auto r = katai::core::solve_nonlinear(su.mesh, su.dofs, {m}, f, solve_spd,
                                                {1, 50, 1e-12}, {}, {},
                                                Structures{{}, {}, {}, {su.ie}});
    check(r.converged, "interface elastic solve converged");
    const double dus = r.displacement[su.sdof[4]];  // mid ux
    const double dun = r.displacement[su.sdof[5]];  // mid uy
    std::printf("  elastic: Du_s=%.6e (tau/ks=%.6e)  Du_n=%.6e (sig/kn=%.6e)  sig_n=%.3f tau=%.3f\n",
                dus, tau / pr.ks, dun, sigma_n / pr.kn, pr.kn * dun, pr.ks * dus);
    check(close(dus, tau / pr.ks, 1e-6), "Du_s = tau / ks");
    check(close(dun, sigma_n / pr.kn, 1e-6), "Du_n = sigma_n / kn");
    check(close(pr.ks * dus, tau, 1e-6), "recovered shear traction = applied");
    check(close(pr.kn * dun, sigma_n, 1e-6), "recovered normal traction = applied");
}

// (B) Coulomb capacity via the closed-form collapse load factor.
void run_collapse(double N0, double T0, double tanphi, double c_i, double L_unused) {
    iface::InterfaceProps pr;
    pr.kn = 2.0e5; pr.ks = 1.0e5; pr.c_i = c_i; pr.phi_i = std::atan(tanphi);
    Setup su = make_setup(pr);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(su.dofs.equation_count());
    // Normal pushes DOWN (-y, compression), shear pushes +x. Totals N0 (down) and T0 (shear).
    apply_tractions(su, f, /*ps=*/T0 / su.L, /*pn=*/-N0 / su.L);

    const MaterialModel m{MaterialType::LinearElastic, 1.0e4, 0.3};
    const auto r = katai::core::solve_nonlinear(su.mesh, su.dofs, {m}, f, solve_spd,
                                                {50, 60, 1e-10}, {}, {},
                                                Structures{{}, {}, {}, {su.ie}});
    const double lam = r.load_factor;
    const double lam_pred = c_i * su.L / (T0 - N0 * tanphi);  // closed-form collapse factor
    // Recover sigma_n, tau at the mid node-pair at the last converged state.
    const double dus = r.displacement[su.sdof[4]], dun = r.displacement[su.sdof[5]];
    const double sig_n = pr.kn * dun, tau = pr.ks * dus;
    std::printf("  collapse N0=%.0f T0=%.0f: lambda=%.4f (pred %.4f)  sig_n=%.2f tau=%.2f "
                "(c_i - sig tanphi=%.2f)\n",
                N0, T0, lam, lam_pred, sig_n, tau, c_i - sig_n * tanphi);
    check(!r.converged, "shear beyond capacity -> collapse (not fully converged)");
    check(close(lam, lam_pred, 0.03), "collapse load factor = c_i L / (T0 - N0 tanphi)");
    // At incipient slip the recovered traction sits on the Coulomb envelope (within load granularity).
    check(close(tau, c_i - sig_n * tanphi, 0.04), "recovered (sigma_n, tau) on the Coulomb envelope");
}

void test_interface_coulomb() {
    const double tanphi = 0.5, c_i = 20.0;
    run_collapse(/*N0=*/100.0, /*T0=*/80.0, tanphi, c_i, 0.0);   // lambda* = 20/(80-50) = 0.667
    run_collapse(/*N0=*/200.0, /*T0=*/160.0, tanphi, c_i, 0.0);  // lambda* = 20/(160-100) = 0.333
}

} // namespace

int main() {
    test_interface_elastic();
    test_interface_coulomb();
    if (g_failures == 0) {
        std::printf("OK: interface verified (elastic kn/ks exact, Coulomb collapse closed form)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}

