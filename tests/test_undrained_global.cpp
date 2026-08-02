// GLOBAL Undrained (A): the pore-fluid bulk stiffness wired into the solver tangent
// and internal force (not just the local D_u identity of test_undrained). A laterally
// confined (oedometer) saturated column is loaded with a vertical total stress q under
// undrained conditions. The constitutive model works on EFFECTIVE stress; the solver
// uses the undrained tangent D_u = D' + (Kw/n) m m^T and total stress sigma = sigma' +
// (Kw/n) eps_v m, and carries the excess pore pressure u = -(Kw/n) eps_v in the Gauss
// state. With only vertical strain (eps_xx = eps_zz = 0):
//   constrained moduli  M' = D'(1,1),  M_u = M' + Kw/n
//   eps_yy = -q / M_u,   u_excess = (Kw/n) q / M_u,   sigma'_yy = -M' q / M_u
// so total sigma_yy = -q (equilibrium) and the share carried by pore water is
//   u / q = (Kw/n) / M_u  ->  1 as nu_u -> 0.5 (incompressible water).
// This is the classic 1D undrained result: immediately after loading, near all of the
// applied stress is taken by excess pore pressure (Skempton, Terzaghi t=0 condition).
// (See docs/references/effective-stress-formulation.md.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
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
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    solver->factorize(k);
    return solver->solve(r);
}

void test_confined_undrained_loading() {
    constexpr double W = 4.0, H = 10.0;
    constexpr double E = 1.0e4, nu = 0.3, nu_u = 0.495;
    constexpr double q = 100.0;  // applied vertical total stress (compression)

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 5);

    // Oedometer constraint: bottom fixed, side rollers (u_x = 0) -> eps_xx = 0.
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0, -q, f);

    // Undrained (A) linear-elastic material: effective E', nu' + nu_u for Kw/n.
    // (Named assignment, not positional aggregate init: MaterialModel grew tension
    // fields between dilatancy_angle and undrained, which would silently shift a
    // positional `true` onto tension_cutoff.)
    MaterialModel m;
    m.type = MaterialType::LinearElastic;
    m.youngs_modulus = E;
    m.poisson_ratio = nu;
    m.undrained = true;
    m.undrained_poisson = nu_u;
    const std::vector<MaterialModel> materials = {m};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, materials, f, solve_spd,
                                                {1, 50, 1e-12});
    check(r.converged, "undrained confined loading converged");

    // Analytic constrained moduli.
    const double Mp = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));  // M'
    const double kwn = m.kw_over_n(nu_u);
    const double Mu = Mp + kwn;                       // undrained constrained modulus
    const double eyy = -q / Mu;                       // uniform vertical strain
    const double u_exact = kwn * q / Mu;              // excess pore pressure
    const double sig_eff_exact = -Mp * q / Mu;        // effective vertical stress
    const double uy_undrained = eyy * H;              // surface settlement
    const double uy_drained = -q * H / Mp;            // drained (reference) settlement

    // (1) Surface settlement matches -q H / M_u (much smaller than drained).
    double uy_top = 0.0;
    for (int n : mesh.top_nodes)
        uy_top = std::min(uy_top, r.displacement[dofs.global_dof(n, 1)]);
    check(close(uy_top, uy_undrained, 1e-6), "undrained settlement = -qH/M_u");

    // (2) Per Gauss point: effective sigma'_yy, excess pore pressure, and equilibrium
    //     of the TOTAL vertical stress (sigma' + (Kw/n) eps_v = -q).
    double max_sig_err = 0.0, max_u_err = 0.0, max_tot_err = 0.0;
    for (const auto& gs : r.gauss_states) {
        const double u_excess = -kwn * gs.eps_vol;
        const double sig_total = gs.stress(1) + kwn * gs.eps_vol;  // = sigma' - u
        max_sig_err = std::fmax(max_sig_err, std::fabs(gs.stress(1) - sig_eff_exact));
        max_u_err = std::fmax(max_u_err, std::fabs(u_excess - u_exact));
        max_tot_err = std::fmax(max_tot_err, std::fabs(sig_total - (-q)));
    }
    check(max_sig_err < 1e-6 * q, "effective sigma'_yy = -M' q / M_u");
    check(max_u_err < 1e-6 * q, "excess pore pressure u = (Kw/n) q / M_u");
    check(max_tot_err < 1e-6 * q, "total sigma_yy = sigma' - u = -q (equilibrium)");

    const double pore_share = u_exact / q;  // = (Kw/n)/M_u = Skempton-like 1D ratio
    std::printf("  confined undrained: u/q=%.4f  sigma'_yy=%.3f  (q=%.0f)\n",
                pore_share, sig_eff_exact, q);
    std::printf("  settlement: undrained=%.3e  drained=%.3e  (%.1fx stiffer)\n",
                uy_undrained, uy_drained, uy_drained / uy_undrained);
    check(pore_share > 0.96 && pore_share < 1.0,
          "near all applied stress carried by pore water (nu_u=0.495)");
}

} // namespace

int main() {
    test_confined_undrained_loading();
    if (g_failures == 0) {
        std::printf("OK: global undrained (A) confined loading verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
