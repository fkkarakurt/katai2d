// Hardening Soil + GLOBAL Undrained-(A) at the BVP level -- a strip footing on undrained NC
// clay. The material-point undrained effective-stress path is validated in test_hs_undrained_clay
// (it terminates on the effective MC envelope, exact for psi=0); this test proves the GLOBAL
// undrained wrapper (pore-fluid bulk Kw/n added to the solver tangent + total stress) works with
// the Hardening Soil constitutive model, i.e. PLAXIS Undrained-(A) with HS.
//
// The defining undrained property is near-INCOMPRESSIBILITY: with nu_u = 0.495 the volumetric
// strain stays ~0 and the load is carried by excess pore pressure (u = -(Kw/n) eps_v), while the
// HS model still works on effective stress. We verify:
//   - convergence of the coupled HS + undrained footing;
//   - the loaded soil is near-incompressible (|eps_vol| << |deviatoric strain|);
//   - excess pore pressure develops under the footing (positive: contraction);
//   - the undrained footing is STIFFER than the same problem drained (immediate vs long-term).
#include <katai/analysis/initial_stress.hpp>
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
using katai::core::GaussState;
using katai::core::K0Options;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

katai::core::HardeningSoilParams clay() {
    katai::core::HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 8e3; p.Eur_ref = 24e3; p.Eoed_ref = 8e3;
    p.m = 1.0; p.nu_ur = 0.2; p.friction = 25 * kDeg;
    p.dilatancy = 0.0; p.cohesion = 10.0; p.Rf = 0.9;   // psi = 0 (undrained clay)
    katai::core::hs_calibrate_cap(p, 1.0 - std::sin(p.friction));
    return p;
}

// Strip footing on HS clay; returns settlement, max |eps_vol| and max excess pore pressure in
// the loaded zone, and convergence. drained=false -> Undrained (A).
struct Res { bool ok; double settle; double evol_max; double u_max; double lf; int iters; };
Res run(bool undrained, double q, double nu_u = 0.495, int steps = 40) {
    constexpr double W = 12.0, H = 6.0, gamma = 18.0, B = 2.0;
    constexpr double x0 = 0.5 * (W - B), x1 = 0.5 * (W + B);
    katai::core::HardeningSoilParams hs = clay();
    const double K0 = 1.0 - std::sin(hs.friction);

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 24, 12);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    std::vector<GaussState> init = katai::core::compute_k0_initial_stress(mesh, K0Options{H, gamma, K0});
    for (auto& gs : init) {
        const Eigen::Vector3d sig_cp(-gs.stress(1), -gs.stress(0), -gs.stress_zz);
        gs.pp = katai::core::hs_initial_pp(hs, sig_cp);
        gs.gamma_p = katai::core::hs_initial_gamma_p(hs, sig_cp);
    }
    Eigen::VectorXd geostatic = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_internal_force(mesh, dofs, init, geostatic);

    std::vector<int> footing;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] >= x0 - 1e-9 && mesh.x[n] <= x1 + 1e-9) footing.push_back(n);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, footing, 0.0, -q, f);

    MaterialModel m;
    m.type = MaterialType::HardeningSoil;
    m.youngs_modulus = hs.Eur_ref; m.poisson_ratio = hs.nu_ur;
    m.friction_angle = hs.friction; m.dilatancy_angle = hs.dilatancy; m.cohesion = hs.cohesion;
    m.hs = hs;
    m.undrained = undrained; m.undrained_poisson = nu_u;
    const std::vector<MaterialModel> mm = {m};

    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_unsym,
                                                {steps, 100, 1e-2}, init, {}, {}, {}, geostatic);
    double settle = 0.0;
    for (int n : footing) settle = std::min(settle, r.displacement[dofs.global_dof(n, 1)]);
    const double kwn = m.kw_over_n(0.495);
    double evol_max = 0.0, u_max = 0.0;
    for (const auto& gs : r.gauss_states) {
        evol_max = std::fmax(evol_max, std::fabs(gs.eps_vol));
        u_max = std::fmax(u_max, -kwn * gs.eps_vol);   // excess pore pressure (tension-pos eps_vol)
    }
    return {r.converged && r.load_factor > 0.999, settle, evol_max, u_max,
            r.load_factor, r.total_iterations};
}

void test_undrained_bvp() {
    // Moderate undrained load. The HS CONTINUUM tangent (linear convergence) limits the
    // attainable undrained load before the global Newton floors (the same limitation that makes
    // the drained footing slow -- the undrained penalty Kw/n amplifies it); a true ALGORITHMIC
    // consistent tangent would lift this ceiling. We verify the mechanism + undrained physics at
    // a converging load.
    const double q = 6.0;
    const Res un = run(true, q);
    const Res dr = run(false, q);
    std::printf("  HS strip footing q=%.0f kPa on NC clay (phi'=25, psi=0):\n", q);
    std::printf("    UNDRAINED: converged=%d lf=%.3f iters=%d settle=%.4e  max|eps_vol|=%.2e  max u=%.1f kPa\n",
                un.ok, un.lf, un.iters, un.settle, un.evol_max, un.u_max);
    std::printf("    DRAINED  : converged=%d lf=%.3f settle=%.4e\n", dr.ok, dr.lf, dr.settle);

    check(un.ok, "HS + Undrained-(A) footing BVP converges");
    check(un.evol_max < 2.0e-3,
          "loaded soil is near-incompressible under undrained loading (|eps_vol| small)");
    check(un.u_max > 5.0, "excess pore pressure develops under the footing (Kw/n carries the load)");
    check(dr.ok && std::fabs(un.settle) < std::fabs(dr.settle),
          "undrained footing is stiffer than drained (immediate vs long-term settlement)");
}

}  // namespace

int main() {
    test_undrained_bvp();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil + Undrained-(A) verified at the BVP level (strip footing on clay)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
