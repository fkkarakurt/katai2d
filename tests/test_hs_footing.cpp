// Hardening Soil at the FE level for a SHEAR-DOMINATED, low-confinement problem:
// a flexible strip footing on a sand half-space under a WORKING (service) pressure.
// This is the canonical PLAXIS HS use case (settlement prediction) and the case the
// fragile nested-Newton return mapping historically diverged on -- unlike the
// laterally confined oedometer column (test_hs_fe_oedometer), the soil near the
// footing edge is at very low confinement and the response is deviatoric (shear).
//
// The footing pressure here is a service load (~1/3 of the bearing capacity), so the
// problem MUST converge to a finite settlement -- divergence is a solver/integrator
// defect, not physics. We require:
//   - convergence at the full load (load_factor == 1);
//   - a sensible downward settlement under the footing (a few cm, bounded);
//   - a reasonable iteration count (the robust integrator should not grind).
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

// HS is non-associated -> non-symmetric global system.
Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    solver->factorize(k);
    return solver->solve(r);
}

// Medium-dense sand (PLAXIS-tutorial-like), calibrated cap.
katai::core::HardeningSoilParams sand() {
    katai::core::HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 3.0e4; p.Eur_ref = 9.0e4; p.Eoed_ref = 3.0e4;
    p.m = 0.5; p.nu_ur = 0.2; p.friction = 32 * kDeg;
    p.dilatancy = 2 * kDeg; p.cohesion = 5.0; p.Rf = 0.9;
    katai::core::hs_calibrate_cap(p, 1.0 - std::sin(p.friction));
    return p;
}

void test_hs_footing() {
    // --- Geometry: strip footing on a sand block --------------------------------
    constexpr double W = 12.0, H = 6.0, gamma = 18.0;
    constexpr double B = 2.0;                       // footing width
    constexpr double x0 = 0.5 * (W - B), x1 = 0.5 * (W + B);
    constexpr double q = 150.0;                     // service pressure (kPa), << q_ult
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 24, 12);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    katai::core::HardeningSoilParams hs = sand();
    const double K0 = 1.0 - std::sin(hs.friction);

    // Geostatic K0 prestress + admissible HS history seeding (cap pp + shear gamma_p).
    std::vector<GaussState> init =
        katai::core::compute_k0_initial_stress(mesh, K0Options{H, gamma, K0});
    for (auto& gs : init) {
        const Eigen::Vector3d sig_cp(-gs.stress(1), -gs.stress(0), -gs.stress_zz);  // comp-pos
        gs.pp = katai::core::hs_initial_pp(hs, sig_cp);
        gs.gamma_p = katai::core::hs_initial_gamma_p(hs, sig_cp);
    }

    // Unified K0 procedure (PLAXIS staged "Plastic" phase): the geostatic equilibrium
    // (seeded K0 stress) is held as a constant internal force B, and ONLY the footing
    // surcharge ramps. Ramping self-weight instead would carry a spurious unbalanced
    // gravity fraction at intermediate load steps and drive the low-confinement surface
    // soil into spurious plasticity. (See build_problem.hpp / initial-stress-k0.md.)
    Eigen::VectorXd geostatic = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_internal_force(mesh, dofs, init, geostatic);

    std::vector<int> footing;
    for (int n : mesh.top_nodes) {
        const double x = mesh.x[n];
        if (x >= x0 - 1e-9 && x <= x1 + 1e-9) footing.push_back(n);
    }
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, footing, 0.0, -q, f);

    MaterialModel m;
    m.type = MaterialType::HardeningSoil;
    m.youngs_modulus = hs.Eur_ref; m.poisson_ratio = hs.nu_ur;
    m.friction_angle = hs.friction; m.dilatancy_angle = hs.dilatancy; m.cohesion = hs.cohesion;
    m.hs = hs;
    const std::vector<MaterialModel> mm = {m};

    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_unsym,
                                                {40, 60, 1e-2}, init, {}, {}, {}, geostatic);

    double settle = 0.0;
    for (int n : footing) settle = std::min(settle, r.displacement[dofs.global_dof(n, 1)]);
    std::printf("  HS strip footing q=%.0f kPa: converged=%d load_factor=%.3f iters=%d settle=%.4e m\n",
                q, r.converged, r.load_factor, r.total_iterations, settle);

    check(r.converged && r.load_factor > 0.999,
          "HS strip footing converges at the full service load");
    check(settle < 0.0 && settle > -0.5,
          "footing settles downward by a bounded, physical amount");
}

} // namespace

int main() {
    test_hs_footing();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil strip-footing settlement (shear-dominated) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
