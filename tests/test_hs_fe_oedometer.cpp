// Hardening Soil at the FE (BVP) level WITH the cap active. The two-surface return mapping
// (shear + cap) is now wired into the FE wrapper (hs_forward -> hs_return_principal); this
// test proves the cap reaches the assembled solver and produces compaction (volumetric
// cap) hardening, not just the earlier shear-only response.
//
// A laterally confined (oedometer) column sits under a geostatic K0 initial stress that
// balances self-weight (zero displacement from gravity alone). The cap pre-consolidation
// pressure pp is initialised to the equivalent isotropic stress of sigma0 (NC: on the cap,
// hs_initial_pp). A surface surcharge q is then applied: the column compacts, the cap yields
// (pp grows) and the lateral stress follows with K0 < 1. We verify:
//   - convergence and downward (compaction) settlement;
//   - the cap is ACTIVE -> pp increases above its geostatic value at loaded Gauss points;
//   - cap-ON settles MORE than cap-OFF (the volumetric cap adds plastic compressibility);
//   - the developed stress ratio is oedometer-like (0 < sigma'_h/sigma'_v < 1).
// (See hardening-soil-formulation.md sec 7; single-point parity is test_hs_berlin / _calibration.)
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
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

// HS is non-associated (asymmetric consistent tangent) -> non-symmetric global system.
Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    solver->factorize(k);
    return solver->solve(r);
}

katai::core::HardeningSoilParams berlin_calibrated() {
    katai::core::HardeningSoilParams p;
    p.p_ref = 100; p.E50_ref = 105e3; p.Eur_ref = 315e3; p.Eoed_ref = 105e3;
    p.m = 0.55; p.nu_ur = 0.2; p.friction = 38 * kPi / 180;
    p.dilatancy = 6 * kPi / 180; p.cohesion = 1.0; p.Rf = 0.9;
    katai::core::hs_calibrate_cap(p, 0.38);
    return p;
}

// Run the confined column with the given cap state; return surface settlement and the
// fraction of Gauss points whose cap hardened (pp grew) + the mean developed K0.
struct Result { bool ok; double settle; double cap_active_frac; double k0_mean; };
Result run(bool cap_on) {
    constexpr double W = 2.0, H = 8.0, gamma = 18.0, q = 250.0;
    katai::core::HardeningSoilParams hs = berlin_calibrated();
    const double K0 = 1.0 - std::sin(hs.friction);  // Jaky
    if (!cap_on) hs.cap_beta = 0.0;  // disable cap (shear-only)

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 2, 8);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();

    // Geostatic K0 prestress (balances gravity); seed cap pp = p_eq(sigma0) (NC).
    std::vector<GaussState> init =
        katai::core::compute_k0_initial_stress(mesh, K0Options{H, gamma, K0});
    for (auto& gs : init) {
        const Eigen::Vector3d sig_cp(-gs.stress(1), -gs.stress(0), -gs.stress_zz);  // comp-pos
        gs.pp = cap_on ? katai::core::hs_initial_pp(hs, sig_cp) : 0.0;
        gs.gamma_p = katai::core::hs_initial_gamma_p(hs, sig_cp);  // K0 state on shear surface
    }
    // Unified K0 procedure: the geostatic equilibrium (K0 prestress) is held as a constant
    // internal force and ONLY the surface surcharge ramps -- ramping self-weight instead
    // double-counts gravity during the increments and yields spuriously once the soil is
    // plastic (see build_problem.hpp / initial-stress-k0.md). A PLAXIS-realistic tolerated
    // error (1%) is used (the continuum HS tangent gives linear, not quadratic, convergence).
    Eigen::VectorXd geostatic = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_internal_force(mesh, dofs, init, geostatic);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, mesh.top_nodes, 0.0, -q, f);

    MaterialModel m;
    m.type = MaterialType::HardeningSoil;
    m.youngs_modulus = hs.Eur_ref; m.poisson_ratio = hs.nu_ur;
    m.friction_angle = hs.friction; m.dilatancy_angle = hs.dilatancy; m.cohesion = hs.cohesion;
    m.hs = hs;
    const std::vector<MaterialModel> mm = {m};

    const auto r = katai::core::solve_nonlinear(mesh, dofs, mm, f, solve_unsym,
                                                {20, 80, 1e-2}, init, {}, {}, {}, geostatic);
    std::printf("    [cap_on=%d converged=%d load_factor=%.3f iters=%d]\n",
                cap_on, r.converged, r.load_factor, r.total_iterations);

    double settle = 0.0;
    for (int n : mesh.top_nodes)
        settle = std::min(settle, r.displacement[dofs.global_dof(n, 1)]);
    int active = 0, n_k0 = 0; double k0_sum = 0.0;
    for (std::size_t i = 0; i < r.gauss_states.size(); ++i) {
        if (cap_on && r.gauss_states[i].pp > init[i].pp * 1.001) ++active;
        const double sv = r.gauss_states[i].stress(1), sh = r.gauss_states[i].stress(0);
        if (sv < -1.0) { k0_sum += sh / sv; ++n_k0; }  // compression (tension-pos negative)
    }
    return {r.converged, settle,
            cap_on ? double(active) / r.gauss_states.size() : 0.0,
            n_k0 ? k0_sum / n_k0 : 0.0};
}

void test_hs_fe_oedometer() {
    const Result on = run(true);
    const Result off = run(false);
    std::printf("  cap-ON : converged=%d settle=%.4e cap-active-frac=%.2f K0_mean=%.3f\n",
                on.ok, on.settle, on.cap_active_frac, on.k0_mean);
    std::printf("  cap-OFF: converged=%d settle=%.4e K0_mean=%.3f\n",
                off.ok, off.settle, off.k0_mean);
    check(on.ok, "HS FE confined column (cap on) converged");
    check(off.ok, "HS FE confined column (cap off) converged");
    check(on.settle < 0.0, "surcharge produces compaction settlement (downward)");
    check(on.cap_active_frac > 0.3, "cap is ACTIVE in the FE (pp hardened at loaded Gauss pts)");
    check(on.settle < off.settle, "cap-on settles more than cap-off (volumetric cap plasticity)");
    check(on.k0_mean > 0.2 && on.k0_mean < 1.0, "oedometer-like stress ratio 0 < sigma_h/sigma_v < 1");
}

} // namespace

int main() {
    test_hs_fe_oedometer();
    if (g_failures == 0) {
        std::printf("OK: Hardening Soil two-surface (cap+shear) verified at the FE level\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
