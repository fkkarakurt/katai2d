// The tolerated band of the hardening soils: an increment Newton cannot close at its tolerance is closed
// when its out-of-balance force is within 1e-2 of the force scale over EVERY degree of freedom, the
// supports with the free nodes.
//
// The force that will not fall sits in the soft modes of a developing mechanism. Measured on the wall of
// case `wall` below, Newton, the full elastic operator and a Krylov accelerator on it all held there,
// at about 2 kN against a thrust of 1000 kN, and that was 1.04e-2 of the force scale of the free equations
// (182 kN) -- while the internal force over every degree of freedom was 1170 kN. So the phase stopped at
// 31% of its movement and was reported as a limit. An embankment built undrained on soft clay stopped the
// same way at 56 to 87% of its load.
//
// verify: KV-NUM-018
//   oracle:   closed_form
//   source:   Prandtl, L. (1920). Ueber die Haerte plastischer Koerper. Nachr. Ges. Wiss. Goettingen, 74-85 -- the limit pressure of a strip on a weightless Tresca half-space, (2 + pi) c; Rankine's at-rest and the classical passive thrust as bounds on the wall; the same model with the band off for the rest
//   locator:  `wall`: a rigid block translated 150 mm into Hardening Soil sand (phi = 32, psi = 2, c = 0, 15-noded, 5442 nodes) past a joint of Rinter = 0.67, in two phases, by the default numerics. `small`: a 12 x 6 m block of the same model family under a small pressure. `strip`: a 16 x 8 m weightless Tresca block (c = 10 kPa), a 2 m strip under twice its limit pressure. `rules`: the band's two decisions on chosen numbers
//   quantity: the load factor reached, the increments kept and the force left in them; the thrust coefficient K = 2 P / (gamma H^2); the quotient of the force scale over every degree of freedom and over the free equations; the displacement; the pressure at collapse [kPa]
//   expected: `wall`: the whole movement, every kept increment within 1e-2 of a scale several times that of the free equations, K between its value at 50 mm and the passive limit; `small`: not one increment kept, the displacement of the same run with the band off; `strip`: the pressure is not carried, the collapse the same with and without the band
//   band:     1e-2 on the force kept, by construction; 3 on the quotient of the two scales; K(150 mm) in 4.2-5.0 (the measured 4.55-4.65 and the plateau 5.05-5.32 on this mesh); 1e-9 relative on the displacement of a run that converges; 1% between the collapses, within +12% of (2 + pi) c on this mesh
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::core::DofMap;
using katai::core::GaussState;
using katai::core::K0Options;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::NewtonOptions;
using katai::core::NewtonResult;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;
namespace linsolve = katai::linsolve;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    std::fflush(stdout);
    if (!ok) ++g_failures;
}

Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto solver = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    solver->factorize(k);
    return solver->solve(r);
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

NewtonOptions numerics(double tol, double band, int steps, int iters) {
    NewtonOptions o{steps, iters, tol};
    o.stagnation_accept = band;
    return o;
}

// ------------------------------------------------------------------------------------- the wall
void wall() {
    std::printf("-- a wall pushed into Hardening Soil sand, default numerics --\n");
    constexpr double kWallH = 5.0, kGamma = 18.0;
    const auto t0 = std::chrono::steady_clock::now();
    m::Project pr;
    m::Material sand; sand.name = "Sand"; sand.model = m::SoilModel::HardeningSoil;
    sand.gamma_unsat = sand.gamma_sat = kGamma;
    sand.E50ref = 4.0e4; sand.Eoedref = 4.0e4; sand.Eurref = 1.2e5;
    sand.m = 0.5; sand.nu_ur = 0.2; sand.p_ref = 100.0; sand.Rf = 0.9;
    sand.c = 0.0; sand.phi = 32.0; sand.psi = 2.0;
    sand.k0nc_auto = false; sand.k0nc = 0.47;
    sand.rinter_rigid = false; sand.Rinter = 0.67;
    pr.materials.push_back(sand);
    m::Material block; block.name = "Rigid wall"; block.model = m::SoilModel::LinearElastic;
    block.E = 1e8; block.nu = 0.2; block.gamma_unsat = block.gamma_sat = kGamma;
    block.k0_auto = false; block.k0 = 0.47;
    pr.materials.push_back(block);
    m::Material joint; joint.name = "Smooth base joint"; joint.model = m::SoilModel::MohrCoulomb;
    joint.E = 2e4; joint.nu = 0.3; joint.c = 0.0; joint.phi = 0.5; joint.psi = 0.0;
    joint.gamma_unsat = joint.gamma_sat = kGamma;
    pr.materials.push_back(joint);

    const int F = (int)m::BCType::Free, H = (int)m::BCType::HorizontallyFixed, X = (int)m::BCType::FullyFixed;
    m::SoilPolygon s; s.name = "Sand"; s.material = 0;
    s.x = {-1, 30, 30, 0, 0, -1, -1}; s.y = {-10, -10, 5, 5, 0, 0, -0.05};
    s.edge_bc = {X, H, F, F, F, F, H};
    m::SoilPolygon b; b.name = "Wall block"; b.material = 1;
    b.x = {-1, 0, 0, -1}; b.y = {0, 0, 5, 5}; b.edge_bc = {F, F, F, F};
    pr.polygons = {s, b};

    m::StructElement itf; itf.kind = m::StructKind::Interface; itf.name = "Wall interface";
    itf.x1 = 0; itf.y1 = 0; itf.x2 = 0; itf.y2 = 5; itf.material = -1; itf.iface_material = 0;
    itf.coarseness = 0.25;
    m::StructElement itb = itf; itb.name = "Base joint";
    itb.x1 = -1; itb.y1 = 0; itb.x2 = 0; itb.y2 = 0; itb.iface_material = 2;
    pr.structs = {itf, itb};

    const auto drive = [](double inc, const char* name) {
        m::PrescribedDisp D; D.name = name;
        D.x1 = -1; D.y1 = 0.025; D.x2 = -1; D.y2 = 5;
        D.set_ux = true; D.ux = inc; D.set_uy = true; D.uy = 0.0;
        D.coarseness = 0.25;
        return D;
    };
    pr.disps = {drive(0.05, "Push 1"), drive(0.10, "Push 2")};
    pr.mesh.elem_size = 2.0;
    pr.mesh.order = 15;
    pr.mesh.auto_refine = true;
    pr.initial.disp_active = {0, 0};
    m::Phase a; a.name = "Push to 50 mm"; a.type = m::PhaseType::Plastic; a.disp_active = {1, 0};
    m::Phase c = a; c.name = "Push to 150 mm"; c.disp_active = {0, 1};
    pr.phases = {a, c};

    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "the model meshes");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    std::printf("      %zu results, %.0f s\n", res.size(), seconds_since(t0));
    check(res.size() == 3, "all the phases ran");
    if (res.size() != 3) return;
    const auto thrust = [&](const katai::core::SolveResult& r) {
        double p = 0.0;
        for (int n = 0; n < r.mesh.node_count && 2 * n < (int)r.reaction.size(); ++n)
            if (std::fabs(r.mesh.x[n] + 1.0) < 1e-6 && r.mesh.y[n] > 0.02) p += r.reaction[2 * n];
        return 2.0 * p / (kGamma * kWallH * kWallH);
    };
    const double k50 = thrust(res[1]), k150 = thrust(res[2]);
    double kept_worst = 0.0;
    bool kept = false;
    for (const auto& d : res[2].diagnostics)
        if (d.code == std::string("K2D-A019")) {
            kept = true;
            const std::string& t = d.message;
            const size_t at = t.find("up to ");
            if (at != std::string::npos) kept_worst = std::atof(t.c_str() + at + 6);
        }
    const double ratio = res[2].convergence.global_error / std::max(res[2].convergence.force_error, 1e-300);
    std::printf("      K(50 mm) = %.4f, K(150 mm) = %.4f, load factor %.3f, kept: %d, worst %.3e, scales %.2f\n",
                k50, k150, res[2].load_factor, (int)kept, kept_worst, ratio);
    check(res[1].ok && res[2].ok && res[2].load_factor > 0.999, "the wall reaches its whole movement");
    check(kept && kept_worst > 0.0 && kept_worst <= 1e-2 + 1e-12, "increments were kept, within 1e-2 of the force scale");
    check(ratio > 3.0, "the internal force on the supports is in the scale");
    check(k50 > 3.0 && k50 < 4.0, "the thrust at 50 mm is where the verified runs have it");
    check(k150 > 4.2 && k150 < 5.0, "the thrust at 150 mm has risen towards the passive limit and is below it");
}

// ------------------------------------------------------------------------------------- a run that converges
void small() {
    std::printf("-- a run that converges: a small pressure, Newton closing every increment --\n");
    const auto t0 = std::chrono::steady_clock::now();
    constexpr double kW = 12.0, kH = 6.0, kGamma = 18.0;
    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 12, 6);
    std::vector<int> strip;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] >= 5.0 - 1e-9 && mesh.x[n] <= 7.0 + 1e-9) strip.push_back(n);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    katai::core::HardeningSoilParams hs;
    hs.p_ref = 100; hs.E50_ref = 3.0e4; hs.Eur_ref = 9.0e4; hs.Eoed_ref = 3.0e4;
    hs.m = 0.5; hs.nu_ur = 0.2; hs.friction = 32 * kDeg;
    hs.dilatancy = 2 * kDeg; hs.cohesion = 5.0; hs.Rf = 0.9;
    katai::core::hs_calibrate_cap(hs, 1.0 - std::sin(hs.friction));
    std::vector<GaussState> init = katai::core::compute_k0_initial_stress(
        mesh, K0Options{kH, kGamma, 1.0 - std::sin(hs.friction)});
    for (auto& gs : init) {
        const Eigen::Vector3d sig_cp(-gs.stress(1), -gs.stress(0), -gs.stress_zz);  // comp-pos
        gs.pp = katai::core::hs_initial_pp(hs, sig_cp);
        gs.gamma_p = katai::core::hs_initial_gamma_p(hs, sig_cp);
    }
    Eigen::VectorXd geostatic = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_internal_force(mesh, dofs, init, geostatic);
    MaterialModel mm;
    mm.type = MaterialType::HardeningSoil;
    mm.youngs_modulus = hs.Eur_ref; mm.poisson_ratio = hs.nu_ur;
    mm.friction_angle = hs.friction; mm.dilatancy_angle = hs.dilatancy; mm.cohesion = hs.cohesion;
    mm.hs = hs;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, strip, 0.0, -30.0, f);
    const NewtonResult a = katai::core::solve_nonlinear(mesh, dofs, {mm}, f, solve_unsym, numerics(1e-3, 0.0, 5, 60),
                                                        init, {}, {}, {}, geostatic);
    const NewtonResult b = katai::core::solve_nonlinear(mesh, dofs, {mm}, f, solve_unsym, numerics(1e-3, 1e-2, 5, 60),
                                                        init, {}, {}, {}, geostatic);
    const double du = (a.displacement - b.displacement).norm() / std::max(a.displacement.norm(), 1e-300);
    std::printf("      %d iterations, %d kept, relative difference %.2e (%.1f s)\n", b.total_iterations,
                b.stagnation_accepted, du, seconds_since(t0));
    check(a.converged && b.converged, "the small pressure is carried either way");
    check(b.stagnation_accepted == 0, "not one increment of a converging run is kept");
    check(du < 1e-9, "and its displacement is what it is without the band");
}

// ------------------------------------------------------------------------------------- past the capacity
NewtonResult tresca_past_limit(double q, const NewtonOptions& o) {
    const RectangularDomain domain{0.0, 0.0, 16.0, 8.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 32, 16);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<int> strip;
    for (int n : mesh.top_nodes)
        if (mesh.x[n] >= 7.0 - 1e-9 && mesh.x[n] <= 9.0 + 1e-9) strip.push_back(n);
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    katai::core::assemble_surface_traction(mesh, dofs, strip, 0.0, -q, f);
    const std::vector<MaterialModel> mats = {{MaterialType::MohrCoulomb, 1.0e4, 0.3, 10.0, 0.0, 0.0}};
    return katai::core::solve_nonlinear(mesh, dofs, mats, f, solve_unsym, o);
}

void strip_past_capacity() {
    std::printf("-- past the capacity: a band is not a way of carrying a load the ground cannot carry --\n");
    const auto t0 = std::chrono::steady_clock::now();
    const double limit = (2.0 + kPi) * 10.0, q = 2.0 * limit;
    const NewtonResult a = tresca_past_limit(q, numerics(1e-6, 0.0, 40, 60));
    const NewtonResult b = tresca_past_limit(q, numerics(1e-6, 1e-2, 40, 60));
    std::printf("      closed form %.2f kPa; collapse %.2f kPa without the band, %.2f with it, %d kept (%.1f s)\n",
                limit, a.load_factor * q, b.load_factor * q, b.stagnation_accepted, seconds_since(t0));
    check(!b.converged && b.load_factor < 0.999, "the pressure past the capacity is not carried");
    // Before kBandCarriedShare two increments were kept here that nothing carried, and the collapse
    // came out 2.5% above the one without the band: on the unsafe side.
    check(std::fabs(b.load_factor - a.load_factor) <= 0.01 * a.load_factor,
          "the collapse with the band is the one without it");
    check(b.load_factor * q <= 1.12 * limit, "and it is not above the closed form by more than this mesh's excess");
}

// ------------------------------------------------------------------------------------- the decisions
void rules() {
    std::printf("-- the band's decisions --\n");
    using katai::core::band_admits;
    using katai::core::band_stalled;
    constexpr double kBand = 1e-2, kScale = 1000.0;   // the band is 10
    constexpr double kOpen = 1e30;                     // no ceiling
    check(band_admits(9.9, kScale, kBand, 0.5, false, kOpen) && !band_admits(10.1, kScale, kBand, 0.5, false, kOpen),
          "the band is read against the scale");
    check(!band_admits(5.0, kScale, kBand, 0.01, false, kOpen), "under a load, not below the stiffness parameter of a collapse");
    check(band_admits(5.0, kScale, kBand, 0.01, true, kOpen), "under an imposed displacement, whatever the stiffness parameter");
    check(band_stalled(8.0, 9.0, kScale, kBand, 0.5, false, kOpen), "a best force not a quarter down is stalled");
    check(!band_stalled(6.0, 9.0, kScale, kBand, 0.5, false, kOpen), "one a third down is not");
    check(!band_stalled(12.0, 13.0, kScale, kBand, 0.5, false, kOpen), "nor one outside the band");
    // What a kept increment left and what the next adds: 2 left, 1 added, the ceiling 2.5.
    const double ceiling = katai::core::band_ceiling(2.0, 1.0);
    check(band_admits(2.2, kScale, kBand, 0.5, false, ceiling), "a stall that leaves about the same force is admitted");
    check(!band_admits(2.9, kScale, kBand, 0.5, false, ceiling), "a load nothing carries, added to the force, is not");
    check(band_admits(2.9, kScale, kBand, 0.5, true, ceiling), "and under an imposed displacement there is no such load");
}

}  // namespace

int main(int argc, char** argv) {
    std::printf("== the tolerated band of the hardening soils (KV-NUM-018) ==\n");
    const std::string which = argc > 1 ? argv[1] : "";
    if (which.empty() || which == "rules") rules();
    if (which.empty() || which == "small") small();
    if (which.empty() || which == "strip") strip_past_capacity();
    if (which.empty() || which == "wall") wall();
    if (g_failures == 0) {
        std::printf("OK: the tolerated band verified\n");
        return 0;
    }
    std::printf("%d check(s) failed\n", g_failures);
    return 1;
}
