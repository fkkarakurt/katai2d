// DIAGNOSTIC (v0.4): why does a fresh model "not converge / collapse" unless the soil is strengthened?
// Reproduces the user's report through the FULL GUI compute path (mesh_from_project -> solve_gravity_le,
// GravityLoading = the plastic analysis) across the models a beginner actually draws: a strip (footing)
// load, a POINT load (the default drawn load), and a slope -- each on the default demo Sand (c'=1 kPa,
// near-cohesionless), with cohesion varied. Goal: find WHICH default action triggers non-convergence,
// and whether it is correct physics (weak soil) or a numerical artefact (a point-load singularity).
//
// Build:  cmake --build build/msvc-rwdi --target study_default_convergence
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cstdio>

namespace m = katai::model;
namespace app = katai::app;

namespace {
m::Material sand(double c) {
    m::Material s;
    s.name = "Sand"; s.model = m::SoilModel::MohrCoulomb;
    s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    s.E = 1.3e4; s.nu = 0.3; s.c = c; s.phi = 30.0; s.psi = 0.0;   // default demo Sand, c varied
    return s;
}
m::SoilPolygon block() {
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    return P;
}
void report(const char* label, const app::SolveResult& R, bool meshed) {
    if (!meshed) { std::printf("  %-34s : MESH FAILED\n", label); return; }
    std::printf("  %-34s : %-9s  equilibrated %5.1f%%\n",
                label, R.ok ? "CONVERGED" : "collapse", R.load_factor * 100.0);
}

void run_strip(const char* label, double c, double q_kPa, double half) {
    m::Project pr; pr.materials.push_back(sand(c)); pr.has_water = false;
    pr.polygons.push_back(block());
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 10.0 - half; L.y1 = 10.0; L.x2 = 10.0 + half; L.y2 = 10.0;
    L.qx1 = L.qx2 = 0.0; L.qy1 = L.qy2 = -q_kPa;
    pr.loads.push_back(L);
    const auto M = app::mesh_from_project(pr, 1.0, 6);
    report(label, M.ok ? app::solve_gravity_le(pr, M.mesh, app::InitialPhase::GravityLoading) : app::SolveResult{}, M.ok);
}
void run_point(const char* label, double c, double F_kN) {
    m::Project pr; pr.materials.push_back(sand(c)); pr.has_water = false;
    pr.polygons.push_back(block());
    m::Load L; L.kind = m::LoadKind::Point; L.x1 = 10.0; L.y1 = 10.0; L.qx1 = 0.0; L.qy1 = -F_kN;
    pr.loads.push_back(L);
    const auto M = app::mesh_from_project(pr, 1.0, 6);
    report(label, M.ok ? app::solve_gravity_le(pr, M.mesh, app::InitialPhase::GravityLoading) : app::SolveResult{}, M.ok);
}
void run_slope(const char* label, double c) {
    m::Project pr; pr.materials.push_back(sand(c)); pr.has_water = false;
    m::SoilPolygon P; P.material = 0;
    // A 1:1 (45 deg) cut slope, steeper than a cohesionless angle of repose (phi=30): base, right,
    // crest, slope face, toe, back. Fix the base; vertical rollers on the two vertical sides.
    P.x = {0, 20, 20, 12, 8, 0}; P.y = {0, 0, 8, 8, 4, 4};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    const auto M = app::mesh_from_project(pr, 1.0, 6);
    report(label, M.ok ? app::solve_gravity_le(pr, M.mesh, app::InitialPhase::GravityLoading) : app::SolveResult{}, M.ok);
}
}  // namespace

int main() {
    std::printf("v0.4 DIAGNOSTIC -- what makes a fresh model 'collapse'? (default demo Sand, c' varied)\n");
    std::printf("  Sand: E=13 MPa, nu=0.3, phi=30 deg, gamma=17/20\n\n");

    std::printf("(A) self-weight only, no load:\n");
    run_strip("  block, no load (q=0)", 1.0, 0.0, 2.0);

    std::printf("\n(B) 4 m strip footing q=50 kPa:\n");
    for (double c : {1.0, 5.0, 25.0}) { char l[48]; std::snprintf(l, sizeof(l), "  strip 4 m, q=50, c'=%g", c); run_strip(l, c, 50.0, 2.0); }

    std::printf("\n(C) narrow 1 m strip, higher q=200 kPa:\n");
    for (double c : {1.0, 25.0}) { char l[48]; std::snprintf(l, sizeof(l), "  strip 1 m, q=200, c'=%g", c); run_strip(l, c, 200.0, 0.5); }

    std::printf("\n(D) POINT load (the default drawn load) on the surface:\n");
    for (double F : {10.0, 50.0, 200.0}) { char l[48]; std::snprintf(l, sizeof(l), "  point F=%g kN, c'=1", F); run_point(l, 1.0, F); }
    run_point("  point F=50 kN, c'=25", 25.0, 50.0);

    std::printf("\n(E) 1:1 cut slope (45 deg > phi=30), self-weight:\n");
    for (double c : {1.0, 10.0, 25.0}) { char l[48]; std::snprintf(l, sizeof(l), "  slope 45deg, c'=%g", c); run_slope(l, c); }

    std::printf("\nWhich rows say 'collapse' pinpoints the real cause (point-load singularity vs slope vs\n"
                "genuine overload). That drives the v0.4 fix: better defaults + clearer UX + input guards.\n");
    return 0;
}
