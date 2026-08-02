// PERFORMANCE BASELINE study (on-demand, EXCLUDE_FROM_ALL):
//   cmake --build build/msvc-rwdi --target study_perf && bin/study_perf
//
// Measures where the wall-clock time of the REAL GUI compute path (build_problem ->
// solve_nonlinear) goes, using NewtonResult::Timings (always-on instrumentation):
//   tan-asm : tangent assembly (K_T + f_int, once per Newton iteration)
//   res-asm : residual-only assembly (line-search trials)
//   csr     : COO -> CSR compile (sort + dedup, every iteration)
//   solve   : LinearSolve callback = symbolic analysis + factorization + back-solve
//             (currently a FRESH solver per iteration -- the prime optimization suspect)
//   other   : everything else inside solve_nonlinear (vector ops, norms, commits)
//
// Scenarios span the three cost regimes:
//   (1) LE large mesh      -- linear-solve dominated (1 step, ~2 iterations)
//   (2) MC footing collapse -- many Newton iterations + cutbacks (limit analysis)
//   (3) HS footing service  -- constitutive integration heavy (substepping, 40 steps)
//
// The numbers feed docs/validation/performance-baseline.md; every optimization step
// must (a) speed these up, (b) keep all ctest results bit-identical / within round-off.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <chrono>
#include <cstdio>
#include <string>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {

struct Row {
    std::string name;
    int elements = 0, nodes = 0, iters = 0;
    int n_tan = 0, n_res = 0, n_solve = 0;
    double load_factor = 0.0;
    katai::core::NewtonResult::Timings t;
    bool ok = false;          // solve ran (collapse runs are "ok" for profiling too)
};

void print_row(const Row& r) {
    const double other = r.t.total - r.t.assemble_tangent - r.t.assemble_residual -
                         r.t.csr_build - r.t.linear_solve;
    std::printf("%-26s %6d %7d %5d | %7.2f %7.2f %6.2f %8.2f %6.2f | %8.2f s\n",
                r.name.c_str(), r.elements, 2 * r.nodes, r.iters,
                r.t.assemble_tangent, r.t.assemble_residual, r.t.csr_build,
                r.t.linear_solve, other, r.t.total);
    std::printf("%-26s   asm x%d + ls x%d, solves x%d, load_factor=%.3f, pardiso %4.1f%%\n",
                "", r.n_tan, r.n_res, r.n_solve, r.load_factor,
                r.t.total > 0 ? 100.0 * r.t.linear_solve / r.t.total : 0.0);
}

Row run(const std::string& name, const m::Project& pr, double max_area, InitialPhase phase) {
    Row row; row.name = name;
    const auto M = katai::app::mesh_from_project(pr, max_area, 6);
    if (!M.ok) { std::printf("%-26s MESH FAILED\n", name.c_str()); return row; }
    row.elements = M.mesh.element_count; row.nodes = M.mesh.node_count;
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, phase);
    row.iters = R.iterations; row.t = R.timings; row.load_factor = R.load_factor;
    row.n_tan = R.timings.n_tangent; row.n_res = R.timings.n_residual;
    row.n_solve = R.timings.n_solve;
    row.ok = true;
    print_row(row);
    return row;
}

m::SoilPolygon rect(double W, double H, int mat) {
    m::SoilPolygon P; P.material = mat;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
    return P;
}

m::Load strip_load(double x0, double x1, double y, double q) {
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = x0; L.y1 = y; L.x2 = x1; L.y2 = y;
    L.qx1 = 0; L.qy1 = -q; L.qx2 = 0; L.qy2 = -q;
    return L;
}

// (1) Linear elastic, self-weight + surcharge on a wide block: the cost is almost
// entirely ONE tangent assembly + ONE linear solve -> isolates solve + assembly cost.
Row scenario_le(double max_area, const char* tag) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 3.0e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);
    pr.polygons.push_back(rect(40.0, 20.0, 0));
    pr.loads.push_back(strip_load(15.0, 25.0, 20.0, 100.0));
    return run(std::string("LE 40x20 gravity ") + tag, pr, max_area, InitialPhase::GravityLoading);
}

// (2) Mohr-Coulomb strip footing ramped PAST collapse (Reissner phi=20): incremental
// limit analysis -> many iterations, line-search activity, cutbacks. Weightless.
Row scenario_mc(double max_area, const char* tag) {
    const double c = 10.0, Nc = 14.83;
    m::Project pr;
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e4; s.nu = 0.3; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    s.c = c; s.phi = 20.0; s.psi = 0.0;
    pr.materials.push_back(s);
    pr.polygons.push_back(rect(6.0, 4.0, 0));
    pr.loads.push_back(strip_load(2.4, 3.6, 4.0, 1.8 * Nc * c));
    return run(std::string("MC footing collapse ") + tag, pr, max_area, InitialPhase::GravityLoading);
}

// (3) Hardening Soil strip footing at a service load (test_hs_footing via the GUI path):
// substepping constitutive integration + 40 load steps -> assembly/material heavy.
Row scenario_hs(double max_area, const char* tag) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::HardeningSoil;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.E50ref = 3.0e4; s.Eoedref = 3.0e4; s.Eurref = 9.0e4;
    s.m = 0.5; s.nu_ur = 0.2; s.p_ref = 100.0; s.Rf = 0.9;
    s.c = 5.0; s.phi = 32.0; s.psi = 2.0;
    pr.materials.push_back(s);
    pr.polygons.push_back(rect(12.0, 6.0, 0));
    pr.loads.push_back(strip_load(5.0, 7.0, 6.0, 150.0));
    return run(std::string("HS footing service ") + tag, pr, max_area, InitialPhase::K0Procedure);
}

} // namespace

int main() {
    std::printf("KATAI 2D performance baseline (GUI compute path)\n");
    std::printf("%-26s %6s %7s %5s | %7s %7s %6s %8s %6s | %8s\n",
                "scenario", "elems", "dofs", "iters",
                "tan-asm", "res-asm", "csr", "pardiso", "other", "total");
    std::printf("--------------------------------------------------------------------------"
                "--------------------\n");
    const auto t0 = std::chrono::steady_clock::now();

    scenario_le(0.50, "coarse");
    scenario_le(0.10, "fine");
    scenario_mc(0.06, "coarse");
    scenario_mc(0.03, "fine");
    scenario_hs(0.30, "coarse");
    scenario_hs(0.15, "fine");

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    std::printf("--------------------------------------------------------------------------"
                "--------------------\n");
    std::printf("total wall (incl. meshing): %.2f s\n", wall);
    return 0;
}
