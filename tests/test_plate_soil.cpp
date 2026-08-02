// Plate-in-soil coupling (P2.4): the 3-node Timoshenko plate (elements/plate.hpp) wired into
// the nonlinear solver -- translational DOFs SHARED with the soil mesh, rotation DOFs appended
// via DofMap::add_extra_dof. A plate strip runs along the top of a soil block.
//
// Validation: with VERY soft soil the plate (ends pinned, central point load) reduces to a
// simply-supported Timoshenko beam -> midspan deflection delta = P L^3/(48 EI) + P L/(4 kGA')
// (bending + shear). Recovering the soil stiffness (stiffer soil) reduces the deflection. The
// plate develops the simply-supported midspan moment M ~ P L/4. This exercises the full
// coupling: shared translational DOFs, extra rotation DOFs, and the solver plate assembly.
// (See docs/references/structural-plate-formulation.md; PLAXIS MMM sec 18.3 Timoshenko plate.)
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::PlateElement;
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

// Build the plate chain along the (x-sorted) top nodes, allocating one rotation DOF per top
// node. Returns the plate elements and the central top node id. nx = soil columns.
struct PlateSetup { std::vector<PlateElement> plates; int mid_node; std::vector<int> top; };
PlateSetup build_top_plate(const Mesh& mesh, DofMap& dofs,
                           const katai::core::plate::PlateProps& pr) {
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    const int nn = static_cast<int>(top.size());      // 2*nx + 1
    std::vector<int> rot(mesh.node_count, -1);
    for (int n : top) rot[n] = dofs.add_extra_dof();
    PlateSetup s; s.top = top; s.mid_node = top[(nn - 1) / 2];
    for (int e = 0; e * 2 + 2 < nn; ++e) {
        const int a = top[2 * e], b = top[2 * e + 2], m = top[2 * e + 1];  // A, B, mid
        s.plates.push_back(PlateElement{{a, b, m}, {rot[a], rot[b], rot[m]}, pr});
    }
    return s;
}

double run(double E_soil, double& mid_moment) {
    constexpr double W = 10.0, H = 2.0, P = 50.0;  // span 10 m, central point load
    const double E = 3.0e7, d = 0.4, nu = 0.15;    // stiff plate (wall)
    katai::core::plate::PlateProps pr;
    pr.EA = E * d; pr.EI = E * d * d * d / 12.0; pr.nu = nu;

    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 8, 2);

    DofMap dofs(mesh.node_count, 2);
    PlateSetup ps = build_top_plate(mesh, dofs, pr);
    // Bottom fixed, side rollers. Pin the plate ends (uy=0 at the two top corners; ux=0 at one).
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.fix_node_component(ps.top.front(), 1);  // left support uy=0
    dofs.fix_node_component(ps.top.back(), 1);   // right support uy=0
    dofs.finalize();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    const int eq = dofs.equation(dofs.global_dof(ps.mid_node, 1));
    f(eq) = -P;  // central point load (downward)

    const MaterialModel m{MaterialType::LinearElastic, E_soil, 0.3};
    const auto r = katai::core::solve_nonlinear(mesh, dofs, {m}, f, solve_spd,
                                                {1, 50, 1e-12}, {}, {},
                                                katai::core::Structures{ps.plates, {}});
    check(r.converged, "plate-in-soil solve converged");

    // Midspan deflection (downward).
    const double delta = -r.displacement[dofs.global_dof(ps.mid_node, 1)];

    // Recover the midspan bending moment from the plate element to the left of midspan
    // (its node B = mid_node, xi=+1).
    mid_moment = 0.0;
    for (const auto& pe : ps.plates) {
        if (pe.nodes[1] != ps.mid_node) continue;
        katai::core::plate::NodeCoords Xe;
        katai::core::plate::Dof u;
        for (int k = 0; k < 3; ++k) {
            Xe(k, 0) = mesh.x[pe.nodes[k]]; Xe(k, 1) = mesh.y[pe.nodes[k]];
            u(3 * k + 0) = r.displacement[dofs.global_dof(pe.nodes[k], 0)];
            u(3 * k + 1) = r.displacement[dofs.global_dof(pe.nodes[k], 1)];
            u(3 * k + 2) = r.displacement[pe.rot_dof[k]];
        }
        mid_moment = std::fabs(katai::core::plate::forces(Xe, pr, u, 1.0).M);
    }
    return delta;
}

void test_plate_soil() {
    constexpr double W = 10.0, P = 50.0;
    const double E = 3.0e7, d = 0.4, nu = 0.15;
    const double EI = E * d * d * d / 12.0;
    const double kGA = (5.0 / 6.0) * (E * d) / (2.0 * (1.0 + nu));
    const double delta_beam = P * W * W * W / (48.0 * EI) + P * W / (4.0 * kGA);
    const double M_beam = P * W / 4.0;  // simply-supported central moment

    double M_soft = 0.0, M_stiff = 0.0;
    const double d_soft = run(1.0e2, M_soft);    // very soft soil -> ~simply-supported beam
    const double d_stiff = run(5.0e4, M_stiff);  // stiffer soil -> less deflection
    std::printf("  soft soil: delta=%.5e (beam=%.5e) M_mid=%.3f (PL/4=%.3f)\n",
                d_soft, delta_beam, M_soft, M_beam);
    std::printf("  stiff soil: delta=%.5e M_mid=%.3f\n", d_stiff, M_stiff);

    check(close(d_soft, delta_beam, 0.06),
          "soft-soil limit = simply-supported Timoshenko beam (PL^3/48EI + PL/4kGA)");
    check(d_stiff < d_soft, "stiffer soil reduces midspan deflection (foundation support)");
    check(close(M_soft, M_beam, 0.10), "plate carries simply-supported midspan moment ~ PL/4");
}

} // namespace

int main() {
    test_plate_soil();
    if (g_failures == 0) {
        std::printf("OK: plate-in-soil coupling verified (shared translational DOFs + extra "
                    "rotation DOFs + solver assembly)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
