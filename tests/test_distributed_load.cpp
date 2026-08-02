// Distributed (line) load -- consistent nodal forces for a linearly-varying surface traction
// (PLAXIS distributed load q1 -> q2 along a segment). Validates assemble_surface_traction_varying:
//   (A) a CONSTANT varying traction reproduces the uniform assemble_surface_traction exactly;
//   (B) the consistent nodal forces sum to the resultant (integral of the traction);
//   (C) a linear ramp q(x) = -p x/W gives the right resultant (-pW/2) AND the right line of
//       action (centroid at x = 2W/3) -- i.e. the load distribution, not just the total, is right.
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/mesh/mesh.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::geometry::RectangularDomain;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol * (1.0 + std::fabs(b)); }

void test_distributed_load() {
    constexpr double W = 6.0, H = 2.0, p = 50.0;
    const RectangularDomain domain{0.0, 0.0, W, H, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 6, 2);

    // Free DOFs everywhere (we only inspect the assembled force vector, no solve).
    DofMap dofs(mesh.node_count, 2);
    dofs.finalize();

    // Top edge chain, ordered along x (corner, mid, corner, ...).
    std::vector<int> top = mesh.top_nodes;  // structured mesh returns it ordered along x
    const int cs = (int)top.size();
    check(cs >= 3 && (cs - 1) % 2 == 0, "top edge chain is a valid quadratic chain");

    // (A) constant traction via the varying API == uniform API.
    Eigen::VectorXd fa = Eigen::VectorXd::Zero(dofs.equation_count());
    Eigen::VectorXd fb = Eigen::VectorXd::Zero(dofs.equation_count());
    std::vector<double> tx(cs, 0.0), ty(cs, -p);
    katai::core::assemble_surface_traction_varying(mesh, dofs, top, tx, ty, fa);
    katai::core::assemble_surface_traction(mesh, dofs, top, 0.0, -p, fb);
    double maxdiff = (fa - fb).cwiseAbs().maxCoeff();
    check(maxdiff < 1e-10 * p, "constant varying-traction == uniform assemble_surface_traction");

    // (B) resultant of the uniform load = -p W.
    double sum_y = 0.0;
    for (int n : top) sum_y += fa[dofs.global_dof(n, 1)];
    check(close(sum_y, -p * W, 1e-12), "uniform distributed load resultant = -p W");

    // (C) linear ramp q(x) = -p x/W: resultant = -p W/2, centroid x = 2W/3.
    Eigen::VectorXd fc = Eigen::VectorXd::Zero(dofs.equation_count());
    std::vector<double> rx(cs, 0.0), ry(cs);
    for (int i = 0; i < cs; ++i) ry[i] = -p * (mesh.x[top[i]] / W);
    katai::core::assemble_surface_traction_varying(mesh, dofs, top, rx, ry, fc);
    double R = 0.0, Mx = 0.0;
    for (int n : top) {
        const double fy = fc[dofs.global_dof(n, 1)];
        R += fy; Mx += fy * mesh.x[n];
    }
    const double centroid = Mx / R;
    std::printf("  distributed load: uniform R=%.3f (exact %.3f) | ramp R=%.3f (exact %.3f) "
                "centroid=%.4f (exact %.4f)\n",
                sum_y, -p * W, R, -0.5 * p * W, centroid, 2.0 * W / 3.0);
    check(close(R, -0.5 * p * W, 1e-12), "linear-ramp distributed load resultant = -p W/2");
    check(close(centroid, 2.0 * W / 3.0, 1e-10),
          "linear-ramp line of action at x = 2W/3 (distribution is correct, not just the total)");
}

}  // namespace

int main() {
    test_distributed_load();
    if (g_failures == 0) {
        std::printf("OK: distributed (line) load -- consistent varying-traction assembly verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
