// Coupled embedded wall + soil excavation (P2.4 supported excavation). A flexural plate
// (Timoshenko beam) is placed along a vertical line in a continuous Mohr-Coulomb soil under
// its K0 geostatic stress; the plate's translational DOFs are shared with the soil (bonded)
// and a rotation DOF is added per wall node. The soil on one side is then excavated to a
// depth, mobilising a net active thrust on the wall, which deflects toward the excavation.
//
// Two checks exercise the full coupled machinery (plate-in-soil + K0 initial stress + staged
// excavation + MC plasticity) end to end:
//   (A) No excavation: K0 is self-equilibrated on the continuous mesh, so the staged-release
//       load (f_int0 -> gravity) is identically zero and the wall does NOT move. This guards
//       against any spurious wall-line thrust from the coupling.
//   (B) Excavation: the wall deflects toward the cut, monotonically increasing with height
//       (a cantilever rotating about its embedded toe), and the solve converges.
//
// Staged excavation uses the SigmaMstage ramp: target(lambda) = f_int0 + lambda*(grav-f_int0),
// i.e. constant_force = f_int0 (initial K0 internal force, residual 0 at lambda=0) and the
// ramped external load f_ext = grav_active - f_int0 (the excavation release). A small cohesion
// regularises the re-entrant excavation corner so the analysis converges to full depth.
// (Reference: PLAXIS staged construction; Terzaghi/Rankine earth pressure; cf. test_excavation.)
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/staged_construction.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

using namespace katai::core;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
Eigen::VectorXd solve_unsym(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// Reduced internal force f_int0 = integral_active B^T sigma0 (u=0) over free DOFs. At u=0 the
// plate contribution vanishes, so this is the soil K0 seed force.
Eigen::VectorXd internal_force0(const Mesh& mesh, const DofMap& dofs,
                                const std::vector<GaussState>& gs,
                                const std::vector<char>& active) {
    namespace tri6 = katai::core::tri6;
    Eigen::VectorXd F = Eigen::VectorXd::Zero(dofs.equation_count());
    const auto gp = tri6::gauss_points();
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active[e]) continue;
        tri6::NodeCoords X;
        for (int k = 0; k < 6; ++k) { X(k, 0) = mesh.x[mesh.node_of(e, k)]; X(k, 1) = mesh.y[mesh.node_of(e, k)]; }
        Eigen::Matrix<double, 12, 1> fe = Eigen::Matrix<double, 12, 1>::Zero();
        for (int g = 0; g < 3; ++g) {
            const auto sg = tri6::strain_displacement(X, gp[g].xi, gp[g].eta);
            fe.noalias() += (gp[g].weight * sg.det_jacobian) * sg.B.transpose() * gs[e * 3 + g].stress;
        }
        for (int k = 0; k < 6; ++k)
            for (int c = 0; c < 2; ++c) {
                const int eq = dofs.equation(dofs.global_dof(mesh.node_of(e, k), c));
                if (eq >= 0) F(eq) += fe(2 * k + c);
            }
    }
    return F;
}

struct WallResult {
    std::vector<double> y, ux;
    bool converged = false;
    // The convergence criteria family at the last accepted iterate. This fixture is here for
    // the MOMENT criterion (Eq. 9-3/9-4): it is the only nonlinear case in the tree that
    // carries rotational freedom, so it is the only place that criterion can be read at all.
    NewtonResult::Convergence conv;
    int iterations = 0;
};

WallResult run(bool excavate, double tol = 1e-6) {
    constexpr double W = 20.0, Htot = 12.0, xw = 10.0, y_exc = 9.0, gamma = 18.0;
    const double pi = std::acos(-1.0);
    const double phi = 30.0 * pi / 180.0, c = 2.0, psi = phi;  // associated; small c regularises corner
    const double K0 = 1.0 - std::sin(phi);
    constexpr double E = 2.0e4, nu = 0.3;
    const double Ep = 2.0e7, d = 0.4;
    const double EA = Ep * d, EI = Ep * d * d * d / 12.0;

    RectangularDomain domain{0.0, 0.0, W, Htot, 0};
    Mesh mesh = katai::mesh::generate_structured_tri6(domain, 10, 6);  // x=10 -> fine column line

    // Bonded plate along the wall line (continuous mesh): translational DOFs shared with the
    // soil, a rotation extra DOF per wall node. Plate triples [corner, corner, mid] up the line.
    std::vector<std::pair<double, int>> wln;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - xw) < 1e-9) wln.push_back({mesh.y[n], n});
    std::sort(wln.begin(), wln.end());

    DofMap dofs(mesh.node_count, 2);
    plate::PlateProps pp; pp.EA = EA; pp.EI = EI; pp.nu = 0.15;
    std::vector<int> wnode, wrot;
    for (auto& p : wln) { wnode.push_back(p.second); wrot.push_back(dofs.add_extra_dof()); }
    std::vector<PlateElement> plates;
    for (size_t e = 0; 2 * e + 2 < wnode.size(); ++e) {
        const size_t a = 2 * e, m = 2 * e + 1, b = 2 * e + 2;
        PlateElement pe;
        pe.nodes = {wnode[a], wnode[b], wnode[m]};
        pe.rot_dof = {wrot[a], wrot[b], wrot[m]};
        pe.props = pp;
        plates.push_back(pe);
    }

    // Excavation: deactivate LEFT elements above the excavation level (centroid based).
    std::vector<char> active(mesh.element_count, 1);
    for (int e = 0; e < mesh.element_count && excavate; ++e) {
        double xc = 0.0, yc = 0.0;
        for (int k = 0; k < 6; ++k) { xc += mesh.x[mesh.node_of(e, k)]; yc += mesh.y[mesh.node_of(e, k)]; }
        if (xc / 6.0 < xw && yc / 6.0 > y_exc) active[e] = 0;
    }

    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes)  dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    fix_inactive_nodes(mesh, active, dofs);
    dofs.fix(wrot.front());  // toe rotation fixed (cantilever base)
    dofs.finalize();

    const std::vector<GaussState> init =
        compute_k0_initial_stress(mesh, K0Options{Htot, gamma, K0});
    Eigen::VectorXd grav = Eigen::VectorXd::Zero(dofs.equation_count());
    assemble_gravity(mesh, dofs, {gamma}, grav, active);

    const Eigen::VectorXd fint0 = internal_force0(mesh, dofs, init, active);
    const Eigen::VectorXd f_ramp = grav - fint0;
    const std::vector<MaterialModel> mm = {{MaterialType::MohrCoulomb, E, nu, c, phi, psi}};
    const auto r = solve_nonlinear(mesh, dofs, mm, f_ramp, solve_unsym, {100, 60, tol}, init,
                                   active, Structures{plates, {}, {}, {}}, {}, fint0);

    WallResult out;
    out.converged = r.converged;
    out.conv = r.convergence;
    out.iterations = r.total_iterations;
    for (size_t i = 0; i < wnode.size(); ++i) {
        out.y.push_back(mesh.y[wnode[i]]);
        out.ux.push_back(r.displacement[dofs.global_dof(wnode[i], 0)]);
    }
    return out;
}

void test_excavation_wall() {
    // (A) No excavation -> K0 self-equilibrium: the wall must not move (no spurious thrust).
    const WallResult a = run(false);
    check(a.converged, "no-excavation phase converged");
    double max_abs_a = 0.0;
    for (double u : a.ux) max_abs_a = std::fmax(max_abs_a, std::fabs(u));
    std::printf("  no-excavation: max|wall ux| = %.3e (K0 self-equilibrium)\n", max_abs_a);
    check(max_abs_a < 1e-9, "no excavation -> wall does not move (K0 self-equilibrium)");

    // (B) Excavation -> wall deflects toward the cut, monotone with height, converged.
    const WallResult b = run(true);
    check(b.converged, "excavation phase converged");
    const double tip = b.ux.back();  // top of wall
    std::printf("  excavation: tip ux = %+.4e (toward cut = negative), profile:\n", tip);
    bool all_toward = true, monotone = true;
    for (size_t i = 0; i < b.y.size(); ++i) {
        if (b.ux[i] > 1e-12) all_toward = false;
        if (i > 0 && std::fabs(b.ux[i]) < std::fabs(b.ux[i - 1]) - 1e-12) monotone = false;
        if (std::fmod(b.y[i], 2.0) < 1e-9)
            std::printf("     y=%5.2f  ux=%+.4e\n", b.y[i], b.ux[i]);
    }
    check(all_toward, "excavation -> all wall displacements toward the cut (ux <= 0)");
    check(monotone, "excavation -> |ux| increases monotonically from toe to top");
    check(tip < -1e-3 && tip > -5e-2, "tip deflection in a physical range for the excavation");
    check(std::fabs(tip) > 5.0 * max_abs_a, "excavation deflection >> no-excavation (excavation-driven)");
}

// The MOMENT criterion (Eq. 9-3/9-4), read on the only kind of case that has rotational
// freedom at all. Until 2026-08-25 the record said this criterion "participates in the gate";
// it did not -- the gate called local_ok(), which does not contain it -- and no case had ever
// been examined through it. Both halves are closed here: it is in the gate now
// (Convergence::enforced_ok), and this is what it reads.
//
// The reading is negative, and that is the point of writing it down. With an ELASTIC plate the
// rotational rows are LINEAR in the plate DOFs, so the linear solve zeroes them exactly and the
// criterion sits at round-off at EVERY tolerance -- including the loose run whose answer is a
// fifth wrong. A criterion that reports "fine" next to a wrong answer is not broken; it is
// answering a different question, and a guard that does not say which question it answers will
// eventually be read as covering the other one.
void test_moment_criterion() {
    std::printf("\n-- moment criterion (Eq. 9-3/9-4) on the wall: what it reads, and what it does not --\n");
    const WallResult tight = run(true, 1e-6);
    const WallResult loose = run(true, 1e-1);
    for (const auto* p : {&tight, &loose}) {
        const auto& c = p->conv;
        std::printf("   tol=%.0e  tip ux=%+.6e  iters=%-4d | moment=%.3e of %.3e tolerated %-4s "
                    "| force=%.3e\n",
                    c.tolerated, p->ux.back(), p->iterations, c.moment_error, c.tolerated,
                    c.moment_ok() ? "ok" : "FAIL", c.force_error);
    }
    check(tight.conv.has_moment && loose.conv.has_moment,
          "the moment criterion is MEASURED wherever a rotational freedom exists (has_moment)");
    // Round-off, not merely "small": four decades below the loosest tolerance the tree ships,
    // so this is a claim about the structure of the equations rather than about a threshold.
    check(tight.conv.moment_error < 1e-10 && loose.conv.moment_error < 1e-10,
          "an ELASTIC plate's rotational rows are linear -> the moment residual is at round-off "
          "at every tolerance (measured 8.5e-15 .. 1.7e-14 over 1e-1 .. 1e-10)");
    // ... and the same loose run is materially wrong, which is what makes the line above a
    // limitation and not a pass mark. What catches this run is the global force gate.
    const double err = std::fabs(loose.ux.back() - tight.ux.back()) / std::fabs(tight.ux.back());
    std::printf("   the loose run is %.1f%% away from the converged tip deflection, and the "
                "moment criterion calls it ok\n", 100.0 * err);
    check(err > 0.05,
          "the tol=1e-1 run IS materially wrong (>5%), so the criterion's 'ok' above is a "
          "statement about the rotational rows and NOT about the answer");
}

} // namespace

int main() {
    test_excavation_wall();
    test_moment_criterion();
    if (g_failures == 0) {
        std::printf("OK: coupled embedded wall + soil excavation verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
