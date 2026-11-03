// Soft Soil Creep in a consolidation phase: the creep of a time step carries the load back into the water.
//
// Time passes in the stress points during a time step, and the step is solved as an increment from the force of the
// state it starts in. That base force was read with the step's time, so the step's relaxation at zero strain was in
// the base as well as in the iterate, and the residual never saw it: the stress relaxed out of equilibrium, with no
// settlement to show for it. The settlement at 1000 days came out -42% at 50 steps and -65% at 800, worse as the step
// was refined, and the combination had to be refused. With no time in the base force the same column converges.
//
// verify: KV-CON-011
//   oracle:   independent_path
//   source:   Vermeer, P.A. & Neher, H.P. (1999). A soft soil model that accounts for creep. Beyond 2000 in Computational Geotechnics, Balkema, 249-261, on Terzaghi's column: Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley, ch. XIII; solved independently of the program by a finite-volume column in space with an implicit integrator in time on the same rate equations
//   locator:  a confined column 1 m high, drained at the top, sigma_v 100 kPa at K0 0.5, POP 30, 100 kPa more placed at t = 0, lam* 0.1, kap* 0.02, mu* 0.005, nu 0.15, c 1, phi 30, K0nc 0.5, k 1e-3 m/day, Kw 2e6, n 0.5; 15-noded, eight layers
//   quantity: the settlement of the top at 1000 days [m]; the effective vertical stress near mid-height [kPa]
//   expected: 0.08305 m; 200 kPa (the excess pore pressure is gone)
//   band:     0.5% on the settlement, at 50 and at 200 steps (measured -0.36% and -0.11%); 0.5% on the stress
#include <katai/analysis/consolidation.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/soft_soil.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace katai::core;
namespace linsolve = katai::linsolve;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

ConsolidationSolveFactory nonsym_factory() {
    return [](const katai::math::CsrMatrix& A) {
        std::shared_ptr<linsolve::DirectSolver> s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
        s->factorize(A);
        return std::function<Eigen::VectorXd(const Eigen::VectorXd&)>([s](const Eigen::VectorXd& b) { return s->solve(b); });
    };
}

struct Run { bool ok = false; double settlement = 0.0, syy_mid = 0.0; };

Run column(int steps) {
    constexpr double H = 1.0, sv0 = 100.0, K0 = 0.5, pop = 30.0, q = 100.0, T = 1000.0;
    constexpr int ny = 8;
    const double W = H / ny;
    const katai::geometry::RectangularDomain dom{0.0, 0.0, W, H, 0};
    katai::mesh::Mesh mesh = katai::mesh::generate_structured_tri15(dom, 1, ny);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.finalize();
    std::vector<char> drained(mesh.node_count, 0);
    for (int n : mesh.top_nodes) drained[n] = 1;

    MaterialModel m;
    m.type = MaterialType::SoftSoilCreep;
    m.ssc.lam_star = 0.1; m.ssc.kap_star = 0.02; m.ssc.mu_star = 0.005; m.ssc.nu_ur = 0.15; m.ssc.c = 1.0;
    m.ssc.phi = 30.0 * kPi / 180.0; m.ssc.psi = 0.0; m.ssc.K0nc = 0.5;
    const softsoil::Params ss = m.ssc.ss();
    const std::vector<MaterialModel> mm = {m};
    const std::vector<Permeability> perm = {{1e-3, 1e-3}};

    const int ngp = Tri15Element::kGaussCount;
    std::vector<GaussState> st((size_t)mesh.element_count * ngp);
    for (auto& s : st) {
        s.stress = Eigen::Vector3d(-K0 * sv0, -sv0, 0.0);
        s.stress_zz = -K0 * sv0;
        const Eigen::Vector3d past(ss.K0nc * (sv0 + pop), sv0 + pop, ss.K0nc * (sv0 + pop));
        const Eigen::Vector3d cur(K0 * sv0, sv0, K0 * sv0);
        s.pp = std::max(softsoil::ss_initial_pp(ss, cur, 1.0), softsoil::ss_initial_pp(ss, past, 1.0));
    }
    // the consistent nodal forces of the pressure on the top, through the Lagrange polynomials along that edge
    Eigen::VectorXd df = Eigen::VectorXd::Zero(dofs.equation_count());
    std::vector<int> top = mesh.top_nodes;
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    const int nn = (int)top.size();
    const double gx[8] = {-0.9602898565, -0.7966664774, -0.5255324099, -0.1834346425, 0.1834346425, 0.5255324099, 0.7966664774, 0.9602898565};
    const double gw[8] = {0.1012285363, 0.2223810345, 0.3137066459, 0.3626837834, 0.3626837834, 0.3137066459, 0.2223810345, 0.1012285363};
    for (int i = 0; i < nn; ++i) {
        double w = 0.0;
        for (int g = 0; g < 8; ++g) {
            const double x = 0.5 * W * (gx[g] + 1.0);
            double L = 1.0;
            for (int j = 0; j < nn; ++j)
                if (j != i) L *= (x - mesh.x[top[j]]) / (mesh.x[top[i]] - mesh.x[top[j]]);
            w += gw[g] * 0.5 * W * L;
        }
        const int eq = dofs.equation(dofs.global_dof(top[i], 1));
        if (eq >= 0) df(eq) += -w * q;
    }
    std::vector<double> pore(mesh.node_count, 0.0);
    const PoreFluidStiffness kw_over_n(2.0e6 / 0.5);
    const auto r = solve_consolidation_plastic(mesh, dofs, mm, perm, 9.81, kw_over_n, drained, st, pore, T / steps, steps, {},
                                               &df, nonsym_factory(), 40, 1e-6);
    Run out;
    out.ok = r.converged && !r.series.displacement.empty();
    if (!out.ok) return out;
    out.settlement = -r.series.displacement.back()(dofs.global_dof(mesh.top_nodes[0], 1));
    int gmid = 0;
    double best = 1e300;
    for (int e = 0; e < mesh.element_count; ++e) {
        double ym = 0.0;
        for (int i = 0; i < 3; ++i) ym += mesh.y[mesh.node_of(e, i)] / 3.0;
        if (std::fabs(ym - 0.5 * H) < best) { best = std::fabs(ym - 0.5 * H); gmid = e * ngp; }
    }
    out.syy_mid = -r.committed[(size_t)gmid].stress(1);
    return out;
}
}  // namespace

int main() {
    std::printf("Soft Soil Creep consolidating: against an independent column, s(1000 d) = 0.08305 m\n\n");
    for (const int steps : {50, 200}) {
        const Run r = column(steps);
        check(r.ok, std::to_string(steps) + " steps: the phase closes");
        std::printf("   %d steps: settlement %.5f m (%+.2f%%), effective vertical stress %.2f kPa\n", steps, r.settlement,
                    100.0 * (r.settlement / 0.08305 - 1.0), r.syy_mid);
        check(std::fabs(r.settlement / 0.08305 - 1.0) < 5e-3, std::to_string(steps) + " steps: the settlement");
        check(std::fabs(r.syy_mid / 200.0 - 1.0) < 5e-3, std::to_string(steps) + " steps: the load has come to the skeleton");
    }
    if (g_failures == 0) {
        std::printf("\nOK: the creep of a time step is in its equilibrium, and the answer converges with the step\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
