// A phase that stops short keeps the state it last equilibrated.
//
// A strip footing on a weightless Tresca soil is loaded to 1.2 x its Prandtl capacity (2 + pi) c B,
// so the phase MUST stop short. Before 2026-09-29 the stopped phase published only its load factor:
// no displacement, no stress, no reactions -- and a displacement-controlled footing that stopped at
// 69% of its settlement reported a footing force of 0.0. The engineer reads the collapse load, and
// the mechanism, off exactly that last equilibrated state, so the phase now carries it (and still
// reports failure; nothing downstream continues from it).
//
// Checked: the phase fails; its displacement, stress and reaction fields are there and sized to the
// mesh; the reactions balance the load that WAS equilibrated (lambda q B, to round-off -- they are
// the discrete B^T sigma of the same committed states); and that load is the Prandtl capacity to
// within the usual limit-analysis band on this mesh.
// verify: KV-FND-016
//   oracle:   closed_form
//   source:   global equilibrium of the kept state (support reactions = the equilibrated share of the footing load) and Prandtl (1921) for a weightless Tresca (phi = 0) half-space
//   locator:  sum of the vertical support reactions = load_factor * q * B; q_ult = (2 + pi) c, N_c = 2 + pi = 5.14159 (stated in full above)
//   quantity: vertical reaction balance of the stopped phase's kept state [-] and its N_c = load_factor * q / c [-]
//   expected: reaction balance 0 (to round-off); N_c 5.142
//   band:     reaction balance 1e-6 relative and N_c 8%, as asserted below -- measured 3.5e-16 and +0.1% on the 0.1 m^2 tri6 mesh
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
}  // namespace

int main() {
    constexpr double W = 6.0, H = 4.0, x0 = 2.4, x1 = 3.6, B = x1 - x0, c = 10.0;
    const double q_applied = 1.2 * (2.0 + kPi) * c;
    m::Project pr;
    m::Material s;
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e4; s.nu = 0.3; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    s.c = c; s.phi = 0.0; s.psi = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = x0; L.y1 = H; L.x2 = x1; L.y2 = H;
    L.qx1 = 0; L.qy1 = -q_applied; L.qx2 = 0; L.qy2 = -q_applied;
    pr.loads.push_back(L);

    const auto M = katai::app::mesh_from_project(pr, 0.1, 6);
    check(M.ok, "the footing meshes");
    if (!M.ok) return 1;
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::GravityLoading);
    const int n = R.mesh.node_count;
    std::printf("  ok=%d load factor %.4f, nodes %d, |disp| %zu, |reaction| %zu, max|u| %.4g m\n",
                (int)R.ok, R.load_factor, n, (size_t)R.disp.size(), (size_t)R.reaction.size(),
                R.max_disp);
    check(!R.ok && R.load_factor > 0.0 && R.load_factor < 1.0,
          "a footing loaded past its capacity stops short and says so");
    check(n > 0 && R.disp.size() == 2 * n && R.reaction.size() == 2 * n &&
              (int)R.stress.stress.size() == n,
          "the stopped phase keeps its last equilibrated displacement, stress and reactions");
    check(R.max_disp > 0.0, "and that displacement is the mechanism, not zeros");
    double Ry = 0.0;
    for (int i = 0; i < n; ++i) Ry += R.reaction[2 * i + 1];
    const double equilibrated = R.load_factor * q_applied * B;
    const double rel = std::fabs(Ry - equilibrated) / equilibrated;
    const double Nc = R.load_factor * q_applied / c;
    std::printf("  support reactions %.4f kN/m vs equilibrated load %.4f (rel %.2e); N_c = %.3f "
                "(2 + pi = %.3f)\n", Ry, equilibrated, rel, Nc, 2.0 + kPi);
    check(rel < 1e-6, "the reactions balance the load that was equilibrated");
    check(std::fabs(Nc - (2.0 + kPi)) / (2.0 + kPi) < 0.08,
          "and that load is the Prandtl capacity (within the limit-analysis band of this mesh)");
    if (g_failures == 0) {
        std::printf("OK: a stopped phase keeps its last equilibrated state\n");
        return 0;
    }
    std::printf("%d check(s) failed\n", g_failures);
    return 1;
}
