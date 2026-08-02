// K0 procedure correctness audit (diagnostic, EXCLUDE_FROM_ALL). Runs the EXACT GUI compute path
// (build_problem.hpp solve_gravity_le) and prints concrete numbers so we can be 1000% sure the
// initial-stress procedure is right -- the user reported "K0 doesn't seem to work in the GUI".
//
// Geotechnical truth being checked (Terzaghi/Lambe&Whitman; PLAXIS Reference Manual, K0 procedure):
//   * K0 procedure: undisturbed ground is in GEOSTATIC EQUILIBRIUM. Displacement ~0 by design,
//     and sigma'_v = -gamma'(H-y), sigma'_h = K0 sigma'_v.
//   * Gravity loading (self-weight from a STRESS-FREE start) reaches the SAME final sigma'_v at
//     equilibrium (vertical equilibrium is path-independent for confined 1D), but WITH settlement.
//     => the two phases differ ONLY in displacement, not in the equilibrium vertical stress.
//   * K0 + surface load q: an added Delta_sigma'_v ~ q under the load, plus a settlement increment.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace model = katai::model;
using katai::app::InitialPhase;

namespace {

// Nearest-node lookup of a stress component (0=sigma'_xx, 1=sigma'_yy, 2=sigma_xy) at (x,y).
double stress_at(const katai::app::SolveResult& R, const katai::mesh::Mesh& mesh,
                 double x, double y, int comp) {
    int best = 0; double bd = 1e30;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double d = std::hypot(mesh.x[n] - x, mesh.y[n] - y);
        if (d < bd) { bd = d; best = n; }
    }
    return R.stress.stress[best](comp);
}

model::Project make_block(double gamma_unsat, double gamma_sat, double E, double nu,
                          double phi_deg, double W, double H) {
    model::Project pr;
    model::Material m;
    m.model = model::SoilModel::LinearElastic;
    m.gamma_unsat = gamma_unsat; m.gamma_sat = gamma_sat;
    m.E = E; m.nu = nu; m.phi = phi_deg; m.c = 1.0;
    m.k0_auto = true;                       // K0 = 1 - sin(phi)
    pr.materials.push_back(m);
    model::SoilPolygon P; P.material = 0;
    P.x = {0.0, W, W, 0.0};
    P.y = {0.0, 0.0, H, H};
    P.edge_bc = {(int)model::BCType::FullyFixed, (int)model::BCType::HorizontallyFixed,
                 (int)model::BCType::Free,       (int)model::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

void audit_k0_vs_gravity() {
    const double W = 20.0, H = 10.0, gamma = 18.0, phi = 30.0;
    const double K0 = 1.0 - std::sin(phi * 3.14159265358979 / 180.0);  // Jaky ~0.5
    model::Project pr = make_block(gamma, 20.0, 2.0e4, 0.3, phi, W, H);
    const auto mr = katai::app::mesh_from_project(pr, 0.5, 6);
    std::printf("=== K0 procedure vs Gravity loading (dry block, gamma=%.0f, K0=%.3f) ===\n", gamma, K0);
    std::printf("    Truth: sigma'_v = -gamma(H-y), sigma'_h = K0 sigma'_v; BOTH phases reach it.\n\n");

    const auto Rk = katai::app::solve_gravity_le(pr, mr.mesh, InitialPhase::K0Procedure);
    const auto Rg = katai::app::solve_gravity_le(pr, mr.mesh, InitialPhase::GravityLoading);
    std::printf("    phase        max|u| [m]\n");
    std::printf("    K0           %.4e   <- geostatic equilibrium (must be ~0)\n", Rk.max_disp);
    std::printf("    Gravity      %.4e   <- self-weight settlement (must be > 0)\n\n", Rg.max_disp);

    std::printf("    depth   y      sigma'_v(K0)   sigma'_v(Grav)   exact      sigma'_h(K0)  K0*sv\n");
    for (double y : {10.0, 7.5, 5.0, 2.5, 0.0}) {
        const double svk = stress_at(Rk, mr.mesh, W / 2, y, 1);
        const double svg = stress_at(Rg, mr.mesh, W / 2, y, 1);
        const double shk = stress_at(Rk, mr.mesh, W / 2, y, 0);
        const double sv_exact = -gamma * (H - y);
        std::printf("    %5.1f  %4.1f   %11.3f   %11.3f   %9.3f   %11.3f  %8.3f\n",
                    H - y, y, svk, svg, sv_exact, shk, K0 * sv_exact);
    }
    std::printf("\n");
}

void audit_k0_with_load() {
    const double W = 20.0, H = 10.0, gamma = 18.0, q = 50.0;  // 50 kPa strip surcharge at the top
    model::Project pr = make_block(gamma, 20.0, 2.0e4, 0.3, 30.0, W, H);
    model::Load L; L.kind = model::LoadKind::Point; L.x1 = W / 2; L.y1 = H; L.qx1 = 0; L.qy1 = -q;
    pr.loads.push_back(L);
    const auto mr = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto R = katai::app::solve_gravity_le(pr, mr.mesh, InitialPhase::K0Procedure);
    std::printf("=== K0 + point surcharge %.0f kN at crest (K0 phase, only the load drives motion) ===\n", q);
    std::printf("    max|u| = %.4e m (load-driven; K0 self-weight alone gives 0).\n", R.max_disp);
    std::printf("    sigma'_v just under the load (x=W/2): %.2f kPa vs geostatic %.2f kPa (load adds compression).\n\n",
                stress_at(R, mr.mesh, W / 2, H - 0.5, 1), -gamma * 0.5);
}

}  // namespace

int main() {
    audit_k0_vs_gravity();
    audit_k0_with_load();
    std::printf("If K0 max|u|~0, both phases match sigma'_v, and the load adds compression -> the\n");
    std::printf("K0 procedure COMPUTES correctly; any 'K0 not working' is a post-processing/display issue.\n");
    return 0;
}
