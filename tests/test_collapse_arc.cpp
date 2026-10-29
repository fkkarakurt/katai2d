// The collapse load of a load-controlled phase, to its own resolution and not to the load step.
//
// A static phase that is given more load than the ground can carry used to stop at the last
// increment that closed: its collapse load came out as a whole number of the smallest increment
// (1/8 of the first), so a comparison with a closed form was as fine as the stepping and no finer
// (KV-CST-024 had to say so: 51.000 = 0.85 x 60, 36.000 = 0.60 x 60). From that increment the
// phase now follows the path by arc length with the load factor an unknown (NewtonOptions::
// limit_by_arc) and reports where the path stops rising (K2D-A040).
//
//   (a) a homogeneous Mohr-Coulomb block, unconfined, its top loaded with three times what it
//       carries: every point fails at once, at the closed-form unconfined strength
//       q_u = 2 c cos(phi) / (1 - sin(phi)); the load factor at collapse is q_u / q.
//   (b) the Prandtl strip footing of the corpus (KV-FND-010) on Tresca clay: (2 + pi) su.
// Each is also run the old way (KATAI_LIMIT_BY_STEPS, the measurement seam), to show what the
// stepping alone reported.
//
// verify: KV-NUM-016
//   oracle:   closed_form
//   source:   the unconfined compressive strength of a Mohr-Coulomb material, q_u = 2 c cos(phi) / (1 - sin(phi)) (Coulomb's criterion with sigma_3 = 0); Prandtl, L. (1920). Ueber die Haerte plastischer Koerper. Nachr. Ges. Wiss. Goettingen, 74-85 -- q_u = (2 + pi) su
//   locator:  (a) a 1 x 1 m plane strain block, tri6 at 0.25 m, smooth base and left side, c = 10 kPa, phi = 30 deg, psi = 0, top loaded with 3 q_u in 10 steps; (b) the corpus file kv-fnd-010-prandtl-strip-footing.k2d as it stands (B = 1.2 m, su = 10 kPa, q = 60 kPa)
//   quantity: the load factor at collapse times the applied load [kPa]
//   expected: (a) q_u = 34.641 kPa within 2e-3 relative, where the stepping alone stops a whole smallest increment (1/80 of the load, 3e-2 relative) below or at it; (b) (2 + pi) su = 51.416 kPa within 1.5% (the mesh), and no lower than what the stepping alone reported
//   band:     (a) 2e-3: the arc's resolution, 1e-3 of the factor, doubled; (b) 1.5%: the mesh's own error on this footing (KV-FND-010)
#include <katai/io/project_io.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }
constexpr double PI = 3.14159265358979323846;

void set_steps_only(bool on) {
#ifdef _WIN32
    _putenv_s("KATAI_LIMIT_BY_STEPS", on ? "1" : "");
#else
    if (on) setenv("KATAI_LIMIT_BY_STEPS", "1", 1); else unsetenv("KATAI_LIMIT_BY_STEPS");
#endif
}

struct Run { bool solved = false; bool collapsed = false; double lf = 0.0; bool said = false; };

Run run(const m::Project& pr, int last) {
    Run r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return r;
    const auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    if ((int)res.size() <= last) return r;
    r.solved = true;
    r.collapsed = !res[(size_t)last].ok;
    r.lf = res[(size_t)last].load_factor;
    for (const auto& d : res[(size_t)last].diagnostics) r.said = r.said || d.code == "K2D-A040";
    return r;
}

}  // namespace

int main() {
    std::printf("== The collapse load by arc length (KV-NUM-016) ==\n");
    {
        // (a) the unconfined block
        const double c = 10.0, phi = 30.0 * PI / 180.0;
        const double qu = 2.0 * c * std::cos(phi) / (1.0 - std::sin(phi));
        const double q = 3.0 * qu;
        m::Project pr;
        pr.x_min = 0; pr.x_max = 1; pr.y_min = 0; pr.y_max = 1;
        pr.has_water = false;
        pr.initial_procedure = m::InitialProcedure::GravityLoading;
        pr.mesh.elem_size = 0.25;
        pr.mesh.auto_refine = false;
        m::Material s;
        s.name = "MC"; s.model = m::SoilModel::MohrCoulomb; s.drainage = m::Drainage::Drained;
        s.gamma_unsat = s.gamma_sat = 0.0; s.E = 1.0e4; s.nu = 0.3; s.c = c; s.phi = 30.0; s.psi = 0.0;
        s.tension_cutoff = false;
        s.k0_auto = false; s.k0 = 1.0;
        pr.materials.push_back(s);
        m::SoilPolygon g;
        g.name = "Block"; g.material = 0; g.x = {0, 1, 1, 0}; g.y = {0, 0, 1, 1};
        g.edge_bc = {(int)m::BCType::VerticallyFixed, 0, 0, (int)m::BCType::HorizontallyFixed};
        pr.polygons.push_back(g);
        m::Load top; top.kind = m::LoadKind::Distributed; top.name = "q";
        top.x1 = 0; top.y1 = 1; top.x2 = 1; top.y2 = 1; top.qx1 = top.qx2 = 0.0; top.qy1 = top.qy2 = -q;
        pr.loads = {top};
        pr.initial.load_active = {0};
        m::Phase ph = pr.initial;
        ph.name = "Load"; ph.type = m::PhaseType::Plastic; ph.load_active = {1}; ph.load_steps = 10;
        pr.phases.push_back(ph);
        set_steps_only(true);
        const Run st = run(pr, 1);
        set_steps_only(false);
        const Run ar = run(pr, 1);
        std::printf("      unconfined block: by steps %.6f kPa, by arc %.6f kPa, closed form %.6f kPa\n",
                    st.lf * q, ar.lf * q, qu);
        check(ar.solved && ar.collapsed && ar.said, "(a) the overloaded block collapses, and says how its collapse load was found (K2D-A040)");
        check(rel(ar.lf * q, qu) < 2e-3, "(a) its collapse load is the unconfined strength within 2e-3");
        check(st.solved && st.lf * q <= qu * (1.0 + 1e-9) && rel(st.lf * q, qu) > rel(ar.lf * q, qu),
              "(a) ...where the steps alone stopped at a whole increment below it");
    }
    {
        // (b) the Prandtl footing
        m::Project pr;
        std::string err;
        const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-fnd-010-prandtl-strip-footing.k2d";
        const bool loaded = m::load_project(path, pr, &err, nullptr);
        check(loaded, "the Prandtl case loads");
        if (loaded) {
            const double q = 60.0, su = pr.materials[0].c, want = (2.0 + PI) * su;
            set_steps_only(true);
            const Run st = run(pr, 1);
            set_steps_only(false);
            const Run ar = run(pr, 1);
            std::printf("      strip footing: by steps %.4f kPa, by arc %.4f kPa, Prandtl %.4f kPa (%+.2f%%)\n",
                        st.lf * q, ar.lf * q, want, 100.0 * (ar.lf * q - want) / want);
            check(ar.solved && ar.said && rel(ar.lf * q, want) < 0.015, "(b) the footing's collapse load is Prandtl's within 1.5%");
            check(st.solved && ar.lf >= st.lf, "(b) ...and no lower than the load the steps alone reached");
        }
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
