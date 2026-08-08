// A prescribed boundary flux: rainfall on a slope, a drainage blanket, a recharge boundary --
// and, until now, an input the engine could compute and the file could not ask for.
//
// The element integral existed (`assemble_seepage_flux`); what was missing was every path from
// a .k2d to it, and a right-hand side in the solvers to receive it. That is the same defect
// class as the anchor lock-off force: a capability the input contract cannot express is a
// capability the user does not have.
//
// WHAT THE MANUAL SAYS, and where the scope comes from. PLAXIS 2D 2025.1 Scientific Manual,
// groundwater flow: the discretised system carries "the prescribed recharges that are given by
// the boundary conditions" (Eq. 3-31), whose boundary term is the surface integral of the
// prescribed flux (Eq. 3-34). The sign here follows the wells of the same chapter -- "the source
// term is positive for a recharge well" (§3.2.7) -- so a positive flux is water entering the
// soil. And the limit is the manual's too: in the CONSOLIDATION formulation (Ch. 4) "it is not
// possible to have boundaries with non-zero prescribed outflow", so this is wired into the flow
// problem only, exactly as PLAXIS does it.
//
// verify: KV-FLW-002
//   oracle:   closed_form
//   source:   Darcy's law for one-dimensional steady flow through a homogeneous column; the boundary term itself is PLAXIS 2D 2025.1 Scientific Manual Eqs. 3-31 and 3-34 (prescribed recharge on the boundary), with the sign convention of §3.2.7 (positive = into the soil) and the scope limit of Ch. 4 (consolidation carries no non-zero prescribed outflow)
//   locator:  a column of height H and width B, head h = h0 prescribed at the base, a uniform inflow q [m/day] prescribed on the top edge, everything else closed: continuity gives the same specific discharge at every level, so Darcy q = k dh/dy integrates to h(y) = h0 + (q/k) y -- a straight line of slope q/k -- and the total inflow is q*B (stated in full)
//   quantity: the nodal head profile of the column and the total discharge recovered from the solution [m; m3/day per metre of wall]
//   expected: h(y) = h0 + (q/k) y at every node, and the recovered inflow equal to q*B
//   band:     1e-9 m on the head profile and 1e-9 on the discharge, as asserted below -- the exact solution is linear in y, so tri6 elements represent it exactly and the residual is the direct solver's

#include <katai/io/validate.hpp>
#include <katai/jobs/flow_driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kW = 4.0, kH = 10.0;    // column [m]
constexpr double kK = 2.0;               // permeability [m/day], isotropic
constexpr double kQ = 0.05;              // prescribed inflow on the top edge [m/day]
constexpr double kH0 = 12.0;             // head prescribed at the base [m]

m::Project column() {
    m::Project pr;
    pr.name = "KV-FLW-002 prescribed flux";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.has_water = true;
    pr.wx = {0.0, kW};
    pr.wy = {kH, kH};

    m::Material s;
    s.name = "Sand";
    s.model = m::SoilModel::LinearElastic;
    s.E = 2.0e4; s.nu = 0.3;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.kx = kK; s.ky = kK;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    // Edges in order: base (head), right (closed), top (prescribed inflow), left (closed).
    P.edge_flow = {(int)m::FlowBCType::Head, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Flux, (int)m::FlowBCType::Closed};
    P.edge_head = {kH0, 0.0, 0.0, 0.0};
    P.edge_flux = {0.0, 0.0, kQ, 0.0};
    pr.polygons.push_back(P);
    return pr;
}

}  // namespace

int main() {
    std::printf("== prescribed boundary flux (KV-FLW-002) ==\n");

    const m::Project pr = column();
    const katai::io::ValidationReport rep = katai::io::validate_project(pr);
    for (const auto& i : rep.issues)
        std::printf("      %s %s: %s\n",
                    i.severity == katai::io::Severity::Error ? "[error]  " : "[warning]",
                    i.path.c_str(), i.message.c_str());
    check(rep.ok(), "a project with a prescribed-flux edge validates");

    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return 1; }

    const auto R = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(R.ok, "the flow solve converges");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return 1; }

    // Darcy: one-dimensional steady flow through a homogeneous column gives a straight head
    // profile whose slope is the specific discharge over the permeability.
    const double slope = kQ / kK;
    double worst = 0.0;
    int worst_node = 0;
    for (int n = 0; n < M.mesh.node_count; ++n) {
        const double want = kH0 + slope * M.mesh.y[n];
        const double err = std::fabs(R.head[n] - want);
        if (err > worst) { worst = err; worst_node = n; }
    }
    std::printf("  q = %.3f m/day, k = %.3f m/day -> dh/dy = %.6f; head at the top should be "
                "%.6f m\n", kQ, kK, slope, kH0 + slope * kH);
    std::printf("  worst nodal deviation from the closed form: %.3e m (node %d at y = %.3f)\n",
                worst, worst_node, M.mesh.y[worst_node]);
    check(worst < 1e-9, "every node lies on the Darcy profile h = h0 + (q/k) y");

    // The solve reports the discharge through its PRESCRIBED-HEAD boundaries. In steady state
    // whatever enters through the flux edge has to leave through the base, so that number is an
    // independent statement of the same quantity -- and it is the one a user reads.
    std::printf("  discharge through the prescribed-head base: %.9f m3/day/m (inflow asked %.9f)\n",
                std::fabs(R.discharge), kQ * kW);
    check(std::fabs(std::fabs(R.discharge) - kQ * kW) < 1e-9,
          "what enters through the flux edge leaves through the base: q x width");
    std::printf("  mass-balance error reported by the solve: %.3e\n", R.balance_err);
    check(R.balance_err < 1e-9, "and the global mass balance closes");

    // A zero flux must be the closed boundary it already was: the natural condition, not a new
    // one. Same model, flux set to zero, and the head must be uniform at the prescribed value.
    m::Project dry = column();
    dry.polygons[0].edge_flux = {0.0, 0.0, 0.0, 0.0};
    const auto Rd = katai::app::solve_groundwater_flow(dry, M.mesh);
    check(Rd.ok, "the zero-flux model solves");
    if (Rd.ok) {
        double spread = 0.0;
        for (int n = 0; n < M.mesh.node_count; ++n)
            spread = std::fmax(spread, std::fabs(Rd.head[n] - kH0));
        std::printf("  zero flux: head spread across the column = %.3e m\n", spread);
        check(spread < 1e-9,
              "a zero prescribed flux is exactly the closed boundary it always was");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
