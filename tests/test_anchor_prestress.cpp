// Anchor prestress: the input that made every anchored excavation unmodellable until now, and
// the first STRUCTURAL case in the verification matrix.
//
// The gap was found the way the plan says to find gaps -- by trying to build a problem the
// industry actually uses as a reference. The DGGT / Schweiger triple-anchored excavation in
// Berlin sand (Schweiger 2002) locks its three anchor
// rows off at 768, 945 and 980 kN. KATAI had no field for that force at all: an anchor could be
// stiff, it could yield, it could not be TENSIONED. A model of an anchored wall built from those
// primitives is a model of a wall with slack anchors, which deflects far more than the real one
// -- the unsafe direction.
//
// The oracle here is an equivalence rather than a formula, and it is exact, in two parts.
//
// IN THE PHASE THAT INSTALLS IT the anchor is the jack: it holds the lock-off force whatever the
// ground does, and adds no stiffness, so
//
//     prestressed anchor  ==  no anchor at all  +  an external force of N0 along its axis,
//
// applied at the node the anchor attaches to, and the anchor reports exactly N0 at the end of the
// phase -- on soft ground as on rigid ground. (Until 2026-09-30 the anchor was a spring from the
// start of that phase, the equivalence was with a SLACK anchor plus the force, and the ground's
// movement towards the anchor unloaded it within its own phase: 15% of the lock-off force was left
// on a tie-back excavation rebuilt from a published tutorial, and only a rigid soil kept it.)
//
// IN THE PHASES AFTER IT the anchor is locked where the ground came to and is an elastic spring
// from the lock-off force: on a linear-elastic block the increment of a later phase is the
// response of the block WITH a slack anchor to that phase's load alone, and the anchor force is
// N0 plus what that slack anchor carries. Both models are built from different primitives -- a
// point load, a slack anchor -- so agreement between them is a statement about the implementation.
//
// verify: KV-STR-001
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, anchors[i].prestress): in the phase that activates it a node-to-node or fixed-end anchor holds its lock-off force and is locked at the end of that phase, after which it behaves as an elastic spring from that state; the equivalences used as the oracle are that a constant force N0 along the anchor axis is, in the residual, an external force of the same magnitude and direction, and that on a linear-elastic model a later phase's increment is superposable -- so the same problem can be built twice from different primitives
//   locator:  installation phase: N = N0 with no stiffness, U_p committed so that N0 + (EA/L)(U - U_p) = N0 at the phase end; later phases: N = N0 + (EA/L)(U - U_p) with U the elongation measured from the installation datum; for a fixed-end anchor the only mesh node carries f_int += N*(-dir), so the N0 term is exactly an external nodal force +N0*dir (stated in full)
//   quantity: the full nodal displacement field of an elastic block anchored by a prestressed fixed-end anchor, against the same block with no anchor and an equivalent point load, and the reported anchor force against the lock-off force; then the increment of a second phase that loads the block, against the same load on the block with a slack anchor, and the anchor force against N0 plus that slack anchor's force [m; kN]
//   expected: the installation-phase fields identical to solver round-off and the anchor at exactly the lock-off force on soft ground; the later phase's increment identical to the slack-anchor response and the anchor force N0 plus the slack anchor's; and the response linear in the lock-off force
//   band:     1e-9 m on both field equivalences and 1e-9 relative on both anchor forces, as asserted below -- the models assemble the same linear systems, so the agreement is exact rather than merely close; 1e-6 on the doubling ratio

#include <katai/analysis/structural_forces.hpp>
#include <katai/jobs/driver.hpp>
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

constexpr double kW = 20.0, kH = 10.0;      // elastic block [m]
constexpr double kN0 = 500.0;               // lock-off force [kN]
constexpr double kEA = 2.0e5;               // anchor axial stiffness [kN]
// The anchor: from a node on the block's surface out to a fixed point beyond the model, so it is
// a fixed-end anchor and exactly one mesh node carries it.
constexpr double kAx = 6.0, kAy = kH;       // soil end
constexpr double kBx = -4.0, kBy = kH + 5.0;   // fixed far end (outside the mesh)

m::Project block(double E) {
    m::Project pr;
    pr.name = "anchor prestress";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Elastic soil";
    s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = 0.3;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;   // weightless: only the anchor acts
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Block";
    P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

void add_anchor(m::Project& pr, double prestress) {
    m::AnchorMaterial am;
    am.name = "Ground anchor";
    am.EA = kEA;
    am.Lspacing = 1.0;
    am.prestress = prestress;
    pr.anchors.push_back(am);
    m::StructElement s;
    s.kind = m::StructKind::Anchor;
    s.name = "Row 1";
    s.x1 = kAx; s.y1 = kAy; s.x2 = kBx; s.y2 = kBy;
    s.material = (int)pr.anchors.size() - 1;
    pr.structs.push_back(s);
}

struct Run {
    bool ok = false;
    std::vector<double> disp;          // the LAST phase's displacement (its own increment)
    double max_disp = 0.0, anchor_N = 0.0;
};

Run solve(const m::Project& pr) {
    Run r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return r; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.empty() || !res.back().ok || res.size() != pr.phases.size() + 1) {
        std::printf("      (solve: %s)\n", res.empty() ? "no phases" : res.back().message.c_str());
        return r;
    }
    r.disp.assign(res.back().disp.data(), res.back().disp.data() + res.back().disp.size());
    r.max_disp = res.back().max_disp;
    for (const auto& sf : res.back().struct_forces)
        if (sf.kind == 1) r.anchor_N = sf.max_N;   // anchors report their axial force
    r.ok = true;
    return r;
}

double field_difference(const Run& a, const Run& b) {
    if (a.disp.size() != b.disp.size()) return 1e300;
    double d = 0.0;
    size_t worst = 0;
    for (size_t i = 0; i < a.disp.size(); ++i) {
        const double e = std::fabs(a.disp[i] - b.disp[i]);
        if (e > d) { d = e; worst = i; }
    }
    std::printf("      worst at dof %zu (node %zu, %s): %.9e vs %.9e\n", worst, worst / 2,
                worst % 2 ? "uy" : "ux", a.disp[worst], b.disp[worst]);
    return d;
}

}  // namespace

int main() {
    std::printf("== anchor prestress (KV-STR-001) ==\n");

    // The equivalent external force: N0 along the anchor axis, from the soil end towards the
    // fixed end -- the direction the anchor pulls the node when it is in tension.
    //
    // The axis is taken from the NODE the anchor attaches to, not from the point it was drawn
    // at, because that is the geometry the element itself uses: an anchor end is carried by the
    // nearest node (the mesher does not insert it as a vertex -- test_anchor_mesh_repro exists
    // to keep it that way). Using the drawn point instead leaves the two models pulling in
    // slightly different directions, which is measurable: it was 0.06% here before this was
    // fixed, and finding that difference is what the equivalence is for.
    const m::Project probe = block(2.0e4);
    const auto PM = katai::app::mesh_from_project(probe);
    check(PM.ok, "the mesh builds");
    if (!PM.ok) return 1;
    int anchor_node = 0;
    double bd = 1e300;
    for (int n = 0; n < PM.mesh.node_count; ++n) {
        const double d = std::hypot(PM.mesh.x[n] - kAx, PM.mesh.y[n] - kAy);
        if (d < bd) { bd = d; anchor_node = n; }
    }
    const double ax = PM.mesh.x[anchor_node], ay = PM.mesh.y[anchor_node];
    std::printf("  anchor drawn at (%.3f, %.3f), carried by node %d at (%.6f, %.6f), %.4f m away\n",
                kAx, kAy, anchor_node, ax, ay, bd);
    const double dx = kBx - ax, dy = kBy - ay, L = std::hypot(dx, dy);
    const double ux = dx / L, uy = dy / L;
    std::printf("  anchor axis (%.4f, %.4f), length %.4f m, lock-off %.1f kN\n", ux, uy, L, kN0);

    // --- the installation phase -------------------------------------------------------------
    m::Project a = block(2.0e4);
    add_anchor(a, kN0);
    m::Project b = block(2.0e4);   // no anchor: the equivalent force alone
    m::Load P;
    P.kind = m::LoadKind::Point;
    P.name = "Equivalent lock-off force";
    P.x1 = ax; P.y1 = ay;   // at the node the anchor is carried by, for the same reason
    P.qx1 = kN0 * ux; P.qy1 = kN0 * uy;
    b.loads.push_back(P);
    m::Project c = block(2.0e4);
    add_anchor(c, 0.0);   // slack, no equivalent force: the control

    const Run ra = solve(a), rb = solve(b), rc = solve(c);
    check(ra.ok && rb.ok && rc.ok, "all three models solve");
    if (!ra.ok || !rb.ok || !rc.ok) { std::printf("\n1 CHECK(S) FAILED\n"); return 1; }
    std::printf("  prestressed      max|u| = %.9e m,  anchor N = %.6f kN\n", ra.max_disp, ra.anchor_N);
    std::printf("  no anchor + load max|u| = %.9e m\n", rb.max_disp);
    std::printf("  slack alone      max|u| = %.9e m\n", rc.max_disp);

    const double diff = field_difference(ra, rb);
    std::printf("  largest nodal difference between the two models = %.3e m\n", diff);
    check(diff < 1e-9,
          "in its installation phase the prestressed anchor IS its lock-off force, with no "
          "stiffness of its own");
    check(rc.max_disp < 1e-12,
          "and a slack anchor on a weightless block does nothing at all -- so the difference "
          "above is the prestress, not the anchor");
    check(std::fabs(ra.anchor_N - kN0) <= 1e-9 * kN0,
          "the anchor ends its installation phase at exactly the lock-off force, on ground that "
          "moved");

    // --- a later phase: locked, and a spring from the lock-off force --------------------------
    // A vertical load on the block's surface next to the anchor, applied in a second phase.
    m::Load Q;
    Q.kind = m::LoadKind::Point;
    Q.name = "Later load";
    Q.x1 = 12.0; Q.y1 = kH;
    Q.qx1 = 0.0; Q.qy1 = -300.0;
    m::Project a2 = a;
    a2.loads.push_back(Q);
    a2.initial.load_active = {0};
    m::Phase later;
    later.name = "Later load";
    later.type = m::PhaseType::Plastic;
    later.load_active = {1};
    a2.phases = {later};
    m::Project c2 = c;             // the slack anchor carrying the later load alone
    c2.loads.push_back(Q);
    const Run ra2 = solve(a2), rc2 = solve(c2);
    check(ra2.ok && rc2.ok, "the two-phase and the slack-anchor models solve");
    if (ra2.ok && rc2.ok) {
        const double d2 = field_difference(ra2, rc2);
        std::printf("  later phase: increment vs slack anchor under the same load: %.3e m; anchor "
                    "N = %.6f kN = N0 + %.6f (slack anchor %.6f)\n", d2, ra2.anchor_N,
                    ra2.anchor_N - kN0, rc2.anchor_N);
        check(d2 < 1e-9,
              "after its installation phase the anchor is an elastic spring: a later phase's "
              "increment is the slack-anchor response to that phase's load");
        check(std::fabs(ra2.anchor_N - (kN0 + rc2.anchor_N)) <= 1e-9 * kN0,
              "and its force is the lock-off force plus what the spring picks up");
        check(std::fabs(rc2.anchor_N) > 1e-3 * kN0,
              "the later load does reach the anchor, so the check above is not trivially met");
    }

    // A tensioned anchor is not a stiffer anchor: doubling the lock-off force doubles the response.
    m::Project dbl = block(2.0e4);
    add_anchor(dbl, 2.0 * kN0);
    const Run rd = solve(dbl);
    check(rd.ok, "the double-prestress model solves");
    if (rd.ok) {
        std::printf("  doubling the lock-off force: max|u| %.9e -> %.9e m (ratio %.6f)\n",
                    ra.max_disp, rd.max_disp, rd.max_disp / ra.max_disp);
        check(std::fabs(rd.max_disp / ra.max_disp - 2.0) < 1e-6,
              "the response is exactly linear in the lock-off force (it is a load, not a "
              "stiffness)");
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
