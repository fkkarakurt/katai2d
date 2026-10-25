// Re-tensioning an anchor in a later phase (`phases[i].sjack`): the jack is put on an anchor that
// is already in place, it holds the new force while the ground moves, and the anchor is locked
// again where the ground comes to. Re-stressing a row of ground anchors after an excavation stage
// is ordinary practice, and until this a model could only install an anchor at its lock-off force
// and leave it there.
//
// The oracle is a path identity on a linear-elastic block, exact rather than close. The jack holds
// a force and adds no stiffness, and the anchor is a spring only from the lock on, so
//
//     install at N1, re-tension to N2   ==   install at N2
//
// in every phase from the re-tension on: the same total displacement field, and the same anchor
// force in a later phase that loads the block -- the ground responds linearly to the forces the
// anchor puts on it, and those are N1 then N2 - N1 on one path and N2 on the other. The two paths
// are built from different inputs (a material prestress against a phase's re-tension). In the
// re-tension phase the jack's force is in the baseline the solve starts from, so the phase begins
// in equilibrium and raises the anchor from N1 to N2 over its steps; an anchor left without the
// jack would stay a spring at N1 and the identity would fail by the whole N2 - N1.
//
// verify: KV-STR-022
//   oracle:   independent_path
//   source:   KATAI 2D input contract (docs/k2d-format.md, anchors[i].prestress and phases[i].sjack): an anchor installed with a lock-off force holds it through its installation phase and is locked at the end; an anchor re-tensioned in a later Plastic phase holds the re-tension force through that phase and is locked again at the end, after which it is an elastic spring from that lock. The identity used as the oracle follows from linear elasticity of the ground and from the jack adding no stiffness: the ground's response depends only on the forces the anchor applies, N1 then N2 - N1 on one path and N2 on the other
//   locator:  a weightless linear-elastic block 20 x 10 m (E 20 MPa, nu 0.3), a fixed-end ground anchor (EA 2e5 kN) from a node on its surface to a point outside the mesh; path A installs it at 300 kN and re-tensions it to 500 kN in the next phase, path B installs it at 500 kN and leaves the next phase empty; both then load the surface with 50 kPa over 4 m. A third run re-tensions an anchor with a tension capacity of 450 kN to 500 kN
//   quantity: the anchor force at the end of the re-tension phase [kN]; the total displacement field after the re-tension phase and after the loading phase, path A against path B [m]; the anchor force in the loading phase, A against B [kN]; the force a re-tension above the capacity leaves [kN]
//   expected: the re-tension phase ends with the anchor at exactly the jack force; the two paths give the same fields and the same anchor force in the loading phase; a jack above the capacity leaves the anchor at its capacity
//   band:     1e-9 relative on the anchor forces and 1e-9 m on the fields -- the two paths assemble the same linear systems, so agreement is exact to solver round-off
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
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

constexpr double kW = 20.0, kH = 10.0;
constexpr double kEA = 2.0e5;
constexpr double kAx = 6.0, kAy = kH;            // soil end, on the surface
constexpr double kBx = -4.0, kBy = kH + 5.0;     // fixed end, outside the mesh

// install = the anchor material's lock-off force; jack = the re-tension of phase 2 (0 = none).
m::Project model(double install, double jack, double capacity = 0.0) {
    m::Project pr;
    pr.name = "anchor re-tension";
    pr.x_min = 0; pr.x_max = kW; pr.y_min = 0; pr.y_max = kH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 2.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Elastic soil"; s.model = m::SoilModel::LinearElastic;
    s.E = 2.0e4; s.nu = 0.3; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Block"; P.material = 0;
    P.x = {0, kW, kW, 0}; P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::AnchorMaterial am; am.name = "Ground anchor"; am.EA = kEA; am.Lspacing = 1.0;
    am.prestress = install;
    if (capacity > 0.0) { am.elastoplastic = true; am.Fmax_tens = capacity; }
    pr.anchors.push_back(am);
    m::StructElement a; a.kind = m::StructKind::Anchor; a.name = "Row 1"; a.material = 0;
    a.x1 = kAx; a.y1 = kAy; a.x2 = kBx; a.y2 = kBy;
    pr.structs.push_back(a);

    m::Load q; q.kind = m::LoadKind::Distributed; q.name = "Surcharge";
    q.x1 = 12.0; q.y1 = kH; q.x2 = 16.0; q.y2 = kH; q.qy1 = q.qy2 = -50.0;
    pr.loads.push_back(q);

    pr.initial.struct_active = {0};
    pr.initial.load_active = {0};
    m::Phase p1; p1.name = "Install"; p1.type = m::PhaseType::Plastic;
    p1.struct_active = {1}; p1.load_active = {0};
    m::Phase p2 = p1; p2.name = "Re-tension";
    if (jack != 0.0) p2.struct_retension = {jack};
    m::Phase p3 = p1; p3.name = "Load"; p3.load_active = {1};
    pr.phases = {p1, p2, p3};
    return pr;
}

struct Phase { std::vector<double> total; double N = 0.0; };
struct Run { bool ok = false; std::vector<Phase> ph; };

Run solve(const m::Project& pr) {
    Run r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { std::printf("      (mesh: %s)\n", M.message.c_str()); return r; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != pr.phases.size() + 1) return r;
    for (size_t k = 1; k < res.size(); ++k) {
        if (!res[k].ok) { std::printf("      (phase %zu: %s)\n", k, res[k].message.c_str()); return r; }
        Phase p;
        p.total.assign(res[k].disp_total.data(), res[k].disp_total.data() + res[k].disp_total.size());
        for (const auto& sf : res[k].struct_forces)
            if (sf.kind == 1) p.N = sf.max_N;
        r.ph.push_back(p);
    }
    r.ok = true;
    return r;
}

double field_difference(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size() || a.empty()) return 1e300;
    double d = 0.0;
    for (size_t i = 0; i < a.size(); ++i) d = std::max(d, std::fabs(a[i] - b[i]));
    return d;
}
double largest(const std::vector<double>& a) {
    double d = 0.0;
    for (double v : a) d = std::max(d, std::fabs(v));
    return d;
}

}  // namespace

int main() {
    std::printf("== anchor re-tension (KV-STR-022) ==\n");
    const double N1 = 300.0, N2 = 500.0;
    const Run A = solve(model(N1, N2));   // install at N1, re-tension to N2
    const Run B = solve(model(N2, 0.0));  // install at N2
    check(A.ok && B.ok, "both paths solve");
    if (!A.ok || !B.ok) { std::fprintf(stderr, "\n%d check(s) failed\n", g_failures + 1); return 1; }

    std::printf("  path A: N = %.9f (installed) -> %.9f (re-tensioned) -> %.9f (loaded) kN\n",
                A.ph[0].N, A.ph[1].N, A.ph[2].N);
    std::printf("  path B: N = %.9f (installed) -> %.9f (empty phase)  -> %.9f (loaded) kN\n",
                B.ph[0].N, B.ph[1].N, B.ph[2].N);
    check(std::fabs(A.ph[0].N - N1) <= 1e-9 * N2, "path A installs the anchor at its lock-off force");
    check(std::fabs(A.ph[1].N - N2) <= 1e-9 * N2, "the re-tension phase ends with the anchor at the jack force");

    const double moved = field_difference(A.ph[1].total, A.ph[0].total);
    const double d2 = field_difference(A.ph[1].total, B.ph[1].total);
    const double d3 = field_difference(A.ph[2].total, B.ph[2].total);
    std::printf("  the re-tension moved the block by up to %.6e m (|u| after it: %.6e m)\n", moved,
                largest(A.ph[1].total));
    std::printf("  fields, A against B: after the re-tension %.3e m, after the loading %.3e m\n", d2, d3);
    check(moved > 1e-4 * largest(A.ph[1].total), "the re-tension changes the ground (its difference is applied)");
    check(d2 < 1e-9, "after the re-tension: the same field as installing at the jack force");
    check(d3 < 1e-9, "after the loading: the same field");
    check(std::fabs(A.ph[2].N - B.ph[2].N) <= 1e-9 * N2,
          "in the loading phase the anchor is a spring from the re-tensioned lock: the same force");
    check(std::fabs(A.ph[2].N - N2) > 1e-6 * N2, "and the loading does change it (the check is not vacuous)");

    // A jack above what the anchor can carry leaves it at its capacity.
    const double cap = 450.0;
    const Run C = solve(model(N1, N2, cap));
    check(C.ok, "the run with a tension capacity solves");
    if (C.ok) {
        std::printf("  capacity %.0f kN, jack %.0f kN: N = %.9f kN\n", cap, N2, C.ph[1].N);
        check(std::fabs(C.ph[1].N - cap) <= 1e-9 * cap, "a jack above the capacity leaves the anchor at its capacity");
    }
    if (g_failures) std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return g_failures ? 1 : 0;
}
