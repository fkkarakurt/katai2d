// The safety property this file exists to pin: KATAI may use an input differently from the
// way it was drawn, but NEVER in silence. Until this gate existed, every geometric object
// that failed to attach to the mesh was skipped by a bare `continue` in the driver's build
// loops: a surcharge drawn five centimetres above the ground surface, a retaining wall drawn
// above the soil, a point load placed outside the model. Each of those produced a converged
// run, exit code 0, "0 warnings" from `validate` -- and the answer to a model the engineer
// never drew. On a weightless benchmark the missing load shows up as a zero displacement; in
// any real model self-weight settlement masks it completely.
//
// That failure mode is worse than a crash, because nothing anywhere says it happened, so the
// fixtures below take one verified corpus case and perturb exactly one field per fixture --
// the same perturbations that were measured on the shipped v0.6.1 binary before the fix.
// Two things are asserted for each:
//
//   (1) the run does NOT succeed silently -- it either refuses (an object that would carry
//       nothing) or warns with a stable code (an object used differently than drawn); and
//   (2) where the run continues, an INDEPENDENT model that draws the object the way the mesh
//       actually took it must give the same answer -- so the warning's account of what was
//       solved is checkable, not a reassurance.
//
// The codes are the contract. Front ends and scripts match on them (`d.code == "K2D-G003"`),
// never on the prose, so a code is never reused or reworded once it ships.
//
// verify: KV-DIA-001
//   oracle:   independent_path
//   source:   KATAI 2D input-safety property (no input may be discarded in silence): the independent path is the same problem with the object drawn where the mesh actually takes it, plus global vertical equilibrium of the support reactions
//   locator:  each fixture changes ONE field of tests/corpus/kv-fnd-008-strip-load.k2d (a weightless elastic half-plane, 40 x 20 m, tri6, 4 m strip at q = 100 kPa) and solves it from that project; sum of the base reactions must equal the load actually applied, and a clipped object must reproduce the explicitly shortened object
//   quantity: the diagnostic severity and code raised by each perturbation [-]; the summed base reaction of the runs that continue [kN/m]; peak displacement of a clipped plate against the explicitly shortened plate [m]
//   expected: refusals K2D-G001 (point load off the mesh), K2D-G003 (line load off the mesh), K2D-G005 (plate off the mesh), K2D-G007 (geogrid off the mesh), K2D-G008 (anchor with no end in the soil); warnings K2D-G002 (point load snapped), K2D-G004 (line load clipped), K2D-G006 (structure clipped), K2D-G009 (wall/interface not on mesh edges -> bonded), K2D-M001 (tension cut-off on a model that ignores it), K2D-A001 (linear dynamic reports zero stress), K2D-A003 (a structure does not receive a prescribed displacement); note K2D-A004 (reactions exclude the structural end force); the unperturbed case raises NOTHING and carries 4 m x 100 kPa = 400 kN/m; the half-outside strip carries 2 m x 100 kPa = 200 kN/m
//   band:     exact on severity and code; 1e-9 relative on the equilibrium sums (the same discrete B^T sigma the supports see, so the residual is round-off, measured ~1e-13); 1e-12 m on the clipped-versus-shortened plate, which is a bit-level identity because the mesher clips structural lines in the arrangement, so both models are the SAME mesh and the same assembly

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/io/project_io.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
namespace core = katai::core;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

// One perturbed run, flattened: what the engine decided, and what it computed.
struct Run {
    bool meshed = false;
    bool ok = false;                          // every phase converged and none refused
    std::vector<core::Diagnostic> diags;      // in phase order
    double max_disp = 0.0;
    double sum_reaction_y = 0.0;              // base reactions [kN/m]: global vertical equilibrium
};

Run solve(const m::Project& pr) {
    Run r;
    const auto M = katai::app::mesh_from_project(pr);
    r.meshed = M.ok;
    if (!M.ok) {
        std::printf("      (mesh: %s)\n", M.message.c_str());
        return r;
    }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    r.ok = !res.empty();
    for (const auto& R : res) {
        r.ok = r.ok && R.ok;
        for (const auto& d : R.diagnostics) r.diags.push_back(d);
    }
    if (!res.empty()) {
        r.max_disp = res.back().max_disp;
        for (int n = 0; n * 2 + 1 < (int)res.back().reaction.size(); ++n)
            r.sum_reaction_y += res.back().reaction[2 * n + 1];
    }
    return r;
}

const char* severity_name(core::DiagnosticSeverity s) {
    return s == core::DiagnosticSeverity::Refusal ? "refusal"
         : s == core::DiagnosticSeverity::Warning ? "warning"
                                                  : "note";
}

// Every diagnostic is printed as it is checked: a fixture that starts raising a SECOND
// diagnostic is a change in behaviour, and the log is where it becomes visible.
void print_diags(const Run& r) {
    for (const auto& d : r.diags)
        std::printf("      [%s %s] %s%s%s\n", severity_name(d.severity), d.code.c_str(),
                    d.subject.c_str(), d.subject.empty() ? "" : ": ", d.message.c_str());
}

bool raised(const Run& r, const char* code, core::DiagnosticSeverity sev) {
    for (const auto& d : r.diags)
        if (d.code == code && d.severity == sev) return !d.message.empty();
    return false;
}

// A refusal is the whole point: the run must stop, name the object, and carry the code. A
// front end that only reads `ok` still cannot mistake it for a successful analysis.
void expect_refusal(const Run& r, const char* code, const std::string& what) {
    print_diags(r);
    check(r.meshed, what + ": the project still meshes (the refusal is about attachment)");
    check(!r.ok, what + ": the run does NOT report success");
    check(raised(r, code, core::DiagnosticSeverity::Refusal), what + ": refuses with " + code);
}

void expect_warning(const Run& r, const char* code, const std::string& what) {
    print_diags(r);
    check(r.ok, what + ": the run completes (the object is usable, just not as drawn)");
    check(raised(r, code, core::DiagnosticSeverity::Warning), what + ": warns with " + code);
}

// ------------------------------------------------------------------- the reference case --
// The corpus strip-load benchmark, loaded from its checked-in file: weightless elastic
// half-plane, so the only vertical force in the system is the strip itself and the summed
// base reaction IS the applied load.
constexpr double kStripQ = 100.0;       // [kPa] downward, as drawn in the file
constexpr double kStripLen = 4.0;       // [m] x = 18 .. 22

m::Project reference() {
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-fnd-008-strip-load.k2d";
    if (!m::load_project(path, pr, &err, nullptr)) {
        std::printf("FAIL: cannot load %s: %s\n", path.c_str(), err.c_str());
        ++g_failures;
    }
    return pr;
}

// The plate material the structural fixtures hang on: a real sheet-pile section, so a run
// that DOES build the structure is a well-posed problem rather than a singular one.
int add_plate_material(m::Project& pr) {
    m::PlateMaterial pm;
    pm.name = "Sheet pile";
    pm.EA = 7.5e6;      // [kN/m]
    pm.EI = 1.0e5;      // [kNm2/m]
    pm.w = 0.0;         // weightless, like the soil in this benchmark
    pm.nu = 0.0;
    pr.plates.push_back(pm);
    return (int)pr.plates.size() - 1;
}

m::StructElement plate_line(int material, double x1, double y1, double x2, double y2) {
    m::StructElement s;
    s.kind = m::StructKind::Plate;
    s.name = "Wall";
    s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2;
    s.material = material;
    return s;
}

// --------------------------------------------------------------------------- fixtures ----

// 0. The control. Nothing is perturbed, so nothing may be said: a tree that warns about a
//    clean model teaches its users to ignore warnings.
void case_reference() {
    std::printf("\n== reference (unperturbed corpus case) ==\n");
    const Run r = solve(reference());
    print_diags(r);
    check(r.ok, "reference: solves");
    check(r.diags.empty(), "reference: raises NO diagnostics");
    const double want = kStripQ * kStripLen;
    check(std::fabs(std::fabs(r.sum_reaction_y) - want) <= 1e-9 * want,
          "reference: base reactions carry the whole 400 kN/m strip");
}

// 1. The measured headline: a surcharge drawn 5 cm above the ground surface. The mesher does
//    embed a load line as a mesh constraint, so this is not an off-node placement problem --
//    the line is simply not on the soil, and every node of it lands in air.
void case_load_above_surface() {
    std::printf("\n== distributed load drawn 5 cm above the surface ==\n");
    m::Project pr = reference();
    pr.loads[0].y1 = pr.loads[0].y2 = 20.05;
    expect_refusal(solve(pr), "K2D-G003", "load above the surface");
}

// 2. Drawn outside the model altogether (x = 41 .. 45 on a 40 m wide domain).
void case_load_outside_model() {
    std::printf("\n== distributed load drawn outside the model ==\n");
    m::Project pr = reference();
    pr.loads[0].x1 = 41.0;
    pr.loads[0].x2 = 45.0;
    expect_refusal(solve(pr), "K2D-G003", "load outside the model");
}

// 3. Half outside: the mesh gives back the stretch it can, so the run is defensible -- and it
//    is a DIFFERENT load from the one drawn. The independent check is equilibrium: exactly the
//    in-soil 2 m must appear in the base reactions, which is also what the warning claims.
void case_load_half_outside() {
    std::printf("\n== distributed load half outside the model ==\n");
    m::Project pr = reference();
    pr.loads[0].x1 = 38.0;
    pr.loads[0].x2 = 45.0;   // 7 m drawn, 2 m of it on the soil
    const Run r = solve(pr);
    expect_warning(r, "K2D-G004", "load half outside");
    const double want = kStripQ * 2.0;
    check(std::fabs(std::fabs(r.sum_reaction_y) - want) <= 1e-9 * want,
          "load half outside: base reactions carry 200 kN/m -- the applied part, as warned");
}

// 4. A point load 20 m beyond the right edge and 15 m above the top. The nearest-node search
//    always succeeds, so this used to be relocated to a corner node and applied there in full;
//    the measured effect on this benchmark was a peak displacement 87 % too large.
void case_point_load_off_mesh() {
    std::printf("\n== point load placed off the mesh ==\n");
    m::Project pr = reference();
    m::Load P;
    P.kind = m::LoadKind::Point;
    P.name = "Stray point load";
    P.x1 = 60.0; P.y1 = 35.0;
    P.qx1 = 0.0; P.qy1 = -100.0;
    pr.loads.push_back(P);
    expect_refusal(solve(pr), "K2D-G001", "point load off the mesh");
}

// 5. A point load INSIDE the soil but between nodes. That is ordinary discretisation, not an
//    error -- the load is carried by the nearest node -- but the engineer should be told how
//    far it moved. The fixture is made deterministic by switching off auto-refinement (a point
//    load with coarseness 1 then contributes nothing to the size field, so the mesh is the
//    same before and after adding it) and by choosing the element centre that is farthest from
//    any node, which is this mesh's own worst case.
void case_point_load_between_nodes() {
    std::printf("\n== point load between nodes (inside the soil) ==\n");
    m::Project pr = reference();
    pr.mesh.auto_refine = false;
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "point load between nodes: reference mesh built");
    if (!M.ok) return;
    double best_x = 0.0, best_y = 0.0, best_d = -1.0;
    for (int e = 0; e < M.mesh.element_count; ++e) {
        double cx = 0.0, cy = 0.0;
        for (int k = 0; k < 3; ++k) { cx += M.mesh.x[M.mesh.node_of(e, k)] / 3.0;
                                      cy += M.mesh.y[M.mesh.node_of(e, k)] / 3.0; }
        double d = 1e300;
        for (int n = 0; n < M.mesh.node_count; ++n)
            d = std::fmin(d, std::hypot(M.mesh.x[n] - cx, M.mesh.y[n] - cy));
        if (d > best_d) { best_d = d; best_x = cx; best_y = cy; }
    }
    std::printf("      worst element centre (%.4f, %.4f) is %.4f m from any node\n",
                best_x, best_y, best_d);
    m::Load P;
    P.kind = m::LoadKind::Point;
    P.name = "Buried point load";
    P.x1 = best_x; P.y1 = best_y;
    P.qx1 = 0.0; P.qy1 = -50.0;
    pr.loads.push_back(P);
    expect_warning(solve(pr), "K2D-G002", "point load between nodes");
}

// 6. A retaining wall drawn entirely above the ground surface. Measured on v0.6.1: the run was
//    bit-identical to the model with NO wall in it -- the most dangerous row in the table,
//    because a wall is the reason such a model is built at all.
void case_plate_above_soil() {
    std::printf("\n== plate drawn entirely above the soil ==\n");
    m::Project pr = reference();
    const int pm = add_plate_material(pr);
    pr.structs.push_back(plate_line(pm, 10.0, 21.0, 30.0, 21.0));
    expect_refusal(solve(pr), "K2D-G005", "plate above the soil");
}

// 7. The same wall drawn 2 m too tall. The mesher clips a structural line to the soil in the
//    planar arrangement, so the run silently models a SHORTER wall. The independent path is
//    the wall drawn where the mesh actually took it: same mesh, same assembly, so the two runs
//    must agree bit for bit -- which is what makes the warning's account checkable.
void case_plate_clipped() {
    std::printf("\n== plate drawn 2 m above the surface (clipped) ==\n");
    m::Project drawn = reference();
    drawn.mesh.auto_refine = false;   // the size field then ignores the line, isolating the clip
    const int pm = add_plate_material(drawn);
    m::Project shortened = drawn;
    drawn.structs.push_back(plate_line(pm, 20.0, 14.0, 20.0, 22.0));       // 8 m drawn
    shortened.structs.push_back(plate_line(pm, 20.0, 14.0, 20.0, 20.0));   // 6 m of soil

    const Run r = solve(drawn);
    expect_warning(r, "K2D-G006", "plate clipped at the surface");
    const Run s = solve(shortened);
    check(s.ok, "plate clipped: the explicitly shortened wall solves");
    check(s.diags.empty(), "plate clipped: the explicitly shortened wall raises nothing");
    std::printf("      clipped max|u| = %.12e m; shortened max|u| = %.12e m\n",
                r.max_disp, s.max_disp);
    check(std::fabs(r.max_disp - s.max_disp) <= 1e-12,
          "plate clipped: the run IS the 6 m wall, to the last bit");
}

// 8. A geogrid drawn off the soil: same rule as the plate, its own code so a script can tell
//    which object it lost.
void case_geogrid_above_soil() {
    std::printf("\n== geogrid drawn entirely above the soil ==\n");
    m::Project pr = reference();
    m::GeogridMaterial gm;
    gm.name = "Geogrid";
    gm.EA = 1.0e4;
    pr.geogrids.push_back(gm);
    m::StructElement s = plate_line((int)pr.geogrids.size() - 1, 10.0, 21.0, 30.0, 21.0);
    s.kind = m::StructKind::Geogrid;
    s.name = "Reinforcement";
    pr.structs.push_back(s);
    expect_refusal(solve(pr), "K2D-G007", "geogrid above the soil");
}

// 9. An anchor with neither end in the soil. A strutted excavation that quietly loses its
//    strut is precisely the run that must never report success.
void case_anchor_outside_soil() {
    std::printf("\n== anchor with neither end in the soil ==\n");
    m::Project pr = reference();
    m::AnchorMaterial am;
    am.name = "Strut";
    am.EA = 2.0e5;
    am.Lspacing = 5.0;
    pr.anchors.push_back(am);
    m::StructElement s = plate_line((int)pr.anchors.size() - 1, 45.0, 25.0, 50.0, 30.0);
    s.kind = m::StructKind::Anchor;
    s.name = "Stray strut";
    pr.structs.push_back(s);
    expect_refusal(solve(pr), "K2D-G008", "anchor outside the soil");
}

// 10. An embedded wall (a plate carrying interfaces) drawn above the soil. Two things are then
//     true at once and both are said: the mesh cannot be split along the line, so the wall
//     would fall back to a plate bonded to the soil (K2D-G009), and that plate does not lie on
//     the mesh either, so there is nothing to build (K2D-G005). The refusal is what stops the
//     run; the warning is what explains it.
void case_wall_above_soil() {
    std::printf("\n== embedded wall (plate + interfaces) drawn above the soil ==\n");
    m::Project pr = reference();
    const int pm = add_plate_material(pr);
    m::StructElement s = plate_line(pm, 20.0, 21.0, 20.0, 27.0);
    s.iface_pos = s.iface_neg = true;
    pr.structs.push_back(s);
    const Run r = solve(pr);
    expect_refusal(r, "K2D-G005", "wall above the soil");
    check(raised(r, "K2D-G009", core::DiagnosticSeverity::Warning),
          "wall above the soil: the bonded fallback is stated too (K2D-G009)");
}

// 11. A parameter the selected model does not read. The Rankine tension cut-off is applied by the
//     Mohr-Coulomb return only; the hardening and soft-soil integrators ignore it, and the schema
//     (like PLAXIS) switches it ON by default -- so the soil takes tension it was told not to,
//     with nothing said. The mesh is coarsened and the load reduced on purpose: this fixture is
//     about what the run declares, and a full Hardening Soil solve on the benchmark's own mesh
//     would cost minutes to assert one string.
void case_tension_cutoff_ignored() {
    std::printf("\n== tension cut-off set on a model that does not read it ==\n");
    m::Project pr = reference();
    pr.materials[0].model = m::SoilModel::HardeningSoil;
    pr.materials[0].tension_cutoff = true;
    pr.materials[0].tensile_strength = 0.0;
    pr.mesh.elem_size = 5.0;
    pr.mesh.auto_refine = false;
    pr.loads[0].qy1 = pr.loads[0].qy2 = -10.0;
    const Run r = solve(pr);
    print_diags(r);
    check(raised(r, "K2D-M001", core::DiagnosticSeverity::Warning),
          "tension cut-off ignored: warns with K2D-M001");
}

// 12. A plate standing on a line that is pushed down. The prescribed-displacement ramp reaches
//     the soil elements only, so the plate reports the internal forces of an undriven element
//     (measured as M ~ 0 when this path was first designed); and the reaction at those now-fixed
//     nodes is the soil's alone. Both facts are stated where they can be acted on.
void case_prescribed_disp_on_structure() {
    std::printf("\n== plate standing on a prescribed-displacement line ==\n");
    m::Project pr = reference();
    const int pm = add_plate_material(pr);
    pr.structs.push_back(plate_line(pm, 18.0, 20.0, 22.0, 20.0));
    pr.loads.clear();                       // the settlement is the action here
    m::PrescribedDisp D;
    D.name = "Footing settlement";
    D.x1 = 18.0; D.y1 = 20.0; D.x2 = 22.0; D.y2 = 20.0;
    D.set_uy = true; D.uy = -0.01;
    pr.disps.push_back(D);
    const Run r = solve(pr);
    print_diags(r);
    check(r.ok, "prescribed displacement on a plate: the run completes");
    check(raised(r, "K2D-A003", core::DiagnosticSeverity::Warning),
          "prescribed displacement on a plate: warns that the plate does not see it (K2D-A003)");
    check(raised(r, "K2D-A004", core::DiagnosticSeverity::Note),
          "structure on a support: the reaction's missing structural share is noted (K2D-A004)");

    // The SAME line prescribed to zero is not a motion, it is a support -- and a support line
    // is the only way this schema can hold a plate at a point (edge_bc reaches model edges
    // only). Nothing is understated there, so the warning must stay quiet; if it did not, it
    // would send the user away from a force diagram that is correct. Pinned as deliberately as
    // the warning above: a later tightening that forgets the distinction would silently put it
    // back. The build's own beam verification (KV-STR-003) reads M off exactly such a model.
    m::Project sup = pr;
    sup.disps[0].name = "Support";
    sup.disps[0].uy = 0.0;
    const Run rs = solve(sup);
    print_diags(rs);
    check(rs.ok, "zero-valued support line on a plate: the run completes");
    check(!raised(rs, "K2D-A003", core::DiagnosticSeverity::Warning),
          "a support line (value 0) drives nothing, so it does NOT raise K2D-A003");
}

// 13. A linear Dynamic phase reports zeros for stress -- not because the soil is unstressed, but
//     because the linear path never recovers the field. Driven from the corpus resonant column,
//     shortened to a few steps: this fixture is about what the phase SAYS, not what it computes.
void case_linear_dynamic_zero_stress() {
    std::printf("\n== linear dynamic phase: stress never recovered ==\n");
    m::Project pr;
    std::string err;
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/kv-dyn-002-resonant-column.k2d";
    check(m::load_project(path, pr, &err, nullptr), "resonant column loads");
    while (pr.phases.size() > 1) pr.phases.pop_back();
    pr.phases[0].duration = 0.1;
    pr.phases[0].time_steps = 10;
    const Run r = solve(pr);
    print_diags(r);
    check(r.ok, "linear dynamic: the run completes");
    check(raised(r, "K2D-A001", core::DiagnosticSeverity::Warning),
          "linear dynamic: the zero stress field is declared (K2D-A001)");
}

}  // namespace

int main() {
    std::printf("== no input may be discarded in silence (KV-DIA-001) ==\n");
    case_reference();
    case_load_above_surface();
    case_load_outside_model();
    case_load_half_outside();
    case_point_load_off_mesh();
    case_point_load_between_nodes();
    case_plate_above_soil();
    case_plate_clipped();
    case_geogrid_above_soil();
    case_anchor_outside_soil();
    case_wall_above_soil();
    case_tension_cutoff_ignored();
    case_prescribed_disp_on_structure();
    case_linear_dynamic_zero_stress();
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
