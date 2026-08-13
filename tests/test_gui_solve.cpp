// End-to-end GUI solve (build_problem.hpp): a drawn soil polygon -> FE mesh -> K0 self-weight
// solve -> effective-stress / pore-pressure post-processing. Validates the wiring that the GUI
// relies on, which earlier ignored water and produced settlement-derived (wrong) stresses.
//
//   (A) Dry block, K0 procedure: undisturbed ground is in geostatic equilibrium, so self-weight
//       gives ~zero displacement and the recovered EFFECTIVE stress is the K0 field
//       sigma'_v = -gamma(H-y), sigma'_h = K0 sigma'_v.
//   (B) Same block with the water table at the surface: buoyant effective stress
//       sigma'_v = -gamma'(H-y), gamma' = gamma_sat - gamma_w, and hydrostatic pore u = gamma_w(H-y).
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace model = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

model::Project make_block() {
    model::Project pr;
    model::Material m;
    m.model = model::SoilModel::LinearElastic;   // these checks verify the LE solve mechanics
    m.gamma_unsat = 17.0; m.gamma_sat = 20.0;
    m.E = 1.0e4; m.nu = 0.3; m.phi = 30.0; m.c = 1.0;
    m.k0_auto = true;                  // K0 = 1 - sin(30) = 0.5
    pr.materials.push_back(m);
    model::SoilPolygon P; P.material = 0;
    P.x = {0.0, 20.0, 20.0, 0.0};
    P.y = {0.0, 0.0, 10.0, 10.0};
    // edges: 0 bottom, 1 right, 2 top, 3 left
    P.edge_bc = {(int)model::BCType::FullyFixed, (int)model::BCType::HorizontallyFixed,
                 (int)model::BCType::Free,       (int)model::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

// Nearest-node stress lookup helper.
void check_k0_field(const katai::app::SolveResult& R, const katai::mesh::Mesh& mesh,
                    double gamma_eff, double K0, double H, double tol, const char* tag) {
    double max_sv = 0.0, max_sh = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double y = mesh.y[n];
        const double sv = R.stress.stress[n](1), sh = R.stress.stress[n](0);
        const double sv_exact = -gamma_eff * (H - y), sh_exact = K0 * sv_exact;
        max_sv = std::fmax(max_sv, std::fabs(sv - sv_exact));
        max_sh = std::fmax(max_sh, std::fabs(sh - sh_exact));
    }
    std::printf("  [%s] max|sigma'_v err|=%.3e  max|sigma'_h err|=%.3e\n", tag, max_sv, max_sh);
    check(max_sv < tol && max_sh < tol, "recovered effective stress matches the K0 field");
}

void test_dry_block() {
    const model::Project pr = make_block();
    const auto mr = katai::app::mesh_from_project(pr, 0.5 * 1.0 * 1.0, 6);   // ~1 m elements
    check(mr.ok, "mesh built");
    const auto R = katai::app::solve_gravity_le(pr, mr.mesh);
    check(R.ok, "dry K0 solve ok");
    std::printf("  dry: max|u|=%.3e\n", R.max_disp);
    check(R.max_disp < 1e-6, "undisturbed K0 ground does not settle (geostatic equilibrium)");
    check_k0_field(R, mr.mesh, 17.0, 0.5, 10.0, 0.5, "dry");
}

void test_submerged_block() {
    model::Project pr = make_block();
    pr.has_water = true;
    pr.wx = {0.0, 20.0}; pr.wy = {10.0, 10.0};   // water table at the surface
    const auto mr = katai::app::mesh_from_project(pr, 0.5 * 1.0 * 1.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, mr.mesh);
    check(R.ok, "submerged K0 solve ok");
    std::printf("  wet: max|u|=%.3e\n", R.max_disp);
    check(R.max_disp < 1e-6, "submerged undisturbed ground does not settle");
    const double gamma_eff = 20.0 - katai::app::kGammaWater;   // buoyant
    check_k0_field(R, mr.mesh, gamma_eff, 0.5, 10.0, 0.5, "wet");
    // Hydrostatic pore pressure: 0 at the surface, gamma_w*H at the base.
    double max_pore_err = 0.0;
    for (int n = 0; n < mr.mesh.node_count; ++n) {
        const double exact = katai::app::kGammaWater * std::fmax(0.0, 10.0 - mr.mesh.y[n]);
        max_pore_err = std::fmax(max_pore_err, std::fabs(R.pore[n] - exact));
    }
    std::printf("  wet: max|pore err|=%.3e\n", max_pore_err);
    check(max_pore_err < 1e-9, "hydrostatic pore pressure u = gamma_w (H - y)");
}

// Embedded wall: a vertical plate with an interface flag becomes a barrier (mesh split + plate on
// independent DOFs + Coulomb interface, K0-seeded). The wished-in-place wall must NOT move under
// self-weight on the unstructured GUI mesh; a surcharge then drives a genuine response.
void test_embedded_wall() {
    auto make = [](bool with_load) {
        model::Project pr = make_block();
        pr.materials[0].rinter_rigid = false; pr.materials[0].Rinter = 0.5;  // reduced soil-wall strength
        model::PlateMaterial pm; pm.EA = 5.0e6; pm.EI = 8.5e3; pm.nu = 0.0;
        pr.plates.push_back(pm);
        model::StructElement s;
        s.kind = model::StructKind::Plate; s.material = 0;
        s.x1 = 10.0; s.y1 = 2.0; s.x2 = 10.0; s.y2 = 10.0;   // vertical wall, toe at y=2
        s.iface_pos = true;
        pr.structs.push_back(s);
        if (with_load) {   // strong horizontal surcharge pushing the left soil AGAINST the wall, so
            // the soil-wall interface (slip / separation) actually engages and differs from bonded.
            model::Load L; L.kind = model::LoadKind::Point; L.x1 = 5; L.y1 = 10; L.qx1 = 1500; L.qy1 = 0;
            pr.loads.push_back(L);
        }
        return pr;
    };
    // Wished in place (self-weight only).
    const model::Project pr0 = make(false);
    const auto mr = katai::app::mesh_from_project(pr0, 0.5 * 1.0 * 1.0, 6);
    check(mr.ok, "wall block meshed");
    const auto R0 = katai::app::solve_gravity_le(pr0, mr.mesh, katai::app::InitialPhase::K0Procedure);
    check(R0.ok, "embedded-wall K0 solve ok");
    if (!R0.ok) { std::printf("  (%s)\n", R0.message.c_str()); return; }
    const bool split = R0.mesh.node_count > mr.mesh.node_count;
    std::printf("  wall wished-in-place: split=%d  max|u|=%.3e\n", (int)split, R0.max_disp);
    check(split, "mesh split along the wall (embedded-wall path taken)");
    check(R0.max_disp < 1e-4, "wished-in-place wall does not move under self-weight (any mesh)");

    // With a surcharge the model responds (genuine, load-driven displacement).
    const model::Project pr1 = make(true);
    const auto mr1 = katai::app::mesh_from_project(pr1, 0.5 * 1.0 * 1.0, 6);
    const auto R1 = katai::app::solve_gravity_le(pr1, mr1.mesh, katai::app::InitialPhase::K0Procedure);
    check(R1.ok, "embedded-wall + load solve ok");
    std::printf("  wall + surcharge: max|u|=%.3e\n", R1.max_disp);
    check(R1.max_disp > 1e-3, "surcharge drives a genuine response from the wished-in-place state");

    // The interface MUST matter: a bonded plate (no interface) gives a different surcharge response
    // than the same plate WITH an interface (which allows soil-wall slip/separation).
    model::Project prb = make(true);
    prb.structs.back().iface_pos = false;   // same plate, but bonded (no interface)
    const auto mrb = katai::app::mesh_from_project(prb, 0.5 * 1.0 * 1.0, 6);
    const auto Rb = katai::app::solve_gravity_le(prb, mrb.mesh, katai::app::InitialPhase::K0Procedure);
    check(Rb.ok, "bonded-plate + load solve ok");
    std::printf("  bonded plate + surcharge: max|u|=%.3e  (interface=%.3e, diff=%.1f%%)\n",
                Rb.max_disp, R1.max_disp, 100.0 * std::fabs(R1.max_disp - Rb.max_disp) / Rb.max_disp);
    check(std::fabs(R1.max_disp - Rb.max_disp) > 0.005 * Rb.max_disp,
          "reduced-strength interface (Rinter<1) changes the result vs a bonded plate");
}

// Embedded beam (pile row): a vertical pile under a surcharge sheds the load to depth (skin + foot),
// so the surface settlement at the load is smaller than without the pile. Verifies the embedded-beam
// wiring (independent DOFs, skin/foot springs, point location in the non-conforming mesh).
void test_embedded_beam() {
    auto make = [](bool with_pile) {
        model::Project pr = make_block();
        model::Load L; L.kind = model::LoadKind::Point; L.x1 = 10; L.y1 = 10; L.qx1 = 0; L.qy1 = -600;
        pr.loads.push_back(L);
        if (with_pile) {
            model::EmbeddedBeamMaterial em;   // defaults: E=3e7, d=0.4, Lspacing=2.5
            pr.embedded.push_back(em);
            model::StructElement s; s.kind = model::StructKind::EmbeddedBeam; s.material = 0;
            s.x1 = 10; s.y1 = 1.5; s.x2 = 10; s.y2 = 9.5;   // vertical pile under the load
            pr.structs.push_back(s);
        }
        return pr;
    };
    const model::Project p0 = make(false), p1 = make(true);
    const auto m0 = katai::app::mesh_from_project(p0, 0.5 * 1.0 * 1.0, 6);
    const auto m1 = katai::app::mesh_from_project(p1, 0.5 * 1.0 * 1.0, 6);
    const auto R0 = katai::app::solve_gravity_le(p0, m0.mesh, katai::app::InitialPhase::K0Procedure);
    const auto R1 = katai::app::solve_gravity_le(p1, m1.mesh, katai::app::InitialPhase::K0Procedure);
    check(R0.ok && R1.ok, "embedded-beam solves ok");
    // mesh must be IDENTICAL with/without the pile (embedded beam is non-conforming, like an anchor).
    // The mesh must be IDENTICAL with and without the pile: the embedded beam is non-conforming
    // end to end -- the shaft crosses elements at any orientation, and even the hinged connection
    // point is tied to the NEAREST EXISTING node rather than inserted as a vertex. Inserting it
    // was measured and rejected: one isolated interior point cost +146 nodes where conforming the
    // whole 8 m shaft as a plate cost +126 (mesh_builder.cpp carries the note).
    check(m0.mesh.node_count == m1.mesh.node_count, "embedded beam does not change the mesh");
    std::printf("  embedded beam: no-pile max|u|=%.4e  pile max|u|=%.4e  ratio=%.3f\n",
                R0.max_disp, R1.max_disp, R1.max_disp / R0.max_disp);
    check(R1.max_disp < 0.98 * R0.max_disp, "embedded pile reduces the settlement under the load");
}

// Robustness: many structural elements in one model (plate wall + interface, anchor strut, geogrid,
// embedded pile) under a load must solve together WITHOUT crashing -- the "complex example" case.
void test_combined_structures() {
    model::Project pr = make_block();
    pr.materials[0].rinter_rigid = false; pr.materials[0].Rinter = 0.6;
    model::PlateMaterial pm; pr.plates.push_back(pm);
    model::AnchorMaterial am; pr.anchors.push_back(am);
    model::GeogridMaterial gm; gm.EA = 1e4; pr.geogrids.push_back(gm);
    model::EmbeddedBeamMaterial em; pr.embedded.push_back(em);
    auto add = [&](model::StructKind k, double x1, double y1, double x2, double y2, bool iface = false) {
        model::StructElement s; s.kind = k; s.material = 0; s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2;
        s.iface_pos = iface; pr.structs.push_back(s);
    };
    add(model::StructKind::Plate, 6, 2, 6, 10, true);          // vertical wall with interface
    add(model::StructKind::Anchor, 6, 9, 14, 9);               // strut from the wall
    add(model::StructKind::Geogrid, 2, 5, 18, 5);              // reinforcement layer
    add(model::StructKind::EmbeddedBeam, 14, 1.5, 14, 9.5);    // pile
    model::Load L; L.kind = model::LoadKind::Point; L.x1 = 3; L.y1 = 10; L.qx1 = 0; L.qy1 = -400;
    pr.loads.push_back(L);

    const auto mr = katai::app::mesh_from_project(pr, 0.5 * 1.0 * 1.0, 6);
    check(mr.ok, "combined model meshed");
    const auto R = katai::app::solve_gravity_le(pr, mr.mesh, katai::app::InitialPhase::K0Procedure);
    check(R.ok, "combined model (plate+interface+anchor+geogrid+pile+load) solved without crashing");
    bool finite = R.ok;
    for (int i = 0; i < (int)R.disp.size(); ++i) if (std::isnan(R.disp[i]) || std::fabs(R.disp[i]) > 1e3) finite = false;
    std::printf("  combined: ok=%d max|u|=%.4e\n", (int)R.ok, R.max_disp);
    check(finite, "combined model field is finite and bounded");
}

} // namespace

int main() {
    test_dry_block();
    test_submerged_block();
    test_embedded_wall();
    test_embedded_beam();
    test_combined_structures();
    if (g_failures == 0) {
        std::printf("OK: GUI solve (K0 + water + effective stress) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
