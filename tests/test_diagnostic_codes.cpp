// Every diagnostic code is SAID when its condition holds and NOT said when it does not.
//
// A diagnostic is the program's statement that a run is not the one its input seems to ask for:
// a well that a consolidation phase cannot pump, a design approach that has nothing to factor,
// a pile with no section. A code that no test triggers is a statement nobody has checked is ever
// made -- the condition guarding it can rot (a renamed field, a branch that is no longer reached)
// and the run then says nothing, which is the failure the code exists to prevent. Each case
// below builds the smallest model that meets the condition and asserts the code, and the same
// model with the condition removed, asserting its absence: a code that is said always is as
// useless as one that is never said.
//
// Codes covered here: K2D-A002, A006, A009, A017, A020, A026, A035, A036, G010, G011, G012,
// G013, G014, G015, M002, M004; and the input and scope refusals I001-I005, S001-S006.
//
// verify: KV-DIA-003
//   oracle:   independent_path
//   source:   KATAI 2D input-safety property (an input may be used differently from the way it was written, but never in silence; docs/diagnostics.md): the independent path is the same model with the condition removed, which must run without the code
//   locator:  one smallest model per code (a 10 x 5 m block, tri6), solved through the staged driver with and without the condition the code names; for K2D-G013 the pile head is placed at the point of a fixed mesh farthest from every node
//   quantity: the diagnostic code raised, and whether the run continued or was refused [-]
//   expected: each code present with its condition and absent without it; the refusals stop the run, the warnings and notes do not
//   band:     exact on the code and the run status
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double W = 10.0, H = 5.0;

// A block of soil on a fixed base with rollers on its sides; dry unless a case adds water.
m::Project block(m::SoilModel model = m::SoilModel::MohrCoulomb) {
    m::Project pr;
    m::Material s; s.name = "Soil"; s.model = model;
    s.E = 2.0e4; s.nu = 0.3; s.c = 5.0; s.phi = 30.0; s.psi = 0.0;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.kx = s.ky = 0.1;
    s.E50ref = 2.0e4; s.Eoedref = 2.0e4; s.Eurref = 6.0e4;
    if (model == m::SoilModel::HoekBrown) {
        s.E = 5.0e6; s.nu = 0.25; s.gamma_unsat = s.gamma_sat = 25.0;
        s.sig_ci = 2.0e4; s.mi = 10.0; s.gsi = 50.0;
    }
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Block"; P.material = 0;
    P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.x_min = 0; pr.x_max = W; pr.y_min = 0; pr.y_max = H;
    pr.has_water = false;
    return pr;
}
void add_water(m::Project& pr, double level) {
    pr.has_water = true; pr.wx = {0.0, W}; pr.wy = {level, level};
}
m::Phase phase(m::PhaseType t, const char* name) {
    m::Phase p; p.type = t; p.name = name;
    return p;
}

std::vector<katai::app::SolveResult> run(const m::Project& pr, double h = 1.0) {
    const auto M = katai::app::mesh_from_project(pr, h, 6);
    if (!M.ok) { std::printf("  (mesh: %s)\n", M.message.c_str()); return {}; }
    return katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
}
bool said(const std::vector<katai::app::SolveResult>& res, const char* code) {
    for (const auto& r : res)
        for (const auto& d : r.diagnostics)
            if (d.code == code) return true;
    return false;
}
bool refused_with(const std::vector<katai::app::SolveResult>& res, const char* code) {
    return !res.empty() && !res.back().ok && said(res, code);
}
bool all_ok(const std::vector<katai::app::SolveResult>& res, size_t n) {
    if (res.size() != n) return false;
    for (const auto& r : res) if (!r.ok) return false;
    return true;
}
void report(const std::vector<katai::app::SolveResult>& res) {
    if (!res.empty() && !res.back().ok) std::printf("  (%s)\n", res.back().message.c_str());
}

// --- materials -------------------------------------------------------------------------------
void rock_on_automatic_k0() {
    std::printf("-- K2D-A017: a Hoek-Brown rock on the automatic K0 --\n");
    m::Project pr = block(m::SoilModel::HoekBrown);
    auto res = run(pr);
    check(all_ok(res, 1) && said(res, "K2D-A017"), "automatic K0 on a rock mass: said (K2D-A017)");
    pr.materials[0].k0_auto = false; pr.materials[0].k0 = 1.0;
    res = run(pr);
    check(all_ok(res, 1) && !said(res, "K2D-A017"), "a K0 entered for the rock: not said");
}

void undrained_hardening_soil() {
    std::printf("-- K2D-M002: an undrained Hardening Soil material --\n");
    m::Project pr = block(m::SoilModel::HardeningSoil);
    add_water(pr, H);
    pr.materials[0].drainage = m::Drainage::Undrained;
    auto res = run(pr);
    check(all_ok(res, 1) && said(res, "K2D-M002"), "undrained Hardening Soil: its pore-fluid stiffness is said (K2D-M002)");
    pr.materials[0].drainage = m::Drainage::Drained;
    res = run(pr);
    check(all_ok(res, 1) && !said(res, "K2D-M002"), "drained: not said");
}

void capped_small_strain_stiffness() {
    std::printf("-- K2D-M004: an HS small G0 above 20 times G_ur --\n");
    m::Project pr = block(m::SoilModel::HSsmall);
    const double Gur = pr.materials[0].Eurref / (2.0 * (1.0 + pr.materials[0].nu_ur));
    pr.materials[0].G0ref = 24.0 * Gur;
    auto res = run(pr);
    check(all_ok(res, 1) && said(res, "K2D-M004"), "G0 = 24 G_ur: capped and said (K2D-M004)");
    pr.materials[0].G0ref = 4.0 * Gur;
    res = run(pr);
    check(all_ok(res, 1) && !said(res, "K2D-M004"), "G0 = 4 G_ur: not said");
}

// --- structures --------------------------------------------------------------------------------
void interface_in_rock() {
    std::printf("-- K2D-G014: an interface that would take its strength from a Hoek-Brown rock --\n");
    m::Project pr = block(m::SoilModel::HoekBrown);
    pr.materials[0].k0_auto = false; pr.materials[0].k0 = 1.0;
    m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Joint";
    j.x1 = 2.0; j.y1 = 1.0; j.x2 = 8.0; j.y2 = 4.0;
    pr.structs.push_back(j);
    auto res = run(pr);
    check(refused_with(res, "K2D-G014"), "joint strength from the rock: refused (K2D-G014)");
    m::Material joint = block().materials[0]; joint.name = "Joint strength"; joint.c = 0.0; joint.phi = 25.0;
    pr.materials.push_back(joint);
    pr.structs[0].iface_material = 1;
    res = run(pr);
    report(res);
    check(all_ok(res, 1) && !said(res, "K2D-G014"), "joint naming a Mohr-Coulomb material: runs, not refused");
}

void wall_toe_outside_the_soil() {
    std::printf("-- K2D-G010: a wall whose toe is drawn below the soil --\n");
    m::Project pr = block();
    m::PlateMaterial pm; pm.name = "Wall"; pm.EA = 7.5e6; pm.EI = 1.0e6;
    pr.plates.push_back(pm);
    m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall"; w.material = 0;
    w.x1 = 5.0; w.y1 = H; w.x2 = 5.0; w.y2 = -2.0;   // toe 2 m below the base
    w.iface_pos = w.iface_neg = true;
    pr.structs.push_back(w);
    auto res = run(pr);
    check(refused_with(res, "K2D-G010"), "toe outside the soil: refused (K2D-G010)");
    pr.structs[0].y2 = 1.5;
    res = run(pr);
    report(res);
    check(all_ok(res, 1) && !said(res, "K2D-G010"), "toe inside the soil: runs");
}

void embedded_beam_without_a_section() {
    std::printf("-- K2D-G012: an embedded beam with no material, or no diameter --\n");
    m::Project pr = block();
    m::StructElement b; b.kind = m::StructKind::EmbeddedBeam; b.name = "Pile";
    b.x1 = 5.0; b.y1 = H; b.x2 = 5.0; b.y2 = 1.0; b.material = -1;
    pr.structs.push_back(b);
    auto res = run(pr);
    check(refused_with(res, "K2D-G012"), "no material: refused (K2D-G012)");
    m::EmbeddedBeamMaterial em; em.name = "Pile"; em.diameter = 0.0;
    pr.embedded.push_back(em);
    pr.structs[0].material = 0;
    res = run(pr);
    check(refused_with(res, "K2D-G012"), "zero diameter: refused (K2D-G012)");
    pr.embedded[0].diameter = 0.4;
    res = run(pr);
    report(res);
    check(all_ok(res, 1) && !said(res, "K2D-G012"), "a section: runs");
}

void pile_head_between_nodes() {
    std::printf("-- K2D-G013: a hinged pile head far from the node it is tied to --\n");
    // The mesh is built once, without the pile, and handed to every run: a pile does not
    // conform the mesh, so the head can be put at a chosen distance from the nodes. The
    // threshold is a quarter of the longest edge of the element, and inside a right-angled
    // 6-noded element no point is farther than exactly that from every node -- a structured
    // mesh of a rectangle can never raise it -- so the block is given a sloping top and the
    // point is searched over the elements and their interiors.
    m::Project pr = block();
    pr.polygons[0].y = {0, 0, H - 0.7, H};
    pr.mesh.auto_refine = false;
    katai::mesh::Mesh mesh;
    int e = -1;
    double cx = 0.0, cy = 0.0, best = 0.0;
    for (const double h : {0.9, 1.1, 1.3, 0.7}) {
        const auto M = katai::app::mesh_from_project(pr, h, 6);
        if (!M.ok) continue;
        const auto& X = M.mesh.x; const auto& Y = M.mesh.y;
        for (int k = 0; k < M.mesh.element_count; ++k) {
            const int a0 = M.mesh.node_of(k, 0), a1 = M.mesh.node_of(k, 1), a2 = M.mesh.node_of(k, 2);
            const double he = std::max({std::hypot(X[a0] - X[a1], Y[a0] - Y[a1]),
                                        std::hypot(X[a1] - X[a2], Y[a1] - Y[a2]),
                                        std::hypot(X[a2] - X[a0], Y[a2] - Y[a0])});
            for (int i = 1; i < 12; ++i)
                for (int j = 1; i + j < 12; ++j) {
                    const double l1 = i / 12.0, l2 = j / 12.0, l0 = 1.0 - l1 - l2;
                    const double x = l0 * X[a0] + l1 * X[a1] + l2 * X[a2], y = l0 * Y[a0] + l1 * Y[a1] + l2 * Y[a2];
                    if (x < 3.0 || x > 7.0 || y < 2.0 || y > 3.8) continue;
                    double d = 1e300;
                    for (int n = 0; n < M.mesh.node_count; ++n) d = std::min(d, std::hypot(X[n] - x, Y[n] - y));
                    if (d / he > best) { best = d / he; e = k; cx = x; cy = y; mesh = M.mesh; }
                }
        }
        if (best > 0.26) break;
    }
    std::printf("  head at (%.3f, %.3f): %.3f of the element size from the nearest node\n", cx, cy, best);
    check(e >= 0 && best > 0.26, "a point more than a quarter of its element from every node");
    if (e < 0 || !(best > 0.26)) return;
    m::EmbeddedBeamMaterial em; em.name = "Pile";
    pr.embedded.push_back(em);
    m::StructElement b; b.kind = m::StructKind::EmbeddedBeam; b.name = "Pile"; b.material = 0;
    b.x1 = cx; b.y1 = cy; b.x2 = cx; b.y2 = cy - 1.5;
    pr.structs.push_back(b);
    auto res = katai::app::solve_phases(pr, mesh, InitialPhase::K0Procedure);
    report(res);
    check(all_ok(res, 1) && said(res, "K2D-G013"), "head between nodes: the snap is said (K2D-G013)");
    const int n0 = mesh.node_of(e, 0);
    pr.structs[0].x1 = pr.structs[0].x2 = mesh.x[n0];
    pr.structs[0].y1 = mesh.y[n0]; pr.structs[0].y2 = mesh.y[n0] - 1.5;
    res = katai::app::solve_phases(pr, mesh, InitialPhase::K0Procedure);
    report(res);
    check(all_ok(res, 1) && !said(res, "K2D-G013"), "head on a node: not said");
}

void displacement_off_the_mesh() {
    std::printf("-- K2D-G011: a prescribed displacement that catches no node --\n");
    m::Project pr = block();
    m::PrescribedDisp d; d.name = "Push"; d.x1 = 3.0; d.y1 = H + 3.0; d.x2 = 7.0; d.y2 = H + 3.0;
    d.set_uy = true; d.uy = -0.01;
    pr.disps.push_back(d);
    pr.initial.disp_active = {0};
    pr.phases.push_back(phase(m::PhaseType::Plastic, "Push"));
    auto res = run(pr);
    check(refused_with(res, "K2D-G011"), "a line 3 m above the ground: refused (K2D-G011)");
    pr.disps[0].y1 = pr.disps[0].y2 = H;
    res = run(pr);
    report(res);
    check(all_ok(res, 2) && !said(res, "K2D-G011"), "on the ground surface: runs");
}

// --- design approach --------------------------------------------------------------------------------
void material_factors_on_rock() {
    std::printf("-- K2D-G015: material partial factors on a Hoek-Brown rock --\n");
    m::Project pr = block(m::SoilModel::HoekBrown);
    pr.materials[0].k0_auto = false; pr.materials[0].k0 = 1.0;
    m::Phase p = phase(m::PhaseType::Plastic, "Design");
    p.design_approach = m::DesignApproach::EC7_DA3;
    pr.phases.push_back(p);
    auto res = run(pr);
    check(refused_with(res, "K2D-G015"), "EC7 DA3 on the rock: refused (K2D-G015)");
    pr.phases[0].design_approach = m::DesignApproach::EC7_DA2;
    res = run(pr);
    report(res);
    check(all_ok(res, 2) && !said(res, "K2D-G015"), "EC7 DA2 (resistance factors): runs");
}

// --- Safety -----------------------------------------------------------------------------------------
void safety_notes() {
    std::printf("-- K2D-A006 / A020 / A026: a loose Safety search on a hardening, cohesive soil --\n");
    m::Project pr = block(m::SoilModel::HardeningSoil);
    // A slope, so the search ends on a mechanism rather than on its step budget.
    pr.polygons[0].x = {0, W, W, 4.0, 0}; pr.polygons[0].y = {0, 0, 2.0, H, H};
    pr.polygons[0].edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                              (int)m::BCType::Free, (int)m::BCType::Free,
                              (int)m::BCType::HorizontallyFixed};
    m::Phase s = phase(m::PhaseType::Safety, "Safety");
    s.tolerance = 1.0e-2;
    pr.phases.push_back(s);
    auto res = run(pr, 1.0);
    report(res);
    check(res.size() == 2 && said(res, "K2D-A006"), "tolerance 1e-2 on a Safety phase: said (K2D-A006)");
    check(res.size() == 2 && said(res, "K2D-A020"), "Hardening Soil in a Safety phase: said (K2D-A020)");
    check(res.size() == 2 && said(res, "K2D-A026"), "cohesion with the tension cut-off on: said (K2D-A026)");
    pr.materials[0].model = m::SoilModel::MohrCoulomb;
    pr.materials[0].tension_cutoff = false;
    pr.phases[0].tolerance = 0.0;
    res = run(pr, 1.0);
    report(res);
    check(res.size() == 2 && !said(res, "K2D-A006") && !said(res, "K2D-A020") && !said(res, "K2D-A026"),
          "Mohr-Coulomb, no cut-off, the default tolerance: none of the three");
}

// --- wells and drains --------------------------------------------------------------------------------
m::Project with_well() {
    m::Project pr = block();
    add_water(pr, H - 0.5);
    m::HydroLine w; w.name = "Well"; w.kind = m::HydroKind::Well;
    w.x1 = 5.0; w.y1 = 4.0; w.x2 = 5.0; w.y2 = 1.0; w.q = 0.5; w.h_min = 0.0;
    pr.hydros.push_back(w);
    return pr;
}

void well_in_consolidation() {
    std::printf("-- K2D-A009: a well in a consolidation phase --\n");
    m::Project pr = with_well();
    m::Phase c = phase(m::PhaseType::Consolidation, "Consolidation");
    c.duration = 1.0; c.time_steps = 4;
    pr.phases.push_back(c);
    auto res = run(pr);
    report(res);
    check(all_ok(res, 2) && said(res, "K2D-A009"), "an active well in consolidation: said (K2D-A009)");
    pr.phases[0].hydro_active = {0};
    res = run(pr);
    check(all_ok(res, 2) && !said(res, "K2D-A009"), "the well switched off: not said");
}

void well_in_a_staged_phase() {
    std::printf("-- K2D-A035 / A036: a well in a Plastic phase --\n");
    m::Project pr = with_well();
    m::Phase p = phase(m::PhaseType::Plastic, "Plastic");
    p.pore_from_flow = true;
    pr.phases.push_back(p);
    auto res = run(pr);
    report(res);
    check(all_ok(res, 2) && said(res, "K2D-A035"), "pore pressures from a flow, a well active: said (K2D-A035)");
    pr.phases[0].pore_from_flow = false;
    res = run(pr);
    check(all_ok(res, 2) && said(res, "K2D-A036") && !said(res, "K2D-A035"),
          "pore pressures from the water table, a well active: said (K2D-A036)");
    pr.initial.hydro_active = {0};   // the initial phase takes the water table too, and says so
    pr.phases[0].hydro_active = {0};
    res = run(pr);
    check(all_ok(res, 2) && !said(res, "K2D-A036") && !said(res, "K2D-A035"), "the well switched off: neither");
}

// --- dynamics -----------------------------------------------------------------------------------------
void linear_dynamics_with_structures() {
    std::printf("-- K2D-A002: a linear Dynamic phase with a geogrid --\n");
    m::Project pr = block(m::SoilModel::LinearElastic);
    m::GeogridMaterial gm; gm.name = "Grid"; gm.EA = 1.0e3;
    pr.geogrids.push_back(gm);
    m::StructElement g; g.kind = m::StructKind::Geogrid; g.name = "Grid"; g.material = 0;
    g.x1 = 1.0; g.y1 = 3.0; g.x2 = 9.0; g.y2 = 3.0;
    pr.structs.push_back(g);
    m::Phase d = phase(m::PhaseType::Dynamic, "Dynamic");
    d.seismic_wave = m::SeismicWave::Harmonic; d.seismic_amp = 0.1; d.seismic_freq = 2.0;
    d.duration = 0.2; d.time_steps = 20;
    pr.phases.push_back(d);
    auto res = run(pr);
    report(res);
    check(all_ok(res, 2) && said(res, "K2D-A002"), "linear dynamics with a geogrid: said (K2D-A002)");
    pr.phases[0].struct_active = {0};
    res = run(pr);
    report(res);
    check(all_ok(res, 2) && !said(res, "K2D-A002"), "the geogrid switched off: not said");
}
// --- input the run cannot be built from (K2D-I...) -------------------------------------------
void input_refusals() {
    std::printf("-- K2D-I001..I005: input the run cannot be built from --\n");
    {   // I001: an elastoplastic plate with no capacity
        m::Project pr = block();
        m::PlateMaterial pm; pm.name = "Wall"; pm.elastoplastic = true;
        pr.plates.push_back(pm);
        m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall"; w.material = 0;
        w.x1 = 5.0; w.y1 = H; w.x2 = 5.0; w.y2 = 2.0;
        pr.structs.push_back(w);
        auto res = run(pr);
        check(refused_with(res, "K2D-I001"), "elastoplastic plate with Mp = Np = 0: refused (K2D-I001)");
        pr.plates[0].Mp = 100.0;
        res = run(pr);
        check(all_ok(res, 1) && !said(res, "K2D-I001"), "with a moment capacity: runs");
    }
    {   // I002: an element without a valid material. The mesher refuses a region drawn without one
        // before it is meshed (it used to give it material 0 in silence); a mesh handed to the
        // run directly with a material the project does not have is refused by the run.
        m::Project pr = block();
        pr.polygons[0].material = -1;
        const auto Mr = katai::app::mesh_from_project(pr, 1.0, 6);
        std::printf("  the mesher: %s\n", Mr.message.c_str());
        check(!Mr.ok && Mr.message.find("has no material") != std::string::npos,
              "a region drawn without a material: the mesher refuses it");
        pr.polygons[0].material = 0;
        auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        check(M.ok, "with its material: meshed");
        if (M.ok) {
            M.mesh.element_material[0] = 7;
            const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
            check(refused_with(res, "K2D-I002"), "an element naming material 8 of 1: refused (K2D-I002)");
        }
    }
    {   // I003: pore pressures from a flow, and no water at all
        m::Project pr = block();
        m::Phase p = phase(m::PhaseType::Plastic, "Flow"); p.pore_from_flow = true;
        pr.phases.push_back(p);
        auto res = run(pr);
        check(refused_with(res, "K2D-I003"), "pore pressures from a flow without water: refused (K2D-I003)");
        add_water(pr, H - 1.0);
        res = run(pr);
        check(all_ok(res, 2) && !said(res, "K2D-I003"), "with a water table: runs");
    }
    {   // I004: an impermeable plate the mesh is not split along
        m::Project pr = block();
        add_water(pr, H - 1.0);
        m::PlateMaterial pm; pm.name = "Screen";
        pr.plates.push_back(pm);
        m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Screen"; w.material = 0;
        w.x1 = 5.0; w.y1 = H; w.x2 = 5.0; w.y2 = 2.0; w.flow_barrier = 1;
        pr.structs.push_back(w);
        m::Phase p = phase(m::PhaseType::Plastic, "Flow"); p.pore_from_flow = true;
        pr.phases.push_back(p);
        auto res = run(pr);
        check(refused_with(res, "K2D-I004"), "an impermeable plate without interfaces: refused (K2D-I004)");
        pr.structs[0].iface_pos = pr.structs[0].iface_neg = true;
        res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-I004"), "with interfaces: runs");
    }
    {   // I005: a phase names a material the project does not have
        m::Project pr = block();
        m::Phase p = phase(m::PhaseType::Plastic, "Improve"); p.poly_material = {3};
        pr.phases.push_back(p);
        auto res = run(pr);
        check(refused_with(res, "K2D-I005"), "a region given material 4 of 1: refused (K2D-I005)");
        pr.phases[0].poly_material = {0};
        res = run(pr);
        check(all_ok(res, 2) && !said(res, "K2D-I005"), "given material 1: runs");
    }
}

// --- combinations this build does not solve (K2D-S...) ------------------------------------------
void scope_refusals() {
    std::printf("-- K2D-S001..S006: combinations this build does not solve --\n");
    {   // S001: a prescribed displacement in a consolidation phase
        m::Project pr = block();
        add_water(pr, H);
        m::PrescribedDisp d; d.name = "Push"; d.x1 = 3.0; d.y1 = H; d.x2 = 7.0; d.y2 = H; d.uy = -0.01;
        pr.disps.push_back(d);
        pr.initial.disp_active = {0};
        m::Phase c = phase(m::PhaseType::Consolidation, "Consolidation"); c.duration = 1.0; c.time_steps = 2;
        c.disp_active = {1};
        pr.phases.push_back(c);
        auto res = run(pr);
        check(refused_with(res, "K2D-S001"), "a prescribed displacement in consolidation: refused (K2D-S001)");
        pr.phases[0].disp_active = {0};
        res = run(pr);
        check(all_ok(res, 2) && !said(res, "K2D-S001"), "switched off there: runs");
    }
    {   // S002: axisymmetry with a plate, and with a consolidation phase
        m::Project pr = block();
        pr.axisymmetric = true;
        pr.polygons[0].edge_bc[3] = (int)m::BCType::HorizontallyFixed;   // the axis
        m::PlateMaterial pm; pm.name = "Wall";
        pr.plates.push_back(pm);
        m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall"; w.material = 0;
        w.x1 = 5.0; w.y1 = H; w.x2 = 5.0; w.y2 = 2.0;
        pr.structs.push_back(w);
        auto res = run(pr);
        check(refused_with(res, "K2D-S002"), "a plate in axisymmetry: refused (K2D-S002)");
        pr.structs.clear();
        add_water(pr, H);
        m::Phase c = phase(m::PhaseType::Consolidation, "Consolidation"); c.duration = 1.0; c.time_steps = 2;
        pr.phases.push_back(c);
        res = run(pr);
        check(refused_with(res, "K2D-S002"), "a consolidation phase in axisymmetry: refused (K2D-S002)");
        pr.phases[0].type = m::PhaseType::Plastic;
        res = run(pr);
        check(all_ok(res, 2) && !said(res, "K2D-S002"), "a Plastic phase in axisymmetry: runs");
    }
    {   // S003: pore pressures from a flow in a consolidation phase
        m::Project pr = block();
        add_water(pr, H - 1.0);
        m::Phase c = phase(m::PhaseType::Consolidation, "Consolidation"); c.duration = 1.0; c.time_steps = 2;
        c.pore_from_flow = true;
        pr.phases.push_back(c);
        auto res = run(pr);
        check(refused_with(res, "K2D-S003"), "pore pressures from a flow in consolidation: refused (K2D-S003)");
        pr.phases[0].pore_from_flow = false;
        res = run(pr);
        check(all_ok(res, 2) && !said(res, "K2D-S003"), "without: runs");
    }
    {   // S004 (retired): a semi-permeable interface in the flow of a staged phase now runs
        m::Project pr = block();
        add_water(pr, H - 1.0);
        m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Screen";
        j.x1 = 5.0; j.y1 = H; j.x2 = 5.0; j.y2 = 2.0; j.flow_barrier = 2; j.hydraulic_resistance = 10.0;
        pr.structs.push_back(j);
        m::Phase p = phase(m::PhaseType::Plastic, "Flow"); p.pore_from_flow = true;
        pr.phases.push_back(p);
        auto res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-S004"),
              "a semi-permeable interface in a staged phase's flow: runs (K2D-S004 retired)");
        pr.structs[0].flow_barrier = 1;
        // An impermeable free interface ending on the ground surface: the flow's end-node split
        // once left a node with no element, and PARDISO's zero pivot was thrown out of the run.
        res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-S004"), "made impermeable: runs");
    }
    {   // S005: an interface activated in a consolidation phase
        m::Project pr = block();
        add_water(pr, H);
        m::StructElement j; j.kind = m::StructKind::Interface; j.name = "Joint";
        j.x1 = 2.0; j.y1 = 1.0; j.x2 = 8.0; j.y2 = 4.0;
        pr.structs.push_back(j);
        pr.initial.struct_active = {0};
        m::Phase c = phase(m::PhaseType::Consolidation, "Consolidation"); c.duration = 1.0; c.time_steps = 2;
        c.struct_active = {1};
        pr.phases.push_back(c);
        auto res = run(pr);
        check(refused_with(res, "K2D-S005"), "an interface activated in consolidation: refused (K2D-S005)");
        pr.phases[0].type = m::PhaseType::Plastic;
        res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-S005"), "activated in a Plastic phase: runs");
    }
    {   // S006: a line contraction on a plate put in place in the same phase
        m::Project pr = block();
        m::PlateMaterial pm; pm.name = "Lining";
        pr.plates.push_back(pm);
        m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Lining"; w.material = 0;
        w.x1 = 3.0; w.y1 = 3.0; w.x2 = 7.0; w.y2 = 3.0;
        pr.structs.push_back(w);
        pr.initial.struct_active = {0};
        m::Phase p = phase(m::PhaseType::Plastic, "Install"); p.struct_active = {1}; p.struct_contraction = {0.01};
        pr.phases.push_back(p);
        auto res = run(pr);
        check(refused_with(res, "K2D-S006"), "a contraction on a plate installed in the same phase: refused (K2D-S006)");
        pr.phases[0].struct_contraction = {0.0};
        res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-S006"), "installed without one: runs");
    }
    {   // A039: the TBDY 2018 slope verdict -- said under that design approach, not otherwise
        m::Project pr = block();
        m::Phase p = phase(m::PhaseType::Plastic, "Shaken"); p.pseudo_kh = 0.05;
        m::Phase s = phase(m::PhaseType::Safety, "Safety");
        s.design_approach = m::DesignApproach::TBDY2018_Seismic;
        pr.phases = {p, s};
        auto res = run(pr);
        check(res.size() == 3 && said(res, "K2D-A039"), "a Safety phase under TBDY 2018 (seismic): the verdict is said (K2D-A039)");
        pr.phases[1].design_approach = m::DesignApproach::None;
        res = run(pr);
        report(res);
        check(res.size() == 3 && res[2].ok && !said(res, "K2D-A039"), "the same without a design approach: no verdict");
    }
    {   // S007: a pseudo-static acceleration in an axisymmetric model
        m::Project pr = block();
        pr.axisymmetric = true;
        m::Phase p = phase(m::PhaseType::Plastic, "Shaken"); p.pseudo_kh = 0.1;
        pr.phases.push_back(p);
        auto res = run(pr);
        check(refused_with(res, "K2D-S007"), "a seismic coefficient in axisymmetry: refused (K2D-S007)");
        pr.axisymmetric = false;
        res = run(pr);
        report(res);
        check(all_ok(res, 2) && !said(res, "K2D-S007"), "the same in plane strain: runs");
    }
}
} // namespace

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);   // a crash must not take the log with it
    rock_on_automatic_k0();
    undrained_hardening_soil();
    interface_in_rock();
    wall_toe_outside_the_soil();
    embedded_beam_without_a_section();
    pile_head_between_nodes();
    capped_small_strain_stiffness();
    displacement_off_the_mesh();
    material_factors_on_rock();
    well_in_consolidation();
    well_in_a_staged_phase();
    linear_dynamics_with_structures();
    safety_notes();
    input_refusals();
    scope_refusals();
    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all diagnostic-code checks passed\n");
    return 0;
}
