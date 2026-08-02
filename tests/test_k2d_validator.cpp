// The .k2d input contract's executable half (validate.hpp) pinned rule by rule.
// A fully-valid project must produce ZERO issues; then every rule is probed by
// mutating one field and asserting the issue lands at the EXACT field path with
// the expected severity. Matched messages are printed, so the engineer-
// readability requirement of the Stage C gate stays reviewable in the test log.
// Also pins the reader's forward-version notes (project_io.hpp): a clamped
// enum is reported by field path, never silently absorbed.
//
// verify: none -- input-contract mechanics (field paths, severities, report
// plumbing); no external oracle exists. The parameter BOUNDS the rules encode
// cite their sources (PLAXIS MMM, van Genuchten 1980) in validate.hpp.
#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>

#include <cstdio>
#include <string>

namespace m = katai::model;
namespace io = katai::io;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// A project that exercises every object kind with defensible values.
m::Project valid_project() {
    m::Project p;
    p.name = "Validator baseline";
    p.has_water = true;
    p.wx = {0.0, 40.0};
    p.wy = {8.0, 8.0};

    m::Material mc;   // Mohr-Coulomb sand (defaults are already valid)
    mc.name = "Sand";
    p.materials.push_back(mc);
    m::Material hs;
    hs.name = "Stiff clay";
    hs.model = m::SoilModel::HardeningSoil;
    p.materials.push_back(hs);
    m::Material ss;
    ss.name = "Soft clay";
    ss.model = m::SoilModel::SoftSoil;
    p.materials.push_back(ss);
    m::Material sc;
    sc.name = "Creep clay";
    sc.model = m::SoilModel::SoftSoilCreep;
    p.materials.push_back(sc);

    p.plates.push_back(m::PlateMaterial{});
    p.anchors.push_back(m::AnchorMaterial{});
    p.geogrids.push_back(m::GeogridMaterial{});
    p.embedded.push_back(m::EmbeddedBeamMaterial{});

    m::SoilPolygon lower;
    lower.name = "Lower";
    lower.material = 0;
    lower.x = {0, 40, 40, 0};
    lower.y = {0, 0, 10, 10};
    lower.edge_bc = {4, 2, 0, 2};
    lower.edge_flow = {1, 0, 2, 0};
    lower.edge_head = {12.0, 0.0, 0.0, 0.0};
    p.polygons.push_back(lower);
    m::SoilPolygon upper;
    upper.name = "Upper";
    upper.material = 1;
    upper.x = {0, 40, 40, 0};
    upper.y = {10, 10, 20, 20};
    p.polygons.push_back(upper);

    m::StructElement wall;
    wall.kind = m::StructKind::Plate;
    wall.name = "Wall";
    wall.x1 = 20; wall.y1 = 4; wall.x2 = 20; wall.y2 = 20;
    wall.material = 0;
    p.structs.push_back(wall);
    m::StructElement tie;
    tie.kind = m::StructKind::Anchor;
    tie.name = "Tie";
    tie.x1 = 20; tie.y1 = 18; tie.x2 = 28; tie.y2 = 16;
    tie.material = 0;
    p.structs.push_back(tie);
    m::StructElement iface;
    iface.kind = m::StructKind::Interface;
    iface.name = "Joint";
    iface.x1 = 10; iface.y1 = 10; iface.x2 = 30; iface.y2 = 10;
    iface.material = -1;
    p.structs.push_back(iface);

    m::Load q;
    q.kind = m::LoadKind::Distributed;
    q.name = "Surcharge";
    q.x1 = 0; q.y1 = 20; q.x2 = 10; q.y2 = 20;
    p.loads.push_back(q);

    m::Phase build;
    build.name = "Build";
    p.phases.push_back(build);
    m::Phase consol;
    consol.name = "Consolidate";
    consol.type = m::PhaseType::Consolidation;
    consol.duration = 30.0;
    p.phases.push_back(consol);
    m::Phase quake;
    quake.name = "Earthquake";
    quake.type = m::PhaseType::Dynamic;
    quake.duration = 10.0;
    quake.time_steps = 500;
    p.phases.push_back(quake);
    m::Phase fos;
    fos.name = "FoS";
    fos.type = m::PhaseType::Safety;
    p.phases.push_back(fos);
    return p;
}

// Mutate one field, expect exactly one rule to fire at `path` with `sev`, and
// print the matched message so the catalogue stays human-reviewable.
template <class Mut>
void probe(const m::Project& base, Mut mutate, const char* path, io::Severity sev,
           const char* what) {
    m::Project p = base;
    mutate(p);
    const io::ValidationReport r = io::validate_project(p);
    const io::Issue* hit = nullptr;
    for (const io::Issue& i : r.issues)
        if (i.path == path && i.severity == sev) { hit = &i; break; }
    check(hit != nullptr, what);
    if (hit)
        std::printf("      %s %s: %s\n",
                    sev == io::Severity::Error ? "[error]  " : "[warning]", hit->path.c_str(),
                    hit->message.c_str());
    if (hit && sev == io::Severity::Error)
        check(!r.ok(), (std::string(what) + " (report not ok)").c_str());
}
}  // namespace

int main() {
    const m::Project base = valid_project();
    {
        const io::ValidationReport r = io::validate_project(base);
        check(r.issues.empty() && r.ok(), "a fully-valid project raises no issues");
        for (const io::Issue& i : r.issues)
            std::printf("      unexpected %s: %s\n", i.path.c_str(), i.message.c_str());
    }

    const auto E = io::Severity::Error;
    const auto W = io::Severity::Warning;

    // -- domain and water ------------------------------------------------------
    probe(base, [](m::Project& p) { p.x_max = p.x_min; }, "x_max", E, "zero-width domain");
    probe(base, [](m::Project& p) { p.y_max = p.y_min - 1.0; }, "y_max", E, "inverted height");
    probe(base, [](m::Project& p) { p.axisymmetric = true; p.x_min = -2.0; }, "x_min", E,
          "negative radius in an axisymmetric model");
    probe(base, [](m::Project& p) { p.wy.pop_back(); }, "wy", E, "water line x/y size mismatch");
    probe(base, [](m::Project& p) { p.wx = {5.0}; p.wy = {8.0}; }, "wx", W,
          "one-point water line is silently dry");

    // -- mesh and initial procedure --------------------------------------------
    probe(base, [](m::Project& p) { p.mesh.elem_size = 0.0; }, "mesh.elem_size", E,
          "zero target element size");
    probe(base, [](m::Project& p) { p.mesh.order = 9; }, "mesh.order", E,
          "element order neither tri6 nor tri15");
    probe(base, [](m::Project& p) { p.initial_procedure = (m::InitialProcedure)7; },
          "initial_procedure", E, "unknown initial procedure");
    probe(base, [](m::Project& p) { p.initial_procedure = m::InitialProcedure::Safety; },
          "initial_procedure", W, "initial Safety with staged phases present");

    // -- soil materials --------------------------------------------------------
    probe(base, [](m::Project& p) { p.materials[0].gamma_sat = -1.0; }, "materials[0].gamma_sat",
          E, "negative unit weight");
    probe(base, [](m::Project& p) { p.materials[0].gamma_sat = 15.0; }, "materials[0].gamma_sat",
          W, "saturated weight below unsaturated");
    probe(base, [](m::Project& p) { p.materials[0].E = 0.0; }, "materials[0].E", E,
          "zero Young's modulus");
    probe(base, [](m::Project& p) { p.materials[0].nu = 0.5; }, "materials[0].nu", E,
          "incompressible Poisson's ratio");
    probe(base, [](m::Project& p) { p.materials[0].phi = 95.0; }, "materials[0].phi", E,
          "friction angle out of range");
    probe(base, [](m::Project& p) { p.materials[0].psi = 40.0; }, "materials[0].psi", E,
          "dilatancy above friction");
    probe(base, [](m::Project& p) { p.materials[0].c = -2.0; }, "materials[0].c", E,
          "negative cohesion");
    probe(base, [](m::Project& p) { p.materials[0].c = 0.0; p.materials[0].phi = 0.0; },
          "materials[0].c", E, "no shear strength at all");
    probe(base, [](m::Project& p) { p.materials[0].tensile_strength = -5.0; },
          "materials[0].tensile_strength", E, "negative tensile strength");
    probe(base,
          [](m::Project& p) {
              p.materials[0].drainage = m::Drainage::UndrainedB;
              p.materials[0].c = 0.0;
          },
          "materials[0].c", E, "Undrained (B) without an undrained strength");
    probe(base, [](m::Project& p) { p.materials[0].drainage = m::Drainage::UndrainedB; },
          "materials[0].phi", W, "Undrained (B) ignores the entered phi");
    probe(base, [](m::Project& p) { p.materials[1].Eurref = p.materials[1].E50ref; },
          "materials[1].Eurref", E, "Eur not above E50");
    probe(base, [](m::Project& p) { p.materials[1].Eurref = 1.5 * p.materials[1].E50ref; },
          "materials[1].Eurref", W, "Eur below the recommended 2 E50");
    probe(base, [](m::Project& p) { p.materials[1].Rf = 1.0; }, "materials[1].Rf", E,
          "failure ratio at the singular limit");
    probe(base, [](m::Project& p) { p.materials[1].m = 1.5; }, "materials[1].m", E,
          "stiffness power out of range");
    probe(base, [](m::Project& p) { p.materials[1].p_ref = 0.0; }, "materials[1].p_ref", E,
          "zero reference pressure");
    probe(base,
          [](m::Project& p) {
              p.materials[1].model = m::SoilModel::HSsmall;
              p.materials[1].G0ref = 1.0e4;   // below Gur = Eur/(2(1+nu_ur)) = 37500
          },
          "materials[1].G0ref", W, "small-strain modulus below Gur");
    probe(base, [](m::Project& p) { p.materials[2].kap_star = p.materials[2].lam_star; },
          "materials[2].kapstar", E, "swelling index not below compression index");
    probe(base, [](m::Project& p) { p.materials[3].mu_star = 0.0; }, "materials[3].mustar", E,
          "zero creep index");
    probe(base, [](m::Project& p) { p.materials[0].kx = -1.0; }, "materials[0].kx", E,
          "negative permeability");
    probe(base, [](m::Project& p) { p.materials[0].gw_gn = 1.0; }, "materials[0].gw_gn", E,
          "van Genuchten n at the undefined limit");
    probe(base, [](m::Project& p) { p.materials[0].gw_Sres = 1.0; }, "materials[0].gw_Sres", E,
          "residual saturation out of range");
    probe(base,
          [](m::Project& p) {
              p.materials[0].rinter_rigid = false;
              p.materials[0].Rinter = 1.2;
          },
          "materials[0].Rinter", E, "interface factor above 1");
    probe(base,
          [](m::Project& p) {
              p.materials[0].k0_auto = false;
              p.materials[0].k0 = 0.0;
          },
          "materials[0].k0", E, "manual K0 not positive");
    probe(base, [](m::Project& p) { p.materials[0].oc_mode = 5; }, "materials[0].oc_mode", E,
          "unknown stress-history mode");
    probe(base,
          [](m::Project& p) {
              p.materials[0].oc_mode = 1;
              p.materials[0].OCR = 0.8;
          },
          "materials[0].OCR", E, "under-consolidated OCR");

    // -- structural material sets ----------------------------------------------
    probe(base, [](m::Project& p) { p.plates[0].EA = 0.0; }, "plates[0].EA", E,
          "plate without axial stiffness");
    probe(base, [](m::Project& p) { p.plates[0].elastoplastic = true; },
          "plates[0].elastoplastic", W, "elastoplastic plate with no capacity");
    probe(base, [](m::Project& p) { p.anchors[0].Lspacing = 0.0; }, "anchors[0].Lspacing", E,
          "anchor without spacing");
    probe(base, [](m::Project& p) { p.geogrids[0].EA = -1.0; }, "geogrids[0].EA", E,
          "geogrid with negative stiffness");
    probe(base, [](m::Project& p) { p.embedded[0].diameter = 0.0; }, "embedded[0].diameter", E,
          "pile without a diameter");

    // -- soil regions ----------------------------------------------------------
    probe(base, [](m::Project& p) { p.polygons.clear(); }, "polygons", E, "nothing to mesh");
    probe(base, [](m::Project& p) { p.polygons[0].y.pop_back(); }, "polygons[0].y", E,
          "polygon x/y size mismatch");
    probe(base,
          [](m::Project& p) {
              p.polygons[1].x = {0, 40};
              p.polygons[1].y = {10, 10};
          },
          "polygons[1].x", E, "two-vertex region");
    probe(base, [](m::Project& p) { p.polygons[0].material = -1; }, "polygons[0].material", E,
          "region without a material");
    probe(base, [](m::Project& p) { p.polygons[0].material = 99; }, "polygons[0].material", E,
          "region with a dangling material index");
    probe(base, [](m::Project& p) { p.polygons[0].coarseness = 0.0; },
          "polygons[0].coarseness", E, "region with zero coarseness");
    probe(base, [](m::Project& p) { p.polygons[0].edge_bc.pop_back(); },
          "polygons[0].edge_bc", E, "edge BC list shorter than the edge count");
    probe(base, [](m::Project& p) { p.polygons[0].edge_bc[1] = 9; },
          "polygons[0].edge_bc[1]", E, "unknown edge BC value");
    probe(base, [](m::Project& p) { p.polygons[0].edge_flow[2] = 7; },
          "polygons[0].edge_flow[2]", E, "unknown flow BC value");
    probe(base, [](m::Project& p) { p.polygons[0].edge_head.pop_back(); },
          "polygons[0].edge_head", E, "edge heads shorter than the edge count");

    // -- structural elements ---------------------------------------------------
    probe(base, [](m::Project& p) { p.structs[0].kind = (m::StructKind)99; },
          "structs[0].kind", E, "unknown element kind");
    probe(base,
          [](m::Project& p) {
              p.structs[0].x2 = p.structs[0].x1;
              p.structs[0].y2 = p.structs[0].y1;
          },
          "structs[0].x2", E, "zero-length element");
    probe(base, [](m::Project& p) { p.structs[0].material = -1; }, "structs[0].material", E,
          "plate without a material set");
    probe(base, [](m::Project& p) { p.structs[1].material = 4; }, "structs[1].material", E,
          "anchor with a dangling material index");
    probe(base, [](m::Project& p) { p.structs[2].material = 0; }, "structs[2].material", W,
          "interface ignores an assigned material");
    probe(base, [](m::Project& p) { p.structs[0].iface_material = 99; },
          "structs[0].iface_material", E, "dangling interface material override");
    probe(base, [](m::Project& p) { p.structs[0].coarseness = -0.5; },
          "structs[0].coarseness", E, "element with negative coarseness");

    // -- loads -----------------------------------------------------------------
    probe(base, [](m::Project& p) { p.loads[0].kind = (m::LoadKind)7; }, "loads[0].kind", E,
          "unknown load kind");
    probe(base,
          [](m::Project& p) {
              p.loads[0].x2 = p.loads[0].x1;
              p.loads[0].y2 = p.loads[0].y1;
          },
          "loads[0].x2", E, "zero-length distributed load");
    probe(base, [](m::Project& p) { p.loads[0].coarseness = 0.0; }, "loads[0].coarseness", E,
          "load with zero coarseness");

    // -- prescribed displacements (schema v2) ----------------------------------
    const auto add_disp = [](m::Project& p) {
        m::PrescribedDisp D;
        D.x1 = 0.0; D.y1 = p.y_max; D.x2 = 1.0; D.y2 = p.y_max;
        D.set_uy = true; D.uy = -0.01;
        p.disps.push_back(D);
    };
    probe(base,
          [&](m::Project& p) {
              add_disp(p);
              p.disps[0].x2 = p.disps[0].x1; p.disps[0].y2 = p.disps[0].y1;
          },
          "disps[0].x2", E, "zero-length prescribed displacement");
    probe(base,
          [&](m::Project& p) {
              add_disp(p);
              p.disps[0].set_ux = false; p.disps[0].set_uy = false;
          },
          "disps[0].set_uy", E, "prescribed displacement that prescribes nothing");
    probe(base,
          [&](m::Project& p) {
              add_disp(p);
              p.disps[0].coarseness = 0.0;
          },
          "disps[0].coarseness", E, "prescribed displacement with zero coarseness");
    probe(base,
          [&](m::Project& p) {
              add_disp(p);
              p.phases[0].type = m::PhaseType::Consolidation;
              p.phases[0].duration = 1.0;   // keep the timed-phase rules quiet
              p.phases[0].disp_active = {1};
          },
          "phases[0].disp", E, "prescribed displacement active in a time-dependent phase");
    probe(base,
          [&](m::Project& p) {
              add_disp(p);   // K0 initial (the base project's procedure) + active by default
          },
          "initial.disp", E, "prescribed displacement active under the K0 initial procedure");

    // -- phases ----------------------------------------------------------------
    probe(base, [](m::Project& p) { p.initial.type = m::PhaseType::Safety; }, "initial.type", W,
          "initial phase type is ignored");
    probe(base, [](m::Project& p) { p.phases[0].type = (m::PhaseType)99; }, "phases[0].type", E,
          "unknown phase type (the driver would silently run Plastic)");
    probe(base, [](m::Project& p) { p.phases[1].duration = 0.0; }, "phases[1].duration", E,
          "consolidation without a time interval");
    probe(base, [](m::Project& p) { p.phases[1].time_steps = 0; }, "phases[1].steps", E,
          "consolidation without time steps");
    probe(base, [](m::Project& p) { p.phases[2].seismic_wave = (m::SeismicWave)9; },
          "phases[2].seiswave", E, "unknown seismic waveform");
    probe(base, [](m::Project& p) { p.phases[2].seismic_wave = m::SeismicWave::Record; },
          "phases[2].rec", E, "record waveform without a record");
    probe(base,
          [](m::Project& p) {
              p.phases[2].seismic_wave = m::SeismicWave::Record;
              p.phases[2].accel_record = {0.0, 1.0, 0.0};
              p.phases[2].record_dt = 0.0;
          },
          "phases[2].recdt", E, "record without a sampling interval");
    probe(base, [](m::Project& p) { p.phases[2].seismic_freq = 0.0; }, "phases[2].seisfreq", E,
          "harmonic input without a frequency");
    probe(base, [](m::Project& p) { p.phases[2].damping_ratio = -0.01; }, "phases[2].damp", E,
          "negative damping");
    probe(base, [](m::Project& p) { p.phases[2].damping_ratio = 1.5; }, "phases[2].damp", W,
          "damping at or above critical");
    probe(base, [](m::Project& p) { p.phases[2].rayleigh_f2 = 0.5; }, "phases[2].rayf1", E,
          "Rayleigh targets out of order");
    probe(base, [](m::Project& p) { p.phases[2].site_class = 9; }, "phases[2].siteclass", E,
          "unknown TBDY site class");
    probe(base, [](m::Project& p) { p.phases[2].tbdy_ss = 0.0; }, "phases[2].tbdyss", E,
          "TBDY coefficients not positive");
    probe(base,
          [](m::Project& p) {
              p.phases[2].ec8_enabled = true;
              p.phases[2].ec8_ground = 7;
          },
          "phases[2].ec8gnd", E, "unknown EC8 ground type");
    probe(base, [](m::Project& p) { p.phases[0].design_approach = (m::DesignApproach)99; },
          "phases[0].design", E, "unknown design approach");
    probe(base, [](m::Project& p) { p.phases[0].poly_active = {1, 1, 1, 1, 1}; },
          "phases[0].poly", W, "more activity flags than soil regions");

    // -- reader forward-version notes (project_io.hpp) -------------------------
    {
        const std::string text =
            "{\"katai2d\":1,\"materials\":[{\"name\":\"Future\",\"model\":9,\"drainage\":7}],"
            "\"polygons\":[]}";
        m::Project p;
        std::string err;
        std::vector<io::Issue> notes;
        const bool loaded = m::project_from_json(text, p, &err, &notes);
        check(loaded, "a forward-version file still loads for display");
        check(notes.size() == 2, "both clamped enums are reported");
        bool model_note = false, drainage_note = false;
        for (const io::Issue& n : notes) {
            if (n.path == "materials[0].model" && n.severity == io::Severity::Error)
                model_note = true;
            if (n.path == "materials[0].drainage" && n.severity == io::Severity::Error)
                drainage_note = true;
            std::printf("      [note]    %s: %s\n", n.path.c_str(), n.message.c_str());
        }
        check(model_note, "the clamped soil model is reported by field path");
        check(drainage_note, "the clamped drainage is reported by field path");
        check(p.materials.size() == 1 && p.materials[0].model == m::SoilModel::LinearElastic,
              "display fallback is Linear elastic, as documented");
    }

    std::printf(g_failures ? "\n%d FAILURE(S)\n" : "\nALL OK\n", g_failures);
    return g_failures ? 1 : 0;
}
