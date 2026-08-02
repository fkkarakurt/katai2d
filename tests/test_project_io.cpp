// Project save/load round-trip (project_io.hpp). EVERY model field is exercised with non-default
// values, serialized to JSON, parsed back, and compared EXACTLY ("%.17g" doubles are IEEE
// round-trip exact). Also: version guard + corrupt-input rejection (no crash, honest error).
#include <katai/io/project_io.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
bool eq(double a, double b) { return a == b; }   // exact: %.17g round-trips IEEE doubles

m::Project build_full() {
    m::Project p;
    p.name = "Round-trip \"quoted\" project\nwith newline";
    p.axisymmetric = true;
    p.initial_procedure = m::InitialProcedure::GravityLoading;
    p.mesh.elem_size = 1.25; p.mesh.order = 15; p.mesh.auto_refine = false;
    p.x_min = -3.5; p.x_max = 41.25; p.y_min = -1.0; p.y_max = 22.125;
    p.has_water = true;
    p.wx = {0.0, 7.3, 20.0}; p.wy = {8.0, 7.1, 5.5};

    m::Material a;   // exercise every field with odd values
    a.name = "Klei"; a.model = m::SoilModel::HardeningSoil; a.drainage = m::Drainage::Undrained;
    a.color[0] = 0.1f; a.color[1] = 0.2f; a.color[2] = 0.3f;
    a.gamma_unsat = 16.123456789; a.gamma_sat = 19.87; a.e_init = 0.71;
    a.E = 12345.6789; a.nu = 0.33; a.c = 4.2; a.phi = 23.456; a.psi = 1.5;
    a.E_inc = 11.0; a.c_inc = 0.5; a.y_ref = 9.0;
    a.tension_cutoff = false; a.tensile_strength = 2.5;
    a.E50ref = 31000; a.Eoedref = 29000; a.Eurref = 93000;
    a.m = 0.55; a.nu_ur = 0.21; a.p_ref = 110; a.Rf = 0.87;
    a.k0nc_auto = false; a.k0nc = 0.43;
    a.G0ref = 1.1e5; a.gamma07 = 1.3e-4;
    a.lam_star = 0.123; a.kap_star = 0.0321;   // Soft Soil (MMM 10) fields
    a.mu_star = 0.0061;                        // Soft Soil Creep (MMM 11) field
    a.kx = 0.86; a.ky = 0.043;
    a.rinter_rigid = false; a.Rinter = 0.67;
    a.k0_auto = false; a.k0 = 0.61; a.oc_mode = 2; a.OCR = 1.8; a.POP = 35.0;
    p.materials.push_back(a);
    m::Material b; b.name = "Sand fill"; p.materials.push_back(b);
    // Enum value 5 (SoftSoilCreep, appended file-stable) must survive the roundtrip.
    m::Material sc; sc.name = "Creep clay"; sc.model = m::SoilModel::SoftSoilCreep;
    sc.lam_star = 0.09; sc.kap_star = 0.018; sc.mu_star = 0.004;
    p.materials.push_back(sc);

    m::PlateMaterial pm; pm.name = "Sheet pile"; pm.elastoplastic = true;
    pm.EA = 4.9e6; pm.EI = 8.4e3; pm.w = 1.2; pm.nu = 0.15; pm.Mp = 120; pm.Np = 900;
    p.plates.push_back(pm);
    m::AnchorMaterial am; am.name = "Tie"; am.elastoplastic = true;
    am.EA = 2.1e5; am.Fmax_tens = 300; am.Fmax_comp = 100; am.Lspacing = 2.5;
    p.anchors.push_back(am);
    m::GeogridMaterial gm; gm.name = "Grid"; gm.elastoplastic = true; gm.EA = 2750; gm.Np = 110;
    p.geogrids.push_back(gm);
    m::EmbeddedBeamMaterial em; em.name = "Pile row"; em.E = 2.9e7; em.gamma = 25;
    em.diameter = 0.55; em.Lspacing = 3.0; em.Tskin_max = 200; em.Fmax_base = 1000;
    p.embedded.push_back(em);

    m::SoilPolygon P1; P1.name = "Lower"; P1.material = 0; P1.coarseness = 0.25;
    P1.x = {0, 20, 20, 0}; P1.y = {0, 0, 6, 6};
    P1.edge_bc = {4, 2, 0, 2};
    P1.edge_flow = {1, 0, 2, 1}; P1.edge_head = {12.0, 0.0, 0.0, 15.0};
    p.polygons.push_back(P1);
    m::SoilPolygon P2; P2.name = "Upper"; P2.material = 1;
    P2.x = {0, 20, 20, 0}; P2.y = {6, 6, 10, 10};
    P2.edge_bc = {0, 2, 0, 2};
    p.polygons.push_back(P2);

    m::StructElement s1; s1.kind = m::StructKind::Plate; s1.name = "Wall";
    s1.x1 = 10; s1.y1 = 2; s1.x2 = 10; s1.y2 = 10; s1.material = 0;
    s1.iface_pos = true; s1.iface_neg = true; s1.iface_material = 1; s1.coarseness = 0.5;
    p.structs.push_back(s1);
    m::StructElement s2; s2.kind = m::StructKind::Anchor; s2.name = "Strut";
    s2.x1 = 10; s2.y1 = 9; s2.x2 = 14; s2.y2 = 13; s2.material = 0;
    p.structs.push_back(s2);

    m::Load L1; L1.kind = m::LoadKind::Distributed; L1.name = "q";
    L1.x1 = 4; L1.y1 = 10; L1.x2 = 16; L1.y2 = 10;
    L1.qx1 = 1.5; L1.qy1 = -50; L1.qx2 = -1.5; L1.qy2 = -75; L1.coarseness = 2.0;
    p.loads.push_back(L1);
    m::Load L2; L2.name = "P"; L2.x1 = 10; L2.y1 = 10; L2.qy1 = -300;
    p.loads.push_back(L2);

    p.initial.name = "Initial"; p.initial.poly_active = {1, 0};
    p.initial.struct_active = {0, 0}; p.initial.load_active = {0, 0};
    m::Phase f1; f1.name = "Fill"; f1.poly_active = {1, 1}; f1.struct_active = {1, 0};
    f1.load_active = {0, 1};
    p.phases.push_back(f1);
    m::Phase f2; f2.name = "FoS"; f2.type = m::PhaseType::Safety; f2.poly_active = {1, 1};
    p.phases.push_back(f2);
    return p;
}

void compare(const m::Project& A, const m::Project& B) {
    check(A.name == B.name, "project name (with quotes + newline) survives");
    check(A.axisymmetric == B.axisymmetric && eq(A.x_min, B.x_min) && eq(A.x_max, B.x_max) &&
          eq(A.y_min, B.y_min) && eq(A.y_max, B.y_max), "settings + extent exact");
    check(A.has_water == B.has_water && A.wx == B.wx && A.wy == B.wy, "water polyline exact");
    check(A.initial_procedure == B.initial_procedure, "initial procedure survives");
    check(eq(A.mesh.elem_size, B.mesh.elem_size) && A.mesh.order == B.mesh.order &&
              A.mesh.auto_refine == B.mesh.auto_refine,
          "mesh settings (elem_size / order / auto_refine) exact");

    bool mats = A.materials.size() == B.materials.size();
    for (size_t i = 0; mats && i < A.materials.size(); ++i) {
        const auto& x = A.materials[i]; const auto& y = B.materials[i];
        mats = x.name == y.name && x.model == y.model && x.drainage == y.drainage &&
               x.color[0] == y.color[0] && x.color[1] == y.color[1] && x.color[2] == y.color[2] &&
               eq(x.gamma_unsat, y.gamma_unsat) && eq(x.gamma_sat, y.gamma_sat) &&
               eq(x.e_init, y.e_init) && eq(x.E, y.E) && eq(x.nu, y.nu) && eq(x.c, y.c) &&
               eq(x.phi, y.phi) && eq(x.psi, y.psi) && eq(x.E_inc, y.E_inc) &&
               eq(x.c_inc, y.c_inc) && eq(x.y_ref, y.y_ref) &&
               x.tension_cutoff == y.tension_cutoff && eq(x.tensile_strength, y.tensile_strength) &&
               eq(x.E50ref, y.E50ref) && eq(x.Eoedref, y.Eoedref) && eq(x.Eurref, y.Eurref) &&
               eq(x.m, y.m) && eq(x.nu_ur, y.nu_ur) && eq(x.p_ref, y.p_ref) && eq(x.Rf, y.Rf) &&
               x.k0nc_auto == y.k0nc_auto && eq(x.k0nc, y.k0nc) &&
               eq(x.G0ref, y.G0ref) && eq(x.gamma07, y.gamma07) &&
               eq(x.lam_star, y.lam_star) && eq(x.kap_star, y.kap_star) &&
               eq(x.mu_star, y.mu_star) &&
               eq(x.kx, y.kx) && eq(x.ky, y.ky) &&
               x.rinter_rigid == y.rinter_rigid && eq(x.Rinter, y.Rinter) &&
               x.k0_auto == y.k0_auto && eq(x.k0, y.k0) &&
               x.oc_mode == y.oc_mode && eq(x.OCR, y.OCR) && eq(x.POP, y.POP);
    }
    check(mats, "soil materials: every field exact");

    check(A.plates.size() == 1 && B.plates.size() == 1 &&
              A.plates[0].name == B.plates[0].name &&
              A.plates[0].elastoplastic == B.plates[0].elastoplastic &&
              eq(A.plates[0].EA, B.plates[0].EA) && eq(A.plates[0].EI, B.plates[0].EI) &&
              eq(A.plates[0].w, B.plates[0].w) && eq(A.plates[0].nu, B.plates[0].nu) &&
              eq(A.plates[0].Mp, B.plates[0].Mp) && eq(A.plates[0].Np, B.plates[0].Np),
          "plate material exact");
    check(B.anchors.size() == 1 && eq(A.anchors[0].Fmax_comp, B.anchors[0].Fmax_comp) &&
              eq(A.anchors[0].Lspacing, B.anchors[0].Lspacing) &&
              A.anchors[0].elastoplastic == B.anchors[0].elastoplastic,
          "anchor material exact");
    check(B.geogrids.size() == 1 && eq(A.geogrids[0].Np, B.geogrids[0].Np), "geogrid material exact");
    check(B.embedded.size() == 1 && eq(A.embedded[0].diameter, B.embedded[0].diameter) &&
              eq(A.embedded[0].Fmax_base, B.embedded[0].Fmax_base),
          "embedded-beam material exact");

    bool polys = A.polygons.size() == B.polygons.size();
    for (size_t i = 0; polys && i < A.polygons.size(); ++i) {
        const auto& x = A.polygons[i]; const auto& y = B.polygons[i];
        polys = x.name == y.name && x.material == y.material && eq(x.coarseness, y.coarseness) &&
                x.x == y.x && x.y == y.y && x.edge_bc == y.edge_bc &&
                x.edge_flow == y.edge_flow && x.edge_head == y.edge_head;
    }
    check(polys, "polygons incl. BC / flow BC / head / coarseness exact");

    bool sts = A.structs.size() == B.structs.size();
    for (size_t i = 0; sts && i < A.structs.size(); ++i) {
        const auto& x = A.structs[i]; const auto& y = B.structs[i];
        sts = x.kind == y.kind && x.name == y.name && eq(x.x1, y.x1) && eq(x.y1, y.y1) &&
              eq(x.x2, y.x2) && eq(x.y2, y.y2) && x.material == y.material &&
              x.iface_pos == y.iface_pos && x.iface_neg == y.iface_neg &&
              x.iface_material == y.iface_material && eq(x.coarseness, y.coarseness);
    }
    check(sts, "structural elements incl. interfaces exact");

    bool lds = A.loads.size() == B.loads.size();
    for (size_t i = 0; lds && i < A.loads.size(); ++i) {
        const auto& x = A.loads[i]; const auto& y = B.loads[i];
        lds = x.kind == y.kind && x.name == y.name && eq(x.qx1, y.qx1) && eq(x.qy1, y.qy1) &&
              eq(x.qx2, y.qx2) && eq(x.qy2, y.qy2) && eq(x.coarseness, y.coarseness);
    }
    check(lds, "loads exact");

    auto same_phase = [](const m::Phase& x, const m::Phase& y) {
        return x.name == y.name && x.type == y.type && x.poly_active == y.poly_active &&
               x.struct_active == y.struct_active && x.load_active == y.load_active;
    };
    bool phs = same_phase(A.initial, B.initial) && A.phases.size() == B.phases.size();
    for (size_t i = 0; phs && i < A.phases.size(); ++i) phs = same_phase(A.phases[i], B.phases[i]);
    check(phs, "initial + staged phases exact");
}

}  // namespace

int main() {
    std::printf("Project file round-trip (JSON, no dependencies)\n\n");
    const m::Project A = build_full();
    const std::string text = m::project_to_json(A);
    std::printf("   file size: %d bytes\n", (int)text.size());

    m::Project B;
    std::string err;
    check(m::project_from_json(text, B, &err), "parse back succeeds");
    if (!err.empty()) std::printf("  (%s)\n", err.c_str());
    compare(A, B);

    // Second generation must be byte-identical (deterministic writer).
    check(m::project_to_json(B) == text, "save(load(save(p))) is byte-identical");

    // Honest failures: corrupt JSON and foreign files are rejected with an error, no crash.
    m::Project C;
    err.clear();
    check(!m::project_from_json(text.substr(0, text.size() / 2), C, &err) && !err.empty(),
          "truncated file rejected with an error message");
    err.clear();
    check(!m::project_from_json("{\"foo\": 1}", C, &err) && !err.empty(),
          "foreign JSON rejected (version guard)");
    err.clear();
    check(!m::project_from_json("not json at all", C, &err) && !err.empty(),
          "garbage rejected");

    // Forward compatibility: an OLD file (no mesh / initial_procedure keys) loads the
    // documented defaults -- the run an old project meant is the run it still gets.
    m::Project D;
    err.clear();
    check(m::project_from_json("{\"katai2d\": 1}", D, &err), "minimal v1 file loads");
    check(D.initial_procedure == m::InitialProcedure::K0Procedure &&
              D.mesh.elem_size == 2.0 && D.mesh.order == 6 && D.mesh.auto_refine,
          "absent mesh / initial_procedure keys load the documented defaults");

    if (g_failures == 0) {
        std::printf("\nOK: project file round-trip exact (all fields), honest failure modes\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
