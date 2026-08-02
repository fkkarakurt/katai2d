// Structural force diagrams (M/Q/N) through the GUI compute path (build_problem.struct_forces).
// The diagram MATH is validated against closed forms in test_struct_forces (cantilever M(s)=P(L-s)
// at 2.7e-9, Barlow shear, anchor/geogrid return mapping); here the WIRING is audited: every drawn
// structural line must come back as a diagram whose stations lie on the line, with physically
// required properties (zero bending at a free plate end, tension-only geogrid, strut force pulling
// its share of the load).
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// 20x10 LE block with a crest point load; structures added per test.
m::Project block(double q) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.3e4; s.nu = 0.3; s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    if (q != 0.0) {
        m::Load L; L.kind = m::LoadKind::Point; L.x1 = 10.0; L.y1 = 10.0; L.qx1 = 0.0; L.qy1 = q;
        pr.loads.push_back(L);
    }
    return pr;
}

void test_plate_diagram() {
    // Horizontal raft on the surface, point load at its CENTRE: real bending (sagging M peak at
    // mid-span) with moment-free free ends -- the closed-form-validated diagram must show both.
    std::printf("-- plate (horizontal raft, centre point load) --\n");
    auto pr = block(-300.0);
    m::PlateMaterial pm; pm.EA = 5e6; pm.EI = 8.5e3;
    pr.plates.push_back(pm);
    m::StructElement e; e.kind = m::StructKind::Plate; e.name = "Plate 1";
    e.x1 = 6.0; e.y1 = 10.0; e.x2 = 14.0; e.y2 = 10.0; e.material = 0;
    pr.structs.push_back(e);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "plate model solved");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.struct_forces.size() == 1 && R.struct_forces[0].kind == 0, "one plate diagram returned");
    if (R.struct_forces.empty()) return;
    const auto& d = R.struct_forces[0];
    std::printf("   %d stations, max|N|=%.4g max|Q|=%.4g max|M|=%.4g\n",
                (int)d.stations.size(), d.max_N, d.max_Q, d.max_M);
    check(d.stations.size() >= 10, "diagram has stations along the wall");
    bool on_line = true, s_mono = true;
    for (size_t i = 0; i < d.stations.size(); ++i) {
        if (std::fabs(d.stations[i].y - 10.0) > 1e-6) on_line = false;
        if (i > 0 && d.stations[i].s < d.stations[i - 1].s - 1e-12) s_mono = false;
    }
    check(on_line, "stations lie on the drawn plate line (y = 10)");
    check(s_mono, "arc length is monotonic");
    check(d.max_M > 1.0, "the centre load bends the raft (M engaged)");
    // Physics: the raft ends are free (no clamp) -> moment-free; the peak is near mid-span.
    const double Mend = std::max(std::fabs(d.stations.front().M), std::fabs(d.stations.back().M));
    double s_peak = 0.0, m_peak = 0.0;
    for (const auto& st : d.stations)
        if (std::fabs(st.M) > m_peak) { m_peak = std::fabs(st.M); s_peak = st.s; }
    const double L = d.stations.back().s;
    std::printf("   |M| ends = %.4g (envelope %.4g), peak at s = %.2f of L = %.2f\n",
                Mend, d.max_M, s_peak, L);
    check(Mend < 0.12 * d.max_M, "bending moment vanishes at the free raft ends");
    check(std::fabs(s_peak - 0.5 * L) < 0.15 * L, "peak moment near mid-span (under the load)");
}

void test_anchor_force() {
    std::printf("-- anchor (fixed-end tie at the crest, downward load) --\n");
    auto pr = block(-100.0);
    m::AnchorMaterial am; am.EA = 1e5; am.Lspacing = 1.0;
    pr.anchors.push_back(am);
    // Fixed end ABOVE the soil, soil end strictly INSIDE (a point exactly on the boundary is not
    // "in soil" for the point-in-polygon test): the settling soil stretches the tie -> tension
    // N > 0, carrying a share of the 100 kN load.
    m::StructElement e; e.kind = m::StructKind::Anchor; e.name = "Anchor 1";
    e.x1 = 10.0; e.y1 = 8.0; e.x2 = 10.0; e.y2 = 14.0; e.material = 0;
    pr.structs.push_back(e);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "anchor model solved");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.struct_forces.size() == 1 && R.struct_forces[0].kind == 1, "one anchor force returned");
    if (R.struct_forces.empty() || R.struct_forces[0].stations.empty()) return;
    const double N = R.struct_forces[0].stations[0].N;
    std::printf("   anchor N = %.4g kN/m (load 100 down)\n", N);
    check(N > 1.0 && N < 100.0, "tie takes a share of the load in tension (0 < N < P)");
}

void test_geogrid_diagram() {
    std::printf("-- geogrid (horizontal, under the crest load) --\n");
    auto pr = block(-300.0);
    m::GeogridMaterial gm; gm.EA = 2e3;
    pr.geogrids.push_back(gm);
    m::StructElement e; e.kind = m::StructKind::Geogrid; e.name = "Geogrid 1";
    e.x1 = 6.0; e.y1 = 8.0; e.x2 = 14.0; e.y2 = 8.0; e.material = 0;
    pr.structs.push_back(e);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "geogrid model solved");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.struct_forces.size() == 1 && R.struct_forces[0].kind == 2, "one geogrid diagram returned");
    if (R.struct_forces.empty()) return;
    const auto& d = R.struct_forces[0];
    bool tension_only = true, on_line = true;
    for (const auto& st : d.stations) {
        if (st.N < -1e-9) tension_only = false;
        if (std::fabs(st.y - 8.0) > 1e-6) on_line = false;
    }
    std::printf("   %d stations, max N = %.4g kN/m\n", (int)d.stations.size(), d.max_N);
    check(d.stations.size() >= 4, "stations along the geogrid");
    check(on_line, "stations lie on the drawn geogrid line (y = 8)");
    check(tension_only, "geogrid axial force is tension-only (N >= 0)");
}

void test_embedded_beam_diagram() {
    // Vertical pile (embedded beam) loaded by a surface surcharge above it: the settling soil sheds
    // load into the pile through the skin, so the wired diagram must come back as N/Q/M stations
    // lying on the pile line with finite, non-trivial axial force (the diagram MATH -- N at head = P,
    // decay to the free toe -- is validated in test_embedded_beam).
    std::printf("-- embedded beam (vertical pile, surface surcharge) --\n");
    auto pr = block(0.0);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 6.0; L.y1 = 10.0; L.x2 = 14.0; L.y2 = 10.0; L.qx1 = L.qx2 = 0.0; L.qy1 = L.qy2 = -200.0;
    pr.loads.push_back(L);
    m::EmbeddedBeamMaterial bm; bm.E = 3.0e7; bm.diameter = 0.4; bm.Lspacing = 2.0;
    pr.embedded.push_back(bm);
    m::StructElement e; e.kind = m::StructKind::EmbeddedBeam; e.name = "Pile 1";
    e.x1 = 10.0; e.y1 = 2.0; e.x2 = 10.0; e.y2 = 10.0; e.material = 0;
    pr.structs.push_back(e);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "embedded beam model solved");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    check(R.struct_forces.size() == 1 && R.struct_forces[0].kind == 0,
          "one embedded beam diagram returned (N/Q/M)");
    if (R.struct_forces.empty()) return;
    const auto& d = R.struct_forces[0];
    bool on_line = true, s_mono = true, finite = true;
    for (size_t i = 0; i < d.stations.size(); ++i) {
        if (std::fabs(d.stations[i].x - 10.0) > 1e-6) on_line = false;
        if (i > 0 && d.stations[i].s < d.stations[i - 1].s - 1e-12) s_mono = false;
        if (!std::isfinite(d.stations[i].N) || !std::isfinite(d.stations[i].M)) finite = false;
    }
    std::printf("   %d stations, max|N|=%.4g max|Q|=%.4g max|M|=%.4g\n",
                (int)d.stations.size(), d.max_N, d.max_Q, d.max_M);
    check(d.stations.size() >= 6, "diagram has stations along the pile");
    check(on_line, "stations lie on the drawn pile line (x = 10)");
    check(s_mono && finite, "arc length monotonic + forces finite");
    check(d.max_N > 1.0, "the pile carries axial load shed from the surcharge (N engaged)");
}

}  // namespace

int main() {
    std::printf("Structural force diagrams through the GUI compute path\n\n");
    test_plate_diagram();
    test_anchor_force();
    test_geogrid_diagram();
    test_embedded_beam_diagram();
    if (g_failures == 0) {
        std::printf("\nOK: GUI structural force output wired (plate M/Q/N, anchor N, geogrid N)\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
