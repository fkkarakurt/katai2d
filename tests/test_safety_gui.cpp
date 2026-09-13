// Safety analysis (phi-c reduction / SRM) through the GUI compute path (build_problem, phase=Safety).
// PLAXIS "Safety" phase -> slope factor of safety + failure mechanism. The core SRM is validated
// directly (test_slope: FoS=1.01 vs ~0.99); here we verify the INTEGRATED GUI path reproduces it and
// returns a non-zero failure mechanism (the slip surface), answering the user's question: a slope DOES
// show collapse -- via the Safety analysis, not the K0 procedure (which is zero-displacement by design).
//
// Rocscience/Slide #1 (Griffiths & Lane): homogeneous 1:2 slope on a foundation, gamma=20.2, c=3 kPa,
// phi=19.6 deg, psi=0 -> FoS ~ 0.99 (Bishop 0.988, Spencer 0.987, Phase2 T6 0.997).
#include <katai/io/validate.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
    std::fflush(stdout);
}

m::Project slope() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3; s.gamma_unsat = 20.2; s.c = 3.0; s.phi = 19.6; s.psi = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    // CCW: base, right, top, slope face, foundation top, back.
    P.x = {20, 70, 70, 50, 30, 20};
    P.y = {20, 20, 35, 35, 25, 25};
    P.edge_bc = {(int)m::BCType::FullyFixed,       // base
                 (int)m::BCType::HorizontallyFixed, // right
                 (int)m::BCType::Free,             // top
                 (int)m::BCType::Free,             // slope face
                 (int)m::BCType::Free,             // foundation top
                 (int)m::BCType::HorizontallyFixed};// back
    pr.polygons.push_back(P);
    pr.has_water = false;
    return pr;
}

void test_safety_le_rejected() {
    // Linear-elastic soil has no strength -> Safety must refuse honestly, not return garbage.
    m::Project pr = slope();
    pr.materials[0].model = m::SoilModel::LinearElastic;
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(!R.ok && R.fos < 0.0, "Safety on linear-elastic soil is rejected (no shear strength)");
}

void test_safety_slope_fos() {
    m::Project pr = slope();
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    check(M.ok, "slope meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(R.ok, "Safety analysis ran");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    const double ref = 0.99;
    const double err = std::fabs(R.fos - ref) / ref * 100.0;
    // Failure mechanism: the slip surface must show NON-zero displacement (unlike the K0 procedure).
    double max_u = R.max_disp;
    std::printf("  GUI Safety: FoS = %.3f  (ref ~0.99, T6 0.997)  err = %.1f%%   mechanism max|u| = %.3e\n",
                R.fos, err, max_u);
    check(err < 8.0, "GUI-path slope factor of safety within 8% of the benchmark");
    check(max_u > 1e-6, "Safety returns a non-zero failure mechanism (slip surface displacement)");
}

void test_safety_hs_gated() {
    // Hardening Soil Safety is GATED honestly: phi-c reduction from a stress-free state is path-
    // unstable with the HS cap (it hangs, or returns a FALSE collapse = a misleading factor of
    // safety). Until a path-stable HS Safety exists, the analysis must REFUSE with guidance (use
    // Mohr-Coulomb c'/phi') rather than return a fake number -- the worst outcome is a wrong FoS the
    // user trusts. (The strength-reduction fix itself -- factoring the HS sub-struct's c/phi -- is in
    // place for when HS Safety is re-enabled.)
    m::Project pr = slope();
    pr.materials[0].model = m::SoilModel::HardeningSoil;
    const auto M = katai::app::mesh_from_project(pr, 25.0, 6);
    check(M.ok, "HS slope meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    std::printf("  HS Safety: ok=%d fos=%.3f  msg=\"%.60s...\"\n", (int)R.ok, R.fos, R.message.c_str());
    check(!R.ok, "HS Safety is refused (gated) -- no misleading factor of safety");
    check(R.fos < 0.0, "HS Safety returns no factor-of-safety number");
}

void test_safety_confined_no_mechanism() {
    // A laterally-confined block under self-weight has NO shear-failure mechanism, so it stays stable
    // up to srf_max. The Safety analysis must report this honestly (fos_lower_bound) instead of
    // presenting the cap (3.0) as a definitive factor of safety.
    m::Project pr;
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3; s.gamma_unsat = 20.0; s.c = 5.0; s.phi = 30.0; s.psi = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 20, 20};            // CCW rectangle
    P.edge_bc = {(int)m::BCType::FullyFixed,        // base
                 (int)m::BCType::HorizontallyFixed, // right (roller)
                 (int)m::BCType::Free,             // top
                 (int)m::BCType::HorizontallyFixed};// left (roller)
    pr.polygons.push_back(P);
    pr.has_water = false;
    const auto M = katai::app::mesh_from_project(pr, 20.0, 6);
    check(M.ok, "confined block meshed");
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety);
    check(R.ok, "Safety ran on the confined block");
    if (!R.ok) { std::printf("  (%s)\n", R.message.c_str()); return; }
    std::printf("  Confined block: FoS = %.3f  lower_bound = %d (expect: no mechanism)\n",
                R.fos, (int)R.fos_lower_bound);
    check(R.fos_lower_bound, "confined block reports NO failure mechanism (FoS is a lower bound, not the cap)");
}

bool raised(const katai::core::SolveResult& R, const char* code) {
    for (const auto& d : R.diagnostics)
        if (d.code == code) return true;
    return false;
}

m::StructElement line(m::StructKind kind, const char* name, double x1, double y1, double x2,
                      double y2) {
    m::StructElement s;
    s.kind = kind; s.name = name; s.material = 0;
    s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2;
    return s;
}

void test_safety_structures_refused() {
    // A Safety run solves the ground alone: the strength-reduction search is handed no structural
    // elements. Measured on this slope before the refusal existed, an active geogrid, anchor, plate
    // carrying 150 kN/m/m and embedded beam each returned the factor of safety of the same mesh
    // with the element deactivated, bit for bit -- stiffness and weight both gone -- and an
    // interface left the two sides of its line unconnected, so a joint as strong as the soil was
    // reported as an unstable slope. Each kind is therefore REFUSED (K2D-G016) with no number,
    // and the remedy the message names -- deactivate the elements -- is checked to RUN.
    m::Project pr = slope();
    m::PlateMaterial slab; slab.name = "Slab"; slab.EA = 1.0e7; slab.EI = 1.0e5; slab.w = 150.0;
    pr.plates.push_back(slab);
    m::GeogridMaterial grid; grid.name = "Grid"; grid.EA = 1.0e5;
    pr.geogrids.push_back(grid);
    m::AnchorMaterial tie; tie.name = "Tie"; tie.EA = 1.0e6;
    pr.anchors.push_back(tie);
    m::EmbeddedBeamMaterial pile; pile.name = "Pile"; pile.diameter = 0.6; pile.Lspacing = 2.0;
    pr.embedded.push_back(pile);
    // The slope face runs (30, 25) -> (50, 35).
    pr.structs.push_back(line(m::StructKind::Plate, "Slab", 52, 35, 68, 35));
    pr.structs.push_back(line(m::StructKind::Geogrid, "Grid", 42.5, 31, 68, 31));
    pr.structs.push_back(line(m::StructKind::Anchor, "Tie", 38.5, 29, 66, 22));
    pr.structs.push_back(line(m::StructKind::EmbeddedBeam, "Pile", 56, 35, 56, 22));
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    check(M.ok, "slope with a slab, a geogrid, an anchor and a pile meshed");
    if (!M.ok) return;

    m::Phase cfg = pr.initial;
    katai::app::PhaseIO io;
    io.config = &cfg;
    for (size_t k = 0; k < pr.structs.size(); ++k) {
        cfg.struct_active.assign(pr.structs.size(), 0);
        cfg.struct_active[k] = 1;
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety,
                                                    nullptr, io);
        const std::string what = "an active " + pr.structs[k].name +
                                 " in a Safety run is REFUSED (K2D-G016), with no factor of safety";
        check(!R.ok && R.fos < 0.0 && raised(R, "K2D-G016"), what.c_str());
    }

    // The remedy runs, and answers the question it states: the ground without the elements.
    cfg.struct_active.assign(pr.structs.size(), 0);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, katai::app::InitialPhase::Safety,
                                                nullptr, io);
    check(R.ok && !raised(R, "K2D-G016"), "with every structural element deactivated it runs");
    if (R.ok) {
        const double err = std::fabs(R.fos - 0.99) / 0.99 * 100.0;
        std::printf("  structures deactivated: FoS = %.3f (ref ~0.99)  err = %.1f%%\n", R.fos, err);
        check(err < 8.0, "and returns the unreinforced slope's factor within 8% of the benchmark");
    }

    // The same refusal in a CHAINED Safety phase, where the elements were active in the phase
    // before it -- and the same remedy there.
    {
        m::Project ch = pr;
        ch.initial.struct_active = {1, 0, 0, 0};
        m::Phase sf; sf.name = "FoS"; sf.type = m::PhaseType::Safety;
        sf.struct_active = {1, 0, 0, 0};
        ch.phases = {sf};
        auto res = katai::app::solve_phases(ch, M.mesh, katai::app::InitialPhase::GravityLoading);
        check(res.size() == 2 && res[0].ok && !res[1].ok && raised(res[1], "K2D-G016"),
              "a chained Safety phase with the slab still active is REFUSED (K2D-G016)");
        ch.phases[0].struct_active = {0, 0, 0, 0};
        res = katai::app::solve_phases(ch, M.mesh, katai::app::InitialPhase::GravityLoading);
        check(res.size() == 2 && res[1].ok && res[1].fos > 0.0,
              "and runs once the slab is deactivated in the Safety phase");
        // And the contract says so before a mesh exists.
        ch.phases[0].struct_active = {1, 0, 0, 0};
        bool schema = false;
        for (const auto& is : katai::io::validate_project(ch).issues)
            schema |= is.severity == katai::io::Severity::Error && is.path == "phases[0].struct";
        check(schema, "and the .k2d contract refuses it at phases[0].struct");
    }

    // An interface is refused the same way, and the remedy is NOT available for it: an interface
    // splits the mesh, which is fixed across phases, so it cannot be deactivated per phase. The
    // message says so; this pins the sentence to the behaviour, so that the day an interface can
    // be switched off per phase this check fails and the message is revisited with it.
    {
        m::Project pj = slope();
        pj.structs.push_back(line(m::StructKind::Interface, "Joint", 41, 30, 68, 30));
        pj.structs.back().material = -1;
        const auto MJ = katai::app::mesh_from_project(pj, 5.0, 6);
        check(MJ.ok, "slope with an interface meshed");
        if (!MJ.ok) return;
        const auto RJ = katai::app::solve_gravity_le(pj, MJ.mesh, katai::app::InitialPhase::Safety);
        check(!RJ.ok && RJ.fos < 0.0 && raised(RJ, "K2D-G016"),
              "an interface in a Safety run is REFUSED (K2D-G016)");
        m::Phase off = pj.initial;
        off.struct_active = {0};
        katai::app::PhaseIO oio;
        oio.config = &off;
        const auto RO = katai::app::solve_gravity_le(pj, MJ.mesh, katai::app::InitialPhase::Safety,
                                                     nullptr, oio);
        check(!RO.ok && !raised(RO, "K2D-G016"),
              "and it cannot be deactivated per phase, which is what the refusal tells the user");
    }
}

}  // namespace

int main() {
    std::printf("Safety analysis (phi-c reduction) through the GUI compute path\n");
    test_safety_le_rejected();
    test_safety_slope_fos();
    test_safety_confined_no_mechanism();
    test_safety_hs_gated();
    test_safety_structures_refused();
    if (g_failures == 0) {
        std::printf("\nOK: GUI Safety analysis -> factor of safety + failure mechanism\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
