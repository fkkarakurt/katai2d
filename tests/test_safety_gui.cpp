// Safety analysis (phi-c reduction / SRM) through the GUI compute path (build_problem, phase=Safety).
// A "Safety" phase -> slope factor of safety + failure mechanism. The core SRM is validated
// directly (test_slope: FoS=1.01 vs 1.00); here we verify the INTEGRATED GUI path reproduces it and
// returns a non-zero failure mechanism (the slip surface), answering the user's question: a slope DOES
// show collapse -- via the Safety analysis, not the K0 procedure (which is zero-displacement by design).
//
// Slope benchmark (after Griffiths & Lane): homogeneous 1:2 slope on a foundation, gamma=20.2,
// c=3 kPa, phi=19.6 deg, psi=0 -> referee FoS 1.00 (Giam & Donald 1989, Monash report 8/1989).
#include <katai/io/validate.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
    std::fflush(stdout);
}
void check(bool ok, const std::string& what) { check(ok, what.c_str()); }

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
    const double ref = 1.00;
    const double err = std::fabs(R.fos - ref) / ref * 100.0;
    // Failure mechanism: the slip surface must show NON-zero displacement (unlike the K0 procedure).
    double max_u = R.max_disp;
    std::printf("  GUI Safety: FoS = %.3f  (ref 1.00)  err = %.1f%%   mechanism max|u| = %.3e\n",
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

m::Project slope_with_structures(double slab_w) {
    m::Project pr = slope();
    m::PlateMaterial slab; slab.name = "Slab"; slab.EA = 1.0e7; slab.EI = 1.0e5; slab.w = slab_w;
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
    pr.structs.push_back(line(m::StructKind::EmbeddedBeam, "Pile", 44, 32, 44, 21));
    return pr;
}

void test_safety_structures_take_part() {
    // The structures take part in the strength-reduction search. Until 2026-09 the search was
    // handed none of them: measured on this slope, an active geogrid, anchor, plate carrying
    // 150 kN/m/m and embedded beam each returned the factor of safety of the same mesh with the
    // element deactivated, bit for bit, and every one was refused (K2D-G016). Now each kind runs,
    // and each is seen to move the factor the way it acts. These are DIRECTION witnesses only --
    // how much is verified against closed forms in test_safety_structures (KV-STR-010). Each is
    // read against the same mesh with the element deactivated: a slope drawn without the element
    // is meshed differently, and the difference between two meshes is not the element's.
    m::Project pr = slope_with_structures(0.0);
    const auto M = katai::app::mesh_from_project(pr, 5.0, 6);
    check(M.ok, "slope with a slab, a geogrid, an anchor and a pile meshed");
    if (!M.ok) return;

    m::Phase cfg = pr.initial;
    katai::app::PhaseIO io;
    io.config = &cfg;
    const auto run = [&](const m::Project& p, std::vector<char> active) {
        cfg.struct_active = std::move(active);
        return katai::app::solve_gravity_le(p, M.mesh, katai::app::InitialPhase::Safety, nullptr,
                                            io);
    };
    const auto none = run(pr, {0, 0, 0, 0});
    check(none.ok && !raised(none, "K2D-A018"),
          "with every element deactivated it runs, and says nothing about structures");
    if (!none.ok) return;
    std::printf("  structures deactivated: FoS = %.5f\n", none.fos);

    // A member across the slip surface holds the slope: the factor rises.
    struct Holder { size_t index; const char* what; };
    for (const Holder& h : {Holder{1, "a geogrid across the slip surface"},
                            Holder{2, "an anchor across the slip surface"},
                            Holder{3, "a pile through the slope face"}}) {
        std::vector<char> on(4, 0);
        on[h.index] = 1;
        const auto R = run(pr, on);
        std::printf("  %-34s FoS = %.5f\n", h.what, R.fos);
        check(R.ok && !raised(R, "K2D-G016") && raised(R, "K2D-A018"),
              std::string(h.what) + " takes part, and the run says so (K2D-A018)");
        check(R.ok && R.fos > none.fos + 0.02, std::string(h.what) + " raises the factor");
    }
    // A weight on the crest loads the slope: the factor falls. The plate's STIFFNESS is present in
    // both runs, so what is compared is its weight alone.
    const auto limp = run(pr, {1, 0, 0, 0});
    const auto heavy = run(slope_with_structures(150.0), {1, 0, 0, 0});
    std::printf("  a crest slab: FoS = %.5f weightless, %.5f at 150 kN/m/m\n", limp.fos, heavy.fos);
    check(limp.ok && heavy.ok && heavy.fos < limp.fos - 0.005,
          "a heavy crest slab lowers the factor -- its weight is in the search");

    // One element still cannot enter: a PRESTRESSED anchor. Its lock-off force belongs to a ground
    // that has already moved, and the search starts from the unstressed one. Refused by the engine
    // and by the contract, and the remedy the message names runs.
    {
        m::Project pre = pr;
        pre.anchors[0].prestress = 100.0;
        const auto R = run(pre, {0, 0, 1, 0});
        check(!R.ok && R.fos < 0.0 && raised(R, "K2D-G016"),
              "a prestressed anchor in a Safety run is REFUSED (K2D-G016), with no factor");
        pre.initial.struct_active = {0, 0, 1, 0};
        pre.initial_procedure = m::InitialProcedure::Safety;
        bool schema = false;
        for (const auto& is : katai::io::validate_project(pre).issues)
            schema |= is.severity == katai::io::Severity::Error && is.path == "initial.struct";
        check(schema, "and the .k2d contract refuses it at initial.struct");
        check(run(pre, {0, 0, 0, 0}).ok, "and deactivated in the Safety run, it runs");
    }

    // A CHAINED Safety phase with the elements still active runs too -- and hands the structural
    // state on unchanged. A Safety phase commits no stress forward; before it passed the structures
    // on as well, the phase after it met a parent with no structural state and re-developed every
    // structural force from zero although nothing had changed.
    {
        m::Project ch = pr;
        ch.initial.struct_active = {0, 0, 1, 0};
        m::Phase sf; sf.name = "FoS"; sf.type = m::PhaseType::Safety;
        sf.struct_active = {0, 0, 1, 0};
        m::Phase nil; nil.name = "Nil"; nil.struct_active = {0, 0, 1, 0};
        ch.phases = {sf, nil};
        const auto res =
            katai::app::solve_phases(ch, M.mesh, katai::app::InitialPhase::GravityLoading);
        check(res.size() == 3 && res[0].ok && res[1].ok && res[1].fos > 0.0 && res[2].ok,
              "gravity -> Safety with the anchor active -> a nil phase: all three run");
        if (res.size() == 3 && res[0].ok && res[2].ok && res[0].struct_forces.size() == 1 &&
            res[2].struct_forces.size() == 1) {
            const double N0 = res[0].struct_forces[0].max_N, N2 = res[2].struct_forces[0].max_N;
            std::printf("  anchor force: %.9f after gravity, %.9f in the nil phase after Safety\n",
                        N0, N2);
            check(N0 > 1.0 && std::fabs(N2 - N0) <= 1e-9 * N0,
                  "the nil phase after Safety carries the anchor force the gravity phase left");
            check(res[2].message.find("supplies no structural state") == std::string::npos,
                  "and does not report a parent without structural state");
        } else {
            check(false, "each static phase reports one anchor force diagram");
        }
    }

    // An interface takes part as well. Its line splits the mesh, and a split mesh with nothing
    // across the seam is two unconnected bodies -- which is what the search used to solve. A joint
    // as strong as the soil now returns a factor of safety near the benchmark's. Deactivated in the
    // Safety run, its seam is tied and the ground is continuous: the factor is the one of the same
    // line drawn as an inactive geogrid, which constrains the mesh alike and splits nothing.
    {
        m::Project pj = slope();
        pj.structs.push_back(line(m::StructKind::Interface, "Joint", 41, 30, 68, 30));
        pj.structs.back().material = -1;
        const auto MJ = katai::app::mesh_from_project(pj, 5.0, 6);
        check(MJ.ok, "slope with an interface meshed");
        if (!MJ.ok) return;
        const auto RJ = katai::app::solve_gravity_le(pj, MJ.mesh, katai::app::InitialPhase::Safety);
        std::printf("  a joint as strong as the soil: FoS = %.5f\n", RJ.fos);
        check(RJ.ok && !RJ.fos_lower_bound && std::fabs(RJ.fos - 1.00) < 0.08,
              "an interface in a Safety run takes part, and the factor is within 8% of the benchmark");
        m::Phase off = pj.initial;
        off.struct_active = {0};
        katai::app::PhaseIO oio;
        oio.config = &off;
        const auto RO = katai::app::solve_gravity_le(pj, MJ.mesh, katai::app::InitialPhase::Safety,
                                                     nullptr, oio);
        m::Project pg = slope();
        m::GeogridMaterial gm; gm.name = "Grid";
        pg.geogrids.push_back(gm);
        pg.structs.push_back(line(m::StructKind::Geogrid, "Joint", 41, 30, 68, 30));
        const auto MG = katai::app::mesh_from_project(pg, 5.0, 6);
        const auto RG = MG.ok ? katai::app::solve_gravity_le(pg, MG.mesh, katai::app::InitialPhase::Safety,
                                                             nullptr, oio)
                              : katai::core::SolveResult{};
        std::printf("  the joint deactivated: FoS = %.10f; the line as an inactive geogrid: %.10f\n",
                    RO.fos, RG.fos);
        check(RO.ok && RG.ok && RO.fos == RG.fos,
              "deactivated, the joint's seam is tied: the factor is the unsplit line's, bit for bit");
    }
}

}  // namespace

int main() {
    std::printf("Safety analysis (phi-c reduction) through the GUI compute path\n");
    test_safety_le_rejected();
    test_safety_slope_fos();
    test_safety_confined_no_mechanism();
    test_safety_hs_gated();
    test_safety_structures_take_part();
    if (g_failures == 0) {
        std::printf("\nOK: GUI Safety analysis -> factor of safety + failure mechanism\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
