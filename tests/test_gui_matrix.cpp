// Systematic integrated-path verification (the "1000% sure" sweep): structures + load + water are
// accounted for and computed correctly across soil models (LE/MC), element orders (tri6/tri15) and
// EVERY structural element. Reference-free but SHARP invariants run through the exact GUI path
// (mesh_from_project -> solve_gravity_le):
//
//   (I1) K0 ADMISSIBILITY: under self-weight the undisturbed geostatic state must give ~zero
//        displacement for every model x order x structure -- a wrong K0 / structure-assembly / pre
//        interaction produces spurious "installation" movement. The sharpest pre+compute check.
//   (I2) Each structural element MATTERS: under a load the result must differ from the no-structure
//        case (the element is actually assembled / accounted for) and stay finite and bounded.
//   (I3) tri6 ~ tri15 AGREEMENT: both element orders give a consistent settlement (mesh-objective;
//        confirms the tri15 GUI path is wired and correct, not just tri6).
//   (I4) WATER: the buoyant effective stress sigma'_v=-gamma'(H-y) is recovered (water accounted for,
//        post-process effective), on both tri6 and tri15.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

namespace m = katai::model;
using katai::app::InitialPhase;
using katai::app::SolveResult;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

enum Struct { None, Plate, Anchor, Geogrid, Embedded };
const char* kStructName[] = {"none", "plate", "anchor", "geogrid", "embedded-beam"};

// 20x10 soil block; optional structure; optional water; optional surcharge point load at the crest.
m::Project block(m::SoilModel model, Struct st, bool water, double q) {
    m::Project pr;
    m::Material s; s.model = model;
    s.E = 2.0e4; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.c = 10.0; s.phi = 25.0; s.psi = 0.0; s.k0_auto = true;
    s.E50ref = 2.0e4; s.Eoedref = 2.0e4; s.Eurref = 6.0e4; s.m = 0.5; s.p_ref = 100.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free,       (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = water;
    if (water) { pr.wx = {0.0, 20.0}; pr.wy = {8.0, 8.0}; }  // water table at y=8
    if (st == Plate) {
        m::PlateMaterial pm; pm.EA = 5.0e6; pm.EI = 8.5e3; pm.nu = 0.0; pr.plates.push_back(pm);
        m::StructElement e; e.kind = m::StructKind::Plate; e.material = 0;
        e.x1 = 10; e.y1 = 2; e.x2 = 10; e.y2 = 10; pr.structs.push_back(e);
    } else if (st == Anchor) {
        // Fixed-end vertical anchor (far end above the soil): resists the settlement of the node it
        // grips, so a vertical crest load engages it (a horizontal node-to-node strut between two
        // symmetric points would see ~zero axial strain under a symmetric vertical load -- correct,
        // but then it would not register as "mattering").
        m::AnchorMaterial am; pr.anchors.push_back(am);
        m::StructElement e; e.kind = m::StructKind::Anchor; e.material = 0;
        e.x1 = 10; e.y1 = 9; e.x2 = 10; e.y2 = 18; pr.structs.push_back(e);
    } else if (st == Geogrid) {
        m::GeogridMaterial gm; gm.EA = 1.0e4; pr.geogrids.push_back(gm);
        m::StructElement e; e.kind = m::StructKind::Geogrid; e.material = 0;
        e.x1 = 2; e.y1 = 5; e.x2 = 18; e.y2 = 5; pr.structs.push_back(e);
    } else if (st == Embedded) {
        // WEIGHTLESS in this matrix, and deliberately so. The embedded beam is the only structure
        // whose default material carries weight (gamma = 24), and a pile's own weight is NOT part
        // of what the K0 procedure seeds -- it is a genuine out-of-balance load, so the pile
        // settles on its skin and foot springs and the phase resolves it (SolveResult.nil_step).
        // That is physics, not a wiring fault, and leaving it in would make invariant (I1) measure
        // it instead of what (I1) is about. The effect is not dropped: test_k0_admissibility pins
        // it explicitly below, because a bound relaxed on purpose has to be pinned as consciously
        // as one that fires.
        m::EmbeddedBeamMaterial em; em.gamma = 0.0; pr.embedded.push_back(em);
        m::StructElement e; e.kind = m::StructKind::EmbeddedBeam; e.material = 0;
        e.x1 = 10; e.y1 = 1.5; e.x2 = 10; e.y2 = 9.5; pr.structs.push_back(e);
    }
    if (q != 0.0) { m::Load L; L.kind = m::LoadKind::Point; L.x1 = 10; L.y1 = 10; L.qx1 = 0; L.qy1 = q; pr.loads.push_back(L); }
    return pr;
}

SolveResult run(m::SoilModel model, int order, Struct st, bool water, double q, InitialPhase ph) {
    const m::Project pr = block(model, st, water, q);
    const auto M = katai::app::mesh_from_project(pr, 1.0, order);
    if (!M.ok) { SolveResult R; R.message = "mesh failed"; return R; }
    return katai::app::solve_gravity_le(pr, M.mesh, ph);
}

void test_k0_admissibility() {
    // (I1) Every model x order x structure: self-weight K0 -> ~zero displacement.
    // The loop said "every model" while covering two of four. Hardening Soil and HSsmall were outside
    // the matrix -- and HSsmall's K0 phase did not converge AT ALL (its seed was computed with E_ur,ref
    // while the model evaluated with E_0,ref, so the pre-stress sat off the yield surface). Nothing
    // caught it because test_hssmall is a material-point test: no GUI path exercised HSsmall in any
    // phase. A model the editor offers must be in the invariant matrix, or "every model" is a claim the
    // test does not keep.
    for (auto model : {m::SoilModel::LinearElastic, m::SoilModel::MohrCoulomb,
                       m::SoilModel::HardeningSoil, m::SoilModel::HSsmall}) {
        for (int order : {6, 15}) {
            for (int s = None; s <= Embedded; ++s) {
                const auto R = run(model, order, (Struct)s, false, 0.0, InitialPhase::K0Procedure);
                char tag[96];
                std::snprintf(tag, sizeof(tag), "K0 admissible: model=%d tri%d struct=%s (max|u|=%.2e)",
                              (int)model, order, kStructName[s], R.max_disp);
                check(R.ok && R.max_disp < 5e-4, tag);
            }
        }
    }
    // The confound removed from the matrix above, pinned here rather than dropped: a pile that
    // DOES carry its own weight settles measurably more in the very same K0 phase, because that
    // weight is not in the geostatic seed. If this ever stopped being true, the reason the matrix
    // runs a weightless pile would have quietly stopped applying.
    m::Project pw = block(m::SoilModel::LinearElastic, Embedded, false, 0.0);
    pw.embedded[0].gamma = 24.0;
    const auto Mw = katai::app::mesh_from_project(pw, 1.0, 6);
    const auto Rw = katai::app::solve_gravity_le(pw, Mw.mesh, InitialPhase::K0Procedure);
    const auto R0 = run(m::SoilModel::LinearElastic, 6, Embedded, false, 0.0, InitialPhase::K0Procedure);
    std::printf("  K0 with a WEIGHTED pile: max|u|=%.3e vs weightless %.3e (ratio %.2f)\n",
                Rw.max_disp, R0.max_disp, Rw.max_disp / std::fmax(R0.max_disp, 1e-30));
    check(Rw.ok && Rw.max_disp > 1.5 * R0.max_disp,
          "a pile's own weight is an out-of-balance load the K0 phase must resolve");
}

void test_structures_matter() {
    // (I2) Under a 300 kN crest load, each structure changes the response vs no-structure (it is
    // assembled / accounted for) and the field is finite & bounded. LE, tri6 (fast).
    std::printf("-- (I2) each structure is accounted for (changes the loaded response) --\n");
    const auto base = run(m::SoilModel::LinearElastic, 6, None, false, -300.0, InitialPhase::K0Procedure);
    check(base.ok && base.max_disp > 0, "loaded no-structure baseline solved");
    for (int s = Plate; s <= Embedded; ++s) {
        const auto R = run(m::SoilModel::LinearElastic, 6, (Struct)s, false, -300.0, InitialPhase::K0Procedure);
        bool finite = R.ok;
        for (int i = 0; i < (int)R.disp.size(); ++i)
            if (std::isnan(R.disp[i]) || std::fabs(R.disp[i]) > 1e3) finite = false;
        const double rel = std::fabs(R.max_disp - base.max_disp) / std::fmax(base.max_disp, 1e-12);
        char tag[96];
        std::snprintf(tag, sizeof(tag), "struct=%s changes the loaded response (%.1f%%, max|u|=%.3e)",
                      kStructName[s], 100.0 * rel, R.max_disp);
        check(finite && rel > 0.01, tag);
    }
}

void test_tri6_tri15_agreement() {
    // (I3) tri6 and tri15 give a consistent settlement under SELF-WEIGHT (a smooth, non-singular field
    // -> mesh-objective). A point load is singular (settlement diverges at the point -> order-dependent),
    // so self-weight is the fair element-order comparison; both must give -gamma H^2/(2 E_oed).
    std::printf("-- (I3) tri6 ~ tri15 agreement (self-weight, non-singular) --\n");
    const auto a = run(m::SoilModel::LinearElastic, 6, None, false, 0.0, InitialPhase::GravityLoading);
    const auto b = run(m::SoilModel::LinearElastic, 15, None, false, 0.0, InitialPhase::GravityLoading);
    const double rel = std::fabs(a.max_disp - b.max_disp) / std::fmax(b.max_disp, 1e-12);
    std::printf("   tri6 max|u|=%.5e  tri15 max|u|=%.5e  diff=%.2f%%\n", a.max_disp, b.max_disp, 100 * rel);
    check(a.ok && b.ok && rel < 0.02, "tri6 and tri15 self-weight settlement agree within 2%");
}

void test_water_buoyancy() {
    // (I4) Water table at y=8: below it sigma'_v = -[gamma_moist*(H-8) + gamma'*(8-y)], buoyant.
    // Check at the base (y=0): sigma'_v = -(gamma_unsat*2 + (gamma_sat-gamma_w)*8). tri6 + tri15.
    std::printf("-- (I4) water -> buoyant effective stress --\n");
    const double gw = katai::app::kGammaWater;
    const double sv_base = -(18.0 * 2.0 + (20.0 - gw) * 8.0);
    for (int order : {6, 15}) {
        const auto pr = block(m::SoilModel::LinearElastic, None, true, 0.0);
        const auto M = katai::app::mesh_from_project(pr, 1.0, order);
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
        double base_sv = 0.0;
        for (int n = 0; n < M.mesh.node_count; ++n)
            if (M.mesh.y[n] < 1e-6) base_sv = std::fmin(base_sv, R.stress.stress[n](1));
        std::printf("   tri%d base sigma'_v=%.3f  exact=%.3f\n", order, base_sv, sv_base);
        char tag[80]; std::snprintf(tag, sizeof(tag), "tri%d buoyant base sigma'_v matches (<3%%)", order);
        check(R.ok && std::fabs(base_sv - sv_base) < 0.03 * std::fabs(sv_base), tag);
    }
}

}  // namespace

int main() {
    std::printf("Systematic GUI-path verification: structures x load x water x model x order\n\n");
    test_k0_admissibility();
    test_structures_matter();
    test_tri6_tri15_agreement();
    test_water_buoyancy();
    if (g_failures == 0) {
        std::printf("\nOK: integrated path verified across models, orders, structures, water\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
