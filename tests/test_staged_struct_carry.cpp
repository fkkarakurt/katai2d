// STATIC-CHAIN STRUCTURAL STATE CARRY -- the nil-phase identity (Track 1a generalized to statics).
//
// History: found 2026-07-19 as a MEASURED defect. solve_phases chained only the committed Gauss
// stresses; solve_nonlinear took no initial structural state, so every chained phase restarted its
// structures at u = 0 with zero plastic memory and the SumMstage imbalance RE-RAMPED the parent's
// structural tractions: in this very model the unchanged nil phase drifted the wall moment by 32%
// (and the K0 -> Surcharge hand-off was ~5x off the K0 phase's own converged wall state). Fixed by
// carrying the parent's displacement datum + committed plastic state into solve_nonlinear and
// putting the full structural internal force into the baseline B.
//
// The oracle is the definition of a nil phase: a chained phase whose configuration is IDENTICAL to
// its parent's changes nothing, so it must be a no-op --
//   (1) its incremental displacement must be ~0 (round-off; nothing new was applied), and
//   (2) its structural forces must equal the parent's, station for station (round-off).
// The tolerances are near-exact ON PURPOSE: any re-development drift is the silent-wrong class (a
// plausible-looking multi-phase result whose forces and plastic memory are wrong), and this test is
// the tripwire that caught it once already.
//
// Model: 20x10 soil, embedded wall (plate + Coulomb interfaces both sides) at x=10 over y=4..10,
// a one-sided surcharge left of the wall. Phases: K0 -> "Surcharge" (Plastic) -> "Nil" (Plastic,
// identical configuration).
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
}  // namespace

int main() {
    std::printf("Static-chain structural state: a NIL phase must be a no-op\n\n");
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
    s.c = 5.0; s.phi = 30.0; s.rinter_rigid = false; s.Rinter = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(P);
    m::PlateMaterial pm;
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = 24.0 * 0.4;
    pr.plates.push_back(pm);
    m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
    st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
    st.iface_pos = true; st.iface_neg = true;
    pr.structs.push_back(st);
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 0; L.y1 = 10; L.x2 = 8; L.y2 = 10;
    L.qx1 = 0; L.qy1 = -50; L.qx2 = 0; L.qy2 = -50;
    pr.loads.push_back(L);
    m::Phase p1; p1.type = m::PhaseType::Plastic; p1.name = "Surcharge";
    pr.phases.push_back(p1);
    m::Phase p2; p2.type = m::PhaseType::Plastic; p2.name = "Nil";   // IDENTICAL configuration
    pr.phases.push_back(p2);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "model meshed");
    if (!M.ok) return 1;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[0].ok && res[1].ok && res[2].ok, "K0 + Surcharge + Nil all solved");
    if (res.size() != 3 || !res[2].ok) return 1;
    const auto& par = res[1];
    const auto& nil = res[2];

    // (1) The nil phase applied nothing new -> its incremental displacement must be ~0.
    std::printf("   |u|: Surcharge phase %.6g m  vs  Nil phase %.6g m  (ratio %.3f)\n",
                par.max_disp, nil.max_disp, nil.max_disp / std::fmax(par.max_disp, 1e-30));
    check(par.max_disp > 1e-4, "the surcharge phase really deforms the model (teeth)");
    check(nil.max_disp < 1e-6 * par.max_disp,
          "NIL phase displacement ~ 0 (round-off; nothing changed, nothing may move)");

    // (2) Structural forces must be unchanged, station for station.
    double difM = 0.0, teethM = 0.0;
    for (size_t i = 0; i < nil.struct_forces.size() && i < par.struct_forces.size(); ++i) {
        const auto& fn = nil.struct_forces[i];
        const auto& fp = par.struct_forces[i];
        if (fn.stations.size() != fp.stations.size()) { check(false, "wall stations align"); continue; }
        for (size_t k = 0; k < fn.stations.size(); ++k) {
            difM = std::fmax(difM, std::fabs(fn.stations[k].M - fp.stations[k].M));
            teethM = std::fmax(teethM, std::fabs(fp.stations[k].M));
        }
    }
    std::printf("   wall:  max |M_nil - M_par| = %.4g kNm/m   (teeth: max |M_par| = %.4g)\n", difM, teethM);
    check(teethM > 1.0, "the surcharge really bends the wall (teeth)");
    check(difM < 1e-6 * teethM, "NIL phase wall forces == parent's (state carried, not re-developed)");

    double difT = 0.0, teethT = 0.0, difS = 0.0;
    for (size_t i = 0; i < nil.interface_forces.size() && i < par.interface_forces.size(); ++i) {
        const auto& jn = nil.interface_forces[i];
        const auto& jp = par.interface_forces[i];
        if (jn.stations.size() != jp.stations.size()) { check(false, "joint stations align"); continue; }
        for (size_t k = 0; k < jn.stations.size(); ++k) {
            difT = std::fmax(difT, std::fabs(jn.stations[k].tau - jp.stations[k].tau));
            difS = std::fmax(difS, std::fabs(jn.stations[k].sigma_n - jp.stations[k].sigma_n));
            teethT = std::fmax(teethT, std::fmax(std::fabs(jp.stations[k].tau),
                                                 std::fabs(jp.stations[k].sigma_n)));
        }
    }
    std::printf("   joint: max |tau_nil - tau_par| = %.4g kPa, |sn_nil - sn_par| = %.4g kPa "
                "(teeth: %.4g kPa)\n", difT, difS, teethT);
    check(teethT > 10.0, "the parent really loads the joint (teeth)");
    check(std::fmax(difT, difS) < 1e-6 * teethT,
          "NIL phase joint stresses == parent's (state carried, not re-developed)");

    if (g_failures == 0) {
        std::printf("\nOK: a nil phase is a no-op -- the static chain carries the structural state\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED -- the static chain RE-DEVELOPS the structural state "
                "(silent-wrong class; see the roadmap audit note)\n", g_failures);
    return 1;
}
