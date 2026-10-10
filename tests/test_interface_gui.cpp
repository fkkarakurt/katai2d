// Standalone interface through the GUI COMPUTE PATH (build_problem). This is the property a user checks
// when they draw an interface and toggle it on/off: the result must CHANGE. Before this work a standalone
// StructKind::Interface was dropped silently (not meshed, not assembled) -> identical results. Now the
// interface line is meshed as a constraint, the mesh is split along it, and a soil-soil Coulomb joint is
// assembled, so a weak (Rinter<1) interface relaxes shear across the seam and changes the settlement.
//   (1) WITHOUT interface vs WITH a vertical interface (Rinter<1) under an off-centre surcharge: the
//       max settlement differs by a clear margin (interface present != absent).
//   (2) The interface actually reaches the solver: SolveResult reports interface elements / the run stays
//       well-posed (ok, finite, no NaN).
//   (3) Holds for tri6 and tri15.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>

// THE START IS THE K0 PROCEDURE, in every model here, with or without an interface. These cases
// used to ask for gravity loading, and a model with an interface in it was then silently started
// by the K0 procedure instead while the one without was not: "the interface changes the result"
// compared a ground that had settled under its own weight with one that had not, and "the wall
// does not deflect under self-weight" held because the weight had never been applied. Gravity
// loading is now gravity loading whatever is in the model, so the statements are made on the
// start they are about -- a geostatic state in balance, the same for the models compared --
// and what they measure is the surcharge and the interface.
namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

// 20x10 LE block, fixed base + laterally fixed sides + free top, with a downward surcharge on the RIGHT
// half of the crest. mode: 0 = plain; 1 = a vertical interface (x=10); 2 = a near-zero-stiffness geogrid
// at x=10 (the CONTROL: it embeds the SAME x=10 constraint line in the mesh but does NOT split it, so it
// isolates the bare mesh-constraint effect from the interface's split+Coulomb joint).
m::Project block(int mode) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.3e4; s.nu = 0.3; s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    s.c = 5.0; s.phi = 25.0; s.Rinter = 0.4; s.rinter_rigid = false;   // weak interface (Rinter < 1)
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 10.0; L.y1 = 10.0; L.x2 = 20.0; L.y2 = 10.0; L.qx1 = 0; L.qy1 = -200.0; L.qx2 = 0; L.qy2 = -200.0;
    pr.loads.push_back(L);
    if (mode == 1) {   // vertical interface at x=10, base to crest (separates loaded/unloaded halves)
        m::StructElement e; e.kind = m::StructKind::Interface; e.name = "iface";
        e.x1 = 10.0; e.y1 = 0.0; e.x2 = 10.0; e.y2 = 10.0; e.material = -1; e.iface_material = -1;
        pr.structs.push_back(e);
    } else if (mode == 2) {   // CONTROL: tiny geogrid embeds the x=10 line without splitting it
        m::GeogridMaterial gm; gm.EA = 1.0e-3; pr.geogrids.push_back(gm);
        m::StructElement e; e.kind = m::StructKind::Geogrid; e.name = "ctrl";
        e.x1 = 10.0; e.y1 = 0.0; e.x2 = 10.0; e.y2 = 10.0; e.material = 0;
        pr.structs.push_back(e);
    }
    return pr;
}

double max_settlement(const Eigen::VectorXd& disp) {
    double s = 0.0;
    for (int i = 1; i < disp.size(); i += 2) s = std::fmax(s, std::fabs(disp[i]));
    return s;
}
// The largest settlement, and the surface settlement 0.6 m either side of the line x = 10 (the
// unloaded side first): the two points an interface there acts on.
struct Settle { double max = 0.0, unloaded = 0.0, loaded = 0.0; };
Settle solve_settlement(const m::Project& pr, int order, bool* ok, bool* finite) {
    const auto M = katai::app::mesh_from_project(pr, 1.0, order);
    const auto R = katai::app::solve_gravity_le(pr, M.mesh, InitialPhase::K0Procedure);
    *ok = R.ok;
    *finite = true;
    for (int i = 0; i < R.disp.size(); ++i) if (!std::isfinite(R.disp[i])) *finite = false;
    Settle out;
    if (!R.ok) return out;
    out.max = max_settlement(R.disp);
    const auto at = [&](double px) {
        int best = 0; double bd = 1e300;
        for (int n = 0; n < M.mesh.node_count; ++n) {
            const double dd = std::hypot(M.mesh.x[n] - px, M.mesh.y[n] - 10.0);
            if (dd < bd) { bd = dd; best = n; }
        }
        return -R.disp[2 * best + 1];
    };
    out.unloaded = at(9.4);
    out.loaded = at(10.6);
    return out;
}

void run(int order) {
    std::printf("-- element order %d --\n", order);
    bool ok_p, ok_c, ok_i, fin_p, fin_c, fin_i;
    const Settle s_plain = solve_settlement(block(0), order, &ok_p, &fin_p);
    const Settle s_ctrl  = solve_settlement(block(2), order, &ok_c, &fin_c);   // line embedded, NOT split
    const Settle s_iface = solve_settlement(block(1), order, &ok_i, &fin_i);   // line embedded + split + joint
    check(ok_p && ok_c && ok_i, "all three models solved (plain / line-control / interface)");
    check(fin_p && fin_c && fin_i, "all solutions finite (well-posed, no NaN)");
    if (!(ok_p && ok_c && ok_i)) return;

    const double line_effect = std::fabs(s_ctrl.max - s_plain.max) / s_plain.max;   // bare constraint-line effect
    std::printf("   max settlement: plain=%.5e  line-control=%.5e  interface=%.5e\n", s_plain.max, s_ctrl.max,
                s_iface.max);
    std::printf("   0.6 m from the line, unloaded / loaded side: line-control %.4e / %.4e   interface %.4e / %.4e\n",
                s_ctrl.unloaded, s_ctrl.loaded, s_iface.unloaded, s_iface.loaded);
    check(s_plain.max > 0.0 && s_iface.max > 0.0, "models settle under the surcharge");
    check(line_effect < 0.03, "the bare x=10 constraint line barely changes the result (mesh effect small)");
    // Decisive: against the SAME-meshed line-control, the interface (split + Coulomb joint) lets the loaded
    // half slide past the unloaded one. Read where that shows -- beside the line: bonded, the ground on the
    // unloaded side is dragged down with the load; with a weak joint it is left behind and the loaded side
    // goes down further. (The LARGEST settlement of the model is far from the line and moves by 0.4%; this
    // used to assert 10% on it, which held only while a model with an interface was started by another
    // procedure than the one it was compared with.) Measured: 39.5 -> 1.9 mm and 74.8 -> 112.9 mm.
    check(s_iface.unloaded < 0.5 * s_ctrl.unloaded,
          "interface (split + joint) markedly changes the result: the unloaded side is no longer dragged down");
    check(s_iface.loaded > 1.2 * s_ctrl.loaded,
          "interface (split + joint) markedly changes the result: the loaded side settles further");
}

// A plate "wall" (with +/- interface flags) must now build as an embedded wall on tri6 AND tri15, at ANY
// orientation. Toggling the interface flag must CHANGE the result -- before this work only a vertical tri6
// wall was built; tri15 or non-vertical plates were silently left bonded (interface dropped).
m::Project wall_model(bool vertical, bool with_iface) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.3e4; s.nu = 0.3; s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    s.c = 5.0; s.phi = 25.0; s.Rinter = 0.4; s.rinter_rigid = false;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    m::Load L; L.kind = m::LoadKind::Distributed;          // surcharge to the LEFT of the wall
    L.x1 = 0.0; L.y1 = 10.0; L.x2 = 8.0; L.y2 = 10.0; L.qx1 = 0; L.qy1 = -200.0; L.qx2 = 0; L.qy2 = -200.0;
    pr.loads.push_back(L);
    m::PlateMaterial pm; pm.EA = 5e6; pm.EI = 8.5e3; pr.plates.push_back(pm);
    m::StructElement e; e.kind = m::StructKind::Plate; e.name = "wall"; e.material = 0;
    if (vertical) { e.x1 = 10; e.y1 = 2; e.x2 = 10; e.y2 = 10; }   // toe embedded at y=2
    else          { e.x1 = 8;  e.y1 = 2; e.x2 = 12; e.y2 = 10; }   // inclined
    if (with_iface) { e.iface_pos = true; e.iface_neg = true; }
    pr.structs.push_back(e);
    return pr;
}
double max_abs_disp(const Eigen::VectorXd& d) {
    double m = 0; for (int i = 0; i < d.size(); ++i) m = std::fmax(m, std::fabs(d[i])); return m;
}

// A wished-in-place wall under SELF-WEIGHT ONLY (no surcharge) must NOT deflect: the geostatic seed
// (soil K0 stress + interface sigma_n0) balances gravity, so residual(0)=0 and the wall installs at
// ~zero displacement. Regression for the tri15 (5-node) interface baseline: build_problem was adding
// sigma_n0 to the staged baseline B only for the 3-node `interfaces`, not the 5-node `interfaces5`,
// so a tri15 embedded wall relaxed the geostatic lateral pressure into a spurious "installation"
// deflection before any excavation -- silently biasing displacements and N/Q/M.
m::Project wall_install_model() {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 1.3e4; s.nu = 0.3; s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    s.c = 5.0; s.phi = 25.0; s.Rinter = 0.4; s.rinter_rigid = false;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;
    m::PlateMaterial pm; pm.EA = 5e6; pm.EI = 8.5e3; pr.plates.push_back(pm);
    m::StructElement e; e.kind = m::StructKind::Plate; e.name = "wall"; e.material = 0;
    e.x1 = 10; e.y1 = 2; e.x2 = 10; e.y2 = 10;   // vertical, toe embedded at y=2
    e.iface_pos = true; e.iface_neg = true;
    pr.structs.push_back(e);
    return pr;   // NO surcharge -> any movement is the installation artifact
}
void run_wall_install(int order) {
    std::printf("-- wall install (self-weight only), order %d --\n", order);
    const auto Pr = wall_install_model();
    const auto M = katai::app::mesh_from_project(Pr, 1.0, order);
    const auto R = katai::app::solve_gravity_le(Pr, M.mesh, InitialPhase::K0Procedure);
    check(R.ok, "self-weight wall install solved");
    if (!R.ok) { std::printf("   (%s)\n", R.message.c_str()); return; }
    const double d = max_abs_disp(R.disp);
    std::printf("   max|disp| at install = %.3e m\n", d);
    // Wished-in-place: displacement is pure numerical residual of the balanced baseline, not a physical
    // deflection. Without the interfaces5 seed the tri15 wall deflects by O(mm-cm); the seed drops it to
    // solver-residual scale. 1e-4 m is far below any real wall movement yet well above the residual floor.
    check(d < 1e-4, "wished-in-place wall does not deflect under self-weight (geostatic seed balances)");
}
void run_wall(const char* label, bool vertical, int order) {
    std::printf("-- wall: %s, order %d --\n", label, order);
    const auto Pb = wall_model(vertical, false); const auto Mb = katai::app::mesh_from_project(Pb, 1.0, order);
    const auto Rb = katai::app::solve_gravity_le(Pb, Mb.mesh, InitialPhase::K0Procedure);
    const auto Pw = wall_model(vertical, true);  const auto Mw = katai::app::mesh_from_project(Pw, 1.0, order);
    const auto Rw = katai::app::solve_gravity_le(Pw, Mw.mesh, InitialPhase::K0Procedure);
    check(Rb.ok && Rw.ok, "bonded-plate and interface-wall both solved");
    if (!Rb.ok || !Rw.ok) { std::printf("   (%s / %s)\n", Rb.message.c_str(), Rw.message.c_str()); return; }
    bool fin = true; for (int i = 0; i < Rw.disp.size(); ++i) if (!std::isfinite(Rw.disp[i])) fin = false;
    check(fin, "interface-wall solution finite");
    const double db = max_abs_disp(Rb.disp), dw = max_abs_disp(Rw.disp);
    std::printf("   max|disp|: bonded-plate=%.5e  interface-wall=%.5e  (%+.1f%%)\n",
                db, dw, 100.0 * (dw - db) / db);
    // What the flag changes is what the WALL carries: bonded, the ground shears with the wall along its
    // whole length; with joints it can slide, and the bending moment is another one. (The largest
    // displacement of the model is under the surcharge, away from the wall, and moves by a fraction of a
    // percent.) Measured: the largest moment changes by 10 to 24%.
    double mb = 0.0, mw = 0.0;
    for (const auto& d : Rb.struct_forces) if (d.kind == 0) mb = std::fmax(mb, std::fabs(d.max_M));
    for (const auto& d : Rw.struct_forces) if (d.kind == 0) mw = std::fmax(mw, std::fabs(d.max_M));
    std::printf("   max|M| of the wall: bonded-plate=%.3f  interface-wall=%.3f kNm/m  (%+.1f%%)\n", mb, mw,
                mb > 0.0 ? 100.0 * (mw - mb) / mb : 0.0);
    check(mb > 0.0 && std::fabs(mw - mb) > 0.05 * mb, "interface flag CHANGES the wall result (embedded wall built)");
    // The wall must also yield an N/Q/M force diagram (structural output), including on tri15.
    bool has_M = false;
    for (const auto& d : Rw.struct_forces)
        if (d.kind == 0 && d.stations.size() >= 3 && d.max_M > 1.0) has_M = true;
    std::printf("   struct_forces: %d diagram(s)%s\n", (int)Rw.struct_forces.size(),
                has_M ? "  (N/Q/M present)" : "");
    check(has_M, "embedded wall produces an N/Q/M force diagram");
    // Interface results (interface output): tau / sigma_n / slip along the wall joint.
    bool iface_ok = false, fin_if = true; double max_tau = 0.0, max_sn = 0.0, max_slip = 0.0;
    for (const auto& ir : Rw.interface_forces) {
        if (ir.stations.size() >= 2) iface_ok = true;
        for (const auto& st : ir.stations) {
            if (!std::isfinite(st.tau) || !std::isfinite(st.sigma_n) || !std::isfinite(st.slip)) fin_if = false;
            max_tau = std::fmax(max_tau, std::fabs(st.tau));
            max_sn = std::fmax(max_sn, std::fabs(st.sigma_n));
            max_slip = std::fmax(max_slip, std::fabs(st.slip));
        }
    }
    std::printf("   interfaces: %d, max|tau|=%.3f  max|sigma_n|=%.3f kPa  max|slip|=%.3e m\n",
                (int)Rw.interface_forces.size(), max_tau, max_sn, max_slip);
    check(iface_ok, "embedded wall produces interface results (tau/sigma_n/slip stations)");
    check(fin_if, "interface results finite");
    check(max_sn > 0.0, "interface carries a normal stress (K0 lateral earth pressure)");
}

}  // namespace

int main() {
    std::printf("Standalone interface through the GUI compute path (build_problem)\n\n");
    run(6);
    run(15);
    run_wall("vertical", true, 15);    // gate removed: tri15 embedded wall
    run_wall("inclined", false, 6);    // gate removed: non-vertical embedded wall
    run_wall_install(6);               // control: tri6 wished-in-place wall installs at ~0
    run_wall_install(15);              // regression: tri15 interface baseline (interfaces5 seed)
    if (g_failures == 0) {
        std::printf("\nOK: a standalone interface is meshed, split, assembled, and changes the result "
                    "(present != absent), tri6 + tri15\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
