// A PHASE THAT IS NOT A PLASTIC PHASE STILL CONTINUES THE STRUCTURES AND HANDS THEM ON.
//
// The phase chain takes a phase's structural state from the result just before it. A Safety phase
// handed its parent's on; nothing else that is not a Plastic phase did. After a consolidation, a
// fully coupled flow-deformation, a transient flow or a dynamic phase, the next phase rebuilt every
// structural force from zero against a ground still in equilibrium with them (a drawdown -- a
// transient flow and the Plastic phase that takes its pore pressures -- moved a wall in that way).
//
// A consolidation or a fully coupled flow-deformation phase solves its plates, anchors, geogrids
// and joints elastically, as a stiffness on top of the state the phase before left them in. Two
// things followed from how that state was treated, and both are checked here against the
// definition of a phase that changes nothing:
//
// (1) What the phase REPORTED. The joints were read on the elastic branch alone, k_n du_n and
//     k_s du_s: without the normal stress the joint was installed with and without the slip it had
//     already made, so a wall at rest against its ground showed a joint stress of a few kPa where
//     it carries tens, and the Coulomb check (K2D-A016) measured the shear against a normal
//     stress that was not there. A plate whose state carries a line contraction, a section change
//     or a hinge, and an anchor or a geogrid that had yielded, were read without that state too.
// (2) What the phase HANDED ON. It published no structural state, so the phase after it rebuilt
//     every structural force from zero against a ground that had been in equilibrium with them --
//     the run said so in a note, and the wall moved.
//
// Oracle: a consolidation (or fully coupled) phase on drained ground with nothing changed applies
// nothing; its reported structural forces and joint stresses are its parent's, and the phase after
// it -- identical again -- does not move and reports the same.
//
// verify: KV-CON-007
//   oracle:   closed_form
//   source:   the definition of a nil phase (a phase whose configuration equals its parent's applies nothing and changes nothing), applied across a consolidation and a fully coupled flow-deformation phase; KATAI 2D input contract (docs/k2d-format.md, phases[i].type): a coupled phase continues the structures of the phase before it
//   locator:  a 20 x 10 m drained linear elastic ground with an embedded wall (plate and interfaces on both faces) loaded by a one-sided 50 kPa surcharge in a Plastic phase; then a consolidation (and, separately, a fully coupled) phase of 10 days with nothing changed; then a Plastic phase identical to it
//   quantity: the coupled phase's and the following phase's wall moment [kNm/m] and joint normal and shear stress [kPa] against the surcharge phase's, and the following phase's displacement [m]
//   expected: equal, station by station; no displacement
//   band:     1e-6 of the parent's values (round-off; nothing is applied)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
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

struct Diff { double dM = 0.0, M = 0.0, dJ = 0.0, J = 0.0; };
Diff diff(const katai::app::SolveResult& a, const katai::app::SolveResult& b) {
    Diff d;
    for (size_t i = 0; i < a.struct_forces.size() && i < b.struct_forces.size(); ++i)
        for (size_t k = 0; k < a.struct_forces[i].stations.size() &&
                           k < b.struct_forces[i].stations.size(); ++k) {
            d.dM = std::fmax(d.dM, std::fabs(a.struct_forces[i].stations[k].M -
                                             b.struct_forces[i].stations[k].M));
            d.M = std::fmax(d.M, std::fabs(b.struct_forces[i].stations[k].M));
        }
    for (size_t i = 0; i < a.interface_forces.size() && i < b.interface_forces.size(); ++i)
        for (size_t k = 0; k < a.interface_forces[i].stations.size() &&
                           k < b.interface_forces[i].stations.size(); ++k) {
            const auto& x = a.interface_forces[i].stations[k];
            const auto& y = b.interface_forces[i].stations[k];
            d.dJ = std::fmax(d.dJ, std::fmax(std::fabs(x.sigma_n - y.sigma_n), std::fabs(x.tau - y.tau)));
            d.J = std::fmax(d.J, std::fabs(y.sigma_n));
        }
    return d;
}

void run(m::PhaseType coupled, const char* label, bool reports_structures) {
    std::printf("-- %s --\n", label);
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
    s.c = 5.0; s.phi = 30.0; s.rinter_rigid = false; s.Rinter = 1.0;
    s.kx = 1e-3; s.ky = 1e-3;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    if (coupled == m::PhaseType::TransientFlow) {
        // A flow needs a head to flow from; at the water table's own level nothing flows.
        pr.has_water = true; pr.wx = {0, 20}; pr.wy = {5, 5};
        P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head,
                       (int)m::FlowBCType::Closed, (int)m::FlowBCType::Head};
        P.edge_head = {0.0, 5.0, 0.0, 5.0};
    }
    pr.polygons.push_back(P);
    m::PlateMaterial pm;
    pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = 0.0;
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
    m::Phase p2 = p1; p2.type = coupled; p2.name = "Between"; p2.duration = 10.0;
    if (coupled == m::PhaseType::Dynamic) { p2.seismic_amp = 0.01; p2.duration = 0.2; }
    m::Phase p3 = p1; p3.name = "Nil";
    pr.phases = {p1, p2, p3};

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    if (!M.ok) { check(false, "model meshed"); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 4 || !res[1].ok || !res[2].ok || !res[3].ok) {
        check(false, "K0 + Surcharge + coupled + Nil all solved");
        for (const auto& r : res) std::printf("  (%s)\n", r.message.c_str());
        return;
    }
    const Diff c = diff(res[2], res[1]), n = diff(res[3], res[1]);
    std::printf("   coupled: |dM| %.3g of %.4g kNm/m, joint |d| %.3g of %.4g kPa\n", c.dM, c.M, c.dJ, c.J);
    std::printf("   after:   |dM| %.3g, joint |d| %.3g kPa, |u| %.3g m (surcharge %.4g m)\n", n.dM, n.dJ,
                res[3].max_disp, res[1].max_disp);
    check(n.M > 1.0 && n.J > 10.0, "the surcharge bends the wall and loads its joints (teeth)");
    if (reports_structures) {
        check(c.dM < 1e-6 * c.M, "the coupled phase reports the parent's wall moment");
        check(c.dJ < 1e-6 * c.J, "the coupled phase reports the parent's joint stresses (installed stress and slip included)");
    }
    check(res[3].max_disp < 1e-6 * res[1].max_disp, "the phase after it does not move");
    check(n.dM < 1e-6 * n.M && n.dJ < 1e-6 * n.J, "the phase after it continues the structures");
    check(res[3].message.find("supplies no structural state") == std::string::npos,
          "the phase after it is handed a structural state");
}
}  // namespace

int main() {
    std::printf("A coupled phase continues the structures and hands them on\n\n");
    run(m::PhaseType::Consolidation, "consolidation", true);
    run(m::PhaseType::FullyCoupled, "fully coupled flow-deformation", true);
    run(m::PhaseType::TransientFlow, "transient groundwater flow", false);
    run(m::PhaseType::Dynamic, "dynamic", false);
    if (g_failures == 0) {
        std::printf("\nOK: a coupled phase reports and hands on the structures' state\n");
        return 0;
    }
    std::printf("\n%d check(s) failed\n", g_failures);
    return 1;
}
