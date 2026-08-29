// Structural elements in a FULLY-COUPLED flow-deformation phase -- the last of the coupled
// families to carry them. The consolidation phase gained plates, anchors and geogrids
// (KV-STR-006) and then interfaces (KV-STR-007); this one refused all of it, which left the
// program saying that a wall may stand in the ground while the pore pressure dissipates but not
// while the water table moves. That is not a distinction the physics makes, and it is not one the
// two solvers make either: they take the same structural stiffness matrix and the same tied seam
// pore equations, and differ only in what this phase adds on its own -- the retention curve and
// Bishop's chi, which the structure neither sees nor is seen by.
//
// THE ORACLE IS THAT SAMENESS, MEASURED. The fully-coupled core reduces to the consolidation core
// in the SATURATED limit (its own header says so, and test_coupled_flow pins it for soil). So the
// same model, with the same wall, run once as a Consolidation phase and once as a FullyCoupled
// one must give the same answer -- settlement AND the structure's own force -- while a third run
// as a drained Plastic phase gives the state both of them end in. Two independent paths and a
// limit, none of which share code beyond the element matrices.
//
// verify: KV-STR-008
//   oracle:   independent_path
//   source:   the saturated limit of the fully-coupled (unsaturated Biot + Bishop) formulation is the consolidation (saturated Biot) formulation -- kernel/analysis/coupled_flow_deformation.hpp states the reduction and test_coupled_flow measures it for soil alone; here the same reduction is asked for with a structure in the system. The drained limit is the one KV-STR-006 uses: once the excess pore pressure is gone, the coupled problem IS the drained problem, which solve_nonlinear solves by a different path
//   locator:  a phase whose excess pore pressure stays non-negative everywhere is fully saturated, so S_eff = 1, Bishop's chi = 1, and the fully-coupled system is term for term the consolidation system. Tv = cv t / H_dr^2 with cv = k Eoed / gamma_w; at Tv = 4 the remaining excess pore pressure is 2e-5 of what the load generated
//   quantity: the maximum settlement [m] and the raft's peak bending moment [kNm/m] at the end of a fully-coupled phase, against the same model as a consolidation phase and as a drained Plastic phase
//   expected: fully-coupled = consolidation (the same system in the saturated limit), and both = the drained Plastic phase (the state the dissipation ends in)
//   band:     5e-4 relative on the saturated-limit pair, measured +0.00000% on the settlement (4.50097547e-02 m from both) and no difference at all in the raft moment (21.621032 kNm/m from both, and from the drained reference too). That is what a REDUCTION should look like: the two solvers are running the same system, so what separates them is their own iteration -- a Picard fixed point against a Newton one -- and here it separates them by nothing measurable. 5e-3 on the drained limit, measured -0.00352%, whose residue is the excess pore pressure still in the ground at Tv = 4 (2e-5 of the load), the same axis KV-STR-006 measures at -0.0025% on the same block without the coupled-flow retention in the way
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
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

constexpr double kW = 20.0, kH = 10.0;
constexpr double kE = 5000.0, kNu = 0.2;
constexpr double kPerm = 1.0e-3;
constexpr double kQ = 50.0;
constexpr double kRaftX1 = 8.0, kRaftX2 = 12.0;

double eoed() { return kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu)); }
double cv() { return kPerm * eoed() / katai::app::kGammaWater; }

// 0 = consolidation, 1 = fully coupled, 2 = drained Plastic.
m::Project block(int kind) {
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = kE; s.nu = kNu; s.gamma_unsat = 17.0; s.gamma_sat = 20.0; s.e_init = 0.7;
    s.kx = kPerm; s.ky = kPerm;
    pr.materials.push_back(s);

    m::SoilPolygon P; P.material = 0;
    P.x = {0, kW, kW, 0};
    P.y = {0, 0, kH, kH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};
    P.edge_head = {0.0, 0.0, kH, 0.0};
    pr.polygons.push_back(P);
    pr.has_water = false;

    m::PlateMaterial pm; pm.EA = 5.0e6; pm.EI = 8.5e3;
    pr.plates.push_back(pm);
    m::StructElement e; e.kind = m::StructKind::Plate; e.name = "Raft";
    e.x1 = kRaftX1; e.y1 = kH; e.x2 = kRaftX2; e.y2 = kH; e.material = 0;
    pr.structs.push_back(e);

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = kRaftX1; L.y1 = kH; L.x2 = kRaftX2; L.y2 = kH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};

    m::Phase ph; ph.load_active = {1};
    if (kind == 2) { ph.name = "Drained"; ph.type = m::PhaseType::Plastic; }
    else {
        ph.name = kind == 0 ? "Consolidation" : "Fully coupled";
        ph.type = kind == 0 ? m::PhaseType::Consolidation : m::PhaseType::FullyCoupled;
        ph.duration = 4.0 * kH * kH / cv();
        ph.time_steps = 60;
    }
    pr.phases.push_back(ph);
    return pr;
}

struct Answer { bool ok = false; double settle = 0.0, moment = 0.0; std::string msg; };

Answer run(const m::Project& pr) {
    Answer a;
    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    if (!M.ok) { a.msg = M.message; return a; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2) { a.msg = "phases did not run"; return a; }
    const auto& R = res[1];
    a.ok = R.ok; a.msg = R.message;
    if (!R.ok) return a;
    for (int n = 0; n < R.mesh.node_count; ++n)
        a.settle = std::fmax(a.settle, std::fabs(R.disp(2 * n + 1)));
    if (!R.struct_forces.empty()) a.moment = R.struct_forces.front().max_M;
    return a;
}

}  // namespace

int main() {
    std::printf("Structural elements in a fully-coupled flow-deformation phase\n\n");
    const Answer fc = run(block(1));
    check(fc.ok, "a fully-coupled phase with a raft solves (it used to be refused)");
    if (!fc.ok) { std::printf("   (%s)\n", fc.msg.c_str()); }

    const Answer co = run(block(0));
    const Answer dr = run(block(2));
    check(co.ok && dr.ok, "the consolidation and drained references solved");
    if (fc.ok && co.ok && dr.ok) {
        std::printf("   settlement:  fully-coupled %.8e  consolidation %.8e  (%+.5f%%)\n",
                    fc.settle, co.settle, 100.0 * (fc.settle / co.settle - 1.0));
        std::printf("   raft |M|max: fully-coupled %.6f  consolidation %.6f  drained %.6f kNm/m\n",
                    fc.moment, co.moment, dr.moment);
        std::printf("   against the drained state:  settlement %+.5f%%,  moment %+.5f%%\n",
                    100.0 * (fc.settle / dr.settle - 1.0), 100.0 * (fc.moment / dr.moment - 1.0));
        check(dr.moment > 1.0, "the reference really bends the raft");
        // The saturated limit: the same system, so the same answer to the tolerance of the two
        // solvers' own iteration (a Picard fixed point against a Newton one).
        check(std::fabs(fc.settle / co.settle - 1.0) < 5e-4,
              "fully coupled = consolidation in the saturated limit (settlement)");
        check(std::fabs(fc.moment / co.moment - 1.0) < 5e-4,
              "and the structure's own force says the same");
        // And both end where the drained problem is.
        check(std::fabs(fc.settle / dr.settle - 1.0) < 5e-3,   // measured -0.00352%
              "and it walks into the drained answer, as the consolidation phase does");
    }

    if (g_failures == 0) {
        std::printf("\nOK: the last coupled family carries the structure too\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
