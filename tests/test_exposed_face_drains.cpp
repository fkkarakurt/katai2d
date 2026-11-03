// The faces an excavation exposes drain in a consolidation phase, whatever was declared on the model's
// own boundary, unless they stand against an impermeable line.
//
// A line inside the model becomes an outer boundary when the ground on one side of it is taken away. The
// flow calculations treated it so; the consolidation phase read its declared conditions on the boundary
// of the WHOLE mesh, excavated regions included. So the floor of a pit was sealed in any model with a flow
// condition anywhere, a head declared on the floor itself included, and nothing said so: a Terzaghi column
// whose top block is removed drained through the two corner nodes of its top edge alone, and reached
// U = 0.060 where the closed form is 0.115. The other way round, a model with no conditions declared
// drained the ground behind an impermeable wall straight through the wall once the pit was dug.
//
// verify: KV-CON-010
//   oracle:   closed_form
//   source:   Terzaghi, K. (1943). Theoretical Soil Mechanics. Wiley -- one-dimensional consolidation, the average degree of consolidation U(Tv) as a series; KATAI 2D input contract (docs/k2d-format.md, edge_flow, flow_barrier) for the conditions on the faces
//   locator:  a 1 x 10 m elastic column (E = 6000 kPa, nu = 0, k = 1.73e-3 m/day) with a removable block on it, single drainage at its top edge once the block is gone, 200 kPa applied at the start of a one-day and a hundred-day phase; a 20 x 10 m block with a wall 6 m deep at x = 10 and a pit 4 m deep beside it; a block with a hole in it, for the mask itself
//   quantity: the degree of consolidation, settlement over its final value; the excess pore pressure on the retained face of the wall [kPa]; the drainage mask around the hole
//   expected: U(1 day) = 0.1149 and U(100 days) = 0.9372 whether the top edge is declared a head or left closed by its default in a model with declared conditions; behind an impermeable wall the pressure stays a sizeable share of the load, behind a permeable one it is zero; every node around the hole drains, the base does not
//   band:     2% on U at one day, 1% at a hundred (the time march; the mesh is fine for the column); 10% of the load as the lower bound of the retained face's pressure behind the impermeable wall; 1e-9 kPa behind the permeable one
#include <katai/analysis/hydraulic_boundary.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kH = 10.0, kQ = 200.0, kE = 6000.0, kBlock = 50.0;
const double kPerm = 1.96e-8 / 9.8 * 86400.0 * 9.81;

double terzaghi_u(double t_day) {
    const double tv = kPerm * kE / katai::app::kGammaWater * t_day / (kH * kH);
    double s = 0.0;
    for (int j = 0; j < 400; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += (2.0 / (M * M)) * std::exp(-M * M * tv);
    }
    return 1.0 - s;
}

enum class TopEdge { Head, DeclaredClosed };

// The degree of consolidation at the end of the phase.
double column_u(TopEdge top, double t_day, int steps, std::string* why) {
    m::Project pr;
    m::Material clay; clay.name = "Clay"; clay.model = m::SoilModel::LinearElastic;
    clay.E = kE; clay.nu = 0.0; clay.gamma_unsat = clay.gamma_sat = 0.0; clay.kx = clay.ky = kPerm;
    pr.materials.push_back(clay);
    m::Material dead; dead.name = "Block"; dead.model = m::SoilModel::LinearElastic;
    dead.E = 1e5; dead.nu = 0.3; dead.gamma_unsat = dead.gamma_sat = kBlock; dead.kx = dead.ky = kPerm;
    pr.materials.push_back(dead);
    const int F = (int)m::BCType::Free, H = (int)m::BCType::HorizontallyFixed, V = (int)m::BCType::VerticallyFixed;
    m::SoilPolygon layer; layer.name = "Layer"; layer.material = 0;
    layer.x = {0, 1, 1, 0}; layer.y = {0, 0, kH, kH}; layer.edge_bc = {V, H, F, H};
    const int C = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
    m::SoilPolygon block; block.name = "Block"; block.material = 1;
    block.x = {0, 1, 1, 0}; block.y = {kH, kH, kH + 1, kH + 1}; block.edge_bc = {F, H, F, H};
    if (top == TopEdge::Head) {
        layer.edge_flow = {C, C, Hd, C};
        layer.edge_head = {0.0, 0.0, kH, 0.0};
    } else {
        // A head stated on the block's top makes this a model with declared conditions; the floor is
        // closed by its edge's default. Taken away, the block's top says nothing, and the floor is an
        // outer boundary now: it drains. What seals a face is an impermeable line.
        layer.edge_flow = {C, C, C, C};
        layer.edge_head = {0.0, 0.0, 0.0, 0.0};
        block.edge_flow = {C, C, Hd, C};
        block.edge_head = {0.0, 0.0, kH + 1.0, 0.0};
    }
    pr.polygons = {layer, block};
    m::Load swap; swap.kind = m::LoadKind::Distributed; swap.name = "In place of the block";
    swap.x1 = 0; swap.y1 = kH; swap.x2 = 1; swap.y2 = kH; swap.qy1 = swap.qy2 = -kBlock;
    m::Load q = swap; q.name = "Load"; q.qy1 = q.qy2 = -kQ;
    pr.loads = {swap, q};
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;
    pr.initial.load_active = {0, 0};
    m::Phase dig; dig.name = "Block off"; dig.type = m::PhaseType::Plastic;
    dig.poly_active = {1, 0}; dig.load_active = {1, 0};
    m::Phase wait = dig; wait.name = "Consolidation"; wait.type = m::PhaseType::Consolidation;
    wait.load_active = {1, 1}; wait.duration = t_day; wait.time_steps = steps; wait.consol_load_at_start = true;
    pr.phases = {dig, wait};
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { *why = "mesh"; return -1.0; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 3 || !res.back().ok || res.back().consol_settlement.empty()) {
        *why = res.empty() ? "no result" : res.back().message;
        return -1.0;
    }
    return res.back().consol_settlement.back() / (kQ * kH / kE);
}

// The largest excess pressure on the retained face of a wall, a pit dug on its other side, a load behind it.
double retained_face_pressure(int barrier, bool declared, std::string* why) {
    constexpr double kLoad = 50.0;
    m::Project pr;
    m::Material s; s.name = "Soil"; s.model = m::SoilModel::LinearElastic;
    s.E = 5000.0; s.nu = 0.2; s.gamma_unsat = s.gamma_sat = 0.0; s.kx = s.ky = 1e-3;
    s.rinter_rigid = true;
    pr.materials.push_back(s);
    const int F = (int)m::BCType::Free, H = (int)m::BCType::HorizontallyFixed, X = (int)m::BCType::FullyFixed;
    m::SoilPolygon behind; behind.name = "Behind"; behind.material = 0;
    behind.x = {0, 10, 10, 10, 0}; behind.y = {0, 0, 6, 10, 10}; behind.edge_bc = {X, F, F, F, H};
    m::SoilPolygon pit; pit.name = "Pit"; pit.material = 0;
    pit.x = {10, 20, 20, 10}; pit.y = {6, 6, 10, 10}; pit.edge_bc = {F, H, F, F};
    m::SoilPolygon below; below.name = "Below"; below.material = 0;
    below.x = {10, 20, 20, 10}; below.y = {0, 0, 6, 6}; below.edge_bc = {X, H, F, F};
    if (declared) {
        const int C = (int)m::FlowBCType::Closed, Hd = (int)m::FlowBCType::Head;
        behind.edge_flow = {C, C, C, Hd, C};
        behind.edge_head = {0.0, 0.0, 0.0, 10.0, 0.0};
    }
    pr.polygons = {behind, pit, below};
    m::PlateMaterial pm; pm.name = "Wall"; pm.EA = 1.0e7; pm.EI = 1.0e5; pm.w = 0.0;
    pr.plates.push_back(pm);
    m::StructElement w; w.kind = m::StructKind::Plate; w.name = "Wall";
    w.x1 = 10; w.y1 = 10; w.x2 = 10; w.y2 = 4; w.material = 0;
    w.iface_pos = true; w.iface_neg = true; w.flow_barrier = barrier;
    pr.structs.push_back(w);
    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Load";
    L.x1 = 4; L.y1 = 10; L.x2 = 10; L.y2 = 10; L.qy1 = L.qy2 = -kLoad;
    pr.loads.push_back(L);
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;
    pr.initial.struct_active = {0};
    pr.initial.load_active = {0};
    m::Phase dig; dig.name = "Dig"; dig.type = m::PhaseType::Plastic;
    dig.struct_active = {1}; dig.poly_active = {1, 0, 1}; dig.load_active = {0};
    m::Phase wait = dig; wait.name = "Wait"; wait.type = m::PhaseType::Consolidation;
    wait.load_active = {1}; wait.duration = 2.0; wait.time_steps = 10; wait.consol_load_at_start = true;
    pr.phases = {dig, wait};
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) { *why = "mesh"; return -1.0; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 3 || !res.back().ok) {
        *why = res.empty() ? "no result" : res.back().message;
        return -1.0;
    }
    const auto& r = res.back();
    double most = 0.0;
    for (int n = 0; n < r.mesh.node_count && n < (int)r.excess_pore.size(); ++n)
        if (std::fabs(r.mesh.x[n] - 10.0) < 1e-6 && r.mesh.y[n] > 6.5 && r.mesh.y[n] < 9.5)
            most = std::max(most, std::fabs(r.excess_pore[n]));
    return most;
}

}  // namespace

int main() {
    std::printf("== the faces an excavation exposes (KV-CON-010) ==\n");
    const double u1 = terzaghi_u(1.0), u100 = terzaghi_u(100.0);
    std::printf("      closed form: U(1 d) = %.4f, U(100 d) = %.4f\n", u1, u100);
    std::fflush(stdout);
    const struct { TopEdge top; const char* name; } cases[] = {
        {TopEdge::Head, "a head declared on the floor"},
        {TopEdge::DeclaredClosed, "the floor closed by default in a model with declared conditions"},
    };
    for (const auto& c : cases) {
        std::string why;
        const double a = column_u(c.top, 1.0, 40, &why);
        const double b = column_u(c.top, 100.0, 40, &why);
        std::printf("      %s: U(1 d) = %.4f, U(100 d) = %.4f %s\n", c.name, a, b, why.c_str());
        std::fflush(stdout);
        check(a > 0.0 && std::fabs(a - u1) <= 0.02 * u1, std::string(c.name) + ": the column drains through its top at one day");
        check(b > 0.0 && std::fabs(b - u100) <= 0.01 * u100, std::string(c.name) + ": and at a hundred days");
    }

    for (const bool declared : {false, true}) {
        std::string why;
        const double closed = retained_face_pressure(1, declared, &why);
        const double open = retained_face_pressure(0, declared, &why);
        std::printf("      %s: retained face behind an impermeable wall %.3f kPa, behind a permeable one %.3e kPa %s\n",
                    declared ? "declared conditions" : "nothing declared", closed, open, why.c_str());
        std::fflush(stdout);
        check(closed > 0.1 * 50.0, "the impermeable wall holds the water behind it");
        check(open >= 0.0 && open < 1e-9, "the permeable one lets it into the pit");
    }

    // The mask itself: a block with a hole. Every face of the hole drains, its crown as well as its floor;
    // the model's base does not.
    {
        const katai::geometry::RectangularDomain domain{0.0, 0.0, 8.0, 8.0, 0};
        const katai::mesh::Mesh mesh = katai::mesh::generate_structured_tri6(domain, 8, 8);
        std::vector<char> act(mesh.element_count, 1);
        for (int e = 0; e < mesh.element_count; ++e) {
            double cx = 0.0, cy = 0.0;
            for (int k = 0; k < 3; ++k) { cx += mesh.x[mesh.node_of(e, k)] / 3.0; cy += mesh.y[mesh.node_of(e, k)] / 3.0; }
            if (cx > 3.0 && cx < 5.0 && cy > 3.0 && cy < 5.0) act[e] = 0;
        }
        for (const bool declared : {false, true}) {
            std::vector<katai::core::FlowEdge> edges;
            if (declared) {
                katai::core::FlowEdge top; top.ax = 0; top.ay = 8; top.bx = 8; top.by = 8;
                top.kind = katai::core::FlowEdgeKind::Head; top.head = 8.0;
                edges.push_back(top);
            }
            const auto drained = katai::core::flow_drained_nodes(edges, declared, mesh, act, 8.0);
            bool hole = true, base = false;
            for (int n = 0; n < mesh.node_count; ++n) {
                const double x = mesh.x[n], y = mesh.y[n];
                const bool on_ring = x >= 3.0 - 1e-9 && x <= 5.0 + 1e-9 && y >= 3.0 - 1e-9 && y <= 5.0 + 1e-9 &&
                                     (std::fabs(x - 3.0) < 1e-9 || std::fabs(x - 5.0) < 1e-9 ||
                                      std::fabs(y - 3.0) < 1e-9 || std::fabs(y - 5.0) < 1e-9);
                if (on_ring) hole = hole && drained[n];
                if (y < 1e-9 && x > 1e-9 && x < 8.0 - 1e-9) base = base || drained[n];
            }
            check(hole, std::string(declared ? "declared: " : "nothing declared: ") + "every node around the hole drains");
            check(!base, std::string(declared ? "declared: " : "nothing declared: ") + "the base does not");
        }
    }

    if (g_failures == 0) {
        std::printf("OK: the faces an excavation exposes verified\n");
        return 0;
    }
    std::printf("%d check(s) failed\n", g_failures);
    return 1;
}
