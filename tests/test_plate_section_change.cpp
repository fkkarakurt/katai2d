// A plate that changes its section in a phase keeps the force it carries; the new section takes
// only what comes after.
//
// A plate lies along a weightless strip of ground a million times softer than itself, held at both
// ends, so it cannot shrink: a line contraction C (half of it a shrinkage strain) puts it in
// tension N = EA C / 2 and nothing else carries anything. A later phase gives it a stiffer
// section and raises the contraction from C1 to C2. Kept, the force the plate already carries is
// the one the first section built, and the new section adds its own share of the change,
//     N = EA1 C1 / 2 + EA2 (C2 - C1) / 2;
// recomputed with the new section over the whole history it would be EA2 C2 / 2. A third phase
// takes the contraction back, and it is taken back through the section that is there: the force
// falls by EA2 C2 / 2 and ends in compression, because part of what it undoes was carried by the
// softer section -- the force is a state, the section a property of each increment.
//
// verify: KV-STR-020
//   oracle:   closed_form
//   source:   a bar restrained at both ends under an imposed free strain eps0 carries N = -EA eps0, and a change of section between two increments of eps0 changes only the stiffness of the later increment (incremental elasticity: the force is a state, the section a property of the increment); KATAI 2D input contract, docs/k2d-format.md, phases[i].smat and scontr
//   locator:  a 10 x 1 m strip of ground with E = 1 kPa, weightless, fully fixed at its two ends, a plate along its top; phase 2 contracts the plate by 1%, phase 3 gives it twice the axial stiffness and contracts it by 2%, phase 4 takes the contraction back
//   quantity: mean axial force of the plate [kN/m]
//   expected: 5000 kN/m after phase 2, 15000 kN/m after phase 3 (not 20000), -5000 kN/m after phase 4 (the 2e6 section undoing 2% from 15000)
//   band:     0.1% (the ground carries the millionth of the force it is that much softer)
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}
double mean_N(const katai::app::SolveResult& r) {
    double s = 0.0; int n = 0;
    for (const auto& f : r.struct_forces)
        for (const auto& q : f.stations) { s += q.N; ++n; }
    return n ? s / n : 0.0;
}
}  // namespace

int main() {
    constexpr double L = 10.0, EA1 = 1.0e6, EA2 = 2.0e6, C1 = 0.01, C2 = 0.02;
    m::Project pr;
    m::Material s; s.name = "Soft"; s.model = m::SoilModel::LinearElastic;
    s.E = 1.0; s.nu = 0.0; s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.name = "Strip"; P.material = 0;
    P.x = {0, L, L, 0}; P.y = {0, 0, 1, 1};
    P.edge_bc = {(int)m::BCType::Free, (int)m::BCType::FullyFixed, (int)m::BCType::Free,
                 (int)m::BCType::FullyFixed};
    pr.polygons.push_back(P);
    pr.x_min = 0; pr.x_max = L; pr.y_min = 0; pr.y_max = 1;
    pr.has_water = false;
    m::PlateMaterial p1; p1.name = "Section 1"; p1.EA = EA1; p1.EI = 1.0e3; p1.nu = 0.0; p1.w = 0.0;
    m::PlateMaterial p2 = p1; p2.name = "Section 2"; p2.EA = EA2;
    pr.plates = {p1, p2};
    m::StructElement b; b.kind = m::StructKind::Plate; b.name = "Bar"; b.material = 0;
    b.x1 = 0.0; b.y1 = 1.0; b.x2 = L; b.y2 = 1.0;
    pr.structs.push_back(b);
    pr.initial.struct_active = {0};
    m::Phase install; install.name = "Install"; install.struct_active = {1};
    m::Phase contract = install; contract.name = "Contract"; contract.struct_contraction = {C1};
    m::Phase stiffer = contract; stiffer.name = "Stiffer and more"; stiffer.struct_material = {1};
    stiffer.struct_contraction = {C2};
    m::Phase undo = stiffer; undo.name = "Undo"; undo.struct_contraction = {0.0};
    pr.phases = {install, contract, stiffer, undo};
    pr.mesh.auto_refine = false;

    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    check(M.ok, "strip meshed");
    if (!M.ok) return 1;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 5 && res[4].ok, "the four phases ran");
    if (res.size() != 5 || !res[4].ok) {
        if (!res.empty()) std::printf("  (%s)\n", res.back().message.c_str());
        return 1;
    }
    const double n2 = mean_N(res[2]), n3 = mean_N(res[3]), n4 = mean_N(res[4]);
    const double w2 = 0.5 * EA1 * C1, w3 = w2 + 0.5 * EA2 * (C2 - C1), w4 = w3 - 0.5 * EA2 * C2;
    std::printf("  after the contraction  N = %.3f kN/m (closed form %.3f)\n", n2, w2);
    std::printf("  after the new section  N = %.3f kN/m (kept: %.3f; recomputed with the new section: %.3f)\n",
                n3, w3, 0.5 * EA2 * C2);
    std::printf("  after taking it back   N = %.3f kN/m (closed form %.3f)\n", n4, w4);
    check(std::fabs(n2 - w2) < 1e-3 * w2, "the first section carries EA1 C1 / 2");
    check(std::fabs(n3 - w3) < 1e-3 * w3, "the new section keeps that force and adds EA2 (C2 - C1) / 2");
    check(std::fabs(n4 - w4) < 1e-3 * std::fabs(w3), "taking the contraction back undoes it through the section there (into compression)");
    if (g_failures) { std::fprintf(stderr, "%d check(s) failed\n", g_failures); return 1; }
    std::printf("all plate section change checks passed\n");
    return 0;
}
