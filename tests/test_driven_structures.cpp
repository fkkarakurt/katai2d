// A structure standing on a node that a prescribed displacement drives, read back from the
// structural force REPORT -- the numbers a user actually takes away.
//
// The analysis has carried a driven node's motion into every structural element since
// 2026-08-13 (internal_forces.hpp, u_at). The reports did not all follow. The plate, geogrid,
// embedded-beam and interface post-processors read the full displacement vector, whose fixed
// entries hold the prescribed values; the anchor force skipped fixed DOFs, and so an anchor with
// an end on a driven line reported N = 0 while the solve carried EA d / L -- under a warning,
// K2D-A003, that said the opposite: that the STRUCTURE did not see the motion. Both are closed
// now, and this case is the evidence, because each of its three problems has a closed form that
// a report ignoring the driven motion cannot meet:
//
//   (a) a weightless two-span continuous plate on pin / roller end supports whose MIDDLE support
//       settles by delta. Statically indeterminate, so the settlement bends it -- a report that
//       read the middle node as standing still would show M = 0. Force method: release the middle
//       support, the simply supported beam of span 2L deflects d11 under a unit mid-span force,
//       and compatibility with the settlement fixes the redundant R = delta / d11, whence
//       M_B = R L / 2. The plate is shear-deformable, so d11 carries both terms.
//   (b) a strut between the two side edges of a weightless elastic block, one edge held at
//       u_x = 0 and the other driven by u_x = D: both ends are driven, so N = EA D / L exactly.
//   (c) a fixed-end anchor from the driven edge to a fixed point outside the block:
//       N = EA (-D . dir) / L_c, i.e. the driven end moving towards the fixed point shortens it.
//
// Before the report was aligned with the solver, (a) already read 8.98816 kNm/m and (b), (c)
// read 0.0 and 0.0.
//
// verify: KV-STR-009
//   oracle:   closed_form
//   source:   the method of consistent deformations (force method) for a statically indeterminate beam with a support settlement, and the linear axial spring of an anchor, both stated in full in the locator; the plate's shear rigidity kGA' = k EA / (2 (1 + nu)), k = 5/6, is PLAXIS 2D Material Models Manual (2025.1) Eq. 18-8, the same expression KV-STR-003 verifies its deflections with
//   locator:  (a) d11 = (2L)^3 / (48 EI) + (2L) / (4 kGA') = L^3 / (6 EI) + L / (2 kGA'), R = delta / d11, M_B = R L / 2, with L = 2 m, EI = 1200 kN m2/m, EA = 1.64e6 kN/m, nu = 0, delta = 0.01 m; (b) N = EA D / L with EA = 1e5 kN, L = 10 m, D = -0.01 m; (c) N = EA (-D) / L_c with L_c = 5 m, the fixed point on the +x side of the driven edge
//   quantity: (a) the peak bending moment of the plate, at the settling support [kNm/m]; (b), (c) the reported axial force of each anchor [kN]
//   expected: (a) M_B = 8.98816193306... kNm/m (Euler-Bernoulli alone would give 3 EI delta / L^2 = 9.0); (b) N = -100 kN (compression); (c) N = +200 kN (tension)
//   band:     1e-9 relative on all three, as asserted below. (a) is exact in the element space -- the moment is linear on each span and the plate element's curvature is linear -- measured 8.988161933076 against 8.988161933064 (1.4e-12). (b) and (c) depend only on the prescribed values at the attached nodes, so the only error is round-off -- measured -100.000000000 and +200.000000000 kN

#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

double rel(double a, double b) { return std::fabs(a - b) / std::fabs(b); }

m::StructElement line(m::StructKind kind, const char* name, double x1, double y1, double x2,
                      double y2) {
    m::StructElement s;
    s.kind = kind; s.name = name; s.material = 0;
    s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2;
    return s;
}

m::PrescribedDisp driven(const char* name, double x1, double y1, double x2, double y2,
                         bool set_ux, double ux, bool set_uy, double uy) {
    m::PrescribedDisp D;
    D.name = name;
    D.x1 = x1; D.y1 = y1; D.x2 = x2; D.y2 = y2;
    D.set_ux = set_ux; D.ux = ux;
    D.set_uy = set_uy; D.uy = uy;
    return D;
}

// The initial phase establishes nothing but the elements; the named phase applies the motion.
std::vector<katai::app::SolveResult> run(m::Project pr, const std::vector<char>& poly) {
    pr.initial.poly_active = poly;
    pr.initial.struct_active.assign(pr.structs.size(), 1);
    pr.initial.disp_active.assign(pr.disps.size(), 0);
    m::Phase ph;
    ph.name = "Drive";
    ph.poly_active = poly;
    ph.struct_active.assign(pr.structs.size(), 1);
    ph.disp_active.assign(pr.disps.size(), 1);
    pr.phases = {ph};
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, pr.name + ": meshed");
    if (!M.ok) return {};
    auto res = katai::app::solve_phases(pr, M.mesh, katai::app::InitialPhase::K0Procedure);
    for (const auto& r : res)
        if (!r.ok) std::printf("      (%s)\n", r.message.c_str());
    check(res.size() == 2 && res[0].ok && res[1].ok, pr.name + ": both phases ran");
    return res;
}

const katai::core::StructForce* force_of(const katai::app::SolveResult& r, const char* name) {
    for (const auto& f : r.struct_forces)
        if (f.name == name) return &f;
    return nullptr;
}

void plate_on_a_settling_support() {
    constexpr double L = 2.0, EI = 1200.0, EA = 1.64e6, nu = 0.0, delta = 0.01;
    m::Project pr;
    pr.name = "two-span plate, middle support settles";
    pr.mesh.elem_size = 0.25;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;
    m::Material soil;             // present only to carry the mesh; never active
    pr.materials.push_back(soil);
    m::SoilPolygon P;
    P.name = "Block"; P.material = 0;
    P.x = {0.0, 2.0 * L, 2.0 * L, 0.0};
    P.y = {0.0, 0.0, 1.0, 1.0};
    P.edge_bc = {(int)m::BCType::Free,             // bottom: the beam lies along it
                 (int)m::BCType::VerticallyFixed,  // right: roller under the beam's right end
                 (int)m::BCType::Free,
                 (int)m::BCType::FullyFixed};      // left: pin under the beam's left end
    pr.polygons.push_back(P);
    m::PlateMaterial pm;
    pm.name = "HEB 200"; pm.EA = EA; pm.EI = EI; pm.nu = nu; pm.w = 0.0;
    pr.plates.push_back(pm);
    pr.structs.push_back(line(m::StructKind::Plate, "Beam", 0.0, 0.0, 2.0 * L, 0.0));
    pr.disps.push_back(driven("Middle support", L, 0.0, L, 1.0, false, 0.0, true, -delta));

    const auto res = run(pr, {0});
    if (res.size() != 2) return;
    const auto* f = force_of(res[1], "Beam");
    check(f != nullptr, "the plate reports a force diagram");
    if (!f) return;
    const double kGA = (5.0 / 6.0) * EA / (2.0 * (1.0 + nu));
    const double d11 = L * L * L / (6.0 * EI) + L / (2.0 * kGA);
    const double MB = (delta / d11) * L / 2.0;
    std::printf("   (a) plate: max|M| = %.12f kNm/m, closed form M_B = %.12f (Euler-Bernoulli %.6f), "
                "deviation %.2e\n", f->max_M, MB, 3.0 * EI * delta / (L * L), rel(f->max_M, MB));
    check(rel(f->max_M, MB) < 1e-9,
          "the settling middle support bends the plate by the force method's M_B, within 1e-9");
}

void anchors_on_driven_edges() {
    constexpr double EA = 1.0e5, D = -0.01, W = 10.0, Lc = 5.0;
    m::Project pr;
    pr.name = "anchors on driven edges";
    pr.mesh.elem_size = 1.0;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;
    m::Material soil;
    soil.model = m::SoilModel::LinearElastic;
    soil.E = 1.0e4; soil.nu = 0.3; soil.gamma_unsat = 0.0; soil.gamma_sat = 0.0;
    pr.materials.push_back(soil);
    pr.has_water = false;
    m::SoilPolygon P;
    P.name = "Block"; P.material = 0;
    // The mid-height points are vertices, so both anchor ends sit on mesh nodes at exactly y = 5
    // and the closed forms need no attachment correction.
    P.x = {0.0, W, W, W, 0.0, 0.0};
    P.y = {0.0, 0.0, 5.0, 10.0, 10.0, 5.0};
    P.edge_bc = {(int)m::BCType::VerticallyFixed, 0, 0, 0, 0, 0};
    pr.polygons.push_back(P);
    m::AnchorMaterial am;
    am.name = "Tie"; am.EA = EA; am.Lspacing = 1.0;
    pr.anchors.push_back(am);
    pr.structs.push_back(line(m::StructKind::Anchor, "Strut", 0.0, 5.0, W, 5.0));
    pr.structs.push_back(line(m::StructKind::Anchor, "Deadman", W, 5.0, W + Lc, 5.0));
    pr.disps.push_back(driven("Left edge", 0.0, 0.0, 0.0, 10.0, true, 0.0, false, 0.0));
    pr.disps.push_back(driven("Right edge", W, 0.0, W, 10.0, true, D, false, 0.0));

    const auto res = run(pr, {1});
    if (res.size() != 2) return;
    const auto* strut = force_of(res[1], "Strut");
    const auto* dead = force_of(res[1], "Deadman");
    check(strut && !strut->stations.empty() && dead && !dead->stations.empty(),
          "both anchors report an axial force");
    if (!strut || strut->stations.empty() || !dead || dead->stations.empty()) return;
    const double n_strut = strut->stations[0].N, n_dead = dead->stations[0].N;
    const double e_strut = EA * D / W, e_dead = EA * (-D) / Lc;
    std::printf("   (b) strut:   N = %.9f kN, closed form %.9f\n", n_strut, e_strut);
    std::printf("   (c) deadman: N = %.9f kN, closed form %.9f\n", n_dead, e_dead);
    check(rel(n_strut, e_strut) < 1e-9,
          "a strut between two driven edges reports N = EA D / L, within 1e-9");
    check(rel(n_dead, e_dead) < 1e-9,
          "a fixed-end anchor from a driven edge reports N = EA (-D) / L_c, within 1e-9");
}

}  // namespace

int main() {
    std::printf("Structures on driven nodes, read from the force report (KV-STR-009)\n\n");
    plate_on_a_settling_support();
    anchors_on_driven_edges();
    if (g_failures == 0) {
        std::printf("\nOK: the force report carries a driven node's motion, for plates and anchors\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
