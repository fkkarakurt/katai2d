// ROCK-SOLID contract: the GUI compute entry points (mesh_from_project, solve_gravity_le,
// solve_phases, solve_groundwater_flow) are called by the GUI WITHOUT a try/catch, so they must
// NEVER throw, NEVER hang, and NEVER hand back a non-finite (NaN/Inf) field. An end user with a
// complex or malformed model must get a clean "ok = false + honest message" -- not a crash.
//
// This battery throws adversarial / degenerate / extreme configurations at the whole pipeline and
// asserts robustness: a clean rejection (ok = false) is fine, a successful finite solve is fine, but
// an escaped exception or a NaN/Inf in the result is a FAILURE. Each case logs (flushed) before it
// runs, so even a hard crash identifies the offending case. Small meshes keep it fast.
#include <katai/jobs/flow_driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <limits>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
const double kNaN = std::numeric_limits<double>::quiet_NaN();
const double kInf = std::numeric_limits<double>::infinity();

bool finite_vec(const Eigen::VectorXd& v) {
    for (int i = 0; i < v.size(); ++i) if (!std::isfinite(v[i])) return false;
    return true;
}

// A baseline well-posed material/polygon a case can start from and then corrupt.
m::Material good_mat() {
    m::Material s; s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e4; s.nu = 0.3; s.c = 5.0; s.phi = 25.0; s.psi = 0.0;
    s.gamma_unsat = 17.0; s.gamma_sat = 20.0; s.kx = 1.0; s.ky = 1.0;
    return s;
}
m::SoilPolygon box(double w = 10.0, double h = 6.0) {
    m::SoilPolygon P; P.material = 0;
    P.x = {0, w, w, 0}; P.y = {0, 0, h, h};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    return P;
}

// Run one project through mesh + single-phase solve. Robust unless it throws or leaks NaN/Inf.
void rc(const char* name, const m::Project& pr,
        InitialPhase phase = InitialPhase::K0Procedure, double area = 1.5) {
    std::printf("RUN  %-42s ", name); std::fflush(stdout);
    try {
        const auto M = katai::app::mesh_from_project(pr, area, 6);
        if (!M.ok) { std::printf("-> mesh rejected (ok)\n"); return; }
        const auto R = katai::app::solve_gravity_le(pr, M.mesh, phase);
        if (!R.ok) { std::printf("-> solve rejected (ok)\n"); return; }
        if (!finite_vec(R.disp) || !std::isfinite(R.max_disp)) {
            std::printf("-> FAIL: NaN/Inf in displacement\n"); ++g_failures; return;
        }
        for (const auto& s : R.stress.stress)
            if (!std::isfinite(s(0)) || !std::isfinite(s(1)) || !std::isfinite(s(2))) {
                std::printf("-> FAIL: NaN/Inf in stress\n"); ++g_failures; return;
            }
        std::printf("-> solved, finite (ok)\n");
    } catch (const std::exception& e) {
        std::printf("-> FAIL: threw std::exception: %s\n", e.what()); ++g_failures;
    } catch (...) {
        std::printf("-> FAIL: threw (unknown)\n"); ++g_failures;
    }
}

// Multi-phase variant (solve_phases must also never throw / leak NaN).
void rc_phases(const char* name, const m::Project& pr, double area = 1.5) {
    std::printf("RUN  %-42s ", name); std::fflush(stdout);
    try {
        const auto M = katai::app::mesh_from_project(pr, area, 6);
        if (!M.ok) { std::printf("-> mesh rejected (ok)\n"); return; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        for (const auto& R : res)
            if (R.ok && !finite_vec(R.disp)) {
                std::printf("-> FAIL: NaN/Inf in a phase displacement\n"); ++g_failures; return;
            }
        std::printf("-> ran %zu phase(s) (ok)\n", res.size());
    } catch (const std::exception& e) {
        std::printf("-> FAIL: threw std::exception: %s\n", e.what()); ++g_failures;
    } catch (...) { std::printf("-> FAIL: threw (unknown)\n"); ++g_failures; }
}

// Flow variant.
void rc_flow(const char* name, const m::Project& pr, double area = 1.5) {
    std::printf("RUN  %-42s ", name); std::fflush(stdout);
    try {
        const auto M = katai::app::mesh_from_project(pr, area, 6);
        if (!M.ok) { std::printf("-> mesh rejected (ok)\n"); return; }
        const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
        if (F.ok && !finite_vec(F.head)) {
            std::printf("-> FAIL: NaN/Inf in flow head\n"); ++g_failures; return;
        }
        std::printf("-> %s (ok)\n", F.ok ? "solved" : "rejected");
    } catch (const std::exception& e) {
        std::printf("-> FAIL: threw std::exception: %s\n", e.what()); ++g_failures;
    } catch (...) { std::printf("-> FAIL: threw (unknown)\n"); ++g_failures; }
}

}  // namespace

int main() {
    std::printf("Robustness battery: the compute pipeline must never crash / leak NaN.\n\n");

    // --- empty / unassigned -----------------------------------------------------------------
    { m::Project pr; rc("empty project (no materials/geometry)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); rc("material but no geometry", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); m::SoilPolygon P = box(); P.material = -1;
      pr.polygons.push_back(P); rc("polygon with no material (-1)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); m::SoilPolygon P = box(); P.material = 7;
      pr.polygons.push_back(P); rc("polygon material out of range", pr); }

    // --- degenerate geometry ----------------------------------------------------------------
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {1, 1, 1, 1}; P.y = {2, 2, 2, 2};
      P.edge_bc.assign(4, (int)m::BCType::FullyFixed); pr.polygons.push_back(P);
      rc("zero-area polygon (coincident vertices)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {0, 5, 10}; P.y = {0, 0, 0};
      P.edge_bc.assign(3, (int)m::BCType::FullyFixed); pr.polygons.push_back(P);
      rc("collinear polygon (degenerate triangle)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {0, 10, 0, 10}; P.y = {0, 10, 10, 0};   // bowtie
      P.edge_bc.assign(4, (int)m::BCType::FullyFixed); pr.polygons.push_back(P);
      rc("self-intersecting polygon (bowtie)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {0, 0, 10, 10, 0}; P.y = {0, 0, 0, 6, 6}; // dup vertex
      P.edge_bc.assign(5, (int)m::BCType::FullyFixed); pr.polygons.push_back(P);
      rc("duplicate consecutive vertices", pr); }
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {0, kNaN, 10, 0}; P.y = {0, 0, 6, 6};
      P.edge_bc.assign(4, (int)m::BCType::FullyFixed); pr.polygons.push_back(P);
      rc("NaN coordinate in geometry", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box(1e-10, 1e-10));
      rc("microscopic domain (1e-10)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box(1e6, 1e6));
      rc("enormous domain (1e6)", pr, InitialPhase::K0Procedure, 1e10); }

    // --- pathological materials --------------------------------------------------------------
    auto with_mat = [](m::Material mm) {
        m::Project pr; pr.materials.push_back(mm); pr.polygons.push_back(box());
        m::Load L; L.kind = m::LoadKind::Distributed; L.x1 = 2; L.y1 = 6; L.x2 = 8; L.y2 = 6;
        L.qy1 = L.qy2 = -50; pr.loads.push_back(L); return pr;
    };
    { auto mm = good_mat(); mm.E = 0.0;        rc("E = 0", with_mat(mm)); }
    { auto mm = good_mat(); mm.E = -1e4;       rc("E negative", with_mat(mm)); }
    { auto mm = good_mat(); mm.nu = 0.5;       rc("nu = 0.5 (incompressible/singular)", with_mat(mm)); }
    { auto mm = good_mat(); mm.nu = 0.7;       rc("nu = 0.7 (> 0.5, invalid)", with_mat(mm)); }
    { auto mm = good_mat(); mm.nu = -3.0;      rc("nu = -3 (< -1, invalid)", with_mat(mm)); }
    { auto mm = good_mat(); mm.E = kNaN;       rc("E = NaN", with_mat(mm)); }
    { auto mm = good_mat(); mm.E = kInf;       rc("E = Inf", with_mat(mm)); }
    { auto mm = good_mat(); mm.E = 1e308;      rc("E = 1e308 (overflow-prone)", with_mat(mm)); }
    { auto mm = good_mat(); mm.c = -10; mm.phi = 120; rc("MC c<0, phi=120 deg", with_mat(mm)); }
    { auto mm = good_mat(); mm.model = m::SoilModel::HardeningSoil; mm.E50ref = 0; mm.Eoedref = 0;
      mm.Eurref = 0; rc("HS with zero stiffness moduli", with_mat(mm)); }

    // --- loads / BCs -------------------------------------------------------------------------
    { m::Project pr; pr.materials.push_back(good_mat()); m::SoilPolygon P = box();
      P.edge_bc.assign(4, (int)m::BCType::Free); pr.polygons.push_back(P);   // no support -> rigid body
      rc("no boundary conditions (rigid body)", pr); }
    { auto pr = with_mat(good_mat()); pr.loads[0].qy1 = kInf; rc("Inf load magnitude", pr); }
    { auto pr = with_mat(good_mat()); pr.loads[0].x1 = 100; pr.loads[0].x2 = 200;
      pr.loads[0].y1 = pr.loads[0].y2 = 100; rc("load far outside the geometry", pr); }
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box());
      m::Load L; L.kind = m::LoadKind::Point; L.x1 = kNaN; L.y1 = 6; L.qy1 = -10;
      pr.loads.push_back(L); rc("point load at NaN location", pr); }

    // --- undrained edge cases ---------------------------------------------------------------
    { auto mm = good_mat(); mm.drainage = m::Drainage::Undrained; mm.nu = 0.5;
      rc("undrained with nu' = 0.5 (Kw/n singular)", with_mat(mm)); }
    { auto mm = good_mat(); mm.drainage = m::Drainage::Undrained; mm.nu = 0.499;
      rc("undrained with nu' = 0.499 (Kw/n -> huge)", with_mat(mm)); }
    { auto mm = good_mat(); mm.drainage = m::Drainage::UndrainedB; mm.c = 10; mm.phi = 45;
      rc("Undrained (B), phi'=45 (forced to 0)", with_mat(mm)); }

    // --- consolidation phase (solve_phases) -------------------------------------------------
    auto consol_proj = [](double duration, int steps, double kx) {
        m::Project pr; auto mm = good_mat(); mm.model = m::SoilModel::LinearElastic; mm.kx = kx; mm.ky = kx;
        pr.materials.push_back(mm); m::SoilPolygon P = box(); P.edge_flow.assign(4, (int)m::FlowBCType::Closed);
        P.edge_flow[2] = (int)m::FlowBCType::Head; P.edge_head.assign(4, 0.0); P.edge_head[2] = 6;
        pr.polygons.push_back(P); pr.has_water = false;
        m::Load L; L.kind = m::LoadKind::Distributed; L.x1 = 0; L.y1 = 6; L.x2 = 10; L.y2 = 6;
        L.qy1 = L.qy2 = -40; pr.loads.push_back(L);
        m::Phase ph; ph.type = m::PhaseType::Consolidation; ph.duration = duration; ph.time_steps = steps;
        ph.load_active = {1}; pr.initial.load_active = {0}; pr.phases.push_back(ph);
        return pr;
    };
    rc_phases("consolidation: duration = 0", consol_proj(0.0, 20, 1.0));
    rc_phases("consolidation: steps = 0", consol_proj(5.0, 0, 1.0));
    rc_phases("consolidation: steps = -5", consol_proj(5.0, -5, 1.0));
    rc_phases("consolidation: kx = 0 (no permeability)", consol_proj(5.0, 20, 0.0));
    rc_phases("consolidation: negative duration", consol_proj(-3.0, 20, 1.0));

    // --- multi-phase activations ------------------------------------------------------------
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box());
      m::Phase ph; ph.type = m::PhaseType::Plastic; ph.poly_active = {0};   // deactivate the only soil
      pr.phases.push_back(ph); rc_phases("phase deactivates the only soil region", pr); }

    // --- flow edge cases --------------------------------------------------------------------
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box());
      rc_flow("flow: no flow BCs", pr); }
    { m::Project pr; auto mm = good_mat(); mm.kx = 0; mm.ky = 0; pr.materials.push_back(mm);
      m::SoilPolygon P = box(); P.edge_flow.assign(4, (int)m::FlowBCType::Closed);
      P.edge_flow[1] = (int)m::FlowBCType::Head; P.edge_flow[3] = (int)m::FlowBCType::Head;
      P.edge_head.assign(4, 0.0); P.edge_head[1] = 2; P.edge_head[3] = 5; pr.polygons.push_back(P);
      rc_flow("flow: zero permeability", pr); }

    // --- structural elements (the real end-user crash surface: split mesh / embedded wall / anchors,
    //     drawn in degenerate or out-of-place positions) -------------------------------------------
    auto struct_proj = [](m::StructKind kind, double x1, double y1, double x2, double y2, bool iface) {
        m::Project pr; pr.materials.push_back(good_mat());
        pr.plates.push_back(m::PlateMaterial{}); pr.anchors.push_back(m::AnchorMaterial{});
        pr.geogrids.push_back(m::GeogridMaterial{}); pr.embedded.push_back(m::EmbeddedBeamMaterial{});
        pr.polygons.push_back(box(10, 8));
        m::StructElement s; s.kind = kind; s.x1 = x1; s.y1 = y1; s.x2 = x2; s.y2 = y2;
        s.material = 0; s.iface_pos = iface; pr.structs.push_back(s);
        return pr;
    };
    using K = m::StructKind;
    rc("plate: zero length", struct_proj(K::Plate, 5, 4, 5, 4, false));
    rc("plate: NaN coordinates", struct_proj(K::Plate, 2, 4, kNaN, 4, false));
    rc("plate: entirely outside the domain", struct_proj(K::Plate, 50, 50, 60, 60, false));
    rc("embedded wall: vertical, mid-domain", struct_proj(K::Plate, 5, 2, 5, 7, true));
    rc("embedded wall: zero length + iface", struct_proj(K::Plate, 5, 4, 5, 4, true));
    rc("embedded wall: outside the domain", struct_proj(K::Plate, 50, 0, 50, 8, true));
    rc("embedded wall: horizontal + iface", struct_proj(K::Plate, 2, 4, 8, 4, true));
    rc("anchor: both ends outside soil", struct_proj(K::Anchor, 50, 50, 60, 60, false));
    rc("anchor: coincident ends", struct_proj(K::Anchor, 5, 4, 5, 4, false));
    rc("anchor: NaN end", struct_proj(K::Anchor, 5, 4, kNaN, 4, false));
    rc("geogrid: zero length", struct_proj(K::Geogrid, 5, 4, 5, 4, false));
    rc("embedded beam: zero length", struct_proj(K::EmbeddedBeam, 5, 4, 5, 4, false));
    { auto pr = struct_proj(K::EmbeddedBeam, 5, 1, 5, 6, false); pr.embedded[0].diameter = 0.0;
      rc("embedded beam: zero diameter", pr); }

    // --- axisymmetric edge cases ------------------------------------------------------------
    { m::Project pr; pr.materials.push_back(good_mat()); pr.polygons.push_back(box()); pr.axisymmetric = true;
      rc("axisymmetric: soil-only (valid)", pr); }
    { m::Project pr; pr.materials.push_back(good_mat());
      m::SoilPolygon P; P.material = 0; P.x = {-5, 5, 5, -5}; P.y = {0, 0, 6, 6};   // x<0 = negative radius
      P.edge_bc.assign(4, (int)m::BCType::FullyFixed); pr.polygons.push_back(P); pr.axisymmetric = true;
      rc("axisymmetric: negative radius (x<0)", pr); }

    // --- a genuinely COMPLEX, well-posed model must solve to a finite result (not just reject) -----
    {
        m::Project pr;
        m::Material clay = good_mat(); clay.name = "Clay"; clay.E = 8000; clay.c = 8; clay.phi = 22;
        clay.drainage = m::Drainage::Undrained;
        m::Material sand = good_mat(); sand.name = "Sand"; sand.E = 25000; sand.c = 1; sand.phi = 33;
        pr.materials = {clay, sand};
        pr.plates.push_back(m::PlateMaterial{}); pr.anchors.push_back(m::AnchorMaterial{});
        m::SoilPolygon lower; lower.material = 0; lower.x = {0, 24, 24, 0}; lower.y = {0, 0, 8, 8};
        lower.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                         (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
        m::SoilPolygon upper; upper.material = 1; upper.x = {0, 24, 24, 0}; upper.y = {8, 8, 14, 14};
        upper.edge_bc = {0, (int)m::BCType::NormallyFixed, (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
        pr.polygons = {lower, upper};
        m::StructElement wall; wall.kind = K::Plate; wall.x1 = 8; wall.y1 = 2; wall.x2 = 8; wall.y2 = 14;
        wall.material = 0; wall.iface_pos = true; wall.iface_neg = true;
        m::StructElement strut; strut.kind = K::Anchor; strut.x1 = 8; strut.y1 = 12; strut.x2 = 4; strut.y2 = 12;
        strut.material = 0;
        pr.structs = {wall, strut};
        m::Load surcharge; surcharge.kind = m::LoadKind::Distributed;
        surcharge.x1 = 12; surcharge.y1 = 14; surcharge.x2 = 20; surcharge.y2 = 14;
        surcharge.qy1 = surcharge.qy2 = -30; pr.loads.push_back(surcharge);
        pr.has_water = true; pr.wx = {0, 24}; pr.wy = {10, 10};
        // initial: soil only; phase 1: excavate the left of the wall above y=8 + install strut + load
        m::Phase ph; ph.type = m::PhaseType::Plastic;
        ph.poly_active = {1, 1}; ph.struct_active = {1, 1}; ph.load_active = {1};
        pr.initial.struct_active = {1, 0}; pr.initial.load_active = {0};
        pr.phases.push_back(ph);
        rc_phases("COMPLEX: 2-layer + embedded wall + strut + surcharge + water + phase", pr, 2.0);
    }

    std::printf("\n");
    if (g_failures == 0) {
        std::printf("OK: the compute pipeline survived every adversarial input (no crash / NaN)\n");
        return 0;
    }
    std::fprintf(stderr, "%d robustness failure(s) -- the pipeline can crash / leak NaN on bad input\n",
                 g_failures);
    return 1;
}
