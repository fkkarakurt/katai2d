// Dynamic (seismic) analysis through the GUI compute path (build_problem PhaseType::Dynamic +
// solve_phases). The dynamic CORE (Newmark + Rayleigh + 2D consistent mass) is validated to machine
// precision in test_dynamics; here the INTEGRATED path -- project + material -> mesh -> Dynamic phase
// (horizontal base acceleration, rigid base, free sides) -- must reproduce 1D site-response physics.
//
// A thin soil column (rigid base, free sides) shaken horizontally responds as a 1D shear column
// (u_x = u_x(y), u_y = 0). Under a base sine a_g = A sin(w t) at the fundamental f_1 = Vs/(4H) the
// steady surface (relative) displacement follows the fundamental-mode SDOF:
//     |u_surf| = Gamma_1 * A / (w_1^2 * 2 xi),   Gamma_1 = 4/pi  (uniform shear column, phi=sin(pi y/2H)).
// Off resonance the response is far smaller -> the driver reproduces the site's resonance + selectivity.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/analysis/response_spectrum.hpp>   // EC8 elastic spectrum oracle (direct engine use)
#include <katai/model/project.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

void test_site_response() {
    std::printf("-- 1D site response through the GUI Dynamic phase path --\n");
    // Column: W x H. E, nu, gamma -> rho = gamma/g, G = E/(2(1+nu)), Vs = sqrt(G/rho), f_1 = Vs/(4H).
    constexpr double W = 2.0, H = 20.0, E = 208000.0, nu = 0.3, gamma = 19.62;  // gamma/9.81 = 2.0
    const double g = 9.81, rho = gamma / g, G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double f1 = Vs / (4 * H), w1 = 2 * kPi * f1, A = 1.0, xi = 0.05;
    const double u_surf = (4.0 / kPi) * A / (w1 * w1 * 2 * xi);  // fundamental-mode steady surface disp
    std::printf("   Vs=%.1f m/s  f_1=%.3f Hz  (expected peak |u_surf| = %.4f m)\n", Vs, f1, u_surf);

    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = E; s.nu = nu; s.gamma_unsat = gamma; s.gamma_sat = gamma; s.e_init = 0.5;
    pr.materials.push_back(s);
    m::SoilPolygon P; P.material = 0;
    P.x = {0, W, W, 0};
    P.y = {0, 0, H, H};
    // bottom = rigid base (fully fixed). Sides = VerticallyFixed (u_y=0, u_x free): on a thin column
    // this suppresses cantilever BENDING (which free sides allow) so the column deforms in pure
    // horizontal SHEAR -> a 1D SH shear column. top = free surface.
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(P);
    pr.has_water = false;

    auto dyn_phase = [&](double freq, double dur, int steps) {
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Dynamic";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = A; p.seismic_freq = freq;
        p.damping_ratio = xi; p.rayleigh_f1 = f1; p.rayleigh_f2 = 3 * f1;
        p.duration = dur; p.time_steps = steps;
        return p;
    };
    pr.phases.push_back(dyn_phase(f1, 8.0, 800));        // resonance (~20 cycles -> steady)
    pr.phases.push_back(dyn_phase(f1 / 3.0, 9.0, 720));  // off-resonance (well below f_1)

    const auto M = katai::app::mesh_from_project(pr, 0.4, 6);
    check(M.ok, "column meshed");
    if (!M.ok) { std::printf("  (%s)\n", M.message.c_str()); return; }

    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3, "initial + 2 dynamic phases ran");
    if (res.size() != 3) return;
    check(res[1].ok && res[2].ok, "both dynamic phases solved");
    if (!res[1].ok) { std::printf("  (%s)\n", res[1].message.c_str()); return; }

    const double res_peak = res[1].max_disp, off_peak = res[2].max_disp;
    check(!res[1].consol_time.empty(), "dynamic phase produced a response-time series");
    std::printf("   resonance  f=%.3f: peak |u|=%.4f m  (theory %.4f, err %.1f%%)\n",
                f1, res_peak, u_surf, 100 * (res_peak - u_surf) / u_surf);
    std::printf("   off-reson. f=%.3f: peak |u|=%.4f m  (ratio to resonance %.2f)\n",
                f1 / 3.0, off_peak, off_peak / res_peak);
    // Resonant surface displacement matches the fundamental-mode SDOF (2D FE + transient + higher modes
    // -> ~15% band around the continuous-column theory).
    check(std::fabs(res_peak - u_surf) < 0.15 * u_surf, "resonant surface displacement = 4/pi A/(w1^2 2xi)");
    // Frequency selectivity: off-resonance is far smaller than at the fundamental.
    check(off_peak < 0.4 * res_peak, "off-resonance response is much smaller (site amplifies at f_1)");

    // --- Surface acceleration + response spectrum + TBDY design spectrum (D4a x D4b) --------------
    const double surf_a = res[1].dyn_peak_surface_a, a_theory = 2.0 / kPi / xi * A;  // 2/(pi xi) A
    check(!res[1].dyn_surface_ax.empty() && res[1].dyn_response_sa.size() == res[1].dyn_period.size(),
          "surface acceleration + response spectrum produced");
    std::printf("   peak surface accel = %.3f m/s^2  (2/(pi xi) A = %.3f, err %.1f%%)\n",
                surf_a, a_theory, 100 * (surf_a - a_theory) / a_theory);
    // Peak surface accel = the transfer-function amplification at resonance (Kramer 7.30; verified in
    // test_dynamics D1). This ties the 2D driver output back to 1D site-response theory.
    check(std::fabs(surf_a - a_theory) < 0.10 * a_theory, "peak surface accel = 2/(pi xi) A (resonance transfer fn)");
    // Response-spectrum PGA (shortest, ~rigid period) equals the peak surface acceleration.
    check(std::fabs(res[1].dyn_response_sa.front() - surf_a) < 0.05 * surf_a,
          "response spectrum PGA (T->0) = peak surface accel");
    // The surface motion is dominated by f_1, so its response spectrum peaks near T = 1/f_1.
    size_t pk = 0; double pv = 0;
    for (size_t i = 0; i < res[1].dyn_response_sa.size(); ++i)
        if (res[1].dyn_response_sa[i] > pv) { pv = res[1].dyn_response_sa[i]; pk = i; }
    std::printf("   response spectrum peak at T=%.3f s  (1/f_1 = %.3f s)\n", res[1].dyn_period[pk], 1.0 / f1);
    check(std::fabs(res[1].dyn_period[pk] - 1.0 / f1) < 0.15 / f1, "response spectrum peaks near T = 1/f_1");
    // TBDY design spectrum (default phase: S_S=1.0, S_1=0.4, ZC): plateau S_ae=S_DS=S_S F_S=1.2 -> 11.77
    // m/s^2 (Ta<=T<=Tb, here Ta=0.1, Tb=0.5). Check the ordinate at T~0.3 (plateau).
    size_t ip = 0; double bd = 1e300;
    for (size_t i = 0; i < res[1].dyn_period.size(); ++i)
        if (std::fabs(res[1].dyn_period[i] - 0.3) < bd) { bd = std::fabs(res[1].dyn_period[i] - 0.3); ip = i; }
    std::printf("   TBDY design spectrum @T~0.3s = %.3f m/s^2  (plateau S_DS*g = %.3f)\n",
                res[1].dyn_design_sa[ip], 1.2 * 9.81);
    check(std::fabs(res[1].dyn_design_sa[ip] - 1.2 * 9.81) < 0.02 * (1.2 * 9.81),
          "TBDY design spectrum plateau = S_DS*g (S_S*F_S, ZC)");
}

// Free-field lateral boundaries through the GUI path (build_problem seismic_free_field). Under a
// uniform horizontal base body force -M r a_g, a column with FREE sides bends as a cantilever (the
// vertical fibres stretch) instead of shearing -- so its surface response departs from the 1D shear-
// column site response. Turning on Lysmer FREE-FIELD sides drives those sides to follow the 1D free
// field, recovering the shear site response (surface amplification ~ (4/pi) A/(w1^2 2xi)), close to
// the VerticallyFixed (pure-SH) reference. This exercises boundary-edge extraction + the 1D free-
// field column + the driving force end-to-end on an unstructured mesh.
void test_freefield_sides() {
    std::printf("-- Free-field lateral boundaries through the GUI Dynamic phase path --\n");
    constexpr double W = 4.0, H = 20.0, E = 208000.0, nu = 0.3, gamma = 19.62;
    const double g = 9.81, rho = gamma / g, G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double f1 = Vs / (4 * H), w1 = 2 * kPi * f1, A = 1.0, xi = 0.05;
    const double u_th = (4.0 / kPi) * A / (w1 * w1 * 2 * xi);
    std::printf("   Vs=%.1f m/s  f_1=%.3f Hz  (1D site-response |u_surf| = %.4f m)\n", Vs, f1, u_th);

    auto run = [&](int side_bc, bool ff) -> double {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = E; s.nu = nu; s.gamma_unsat = gamma; s.gamma_sat = gamma; s.e_init = 0.5;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, side_bc, (int)m::BCType::Free, side_bc};
        pr.polygons.push_back(P); pr.has_water = false;
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Dyn";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = A; p.seismic_freq = f1;
        p.damping_ratio = xi; p.rayleigh_f1 = f1; p.rayleigh_f2 = 3 * f1;
        p.duration = 8.0; p.time_steps = 800; p.seismic_free_field = ff;
        pr.phases.push_back(p);
        const auto M = katai::app::mesh_from_project(pr, 0.4, 6);
        if (!M.ok) { std::printf("  (mesh: %s)\n", M.message.c_str()); return -1.0; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.size() != 2 || !res[1].ok) return -1.0;
        return res[1].max_disp;
    };

    const double sh    = run((int)m::BCType::VerticallyFixed, false);  // pure-SH reference (uy=0 sides)
    const double freeS = run((int)m::BCType::Free, false);            // free sides, NO free-field
    const double ff    = run((int)m::BCType::Free, true);             // free sides + free-field boundary
    check(sh > 0 && freeS > 0 && ff > 0, "all three dynamic configurations solved");
    std::printf("   SH (uy=0 sides)      peak |u|=%.4f m  (err vs 1D %.1f%%)\n", sh, 100 * (sh - u_th) / u_th);
    std::printf("   free sides, no FF    peak |u|=%.4f m  (err vs 1D %.1f%%)\n", freeS, 100 * (freeS - u_th) / u_th);
    std::printf("   free sides + FF      peak |u|=%.4f m  (err vs 1D %.1f%%)\n", ff, 100 * (ff - u_th) / u_th);
    // The SH reference reproduces the 1D site response (regression of the tested path).
    check(std::fabs(sh - u_th) < 0.15 * u_th, "SH sides reproduce the 1D site response");
    // Free-field sides recover the shear site response from FREE sides (measured err ~0.1%, vs
    // plain free sides ~-34%): the sides are driven to follow the 1D free-field column.
    check(std::fabs(ff - u_th) < 0.10 * u_th, "free-field sides recover the 1D site response");
    // ... and do so far better than plain free sides (which bend / reflect).
    check(std::fabs(ff - u_th) < std::fabs(freeS - u_th),
          "free-field sides are closer to the 1D site response than plain free sides");
}

// Robustness: a Dynamic phase must NEVER crash the process, whatever the model has. Soil-only runs
// (plain / free-field / tri15 / with water) must solve; so must a model with a wall or an anchor, now
// that the structures take part in the dynamic system (soil-structure interaction). An embedded beam
// -- the one element still outside the dynamic assembly -- must be REJECTED gracefully (ok = false +
// a clear message), never crash. Each label is flushed before the solve so a crash (an access
// violation, which /EHsc catch(...) does NOT catch) still names the culprit.
void test_dynamic_robustness() {
    std::printf("-- Dynamic phase robustness (no crash on any model) --\n");
    auto base = []() {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free,
                     (int)m::BCType::Free, (int)m::BCType::Free};
        pr.polygons.push_back(P);
        return pr;
    };
    auto dyn = [](bool ff) {
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Ricker; p.seismic_amp = 1.0; p.seismic_freq = 3.0;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 2.0; p.time_steps = 200; p.seismic_free_field = ff;
        return p;
    };
    auto run = [&](const char* label, m::Project pr, int order) -> katai::app::SolveResult {
        std::printf("   [run] %-24s ... ", label); std::fflush(stdout);   // flush -> a crash names it
        const auto M = katai::app::mesh_from_project(pr, 1.0, order);
        if (!M.ok) { std::printf("MESH FAILED\n"); return {}; }
        std::fflush(stdout);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        katai::app::SolveResult last = res.empty() ? katai::app::SolveResult{} : res.back();
        std::printf("%s\n", last.ok ? "solved" : ("rejected: " + last.message).c_str());
        std::fflush(stdout);
        return last;
    };
    auto with_plate = [&]() {
        auto pr = base();
        m::PlateMaterial pm; pr.plates.push_back(pm);
        m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
        st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
        pr.structs.push_back(st); pr.phases.push_back(dyn(false)); return pr;
    };
    auto with_anchor = [&]() {
        auto pr = base();
        m::AnchorMaterial am; pr.anchors.push_back(am);
        m::StructElement st; st.kind = m::StructKind::Anchor; st.name = "Anchor";
        st.x1 = 5; st.y1 = 8; st.x2 = 9; st.y2 = 6; st.material = 0;
        pr.structs.push_back(st); pr.phases.push_back(dyn(false)); return pr;
    };
    auto with_wall = [&]() {           // plate + interface flag = embedded wall (own DOFs + split mesh)
        auto pr = base();
        m::PlateMaterial pm; pm.w = 5.0; pr.plates.push_back(pm);
        m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
        st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
        st.iface_pos = true; st.iface_neg = true;
        pr.structs.push_back(st); pr.phases.push_back(dyn(false)); return pr;
    };
    auto with_geogrid = [&]() {
        auto pr = base();
        m::GeogridMaterial gm; pr.geogrids.push_back(gm);
        m::StructElement st; st.kind = m::StructKind::Geogrid; st.name = "Geogrid";
        st.x1 = 4; st.y1 = 6; st.x2 = 16; st.y2 = 6; st.material = 0;
        pr.structs.push_back(st); pr.phases.push_back(dyn(false)); return pr;
    };
    auto with_embedded_beam = [&]() {  // still outside the dynamic assembly -> must reject, not crash
        auto pr = base();
        m::EmbeddedBeamMaterial em; pr.embedded.push_back(em);
        m::StructElement st; st.kind = m::StructKind::EmbeddedBeam; st.name = "Pile";
        st.x1 = 10; st.y1 = 2; st.x2 = 10; st.y2 = 10; st.material = 0;
        pr.structs.push_back(st); pr.phases.push_back(dyn(false)); return pr;
    };
    m::Project p1 = base(); p1.phases.push_back(dyn(false));
    m::Project p2 = base(); p2.phases.push_back(dyn(true));
    m::Project p3 = base(); p3.phases.push_back(dyn(false));
    m::Project p4 = base(); p4.wx = {0, 20}; p4.wy = {6, 6}; p4.has_water = true; p4.phases.push_back(dyn(false));

    check(run("plain soil tri6", p1, 6).ok, "soil-only dynamic solves (tri6)");
    check(run("free-field tri6", p2, 6).ok, "soil-only dynamic with free-field solves");
    check(run("plain soil tri15", p3, 15).ok, "soil-only dynamic solves (tri15)");
    check(run("water tri6", p4, 6).ok, "soil-only dynamic with a water table solves");
    // Structural elements now take part in the dynamic system (SSI) -> they must SOLVE. Each of these
    // used to corrupt the heap (unassembled structural DOFs -> singular K_eff -> access violation).
    check(run("bonded plate tri6", with_plate(), 6).ok, "dynamic phase with a bonded plate solves");
    check(run("embedded wall tri6", with_wall(), 6).ok,
          "dynamic phase with an embedded wall (plate + interface, own DOFs) solves");
    check(run("anchor tri6", with_anchor(), 6).ok, "dynamic phase with an anchor solves");
    check(run("geogrid tri6", with_geogrid(), 6).ok, "dynamic phase with a geogrid solves");
    check(run("embedded wall tri15", with_wall(), 15).ok,
          "dynamic phase with an embedded wall solves (tri15: plate5 + interface5)");
    // Embedded beams (pile rows) are in the dynamic system too now: their skin/foot coupling is
    // mesh-NONCONFORMING, so assemble_structural_stiffness dispatches on the soil element type.
    check(run("embedded beam tri6", with_embedded_beam(), 6).ok,
          "dynamic phase with an embedded beam (pile row) solves");
    check(run("embedded beam tri15", with_embedded_beam(), 15).ok,
          "dynamic phase with an embedded beam solves (tri15 skin coupling)");
}

// Seismic structural forces through the GUI Dynamic path: the wall's N/Q/M ENVELOPE over the shaking
// (the numbers a section is designed with). The physics oracle -- the envelope equals the static
// moment under the same body force in the quasi-static limit -- lives in test_ssi_dynamics (g); here
// the wiring must hold: the envelope reaches the result, is flagged as an envelope, is a real force,
// and scales the way the linear system must.
void test_seismic_structural_forces() {
    std::printf("-- Seismic wall force envelope through the GUI Dynamic phase path --\n");
    auto model = [&](double amp) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P);
        m::PlateMaterial pm;                       // 0.4 m concrete diaphragm wall
        // w = 0 DELIBERATELY: this sub-test pins the LINEARITY of the dynamic increment
        // (2x a_g -> 2x envelope). Since the self-weight audit fix the static parent
        // carries a real weight-induced M/N/tau baseline, and the dynamic envelope
        // reports the TOTAL (static + dynamic) -- with w > 0 the total no longer
        // scales, which is correct physics, not linearity loss. The weight statics
        // are pinned in test_input_audit (f)/(g); the mass effect and the static+
        // dynamic superposition are pinned by the neighbouring sub-tests with w > 0.
        pm.EA = 3.0e7 * 0.4; pm.EI = 3.0e7 * 0.064 / 12.0; pm.w = 0.0;
        pr.plates.push_back(pm);
        m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
        st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
        st.iface_pos = true; st.iface_neg = true;   // embedded wall -> Coulomb joints to report too
        pr.structs.push_back(st);
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = amp; p.seismic_freq = 2.5;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 4.0; p.time_steps = 400;
        pr.phases.push_back(p);
        return pr;
    };
    auto run = [&](double amp, katai::app::SolveResult& out) {
        m::Project pr = model(amp);
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) return false;
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.empty() || !res.back().ok) return false;
        out = res.back();
        return true;
    };
    katai::app::SolveResult r1, r2;
    check(run(1.0, r1) && run(2.0, r2), "dynamic runs with a wall solved");
    if (r1.struct_forces.empty()) { check(false, "dynamic phase reports the wall's forces"); return; }
    const auto& w1 = r1.struct_forces.front();
    check(r1.struct_forces.size() == 1 && w1.name == "Wall", "the wall's force diagram is reported");
    check(w1.envelope, "the dynamic phase flags its forces as an ENVELOPE (not one instant)");
    check(!w1.stations.empty(), "the envelope has stations along the wall");
    std::printf("   wall envelope: max|M| = %.4g kNm/m, max|Q| = %.4g kN/m, max|N| = %.4g kN/m over %d stations\n",
                w1.max_M, w1.max_Q, w1.max_N, (int)w1.stations.size());
    check(w1.max_M > 0.0 && w1.max_Q > 0.0, "the shaken wall carries a real moment and shear");
    // Every station of an envelope is a max of |.| -> non-negative by construction.
    bool nonneg = true;
    for (const auto& st : w1.stations) if (st.M < 0.0 || st.Q < 0.0 || st.N < 0.0) nonneg = false;
    check(nonneg, "envelope stations are max|.| (non-negative)");
    // Linear system: doubling the base amplitude doubles the forces.
    const double ratio = r2.struct_forces.front().max_M / w1.max_M;
    std::printf("   amplitude x2 -> max|M| ratio = %.6f  (linear system: 2)\n", ratio);
    check(std::fabs(ratio - 2.0) < 1e-6, "doubling a_g doubles the wall moment envelope (linearity)");

    // The wall<->soil Coulomb joints report their seismic envelope too (elastic branch).
    if (r1.interface_forces.empty()) { check(false, "the wall's joints report a seismic envelope"); return; }
    const auto& j1 = r1.interface_forces.front();
    check(j1.envelope, "the joint envelope is flagged as an envelope");
    check(!j1.any_slip, "the dynamic joint reports no slip (it is elastic by construction)");
    std::printf("   joint envelope: max|tau| = %.4g kPa, max|sigma_n| = %.4g kPa over %d stations\n",
                j1.max_abs_tau, j1.max_abs_sigma_n, (int)j1.stations.size());
    check(j1.max_abs_tau > 0.0, "the shaken joint carries a real shear stress");
    const double jratio = r2.interface_forces.front().max_abs_tau / j1.max_abs_tau;
    std::printf("   amplitude x2 -> joint max|tau| ratio = %.6f  (linear system: 2)\n", jratio);
    check(std::fabs(jratio - 2.0) < 1e-6, "doubling a_g doubles the joint shear envelope (linearity)");
}

// A Dynamic phase reports the TOTAL design action: its own linear increment superposed on the parent
// phase's static state (PLAXIS reports the total too -- Ref 9.4.5 -- and continues the parent's state,
// Sci 6.4). Two exact statements pin the machinery:
//
//   (1) ZERO-AMPLITUDE IDENTITY. Shake with a_g = 0 and the dynamic increment vanishes, so the
//       reported total MUST equal the parent phase's static forces -- EXACTLY, station for station.
//       Independent oracle: the static phase's own, already-validated output. This is what catches a
//       mismatched station mapping, a dropped static term, or a sign error in the superposition; none
//       of them survive an exact identity.
//   (2) The joint's Coulomb DEMAND/CAPACITY becomes computable only because the static state supplies
//       the total normal stress. A weak joint under strong shaking must be FLAGGED as unconservative,
//       because the elastic dynamic solve can never discover the slip itself.
void test_total_action_superposition() {
    std::printf("-- Total design action = parent static state + dynamic increment --\n");
    auto model = [&](bool with_dynamic, double amp, double rinter) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        // rinter_rigid defaults to TRUE and overrides Rinter (build_problem: Rinter = rigid ? 1.0 : R)
        // -- without clearing it, the joint is rigid whatever Rinter says.
        s.c = 5.0; s.phi = 30.0; s.rinter_rigid = false; s.Rinter = rinter;
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
        if (with_dynamic) {
            m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
            p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = amp; p.seismic_freq = 2.5;
            p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
            p.duration = 4.0; p.time_steps = 400;
            pr.phases.push_back(p);
        }
        return pr;
    };
    auto run = [&](m::Project pr, std::vector<katai::app::SolveResult>& out) {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) return false;
        out = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        return !out.empty() && out.back().ok;
    };

    // (1) a_g = 0 -> the dynamic increment is identically zero -> total == the parent's static action.
    std::vector<katai::app::SolveResult> rs, rd;
    check(run(model(false, 0.0, 1.0), rs), "static-only run solved");
    check(run(model(true, 0.0, 1.0), rd), "zero-amplitude dynamic run solved");
    if (rs.empty() || rd.empty() || !rs.back().ok || !rd.back().ok) return;
    const auto& stat = rs.back();          // the initial (K0) phase = the dynamic phase's parent
    const auto& dyn = rd.back();
    check(!dyn.struct_forces.empty() && dyn.struct_forces.front().superposed,
          "the dynamic phase superposed its parent's static state");
    check(!dyn.interface_forces.empty() && dyn.interface_forces.front().superposed,
          "the joint report superposed its parent's static state");
    if (dyn.struct_forces.empty() || stat.struct_forces.empty()) return;

    double dM = 0.0, mM = 0.0;
    const auto& ws = stat.struct_forces.front();
    const auto& wd = dyn.struct_forces.front();
    check(ws.stations.size() == wd.stations.size(), "the wall diagram lines up station for station");
    for (size_t k = 0; k < std::min(ws.stations.size(), wd.stations.size()); ++k) {
        // The static station is signed; the dynamic total is the design |.| of (static + 0).
        dM = std::fmax(dM, std::fabs(wd.stations[k].M - std::fabs(ws.stations[k].M)));
        mM = std::fmax(mM, std::fabs(ws.stations[k].M));
    }
    std::printf("   a_g = 0: max |M_total - |M_static|| = %.3e  (max |M_static| = %.3e)\n", dM, mM);
    check(dM < 1e-9 * std::fmax(mM, 1e-12),
          "with no shaking the reported total EQUALS the parent's static wall moment");

    double dT = 0.0, mT = 0.0;
    if (!stat.interface_forces.empty() && !dyn.interface_forces.empty()) {
        const auto& js = stat.interface_forces.front();
        const auto& jd = dyn.interface_forces.front();
        for (size_t k = 0; k < std::min(js.stations.size(), jd.stations.size()); ++k) {
            dT = std::fmax(dT, std::fabs(jd.stations[k].sigma_n - std::fabs(js.stations[k].sigma_n)));
            mT = std::fmax(mT, std::fabs(js.stations[k].sigma_n));
        }
        std::printf("   a_g = 0: max |sigma_n_total - |sigma_n_static|| = %.3e  (max = %.3e)\n", dT, mT);
        check(mT > 1e-6, "the K0 parent really does put normal stress on the joint (the check has teeth)");
        check(dT < 1e-9 * std::fmax(mT, 1e-12),
              "with no shaking the reported total EQUALS the parent's static joint stress");
    }

    // (2) Demand/capacity: computable ONLY because the parent supplies the total normal stress.
    // It is reported as a measure, not a verdict: sigma_n -> 0 at the ground surface, so tau_max -> c_i
    // and a short zone at the top of ANY wall exceeds capacity under real shaking (genuine local slip).
    // So do not assert "strong joint => never over"; assert that the measure DISCRIMINATES -- a weaker
    // joint, or harder shaking, must come out worse. That is a property of the physics, not a threshold
    // of my choosing.
    auto util_of = [&](double amp, double rinter, double& over) {
        std::vector<katai::app::SolveResult> r;
        if (!run(model(true, amp, rinter), r) || r.back().interface_forces.empty()) return -1.0;
        over = r.back().interface_forces.front().over_fraction;
        return r.back().interface_forces.front().max_utilisation;
    };
    // Discriminate by AMPLITUDE, at a fixed joint. Two traps make that the only clean axis:
    //  - hard enough shaking separates the joint somewhere (tau_max = 0 -> the ratio is unbounded and
    //    reports the cap), and a capped peak cannot rank anything;
    //  - ranking by Rinter does NOT work, and that is physics rather than a defect: PLAXIS's virtual-
    //    thickness formulation reduces interface STIFFNESS as Rinter^2 (G_i = Rinter^2 G -> k_s = G_i/t_i)
    //    but strength only as Rinter^1, so a weaker joint is also a SOFTER joint, its dynamic
    //    tau = k_s du_s falls faster than its capacity, and utilisation ~ Rinter -- a weaker joint can
    //    report a LOWER ratio. Measured here: Rinter 0.05 -> 0.99x vs Rinter 1.0 -> 3.42x.
    double o_hard = 0, o_mid = 0, o_soft = 0;
    const double u_hard = util_of(4.0, 1.00, o_hard);
    const double u_mid  = util_of(1.0, 1.00, o_mid);
    const double u_soft = util_of(0.2, 1.00, o_soft);
    std::printf("   Rinter=1.0:  a_g 4.0 -> %.2fx (%.0f%% over) | 1.0 -> %.2fx (%.0f%% over) | "
                "0.2 -> %.2fx (%.0f%% over)\n",
                u_hard, 100 * o_hard, u_mid, 100 * o_mid, u_soft, 100 * o_soft);
    check(u_hard > 0 && u_mid > 0 && u_soft > 0, "all utilisation runs solved");
    check(o_hard >= o_mid && o_mid >= o_soft, "harder shaking drives more of the joint past capacity");
    check(u_hard >= u_mid && u_mid >= u_soft, "harder shaking never reports a lower peak demand/capacity");
    check(o_soft < 1.0, "mild shaking leaves most of the joint within capacity (the measure discriminates)");
    // The message must carry the finding, and must NOT cry wolf when nothing is over.
    std::vector<katai::app::SolveResult> rw;
    check(run(model(true, 4.0, 0.05), rw), "weak-joint + hard-shaking run solved");
    if (!rw.empty() && rw.back().ok) {
        const auto& j = rw.back().interface_forces.front();
        std::printf("   weak + hard: %.2fx over %.0f%% -> message: %s\n", j.max_utilisation,
                    100 * j.over_fraction,
                    rw.back().message.find("unconservative") != std::string::npos ? "warns" : "SILENT");
        check(j.max_utilisation > 1.0, "the weak joint is driven past capacity");
        check(rw.back().message.find("unconservative") != std::string::npos,
              "the phase message reports the joint is unconservative there");
        check(rw.back().message.find("% of the joint") != std::string::npos,
              "the message says HOW MUCH of the joint is over (a measure, not a verdict)");
    }
}

// Soil-structure interaction through the GUI Dynamic path. Two independent statements, neither of
// which a wrong assembly could satisfy by accident:
//   (1) VANISHING-STRUCTURE LIMIT: shrink the wall's stiffness AND weight towards zero and the site
//       response must converge back onto the verified soil-only answer. This is the continuity check
//       -- it says the SSI code adds nothing spurious when there is nothing to add.
//   (2) The wall MATTERS: a real (stiff, heavy) wall must change the response measurably. Without
//       this, (1) would also pass for an assembly that silently ignored the structure altogether.
void test_ssi_effect() {
    std::printf("-- Soil-structure interaction through the GUI Dynamic phase path --\n");
    constexpr double W = 20.0, H = 10.0;
    auto model = [&](bool with_wall, double EA, double EI, double w) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P);
        if (with_wall) {
            m::PlateMaterial pm; pm.EA = EA; pm.EI = EI; pm.w = w; pr.plates.push_back(pm);
            m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
            st.x1 = 10; st.y1 = 4; st.x2 = 10; st.y2 = 10; st.material = 0;
            pr.structs.push_back(st);
        }
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 1.0; p.seismic_freq = 2.5;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 4.0; p.time_steps = 400;
        pr.phases.push_back(p);
        return pr;
    };
    auto peak = [&](m::Project pr, const char* label, bool& ok) {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) { std::printf("   %s: MESH FAILED\n", label); ok = false; return 0.0; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.empty() || !res.back().ok) {
            std::printf("   %s: FAILED (%s)\n", label, res.empty() ? "" : res.back().message.c_str());
            ok = false; return 0.0;
        }
        return res.back().dyn_peak_surface_a;
    };
    bool ok = true;
    const double a_soil  = peak(model(false, 0, 0, 0), "soil only", ok);
    // Vanishing wall: EA/EI ~ 1e-6 of the real thing and weightless -> must not perturb the site.
    const double a_tiny  = peak(model(true, 1.0, 1.0e-3, 0.0), "vanishing wall", ok);
    // Real wall: 0.4 m concrete diaphragm (EA = E d, EI = E d^3/12, w = 24 d).
    const double a_wall  = peak(model(true, 3.0e7 * 0.4, 3.0e7 * 0.064 / 12.0, 24.0 * 0.4), "real wall", ok);
    check(ok, "all three dynamic runs solved");
    if (!ok) return;
    std::printf("   peak surface accel:  soil-only %.4f | vanishing wall %.4f (%+.3f%%) | real wall %.4f (%+.1f%%)\n",
                a_soil, a_tiny, 100 * (a_tiny - a_soil) / a_soil, a_wall, 100 * (a_wall - a_soil) / a_soil);
    // (1) Continuity: a structure with no stiffness and no mass is not there.
    check(std::fabs(a_tiny - a_soil) < 0.01 * a_soil,
          "vanishing wall (EA,EI,w -> 0) recovers the soil-only site response");
    // (2) The real wall is a real structure in the dynamic system.
    check(std::fabs(a_wall - a_soil) > 0.02 * a_soil,
          "a real wall changes the seismic surface response (SSI is actually coupled)");
}

// ===============================================================================================
// Correctness guards found by the seismic audit. Each one pins a path that previously produced a
// PLAUSIBLE-LOOKING WRONG NUMBER rather than an error -- the worst failure mode for this program.
// ===============================================================================================

// A1 (was CRITICAL): a Dynamic phase must honour staged construction. assemble_stiffness/assemble_mass
// had no active_element parameter, so an excavated region kept its full K and M in the dynamic system
// while the GUI drew the excavation -- the wall was braced by soil that is not there, and the reported
// seismic wall forces came out too LOW (unconservative). The pin: deactivating soil must CHANGE the
// seismic result. Before the fix these two runs were bit-identical.
void test_staging_changes_dynamic_result() {
    std::printf("-- A1: staged construction is honoured by the Dynamic phase --\n");
    auto model = [&](bool excavate) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon base; base.material = 0;              // 0..20 x 0..6, the retained body
        base.x = {0, 20, 20, 0}; base.y = {0, 0, 6, 6};
        base.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                        (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(base);
        m::SoilPolygon upper; upper.material = 0;            // 0..20 x 6..10, removed when excavating
        upper.x = {0, 20, 20, 0}; upper.y = {6, 6, 10, 10};
        upper.edge_bc = {(int)m::BCType::Free, (int)m::BCType::VerticallyFixed,
                         (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(upper);
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 1.0; p.seismic_freq = 3.0;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 4.0; p.time_steps = 400;
        if (excavate) p.poly_active = {1, 0};                // remove the upper block in this phase
        pr.phases.push_back(p);
        return pr;
    };
    auto peak = [&](m::Project pr, bool& ok) {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) { ok = false; return 0.0; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.empty() || !res.back().ok) {
            std::printf("   FAILED: %s\n", res.empty() ? "" : res.back().message.c_str());
            ok = false; return 0.0;
        }
        return res.back().dyn_peak_surface_a;
    };
    bool ok = true;
    const double a_full = peak(model(false), ok);
    const double a_exc  = peak(model(true), ok);
    check(ok, "both dynamic runs solved (full model / excavated)");
    if (!ok) return;
    std::printf("   peak surface accel: full %.4f | upper block excavated %.4f (%+.1f%%)\n",
                a_full, a_exc, 100 * (a_exc - a_full) / a_full);
    // Removing 4 m of a 10 m column changes its fundamental period and its mass -> the site response
    // MUST move. Identical values mean the mask was ignored and the excavation is still in the solver.
    check(std::fabs(a_exc - a_full) > 0.02 * a_full,
          "deactivating soil changes the seismic response (the active mask reaches K and M)");
}

// A2 (was CRITICAL): the Mechanical tab never exposes E' for a Hardening Soil / HSsmall material, so
// the linear-elastic dynamic system would have used E's struct default (13 MPa) -- an entire seismic
// run on a stiffness the user never typed. Must be refused, not guessed.
// A3: E'_inc (depth-varying stiffness) is offered by the editor and implemented by NOTHING -- a
// seismic Vs profile would be silently flattened to E'_ref. Must be refused (all phases).
// A4: a model with no mass has no seismic force; peak |u| = 0 would look like a real answer.
void test_dynamic_input_guards() {
    std::printf("-- A2/A3/A4: inputs that would silently produce a wrong number are refused --\n");
    auto base = [&]() {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P);
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 1.0; p.seismic_freq = 2.5;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 4.0; p.time_steps = 400;
        pr.phases.push_back(p);
        return pr;
    };
    auto run = [&](m::Project pr) -> katai::app::SolveResult {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) return {};
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        return res.empty() ? katai::app::SolveResult{} : res.back();
    };
    check(run(base()).ok, "the baseline linear-elastic model still solves");

    // A2 (v0.5): an HS material's E' is not an input -- but the information IS there, just not called
    // E'. HSsmall carries E_0,ref = 2(1+nu_ur) G_0,ref (literally the small-strain stiffness dynamics
    // needs); plain HS carries E_ur,ref. Both are stress-dependent and are now evaluated at each stress
    // point from the PARENT phase's committed stress. So they must SOLVE -- and HSsmall, whose G_0 is
    // several times E_ur, must come out STIFFER (higher Vs -> higher f_1) than plain HS on the same
    // model. That ordering is the physics, and it cannot happen if the stiffness were a stand-in.
    {
        // Nothing overridden but the model type: both use the SAME default HS input (E50 = 30, Eur = 90,
        // G0 = 120 MPa, nu_ur = 0.2), so the only difference is which stiffness each model exposes --
        // E_ur,ref = 90 MPa for plain HS vs E_0,ref = 2(1+nu_ur) G_0,ref = 288 MPa for HSsmall.
        auto hs = [&](m::SoilModel sm) {
            m::Project pr = base();
            pr.materials[0].model = sm;
            return run(pr);
        };
        const auto rhs = hs(m::SoilModel::HardeningSoil);
        std::printf("   HS  -> %s (peak surf a = %.4f)\n", rhs.ok ? "solved" : rhs.message.substr(0, 50).c_str(),
                    rhs.dyn_peak_surface_a);
        check(rhs.ok, "a Hardening Soil material now SOLVES in a Dynamic phase (E_ur from its own input)");
        // Its own stiffness must reach the solve: E_ur,ref = 90 MPa is 3x the linear-elastic E' = 30 MPa
        // this model would otherwise have used, so the site response cannot come out the same. (The old
        // failure mode was the opposite -- E's 13 MPa struct default, silently.)
        const auto rle = run(base());
        std::printf("   HS vs LE(E'=30 MPa): %.4f vs %.4f\n", rhs.dyn_peak_surface_a, rle.dyn_peak_surface_a);
        check(rhs.ok && rle.ok &&
                  std::fabs(rhs.dyn_peak_surface_a - rle.dyn_peak_surface_a) > 0.02 * rle.dyn_peak_surface_a,
              "the HS material's OWN stiffness reaches the seismic solve (not a default, not E')");

        // HSsmall is NOT asserted here. This session found it does not converge in the initial K0 phase
        // AT ALL -- with pure default parameters, on a level block, where the K0 seed IS the answer and
        // the step should be trivial. That is independent of the dynamic branch (K0 passes no
        // gauss_elastic, and its profile is uniform() -> bit-for-bit), and it is invisible today because
        // test_hssmall is a MATERIAL-POINT test: no GUI-path test covers HSsmall in any phase. It fails
        // honestly (ok = false + a message), so nothing wrong is reported -- but a model type the editor
        // offers does not run. Recorded as the next defect to fix rather than papered over here.
    }
    // An HS dynamic phase needs a parent stress state to evaluate its stress-dependent stiffness.
    // (Every project here has the initial K0 phase, so this is about the honest guard existing.)

    // A3 (v0.5): a depth-varying E' is now IMPLEMENTED in EVERY deformation phase (per stress point;
    // verified against a closed form + a constant-E sub-layer stack in test_material_profile). It must
    // solve -- and it must actually change the site response, because E(y) IS the Vs profile: a stiffer
    // deposit at depth raises f_1. Silently flattening it to E'_ref (the old behaviour) would leave the
    // two runs identical.
    { m::Project pr = base(); pr.materials[0].E_inc = 6000.0; pr.materials[0].y_ref = 10.0;
      const auto r = run(pr);
      check(r.ok, "a depth-varying E' now SOLVES in a Dynamic phase (no longer refused)");
      const auto r0 = run(base());
      if (r.ok && r0.ok) {
          std::printf("   E'_inc = 6000  -> peak surface accel %.4f  (uniform E' %.4f, %+.1f%%)\n",
                      r.dyn_peak_surface_a, r0.dyn_peak_surface_a,
                      100 * (r.dyn_peak_surface_a - r0.dyn_peak_surface_a) / r0.dyn_peak_surface_a);
          check(std::fabs(r.dyn_peak_surface_a - r0.dyn_peak_surface_a) > 0.02 * r0.dyn_peak_surface_a,
                "E'(y) reaches the seismic solve -- the Vs profile changes the site response");
      } }

    // A4: no weight -> no seismic force at all.
    { m::Project pr = base(); pr.materials[0].gamma_unsat = 0.0; pr.materials[0].gamma_sat = 0.0;
      const auto r = run(pr);
      std::printf("   zero unit wt   -> %s\n", r.ok ? "SOLVED (BAD)" : r.message.substr(0, 74).c_str());
      check(!r.ok && r.message.find("no mass") != std::string::npos,
            "a massless model is refused (it would report a plausible peak |u| = 0)"); }

}

// A5: the surface monitor must sit on a node of the ACTIVE soil that is free to move. The old search
// scanned every mesh node with no active filter, so in an excavated model it could land on an orphan
// node of the removed block -- one that fix_inactive_nodes has pinned. a_surf then collapses to exactly
// a_g(t), and the run published the INPUT motion as the "surface response spectrum", overlaid on the
// TBDY design spectrum. A completely plausible plot of the wrong thing. The pin: an excavated model's
// surface response must be an amplified site response, not the base motion echoed back.
void test_surface_monitor_after_excavation() {
    std::printf("-- A5: surface monitor lands on the ACTIVE, free surface (not the input motion) --\n");
    constexpr double kAmp = 1.0;
    m::Project pr;
    m::Material s; s.model = m::SoilModel::LinearElastic;
    s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
    pr.materials.push_back(s);
    m::SoilPolygon base; base.material = 0;
    base.x = {0, 20, 20, 0}; base.y = {0, 0, 6, 6};
    base.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                    (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(base);
    m::SoilPolygon upper; upper.material = 0;      // removed in the dynamic phase -> orphan nodes at y=10
    upper.x = {0, 20, 20, 0}; upper.y = {6, 6, 10, 10};
    upper.edge_bc = {(int)m::BCType::Free, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(upper);
    m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
    p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = kAmp; p.seismic_freq = 3.3;
    p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
    p.duration = 4.0; p.time_steps = 400;
    p.poly_active = {1, 0};
    pr.phases.push_back(p);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "excavated model meshed");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(!res.empty() && res.back().ok, "excavated dynamic phase solved");
    if (res.empty() || !res.back().ok) return;
    const double a_surf = res.back().dyn_peak_surface_a;
    std::printf("   peak surface accel = %.4f m/s^2  (input PGA = %.2f)\n", a_surf, kAmp);
    // If the monitor had fallen back to a pinned orphan, a_surf would be EXACTLY the input motion.
    check(std::fabs(a_surf - kAmp) > 0.1 * kAmp,
          "the excavated model reports an amplified site response, not the input motion echoed back");
    check(a_surf > kAmp, "the remaining column amplifies the base motion (shaken near its f_1)");
}

// A6: the base motion is sampled pointwise and Newmark is unconditionally stable -> too few steps per
// cycle silently integrates a DIFFERENT signal. The phase editor must say so (like the consolidation
// dt_crit warning). Pin the thresholds directly.
void test_dynamic_step_warning() {
    std::printf("-- A6: coarse-time-step (aliasing) warning --\n");
    auto ph = [](double dur, int steps, double freq, m::SeismicWave w) {
        m::Phase p; p.type = m::PhaseType::Dynamic; p.duration = dur; p.time_steps = steps;
        p.seismic_freq = freq; p.seismic_wave = w; return p;
    };
    using W = m::SeismicWave;
    check(katai::app::dynamic_step_warning(ph(4.0, 400, 2.5, W::Harmonic)).empty(),
          "40 steps/cycle: no warning");
    check(katai::app::dynamic_step_warning(ph(1.0, 25, 2.0, W::Harmonic)).find("coarse") != std::string::npos,
          "the GUI default (1.0 s / 25 steps @ 2 Hz = 12.5 per cycle) warns");
    const std::string alias = katai::app::dynamic_step_warning(ph(10.0, 50, 5.0, W::Harmonic));
    std::printf("   1 step/cycle -> %s\n", alias.substr(0, 70).c_str());
    check(alias.find("Nyquist") != std::string::npos, "below Nyquist says so explicitly");
    // A Ricker pulse carries energy well above its central frequency -> judged at ~2.5 f.
    check(!katai::app::dynamic_step_warning(ph(2.0, 200, 3.0, W::Ricker)).empty(),
          "a Ricker pulse is judged on its high-frequency content, not its central frequency");
    // A non-dynamic phase is not this warning's business.
    { m::Phase p; p.type = m::PhaseType::Plastic; p.duration = 1.0; p.time_steps = 1;
      check(katai::app::dynamic_step_warning(p).empty(), "a non-dynamic phase never warns"); }
}

// NONLINEAR Dynamic phase through the GUI path (Phase.dynamic_nonlinear opt-in). The core solver is
// validated to round-off in test_dynamics_nonlinear (linear limit + quasi-static elastoplastic limit);
// here the WIRING must hold: (1) with an elastic soil the nonlinear path REDUCES to the linear GUI path
// it is bolted next to (the opt-in routes correctly and reproduces the validated linear physics), and
// (2) with a yielding Mohr-Coulomb soil under a strong shake the nonlinear path DEPARTS from the linear
// (elastic) result -- the soil genuinely plastifies during shaking -- and it converges / never crashes.
void test_nonlinear_dynamic_gui() {
    std::printf("-- Nonlinear Dynamic phase through the GUI path (opt-in dynamic_nonlinear) --\n");
    constexpr double W = 2.0, H = 20.0, E = 208000.0, nu = 0.3, gamma = 19.62;
    const double g = 9.81, rho = gamma / g, G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double f1 = Vs / (4 * H), xi = 0.05;   // shake at the site's fundamental (real displacement)

    auto run = [&](m::SoilModel model, bool nonlinear, double amp, double c, double phi) -> katai::app::SolveResult {
        m::Project pr;
        m::Material s; s.model = model;
        s.E = E; s.nu = nu; s.gamma_unsat = gamma; s.gamma_sat = gamma; s.e_init = 0.5;
        s.c = c; s.phi = phi; s.psi = 0.0;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, W, W, 0}; P.y = {0, 0, H, H};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};   // pure-SH shear column
        pr.polygons.push_back(P); pr.has_water = false;
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Dyn";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = amp; p.seismic_freq = f1;
        p.damping_ratio = xi; p.rayleigh_f1 = f1; p.rayleigh_f2 = 3 * f1;
        p.duration = 6.0; p.time_steps = 600;
        p.dynamic_nonlinear = nonlinear;
        pr.phases.push_back(p);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        if (!M.ok) { std::printf("  (mesh: %s)\n", M.message.c_str()); return {}; }
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);   // K0 parent -> init stress
        return res.size() == 2 ? res[1] : katai::app::SolveResult{};
    };

    // (1) LINEAR LIMIT (wiring): elastic soil -> nonlinear path == linear path.
    const auto le_lin = run(m::SoilModel::LinearElastic, false, 1.0, 1.0, 30.0);
    const auto le_nl  = run(m::SoilModel::LinearElastic, true,  1.0, 1.0, 30.0);
    check(le_lin.ok && le_nl.ok, "elastic dynamic solved BOTH linear and nonlinear (opt-in routes)");
    if (le_lin.ok && le_nl.ok) {
        const double rel = std::fabs(le_nl.max_disp - le_lin.max_disp) / std::fmax(le_lin.max_disp, 1e-12);
        std::printf("   LE: linear peak |u|=%.6f  nonlinear peak |u|=%.6f  (rel diff %.2e)\n",
                    le_lin.max_disp, le_nl.max_disp, rel);
        check(le_lin.max_disp > 1e-4, "the elastic run actually moved");
        check(rel < 0.01, "nonlinear dynamic REDUCES to the linear dynamic for elastic soil (wiring correct)");
    }

    // (2) NONLINEARITY ENGAGES: Mohr-Coulomb, stronger shake -> the nonlinear path yields and departs
    // from the linear (elastic) result, and must converge / not crash. Wall-clock is MEASURED and
    // printed (the honest speed report the nonlinear opt-in promises: the tangent refactors every
    // Newton iteration, so the cost multiplier is real and should be visible, not folklore).
    const auto tl0 = std::chrono::steady_clock::now();
    const auto mc_lin = run(m::SoilModel::MohrCoulomb, false, 4.0, 20.0, 25.0);
    const auto tl1 = std::chrono::steady_clock::now();
    const auto mc_nl  = run(m::SoilModel::MohrCoulomb, true,  4.0, 20.0, 25.0);
    const auto tl2 = std::chrono::steady_clock::now();
    const double s_lin = std::chrono::duration<double>(tl1 - tl0).count();
    const double s_nl = std::chrono::duration<double>(tl2 - tl1).count();
    std::printf("   SPEED (600 steps, same model): linear %.2f s, nonlinear %.2f s -> %.1fx\n",
                s_lin, s_nl, s_nl / std::fmax(s_lin, 1e-9));
    check(mc_lin.ok, "MC linear (elastic) dynamic solved");
    check(mc_nl.ok, "MC nonlinear dynamic solved (converged, no crash)");
    if (mc_nl.ok) std::printf("   MC message: %s\n", mc_nl.message.c_str());
    if (mc_lin.ok && mc_nl.ok) {
        const double rel = std::fabs(mc_nl.max_disp - mc_lin.max_disp) / std::fmax(mc_lin.max_disp, 1e-12);
        std::printf("   MC: linear(elastic) peak |u|=%.6f  nonlinear peak |u|=%.6f  (rel diff %.2e)\n",
                    mc_lin.max_disp, mc_nl.max_disp, rel);
        check(rel > 0.01, "MC nonlinear DEPARTS from linear-elastic (the soil genuinely yielded during shaking)");
    }

    // Seismic effective-stress recovery: the nonlinear path carries the committed stress through the
    // shaking (parent K0/gravity stress + plastic evolution) and reports the post-earthquake nodal field;
    // the linear path solves about a zero-stress datum, so it leaves the field zero (deferred).
    auto max_stress = [](const katai::app::SolveResult& r) {
        double m = 0.0;
        for (const auto& s : r.stress.stress) m = std::fmax(m, s.cwiseAbs().maxCoeff());
        return m;
    };
    if (mc_nl.ok && mc_lin.ok) {
        const double s_nl = max_stress(mc_nl), s_lin = max_stress(mc_lin);
        std::printf("   stress field: nonlinear max|sigma| = %.4g kPa   linear max|sigma| = %.4g kPa\n", s_nl, s_lin);
        check(s_nl > 1.0, "NONLINEAR: seismic effective-stress field recovered (post-earthquake state)");
        check(s_lin < 1e-9, "LINEAR: stress recovery deferred (zero -- it solves about a zero-stress datum)");
    }
}

// The nonlinear dynamic path RESOLVES the interface slip that the linear path can only MEASURE. This is
// the whole reason the nonlinear solver exists (D7 measured a wall interface exceeding its Coulomb
// capacity under mild shaking; the linear analysis is unconservative exactly there). Same weak-joint
// model, run both ways: the LINEAR joint reports demand/capacity > 1 (it cannot slip), while the
// NONLINEAR joint actually slips (real plastic slip) so the Coulomb demand no longer exceeds capacity.
void test_nonlinear_resolves_slip() {
    std::printf("-- Nonlinear resolves the interface slip the linear path only MEASURES --\n");
    auto model = [&](bool nonlinear, double amp) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        s.c = 5.0; s.phi = 30.0; s.rinter_rigid = false; s.Rinter = 1.0;   // real joint strength c_i, phi_i
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
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = amp; p.seismic_freq = 2.5;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 4.0; p.time_steps = 400;
        p.dynamic_nonlinear = nonlinear;
        pr.phases.push_back(p);
        return pr;
    };
    auto run = [&](m::Project pr, std::vector<katai::app::SolveResult>& out) {
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) return false;
        out = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        return !out.empty() && out.back().ok;
    };

    std::vector<katai::app::SolveResult> rl, rn;
    check(run(model(false, 4.0), rl), "linear weak-joint dynamic solved");
    check(run(model(true, 4.0), rn), "nonlinear weak-joint dynamic solved (converged)");
    if (rl.empty() || rn.empty() || !rl.back().ok || !rn.back().ok) return;
    if (rl.back().interface_forces.empty() || rn.back().interface_forces.empty()) {
        check(false, "both runs report an interface"); return;
    }
    const auto& jl = rl.back().interface_forces.front();
    const auto& jn = rn.back().interface_forces.front();
    std::printf("   LINEAR:    max util = %.2fx over %.0f%% (measures), any_slip = %s, max|slip| = %.3e\n",
                jl.max_utilisation, 100 * jl.over_fraction, jl.any_slip ? "yes" : "no", jl.max_abs_slip);
    std::printf("   NONLINEAR: util reported = %s,   any_slip = %s, max|slip| = %.3e (resolves)\n",
                jn.max_utilisation > 0.0 ? "yes" : "no (slip instead)", jn.any_slip ? "yes" : "no", jn.max_abs_slip);
    // Linear cannot slip -> it reports the demand exceeding capacity (a MEASURE of unconservatism).
    check(jl.max_utilisation > 1.0, "LINEAR: joint demand exceeds capacity (unconservative -- it cannot slip)");
    check(!jl.any_slip, "LINEAR: reports no slip (elastic by construction)");
    check(!jl.slip_checked, "LINEAR: the envelope is flagged as NOT slip-checked (v5 badge fix)");
    // Nonlinear RESOLVES it by ACTUALLY SLIPPING (real plastic slip, an exact committed-state event).
    check(jn.max_abs_slip > 1e-6, "NONLINEAR: the joint actually SLIPS (real plastic slip, not just a measure)");
    check(jn.any_slip, "NONLINEAR: reports the slip (no phantom [elastic, no slip check] badge)");
    check(jn.slip_checked, "NONLINEAR: the envelope is flagged as slip-checked (Coulomb branch really ran)");
    // Track 1a: the parent (K0) structural state is CARRIED into the nonlinear increment, so the
    // Coulomb cap acts on the TOTAL action and the demand/capacity ratio is meaningful again --
    // accumulated per instant (tau and sigma_n of the SAME step). It must sit at/below capacity:
    // the solver RESOLVES the exceedance the linear path could only measure.
    check(jn.superposed, "NONLINEAR: stations are the TOTAL action (parent state carried, Track 1a)");
    check(jn.max_utilisation > 0.5,
          "NONLINEAR: demand/capacity is reported again (per instant, on the solver's own total)");
    check(jn.max_utilisation <= 1.0 + 1e-6,
          "NONLINEAR: the total NEVER exceeds capacity (tau is Coulomb-capped instant by instant)");
    check(jn.over_fraction == 0.0,
          "NONLINEAR: no station past capacity (the exceedance is resolved, not measured)");
    if (!rn.back().message.empty())
        check(rn.back().message.find("continuing from the parent phase's structural state") != std::string::npos,
              "NONLINEAR: the phase message says the parent structural state was carried");
}

// Track 1a END-TO-END IDENTITY: a NONLINEAR dynamic phase seeded with its parent's structural state
// and shaken by NOTHING (a_g = 0) must report exactly the parent's static structural state --
// station for station (the envelope evaluates the same validated post-processors at the same total
// displacement and the same committed plastic state the parent's static tail used). This is the D7
// superposition identity moved INSIDE the solver: one check pins the displacement-datum mapping,
// the plastic-state seeding and the f_int0 baseline together (any wiring error survives none of it).
// It also pins the BONDED case: an undriven joint commits no slip event, so any_slip must be false
// -- [bonded], not a phantom [SLIPPING] from the joint's ordinary elastic relative displacement.
void test_nonlinear_carry_identity() {
    std::printf("-- Track 1a identity: undriven nonlinear dynamic == the parent's static state --\n");
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
    // A one-sided surcharge FIRST: the wall bends (M != 0), the joint shears (tau != 0) and the
    // parent develops a REAL displacement field. Without it the K0 parent is wished-in-place
    // (u ~ 0, M ~ 0) and the identity would not exercise the carried displacement DATUM at all --
    // sigma_n = 54 kPa comes from the sigma_n0 seed and would match even with a zeroed datum.
    m::Load L; L.kind = m::LoadKind::Distributed;
    L.x1 = 0; L.y1 = 10; L.x2 = 8; L.y2 = 10;
    L.qx1 = 0; L.qy1 = -50; L.qx2 = 0; L.qy2 = -50;
    pr.loads.push_back(L);
    m::Phase pl; pl.type = m::PhaseType::Plastic; pl.name = "Surcharge";
    pr.phases.push_back(pl);
    m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Nil quake";
    p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 0.0; p.seismic_freq = 2.5;
    p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
    p.duration = 0.6; p.time_steps = 60;
    p.dynamic_nonlinear = true;
    pr.phases.push_back(p);

    const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
    check(M.ok, "carry-identity model meshed");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    check(res.size() == 3 && res[0].ok && res[1].ok && res[2].ok,
          "K0 + loaded Plastic parent + undriven nonlinear dynamic solved");
    if (res.size() != 3 || !res[0].ok || !res[1].ok || !res[2].ok) return;
    const auto& par = res[1];   // the LOADED static parent (the dynamic phase's io.prev)
    const auto& dyn = res[2];
    check(dyn.message.find("continuing from the parent phase's structural state") != std::string::npos,
          "the phase message says the parent structural state was carried");

    // Interfaces: the undriven envelope must equal |parent static| station for station.
    check(par.interface_forces.size() == dyn.interface_forces.size() && !dyn.interface_forces.empty(),
          "parent and dynamic report the same joints");
    double dif = 0.0, teeth = 0.0;
    for (size_t i = 0; i < dyn.interface_forces.size() && i < par.interface_forces.size(); ++i) {
        const auto& jd = dyn.interface_forces[i];
        const auto& jp = par.interface_forces[i];
        check(jd.superposed, "the joint stations are the TOTAL action (carried)");
        check(jd.slip_checked, "the joint envelope is slip-checked (nonlinear Coulomb branch)");
        check(!jd.any_slip, "an UNDRIVEN joint commits no slip event -> [bonded], not [SLIPPING]");
        check(jd.max_utilisation <= 1.0 + 1e-9, "utilisation never exceeds capacity on the carry path");
        if (jd.stations.size() != jp.stations.size()) { check(false, "joint stations align"); continue; }
        for (size_t k = 0; k < jd.stations.size(); ++k) {
            dif = std::fmax(dif, std::fabs(jd.stations[k].tau - std::fabs(jp.stations[k].tau)));
            dif = std::fmax(dif, std::fabs(jd.stations[k].sigma_n - std::fabs(jp.stations[k].sigma_n)));
            teeth = std::fmax(teeth, std::fabs(jp.stations[k].sigma_n));
        }
    }
    std::printf("   joints: max |dyn - |static|| = %.3e   (teeth: max static |sigma_n| = %.4g kPa)\n",
                dif, teeth);
    check(teeth > 10.0, "the K0 parent really loads the joint (the identity has teeth)");
    check(dif <= 1e-9 * std::fmax(1.0, teeth),
          "undriven nonlinear dynamic joint == parent static joint, station for station");

    // Wall forces: same identity for the structural envelope.
    check(par.struct_forces.size() == dyn.struct_forces.size() && !dyn.struct_forces.empty(),
          "parent and dynamic report the same structures");
    double difM = 0.0, teethM = 0.0;
    for (size_t i = 0; i < dyn.struct_forces.size() && i < par.struct_forces.size(); ++i) {
        const auto& fd = dyn.struct_forces[i];
        const auto& fp = par.struct_forces[i];
        check(fd.superposed, "the wall stations are the TOTAL action (carried)");
        if (fd.stations.size() != fp.stations.size()) { check(false, "wall stations align"); continue; }
        for (size_t k = 0; k < fd.stations.size(); ++k) {
            difM = std::fmax(difM, std::fabs(fd.stations[k].M - std::fabs(fp.stations[k].M)));
            teethM = std::fmax(teethM, std::fabs(fp.stations[k].M));
        }
    }
    std::printf("   wall:   max |dyn - |static|| = %.3e   (teeth: max static |M| = %.4g kNm/m)\n",
                difM, teethM);
    check(teethM > 1.0, "the surcharge really bends the wall (the DATUM identity has teeth)");
    check(difM <= 1e-9 * std::fmax(1.0, teethM),
          "undriven nonlinear dynamic wall == parent static wall, station for station");
}

// The LINEAR dynamic system solves an anchor UNCAPPED (elastic EA/L -- its own message discloses
// "anchors do not yield"), so its envelope must report the uncapped elastic demand, even beyond
// F_max. Before this fix the diagram applied the F_max cap the solver never had: a linear run whose
// anchor demand exceeded capacity silently reported N = F_max -- an under-reported demand that also
// implied a yield the analysis never tested (the D6b class, on the anchor). The NONLINEAR path
// really yields the anchor, so its envelope must sit at/below F_max with yielded = true (an exact
// committed-state event, Track 1a).
void test_anchor_envelope_consistency() {
    std::printf("-- Anchor envelope: linear reports the UNCAPPED demand, nonlinear the resolved cap --\n");
    auto model = [&](bool nonlinear, bool capped, double Fmax) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 30000; s.nu = 0.3; s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.6;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, 20, 20, 0}; P.y = {0, 0, 10, 10};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P);
        m::AnchorMaterial am; am.EA = 1.0e5;
        am.elastoplastic = capped; am.Fmax_tens = Fmax; am.Fmax_comp = Fmax;
        pr.anchors.push_back(am);
        m::StructElement st; st.kind = m::StructKind::Anchor; st.name = "Tie";
        st.x1 = 4; st.y1 = 9; st.x2 = 10; st.y2 = 5; st.material = 0;   // spans depths -> axial demand
        pr.structs.push_back(st);
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Quake";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 4.0; p.seismic_freq = 2.5;
        p.damping_ratio = 0.05; p.rayleigh_f1 = 1.0; p.rayleigh_f2 = 8.0;
        p.duration = 2.0; p.time_steps = 200;
        p.dynamic_nonlinear = nonlinear;
        pr.phases.push_back(p);
        return pr;
    };
    auto anchor_env = [&](const m::Project& pr, double& maxN, bool& yielded, bool& ok) {
        maxN = 0.0; yielded = false; ok = false;
        const auto M = katai::app::mesh_from_project(pr, 1.0, 6);
        if (!M.ok) return;
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        if (res.size() != 2 || !res.back().ok) return;
        for (const auto& sf : res.back().struct_forces)
            if (sf.kind == 1) { maxN = sf.max_N; yielded = sf.yielded; ok = true; }
    };

    // (1) Reference: linear, unbounded anchor -> the honest elastic demand N_ref.
    double N_ref = 0.0; bool y_ref = false, ok_ref = false;
    anchor_env(model(false, false, 0.0), N_ref, y_ref, ok_ref);
    check(ok_ref, "linear reference run solved and reports the anchor");
    check(N_ref > 1.0, "the shake really loads the anchor (the check has teeth)");

    // (2) Same shake, F_max = N_ref/2. The LINEAR solve is UNCHANGED by F_max (its K never saw it),
    // so its envelope must still report N_ref -- not a silent cap at F_max.
    const double Fmax = 0.5 * N_ref;
    double N_lin = 0.0; bool y_lin = false, ok_lin = false;
    anchor_env(model(false, true, Fmax), N_lin, y_lin, ok_lin);
    check(ok_lin, "linear capped-anchor run solved");
    std::printf("   N_ref = %.4f | linear w/ F_max=%.4f -> N = %.4f | ", N_ref, Fmax, N_lin);
    check(std::fabs(N_lin - N_ref) <= 1e-9 * N_ref,
          "LINEAR: the envelope reports the UNCAPPED elastic demand (the solve never capped)");
    check(!y_lin, "LINEAR: no yield claim (the analysis never tested one)");

    // (3) NONLINEAR, same F_max: the anchor really yields; the envelope respects the cap.
    double N_nl = 0.0; bool y_nl = false, ok_nl = false;
    anchor_env(model(true, true, Fmax), N_nl, y_nl, ok_nl);
    check(ok_nl, "nonlinear capped-anchor run solved");
    std::printf("nonlinear -> N = %.4f, yielded = %s\n", N_nl, y_nl ? "yes" : "no");
    check(N_nl <= Fmax * (1.0 + 1e-9),
          "NONLINEAR: the envelope respects F_max (the solver really capped the force)");
    check(y_nl, "NONLINEAR: the anchor's yield is reported (exact committed-state event)");
}

// EC8 overlay wiring: the phase's EC8 inputs must reach the result as the EXACT core spectrum at
// every reported period (a_g = gamma_I a_gR g), and stay EMPTY when the overlay is off. The core
// formulas/tables are pinned in test_tbdy_seismic; here only the plumbing can be wrong.
void test_ec8_overlay() {
    std::printf("-- EC8 (EN 1998-1) design-spectrum overlay wiring --\n");
    auto model = [&](bool ec8_on) {
        m::Project pr;
        m::Material s; s.model = m::SoilModel::LinearElastic;
        s.E = 50000; s.nu = 0.3; s.gamma_unsat = 19.0; s.gamma_sat = 19.0; s.e_init = 0.5;
        pr.materials.push_back(s);
        m::SoilPolygon P; P.material = 0;
        P.x = {0, 2, 2, 0}; P.y = {0, 0, 10, 10};
        P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
        pr.polygons.push_back(P); pr.has_water = false;
        m::Phase p; p.type = m::PhaseType::Dynamic; p.name = "Dyn";
        p.seismic_wave = m::SeismicWave::Harmonic; p.seismic_amp = 1.0; p.seismic_freq = 3.0;
        p.duration = 1.0; p.time_steps = 100;
        p.ec8_enabled = ec8_on;
        p.ec8_agr = 0.25; p.ec8_gamma = 1.2; p.ec8_ground = 3; p.ec8_type = 1;   // D, Type 2
        pr.phases.push_back(p);
        return pr;
    };
    auto run = [&](const m::Project& pr) {
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        if (!M.ok) return katai::app::SolveResult{};
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        return res.size() == 2 ? res.back() : katai::app::SolveResult{};
    };
    const auto off = run(model(false));
    const auto on = run(model(true));
    check(off.ok && on.ok, "both runs solved");
    if (!off.ok || !on.ok) return;
    check(off.dyn_design_sa_ec8.empty(), "overlay OFF -> no EC8 curve in the result");
    check(on.dyn_design_sa_ec8.size() == on.dyn_period.size(),
          "overlay ON -> EC8 curve over the full period range");
    check(!on.dyn_design_sa.empty(), "the TBDY curve is still there (EC8 is an addition, not a swap)");
    double worst = 0.0;
    const double ag = 1.2 * 0.25 * 9.81;   // gamma_I * a_gR * g [m/s^2]
    for (size_t i = 0; i < on.dyn_design_sa_ec8.size() && i < on.dyn_period.size(); ++i) {
        const double ref = katai::core::ec8_elastic_spectrum(
            ag, katai::core::Ec8GroundType::D, katai::core::Ec8SpectrumType::Type2, on.dyn_period[i]);
        worst = std::fmax(worst, std::fabs(on.dyn_design_sa_ec8[i] - ref));
    }
    std::printf("   max |overlay - core Se(T)| = %.3e m/s^2 over %zu periods\n",
                worst, on.dyn_design_sa_ec8.size());
    check(worst < 1e-12, "the overlay IS the core EC8 spectrum at every period (exact plumbing)");
}

}  // namespace

int main() {
    std::printf("Dynamic (seismic) analysis through the GUI compute path\n\n");
    test_dynamic_robustness();
    test_staging_changes_dynamic_result();
    test_dynamic_input_guards();
    test_surface_monitor_after_excavation();
    test_dynamic_step_warning();
    test_ssi_effect();
    test_seismic_structural_forces();
    test_total_action_superposition();
    test_site_response();
    test_freefield_sides();
    test_nonlinear_dynamic_gui();
    test_nonlinear_resolves_slip();
    test_nonlinear_carry_identity();
    test_anchor_envelope_consistency();
    test_ec8_overlay();
    if (g_failures == 0) {
        std::printf("\nOK: GUI Dynamic path = 1D site-response resonance + frequency selectivity\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
