// INPUT-AUDIT pins (2026-07 audit): the fixes for silently-ignored or silently-wrong user
// inputs, verified through the FULL app compute path.
//   (a) OCR raises the automatic K0 (PLAXIS Ref: K0 = K0nc*OCR - nu/(1-nu)*(OCR-1)); it was
//       silently ignored for MC/LE geostatics before -- pinned against the hand formula;
//   (b) NC control: OCR = 1 leaves K0 = 1 - sin(phi) untouched;
//   (c) Hardening Soil + Undrained (B) is REFUSED (it silently behaved like Undrained (A));
//   (d) plate Mp/Np is a FEATURE (M-N hinge): a real capacity SOLVES; the GUI-default trap
//       (elastoplastic checked, both capacities 0 = unlimited) stays an honest refusal.
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/analysis/structural_dynamics.hpp>   // assemble_structural_weight (direct engine use)
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

constexpr double kH = 8.0, kW = 4.0;

m::Project mc_column(double ocr) {
    m::Project pr;
    m::Material s;
    s.name = "MC";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 2.0e4; s.nu = 0.3; s.c = 5.0; s.phi = 30.0; s.psi = 0.0;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    s.k0_auto = true;
    if (ocr > 1.0) { s.oc_mode = 1; s.OCR = ocr; }
    pr.materials.push_back(s);
    m::SoilPolygon L; L.material = 0;
    L.x = {0, kW, kW, 0}; L.y = {0, 0, kH, kH};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(L);
    pr.has_water = false;
    return pr;
}

// K0 of the geostatic field at mid-depth, from the recovered nodal stresses.
double measured_k0(const katai::app::SolveResult& R) {
    double sxx = 0.0, syy = 0.0; int n = 0;
    for (int i = 0; i < R.mesh.node_count; ++i)
        if (std::fabs(R.mesh.y[i] - kH / 2) < 0.26 && R.mesh.x[i] > 0.9 && R.mesh.x[i] < kW - 0.9) {
            sxx += R.stress.stress[i](0); syy += R.stress.stress[i](1); ++n;
        }
    return n > 0 && syy != 0.0 ? sxx / syy : 0.0;
}

void test_ocr_k0() {
    std::printf("-- (a,b) automatic K0 with OCR (PLAXIS formula), NC control --\n");
    for (double ocr : {1.0, 2.0}) {
        m::Project pr = mc_column(ocr);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 1 && res[0].ok;
        check(ok, "K0 phase ran");
        if (!ok) { if (!res.empty()) std::printf("   (%s)\n", res[0].message.c_str()); continue; }
        const double k0nc = 1.0 - std::sin(30.0 * 3.14159265358979323846 / 180.0);   // 0.5
        const double k0_ex = ocr > 1.0 ? k0nc * ocr - 0.3 / 0.7 * (ocr - 1.0) : k0nc;
        const double k0_fe = measured_k0(res[0]);
        std::printf("   OCR = %.0f: K0 = %.4f (formula %.4f, err %+.2f%%)\n", ocr, k0_fe, k0_ex,
                    100.0 * (k0_fe - k0_ex) / k0_ex);
        check(std::fabs(k0_fe - k0_ex) < 0.02 * k0_ex,
              ocr > 1.0 ? "OCR = 2 raises K0 per K0nc*OCR - nu/(1-nu)*(OCR-1) (2%)"
                        : "OCR = 1 leaves K0 = 1 - sin(phi) (2%)");
    }
}

void test_hs_undrained_b_refusal() {
    std::printf("-- (c) Hardening Soil + Undrained (B): honest refusal --\n");
    m::Project pr = mc_column(1.0);
    pr.materials[0].model = m::SoilModel::HardeningSoil;
    pr.materials[0].drainage = m::Drainage::UndrainedB;
    const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    bool refused = false;
    for (const auto& r : res)
        if (!r.ok && r.message.find("Undrained (B)") != std::string::npos) refused = true;
    check(refused, "refused with an explicit Hardening Soil + Undrained (B) message");
}

// (d) Plate Mp/Np is a FEATURE now (the M-N hinge, structural-plate-formulation.md §10) --
// an elastoplastic plate with a real capacity SOLVES; the remaining honest gate is the
// GUI-default trap: elastoplastic checked with BOTH capacities 0 promises a hinge the run
// would never form (unlimited-elastic), so that exact combination is refused.
void test_plate_elastoplastic_gate() {
    std::printf("-- (d) plate Mp/Np: capacity solves, zero-capacity elastoplastic refused --\n");
    auto with_plate = [](double Mp) {
        m::Project pr = mc_column(1.0);
        m::PlateMaterial pm; pm.name = "Wall"; pm.elastoplastic = true; pm.Mp = Mp;
        pr.plates.push_back(pm);
        m::StructElement st; st.kind = m::StructKind::Plate; st.name = "Wall";
        st.material = 0; st.x1 = 2.0; st.y1 = kH; st.x2 = 2.0; st.y2 = 3.0;
        pr.structs.push_back(st);
        return pr;
    };
    {   // real capacity: the phase runs (the hinge machinery is live, not a refusal)
        m::Project pr = with_plate(100.0);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool solved = !res.empty() && res.back().ok;
        if (!solved && !res.empty()) std::printf("   (%s)\n", res.back().message.c_str());
        check(solved, "elastoplastic plate with Mp > 0 solves (plate hinge is a feature now)");
    }
    {   // GUI-default trap: elastoplastic + Mp = Np = 0 -> honest refusal
        m::Project pr = with_plate(0.0);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool refused = false;
        for (const auto& r : res)
            if (!r.ok && r.message.find("elastoplastic") != std::string::npos) refused = true;
        check(refused, "elastoplastic with Mp = Np = 0 refused (would silently run elastic)");
    }
}

// (e) MC tension cut-off IS applied to soil elements now (audit item #1 -- the last
// big silent-wrong: the box was checkable but soil carried tension unclipped). The
// canonical scenario is the one PLAXIS itself cites (MMM 3.3.10, tensile cracks near
// a trench in clay): excavating a deep vertical cut in cohesive soil relieves the
// horizontal stress behind the crest into TENSION. With the cut-off ON (default,
// sigma_t = 0) the recovered major principal stress stays at ~0 there; with it OFF
// the same model carries real tension (MC alone allows sigma_1 up to
// 2c cos(phi)/(1+sin(phi)) = 28.9 kPa at zero confinement). Both runs go through
// the FULL app path (build_problem -> solver -> recovery): this pins the wiring.
void test_soil_tension_cutoff_applied() {
    std::printf("-- (e) MC tension cut-off active in soil elements (full path) --\n");
    auto max_s1_at_crest = [&](bool cutoff, bool& ok) {
        m::Project pr;
        pr.x_min = 0.0; pr.x_max = 8.0; pr.y_min = 0.0; pr.y_max = 8.0;
        pr.has_water = false;
        m::Material s;
        s.name = "clay";
        s.model = m::SoilModel::MohrCoulomb;
        s.E = 2.0e4; s.nu = 0.3; s.c = 25.0; s.phi = 30.0; s.psi = 0.0;
        s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
        s.k0_auto = true;
        s.tension_cutoff = cutoff;
        s.tensile_strength = 0.0;
        pr.materials.push_back(s);
        const int F = (int)m::BCType::FullyFixed, H = (int)m::BCType::HorizontallyFixed,
                  N = (int)m::BCType::Free;
        m::SoilPolygon slab;   // lower slab, full width
        slab.material = 0;
        slab.x = {0, 8, 8, 0}; slab.y = {0, 0, 3, 3};
        slab.edge_bc = {F, H, N, H};
        m::SoilPolygon cut;    // upper-left block: excavated in the phase (5 m cut)
        cut.material = 0;
        cut.x = {0, 3, 3, 0}; cut.y = {3, 3, 8, 8};
        cut.edge_bc = {N, N, N, H};
        m::SoilPolygon keep;   // upper-right block: stays; its crest goes tensile
        keep.material = 0;
        keep.x = {3, 8, 8, 3}; keep.y = {3, 3, 8, 8};
        keep.edge_bc = {N, H, N, N};
        pr.polygons.push_back(slab);
        pr.polygons.push_back(cut);
        pr.polygons.push_back(keep);
        m::Phase ph;
        ph.name = "excavate";
        ph.type = m::PhaseType::Plastic;
        ph.poly_active = {1, 0, 1};
        pr.phases.push_back(ph);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        ok = res.size() == 2 && res[0].ok && res[1].ok;
        if (!ok) {
            for (const auto& r : res)
                if (!r.ok) std::printf("   phase failed: %s\n", r.message.c_str());
            return 0.0;
        }
        const auto& R = res[1];
        double s1max = -1e30;
        for (int i = 0; i < R.mesh.node_count; ++i) {
            // Crest region behind the cut face (the retained block's top corner).
            if (R.mesh.y[i] < 6.0 || R.mesh.x[i] < 3.0 || R.mesh.x[i] > 6.5) continue;
            const auto& s = R.stress.stress[i];
            const double mean = 0.5 * (s(0) + s(1));
            const double rad =
                std::sqrt(0.25 * (s(0) - s(1)) * (s(0) - s(1)) + s(2) * s(2));
            s1max = std::max(s1max, mean + rad);
        }
        return s1max;
    };
    bool ok_on = false, ok_off = false;
    const double s1_on = max_s1_at_crest(true, ok_on);
    const double s1_off = max_s1_at_crest(false, ok_off);
    check(ok_on && ok_off, "both excavation runs solved");
    if (ok_on && ok_off) {
        std::printf("   max sigma_1 behind the crest: cut-off ON %.3f kPa, OFF %.3f kPa\n",
                    s1_on, s1_off);
        check(s1_on <= 1.5,
              "cut-off ON: recovered sigma_1 held at ~sigma_t = 0 (recovery slack)");
        check(s1_off >= s1_on + 2.0 && s1_off >= 3.0,
              "cut-off OFF: the same excavation carries real tension (witness)");
    }
}

// (f) Structural self-weight consistent forces (audit item #2: plate w / pile gamma*A
// fed only the dynamic mass -- statics silently carried weightless walls and piles).
// Direct checks of assemble_structural_weight against hand integrals: on a straight
// uniform 3-node element the consistent load is EXACTLY Simpson (L/6, L/6, 2L/3); on
// the 5-node quartic element EXACTLY Boole (7,32,12,32,7)/90 * L; obliquity leaves the
// force vertical and the total = -w*L; independent (embedded-wall) translation DOFs
// receive the load instead of the mesh nodes; the embedded beam loads its OWN DOFs.
void test_struct_weight_consistent_forces() {
    std::printf("-- (f) structural self-weight: consistent nodal forces (hand integrals) --\n");
    namespace kc = katai::core;
    const double g = 9.81, w = 8.3;

    auto close = [](double a, double b) {
        return std::fabs(a - b) <= 1e-11 * (1.0 + std::fabs(b));
    };

    // 3-node vertical plate, L = 4: Simpson distribution on mesh uy DOFs.
    {
        katai::mesh::Mesh msh;
        msh.node_count = 3;
        msh.x = {0.0, 0.0, 0.0};
        msh.y = {0.0, 4.0, 2.0};          // [A, B, mid]
        kc::DofMap dofs(3);
        kc::PlateElement pe;
        pe.nodes = {0, 1, 2};
        for (int k = 0; k < 3; ++k) pe.rot_dof[k] = dofs.add_extra_dof();
        pe.props.rho_A = w / g;
        dofs.finalize();
        kc::Structures S;
        S.plates.push_back(pe);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        kc::assemble_structural_weight(msh, dofs, S, g, f);
        const double L = 4.0;
        check(close(f(dofs.equation(dofs.global_dof(0, 1))), -w * L / 6.0) &&
                  close(f(dofs.equation(dofs.global_dof(1, 1))), -w * L / 6.0) &&
                  close(f(dofs.equation(dofs.global_dof(2, 1))), -w * 2.0 * L / 3.0),
              "3-node plate: exact Simpson (L/6, L/6, 2L/3) on uy");
        check(f(dofs.equation(dofs.global_dof(0, 0))) == 0.0 &&
                  f(dofs.equation(pe.rot_dof[0])) == 0.0,
              "3-node plate: no ux / rotation load (gravity is vertical)");
    }
    // Oblique 3-node plate (3-4-5), L = 5: still vertical, total = -w*L, Simpson holds.
    {
        katai::mesh::Mesh msh;
        msh.node_count = 3;
        msh.x = {0.0, 3.0, 1.5};
        msh.y = {0.0, 4.0, 2.0};
        kc::DofMap dofs(3);
        kc::PlateElement pe;
        pe.nodes = {0, 1, 2};
        for (int k = 0; k < 3; ++k) pe.rot_dof[k] = dofs.add_extra_dof();
        pe.props.rho_A = w / g;
        dofs.finalize();
        kc::Structures S;
        S.plates.push_back(pe);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        kc::assemble_structural_weight(msh, dofs, S, g, f);
        double sum_y = 0.0, sum_x = 0.0;
        for (int n = 0; n < 3; ++n) {
            sum_y += f(dofs.equation(dofs.global_dof(n, 1)));
            sum_x += f(dofs.equation(dofs.global_dof(n, 0)));
        }
        check(close(sum_y, -w * 5.0) && sum_x == 0.0,
              "oblique plate: total = -w*L, purely vertical");
        check(close(f(dofs.equation(dofs.global_dof(2, 1))), -w * 2.0 * 5.0 / 3.0),
              "oblique plate: Simpson mid-node share holds (straight uniform)");
    }
    // Independent translation DOFs (embedded-wall barrier): the load lands THERE.
    {
        katai::mesh::Mesh msh;
        msh.node_count = 3;
        msh.x = {0.0, 0.0, 0.0};
        msh.y = {0.0, 4.0, 2.0};
        kc::DofMap dofs(3);
        kc::PlateElement pe;
        pe.nodes = {0, 1, 2};
        for (int k = 0; k < 3; ++k) pe.rot_dof[k] = dofs.add_extra_dof();
        for (int k = 0; k < 6; ++k) pe.trans_dof[k] = dofs.add_extra_dof();
        pe.props.rho_A = w / g;
        dofs.finalize();
        kc::Structures S;
        S.plates.push_back(pe);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        kc::assemble_structural_weight(msh, dofs, S, g, f);
        check(close(f(dofs.equation(pe.trans_dof[1])), -w * 4.0 / 6.0) &&
                  close(f(dofs.equation(pe.trans_dof[5])), -w * 2.0 * 4.0 / 3.0) &&
                  f(dofs.equation(dofs.global_dof(0, 1))) == 0.0,
              "independent trans DOFs (embedded wall): load on the wall, not the mesh");
    }
    // 5-node quartic plate, horizontal, L = 2: exact Boole (7,32,12,32,7)/90 * L.
    {
        katai::mesh::Mesh msh;
        msh.node_count = 5;
        msh.x = {0.0, 0.5, 1.0, 1.5, 2.0};   // kXi5 order: xi = -1,-0.5,0,0.5,1
        msh.y = {0.0, 0.0, 0.0, 0.0, 0.0};
        kc::DofMap dofs(5);
        kc::PlateElement5 pe;
        pe.nodes = {0, 1, 2, 3, 4};
        for (int k = 0; k < 5; ++k) pe.rot_dof[k] = dofs.add_extra_dof();
        pe.props.rho_A = w / g;
        dofs.finalize();
        kc::Structures S;
        S.plates5.push_back(pe);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        kc::assemble_structural_weight(msh, dofs, S, g, f);
        const double L = 2.0, boole[5] = {7.0, 32.0, 12.0, 32.0, 7.0};
        bool ok = true;
        for (int k = 0; k < 5; ++k)
            ok = ok && close(f(dofs.equation(dofs.global_dof(k, 1))),
                             -w * L * boole[k] / 90.0);
        check(ok, "5-node plate: exact Boole (7,32,12,32,7)/90 * L on uy");
    }
    // Embedded beam: its OWN nodes and extra DOFs, Simpson again, total = -w_row*L.
    {
        katai::mesh::Mesh msh;
        msh.node_count = 1;                  // the beam never touches mesh nodes here
        msh.x = {0.0};
        msh.y = {0.0};
        kc::DofMap dofs(1);
        katai::core::ebeam::EmbeddedBeam eb;
        eb.node_x = {0.0, 0.0, 0.0};
        eb.node_y = {0.0, 1.5, 3.0};         // bottom -> top, mid at node 1
        for (int i = 0; i < 3; ++i) {
            eb.dof_x.push_back(dofs.add_extra_dof());
            eb.dof_y.push_back(dofs.add_extra_dof());
            eb.dof_phi.push_back(dofs.add_extra_dof());
        }
        eb.elements.push_back({0, 2, 1});    // [A, B, mid] convention
        const double w_row = 2.0;
        eb.props.rho_A = w_row / g;
        dofs.finalize();
        kc::Structures S;
        S.embedded_beams.push_back(eb);
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        kc::assemble_structural_weight(msh, dofs, S, g, f);
        const double L = 3.0;
        check(close(f(dofs.equation(eb.dof_y[0])), -w_row * L / 6.0) &&
                  close(f(dofs.equation(eb.dof_y[2])), -w_row * L / 6.0) &&
                  close(f(dofs.equation(eb.dof_y[1])), -w_row * 2.0 * L / 3.0),
              "embedded beam: Simpson on the beam's own uy DOFs");
    }
}

// (g) Wall self-weight through the FULL app path: a wished-in-place wall in a level
// K0 initial phase. w = 0 keeps the undisturbed identity (max|u| = 0 class); w = 8.3
// makes the wall settle under its own weight (witness that the load is live); a
// chained NIL phase stays a true no-op (the parent's committed state balances the
// weight -- the f/B single-counting contract of the static chain).
void test_wall_weight_static() {
    std::printf("-- (g) wall self-weight in statics: witness + nil identity (full path) --\n");
    auto run = [&](double w_plate) {
        m::Project pr = mc_column(1.0);
        m::PlateMaterial pm;
        pm.name = "Wall";
        pm.EA = 5.0e6;
        pm.EI = 8.5e3;
        pm.w = w_plate;
        pr.plates.push_back(pm);
        m::StructElement st;
        st.kind = m::StructKind::Plate;
        st.name = "Wall";
        st.material = 0;
        st.x1 = 2.0; st.y1 = kH; st.x2 = 2.0; st.y2 = 3.0;
        pr.structs.push_back(st);
        m::Phase nil;
        nil.name = "nil";
        nil.type = m::PhaseType::Plastic;
        pr.phases.push_back(nil);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        return katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    };
    auto max_u = [](const katai::app::SolveResult& R) {
        double m = 0.0;
        for (Eigen::Index i = 0; i < R.disp.size(); ++i) m = std::max(m, std::fabs(R.disp(i)));
        return m;
    };
    const auto r0 = run(0.0);
    const auto rw = run(8.3);
    const bool ok = r0.size() == 2 && r0[0].ok && r0[1].ok && rw.size() == 2 && rw[0].ok &&
                    rw[1].ok;
    check(ok, "all wall phases solved");
    if (ok) {
        const double u0 = max_u(r0[0]), uw = max_u(rw[0]), unil = max_u(rw[1]);
        std::printf("   K0 max|u|: w=0 -> %.3e m, w=8.3 -> %.3e m; nil phase %.3e m\n",
                    u0, uw, unil);
        check(u0 <= 1e-10, "w = 0: undisturbed K0 identity preserved (bit-identical class)");
        check(uw > 1e-6, "w = 8.3: the wall genuinely settles under its own weight");
        check(unil <= 1e-7, "nil phase after the weighted wall stays a true no-op");
    }
}

// (h) NonPorous drainage is now real (audit item #3 -- it was accepted by the GUI and
// silently treated as Drained: a submerged concrete block got buoyancy AND pore
// pressure). PLAXIS rule: a non-porous material holds no water at all -- gamma_unsat
// everywhere, no pore pressure, total-stress equilibrium in its region. Pins on a
// fully submerged column (water table at the surface):
//   NonPorous:  sigma_yy(z) = -gamma_unsat * z (TOTAL; no buoyancy, no pore) and the
//               K0 identity max|u| ~ 0 still holds (seed, gravity and the pore-load
//               exclusion are mutually consistent through the full path);
//   Drained:    sigma'_yy(z) = -(gamma_sat - gamma_w) * z (the effective/buoyant
//               witness -- the two columns differ by design);
//   Consolidation with NonPorous is refused honestly (pore DOF + single Kw would
//   make concrete silently water-filled).
void test_nonporous_drainage() {
    std::printf("-- (h) NonPorous: total stress, no buoyancy/pore, honest guards --\n");
    auto column = [&](m::Drainage dr) {
        m::Project pr;
        m::Material s;
        s.name = "Concrete";
        s.model = m::SoilModel::LinearElastic;
        s.E = 3.0e7; s.nu = 0.2;
        s.gamma_unsat = 24.0; s.gamma_sat = 25.0;
        s.k0_auto = false; s.k0 = 0.5;
        s.drainage = dr;
        pr.materials.push_back(s);
        m::SoilPolygon L; L.material = 0;
        L.x = {0, kW, kW, 0}; L.y = {0, 0, kH, kH};
        L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
        pr.polygons.push_back(L);
        pr.has_water = true;                      // water table AT the surface: all submerged
        pr.wx = {-1.0, kW + 1.0};
        pr.wy = {kH, kH};
        return pr;
    };
    auto sigma_yy_mid = [&](const katai::app::SolveResult& R, double& u_max) {
        double s = 0.0; int n = 0;
        u_max = R.max_disp;
        for (int i = 0; i < R.mesh.node_count; ++i)
            if (std::fabs(R.mesh.y[i] - kH / 2) < 0.26) { s += R.stress.stress[i](1); ++n; }
        return n > 0 ? s / n : 0.0;
    };
    // NonPorous column: total stress, no buoyancy, K0 identity.
    {
        m::Project pr = column(m::Drainage::NonPorous);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 1 && res[0].ok;
        check(ok, "submerged NonPorous column solved (K0)");
        if (ok) {
            double umax = 0.0;
            const double sv = sigma_yy_mid(res[0], umax);
            const double sv_ex = -24.0 * kH / 2;              // -96 kPa: TOTAL, gamma_unsat
            std::printf("   mid-depth sigma_yy = %.3f kPa (total-form %.3f), K0 max|u| = %.2e m\n",
                        sv, sv_ex, umax);
            check(std::fabs(sv - sv_ex) < 0.01 * std::fabs(sv_ex),
                  "NonPorous: sigma_v = -gamma_unsat*z (no buoyancy, no pore; 1%)");
            check(umax < 1e-8, "NonPorous: the submerged K0 identity still holds (max|u| ~ 0)");
        }
    }
    // Drained witness: the same column responds with the buoyant EFFECTIVE stress.
    {
        m::Project pr = column(m::Drainage::Drained);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        const bool ok = res.size() == 1 && res[0].ok;
        check(ok, "submerged Drained column solved (K0)");
        if (ok) {
            double umax = 0.0;
            const double sv = sigma_yy_mid(res[0], umax);
            const double sv_ex = -(25.0 - 9.81) * kH / 2;     // -60.76 kPa: effective/buoyant
            std::printf("   mid-depth sigma'_yy = %.3f kPa (buoyant-form %.3f)\n", sv, sv_ex);
            check(std::fabs(sv - sv_ex) < 0.01 * std::fabs(sv_ex),
                  "Drained witness: sigma'_v = -(gamma_sat - gamma_w)*z (the columns differ)");
        }
    }
    // Consolidation with NonPorous: honest refusal.
    {
        m::Project pr = column(m::Drainage::NonPorous);
        m::Phase ph; ph.name = "consol"; ph.type = m::PhaseType::Consolidation;
        ph.duration = 10.0; ph.time_steps = 10;
        pr.phases.push_back(ph);
        const auto M = katai::app::mesh_from_project(pr, 0.5, 6);
        const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
        bool refused = false;
        for (const auto& r : res)
            if (!r.ok && r.message.find("Non-porous") != std::string::npos) refused = true;
        check(refused, "NonPorous + consolidation refused with an explicit message");
    }
}

}  // namespace

int main() {
    std::printf("INPUT AUDIT pins -- silently-ignored inputs now honest (fixed or refused)\n\n");
    test_ocr_k0();
    std::printf("\n");
    test_hs_undrained_b_refusal();
    std::printf("\n");
    test_plate_elastoplastic_gate();
    std::printf("\n");
    test_soil_tension_cutoff_applied();
    std::printf("\n");
    test_struct_weight_consistent_forces();
    std::printf("\n");
    test_wall_weight_static();
    std::printf("\n");
    test_nonporous_drainage();
    if (g_failures == 0) {
        std::printf("\nOK: OCR->K0 wired (PLAXIS formula), HS+Undrained(B) refused, plate Mp/Np "
                    "hinge live (zero-capacity elastoplastic refused)\n");
        return 0;
    }
    std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return 1;
}
