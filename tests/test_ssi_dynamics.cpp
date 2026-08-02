// Soil-structure interaction in a DYNAMIC (seismic) analysis -- the core assembly
// (analysis/structural_dynamics.hpp) that puts the structural stiffness into K, the plate inertia
// into M and builds the seismic influence vector r. Until now the Dynamic phase was soil-only and
// rejected any model with a wall; this is the assembly that lets soil and structure shake together.
//
// The user directive for this project is that self-consistency is NOT enough: every result must be
// cross-checked by an INDEPENDENT route. So this file pins the new assembly on four oracles that do
// not share a code path with it:
//
//   (a) ELASTIC-ASSEMBLY IDENTITY. The structural stiffness assembled here must BE the u=0 elastic
//       tangent that the long-validated static solver (solve_nonlinear, verified in test_plate_soil /
//       test_interface / test_anchor / test_geogrid) assembles internally. Independent oracle: solve
//       the same soil+wall+interface+anchor+geogrid problem BOTH ways and demand agreement to machine
//       precision. This is what stops the two DOF-mapping conventions from silently drifting apart.
//   (b) RIGID-BODY MODE. K r = 0 for the unsupported system: a rigid horizontal translation stores no
//       energy, in the soil OR in any structural element. Independent oracle: a theorem, not a number.
//       It catches a mis-mapped wall translation DOF, which would otherwise quietly produce a wall
//       that resists (or ignores) the base motion.
//   (c) TOTAL MASS. r^T M r = total translational mass = soil mass + sum(plate rho_A * L), exactly
//       (partition of unity). Independent oracle: hand arithmetic. It checks that the plate mass lands
//       on the DOFs the influence vector actually drives -- (b) and (c) together are what make
//       F = -M r a_g the right seismic force for a wall.
//   (d) PLATE CANTILEVER FREQUENCY. K phi = w^2 M phi on a plate-only system reproduces the
//       Euler-Bernoulli cantilever f_n = (beta_n L)^2/(2 pi L^2) sqrt(EI/rho_A) (beta_1 L = 1.8751,
//       beta_2 L = 4.6941; Blevins, Formulas for Natural Frequency and Mode Shape, Table 8-1).
//       Independent oracle: a closed-form textbook frequency. This is the check that rho_A is a real
//       mass in the right place -- (c) alone would pass even with the mass smeared wrongly along the
//       plate, because it only constrains the SUM of the mass matrix.
//
// Math: docs/references/dynamic-seismic-formulation.md sec 9. GUI-path SSI (a wall shaken through
// build_problem's Dynamic phase) is verified in test_dynamic_gui.
#include <katai/analysis/dynamics.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/structural_dynamics.hpp>
#include <katai/analysis/structural_forces.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::AnchorElement;
using katai::core::DofMap;
using katai::core::GeogridElement;
using katai::core::InterfaceElement;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::PlateElement;
using katai::core::Structures;
using katai::geometry::RectangularDomain;
using katai::math::CsrMatrix;
using katai::math::SparseMatrixBuilder;
using katai::mesh::Mesh;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

Eigen::MatrixXd dense(const CsrMatrix& A) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows, A.cols);
    for (int r = 0; r < A.rows; ++r)
        for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p) D(r, A.col_indices[p]) = A.values[p];
    return D;
}
// Dense LDLT -- the systems here are small, and it keeps this a pure-Eigen (MKL-free) core test.
Eigen::VectorXd solve_dense(const CsrMatrix& A, const Eigen::VectorXd& b) {
    return dense(A).ldlt().solve(b);
}

// A soil block with a WALL along an interior vertical line, built the way build_problem builds an
// embedded wall: the plate sits on its OWN translation + rotation DOFs (add_extra_dof) and is tied to
// the soil column by a zero-thickness interface. Plus a geogrid strip and a node-to-node anchor on the
// top surface, so a single model exercises every element the dynamic assembly claims to support.
// Soil E, nu -- shared by BOTH routes of check (a): solve_nonlinear takes MaterialModel, the direct
// assembly takes LinearElastic. Deriving both from the same pair keeps the comparison honest.
constexpr double kSoilE = 3.0e4, kSoilNu = 0.3;

struct Model {
    Mesh mesh;
    Structures structures;
    std::vector<MaterialModel> materials;               // solve_nonlinear (route 1)
    std::vector<katai::core::LinearElastic> elastic;    // assemble_stiffness (route 2)
    std::vector<int> wall_nodes;      // soil nodes on the wall line (y-sorted)
    double plate_len = 0.0;
    double beam_len = 0.0, beam_rho_A = 0.0;   // embedded beam (pile row), when built
};

// Soil nodes on the x = xw line, sorted by y.
std::vector<int> nodes_on_line(const Mesh& mesh, double xw) {
    std::vector<int> ns;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - xw) < 1e-9) ns.push_back(n);
    std::sort(ns.begin(), ns.end(), [&](int a, int b) { return mesh.y[a] < mesh.y[b]; });
    return ns;
}

// Build the model. `dofs` gains the plate's independent translation/rotation DOFs.
// with_mass:     give the plate inertia (only the dynamic checks need it).
// with_geogrid:  include the top geogrid strip. Check (a) leaves it OUT and tests it separately in
//                test_geogrid_linearization, because the geogrid is the one element whose static
//                branch is state-dependent (tension-only) -- see that test for the full story.
Model build_model(DofMap& dofs, const Mesh& mesh, bool with_mass, bool with_geogrid,
                  bool with_beam = false) {
    Model m;
    m.mesh = mesh;

    MaterialModel soil;
    soil.type = MaterialType::LinearElastic;
    soil.youngs_modulus = kSoilE; soil.poisson_ratio = kSoilNu;
    m.materials.push_back(soil);
    m.elastic.push_back(katai::core::LinearElastic{kSoilE, kSoilNu});

    // --- Wall: plate on independent DOFs + interface to the soil column at x = xw.
    const double xw = 2.0;
    m.wall_nodes = nodes_on_line(mesh, xw);
    katai::core::plate::PlateProps pp;
    const double Ep = 3.0e7, d = 0.3;
    pp.EA = Ep * d; pp.EI = Ep * d * d * d / 12.0; pp.nu = 0.15;
    if (with_mass) { pp.rho_A = 24.0 * d / 9.81; pp.rho_I = pp.rho_A * d * d / 12.0; }
    m.plate_len = mesh.y[m.wall_nodes.back()] - mesh.y[m.wall_nodes.front()];

    const int nw = (int)m.wall_nodes.size();
    std::vector<int> wx(nw), wy(nw), wphi(nw);
    for (int i = 0; i < nw; ++i) {
        wx[i] = dofs.add_extra_dof(); wy[i] = dofs.add_extra_dof(); wphi[i] = dofs.add_extra_dof();
    }
    katai::core::iface::InterfaceProps ip;
    ip.kn = 1.0e6; ip.ks = 1.0e6;
    // Strength far above anything the test loads reach -> the Coulomb return stays on its ELASTIC
    // branch, which is exactly the branch (a) compares against. (v1 dynamics is elastic by design.)
    ip.c_i = 1.0e9; ip.phi_i = 0.0; ip.sigma_t = 1.0e9;
    for (int e = 0; 2 * e + 2 < nw; ++e) {
        const int A = 2 * e, B = 2 * e + 2, M = 2 * e + 1;   // [A, B, mid] on the wall line
        m.structures.plates.push_back(PlateElement{
            {m.wall_nodes[A], m.wall_nodes[B], m.wall_nodes[M]},
            {wphi[A], wphi[B], wphi[M]}, pp,
            {wx[A], wy[A], wx[B], wy[B], wx[M], wy[M]}});
        m.structures.interfaces.push_back(InterfaceElement{
            {m.wall_nodes[A], m.wall_nodes[B], m.wall_nodes[M]},
            {wx[A], wy[A], wx[B], wy[B], wx[M], wy[M]}, ip, {0.0, 0.0, 0.0}});
    }

    // --- Geogrid strip along the top (shared soil DOFs, no rotation).
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });
    if (with_geogrid) {
        katai::core::geogrid::GeogridProps gp;
        gp.EA = 2.0e4; gp.Np = -1.0;   // no yield cap
        for (int e = 0; 2 * e + 2 < (int)top.size(); ++e)
            m.structures.geogrids.push_back(GeogridElement{
                {top[2 * e], top[2 * e + 2], top[2 * e + 1]}, gp});
    }

    // --- Node-to-node anchor between two top nodes (a strut). Elastic: no force cap.
    AnchorElement an;
    an.node_a = top.front(); an.node_b = top.back();
    an.EA = 1.0e5; an.L = mesh.x[top.back()] - mesh.x[top.front()];
    an.Fmax_tens = -1.0; an.Fmax_comp = -1.0;
    m.structures.anchors.push_back(an);

    // --- Embedded beam (pile row): its skin/foot coupling is mesh-NONCONFORMING -- the pile's nodes sit
    // at arbitrary (xi, eta) INSIDE soil elements, coupled through the soil shape functions. That makes
    // it the element most likely to get its DOF mapping wrong, and the rigid-mode / total-mass
    // invariants below are exactly what would catch it: under a rigid translation the soil node and the
    // pile node at the same place move together, so the skin springs must store no energy.
    if (with_beam) {
        const double bx0 = 1.0, by0 = 0.5, by1 = 2.5;   // a short vertical pile, inside the block
        std::vector<double> px, py;
        for (int i = 0; i < 5; ++i) { px.push_back(bx0); py.push_back(by0 + (by1 - by0) * i / 4.0); }
        katai::core::plate::PlateProps bp;
        bp.EA = 3.0e6; bp.EI = 1.0e4; bp.nu = 0.2;
        if (with_mass) { bp.rho_A = 24.0 * 0.13 / 9.81; bp.rho_I = bp.rho_A * 0.01 / 12.0; }
        m.structures.embedded_beams.push_back(katai::core::ebeam::build_embedded_beam(
            mesh, dofs, px, py, bp, /*k_axial=*/5.0e3, /*k_lateral=*/5.0e3));
        m.beam_len = by1 - by0;
        m.beam_rho_A = bp.rho_A;
    }
    return m;
}

// K = soil + structures, through the NEW assembly.
CsrMatrix assemble_full_K(const Mesh& mesh, const DofMap& dofs, const Model& m) {
    SparseMatrixBuilder bK(dofs.equation_count());
    katai::core::assemble_stiffness(mesh, dofs, m.elastic, bK);
    katai::core::assemble_structural_stiffness(mesh, dofs, m.structures, bK);
    return bK.build();
}

// ============================================================================================
// (a) The new elastic assembly IS the static solver's tangent -- to machine precision.
// ============================================================================================
void test_elastic_assembly_identity() {
    std::printf("-- (a) structural assembly == solve_nonlinear elastic tangent --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 3.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);

    DofMap dofs(mesh.node_count, 2);
    const Model m = build_model(dofs, mesh, false, false);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();
    check(neq > 0 && !m.structures.plates.empty() && !m.structures.interfaces.empty(),
          "soil + wall (plate on own DOFs + interface) + anchor built");

    // A horizontal load on the top-left corner: bends the wall, shears the interface and works the
    // strut -> every element carries force (a zero-load identity would prove nothing).
    Eigen::VectorXd f = Eigen::VectorXd::Zero(neq);
    int corner = -1; double best = 1e300;
    for (int n : mesh.top_nodes) if (mesh.x[n] < best) { best = mesh.x[n]; corner = n; }
    const int ceq = dofs.equation(dofs.global_dof(corner, 0));
    check(ceq >= 0, "load DOF is free");
    f[ceq] = 25.0;

    // Route 1 -- the validated static path (Newton; LE converges to the linear answer at round-off).
    katai::core::NewtonOptions opt{1, 30, 1e-12};
    const katai::core::NewtonResult nr = katai::core::solve_nonlinear(
        mesh, dofs, m.materials, f,
        [](const CsrMatrix& A, const Eigen::VectorXd& b) { return solve_dense(A, b); },
        opt, {}, {}, m.structures);
    check(nr.converged, "static reference solve converged");

    // Route 2 -- assemble K with the NEW code and solve the linear system once.
    const Eigen::VectorXd u_new = solve_dense(assemble_full_K(mesh, dofs, m), f);

    // Compare on the free DOFs (nr.displacement is total_dofs, indexed by global DOF).
    double dmax = 0.0, umax = 0.0;
    for (int g = 0; g < dofs.total_dofs(); ++g) {
        const int eq = dofs.equation(g);
        if (eq < 0) continue;
        dmax = std::fmax(dmax, std::fabs(u_new[eq] - nr.displacement[g]));
        umax = std::fmax(umax, std::fabs(nr.displacement[g]));
    }
    std::printf("   max |u_assembled - u_static| = %.3e   (max |u| = %.3e, relative %.2e)\n",
                dmax, umax, dmax / umax);
    check(umax > 1e-9, "reference solution is nontrivial (elements actually loaded)");
    // Same matrix, same DOF map => the same linear system. Round-off only (Newton residual 1e-12).
    check(dmax < 1e-9 * umax, "dynamic structural assembly reproduces the static tangent exactly");
}

// ============================================================================================
// (b) K r = 0 and (c) r^T M r = total mass -- on the UNSUPPORTED system.
// ============================================================================================
void test_rigid_mode_and_total_mass() {
    std::printf("-- (b) rigid-body mode K r = 0  /  (c) r^T M r = total mass --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 3.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);

    DofMap dofs(mesh.node_count, 2);
    // Everything at once -- geogrid AND the mesh-nonconforming embedded beam: neither may break the
    // rigid-body mode or the mass balance.
    const Model m = build_model(dofs, mesh, true, true, /*with_beam=*/true);
    dofs.finalize();   // NO supports: a rigid translation must be a genuine zero-energy mode
    const int neq = dofs.equation_count();

    const CsrMatrix K = assemble_full_K(mesh, dofs, m);
    const Eigen::VectorXd r = katai::core::seismic_influence_x(mesh, dofs, m.structures);

    // The influence vector must reach every INDEPENDENT translation DOF -- the wall's and the pile's --
    // not just the soil nodes. A pile left out would take no inertia force from the base: silent.
    const int n_beam = (int)m.structures.embedded_beams.front().dof_x.size();
    const int expect_r = mesh.node_count + (int)m.wall_nodes.size() + n_beam;
    std::printf("   r marks %d DOFs (soil ux %d + wall ux %d + pile ux %d = %d)\n",
                (int)r.sum(), mesh.node_count, (int)m.wall_nodes.size(), n_beam, expect_r);
    check((int)std::lround(r.sum()) == expect_r,
          "influence vector covers soil AND wall AND pile translation DOFs");

    // K r = 0: no strain anywhere under a rigid horizontal translation.
    const Eigen::VectorXd Kr = dense(K) * r;
    const double kscale = dense(K).cwiseAbs().maxCoeff();
    std::printf("   |K r|_inf = %.3e   (|K|_inf = %.3e, ratio %.2e)\n",
                Kr.cwiseAbs().maxCoeff(), kscale, Kr.cwiseAbs().maxCoeff() / kscale);
    check(Kr.cwiseAbs().maxCoeff() < 1e-9 * kscale, "rigid horizontal translation is a zero-energy mode");

    // r^T M r = soil mass + plate mass, exactly (partition of unity in each mass matrix).
    SparseMatrixBuilder bM(neq);
    const double rho_soil = 20.0 / 9.81;
    katai::core::assemble_mass(mesh, dofs, {rho_soil}, bM);
    katai::core::assemble_structural_mass(mesh, dofs, m.structures, bM);
    const CsrMatrix M = bM.build();

    const double m_soil = rho_soil * 4.0 * 3.0;                             // rho * area
    const double m_plate = m.structures.plates.front().props.rho_A * m.plate_len;
    const double m_beam = m.beam_rho_A * m.beam_len;
    const double m_tot = m_soil + m_plate + m_beam;
    const double m_num = r.transpose() * (dense(M) * r);
    std::printf("   r^T M r = %.9f   (soil %.6f + plate %.6f + pile %.6f = %.9f)\n",
                m_num, m_soil, m_plate, m_beam, m_tot);
    check(std::fabs(m_num - m_tot) < 1e-10 * m_tot,
          "total translational mass = soil + plate + pile (structural mass on the driven DOFs)");
    check(m_plate > 0.01 * m_soil, "plate mass is a meaningful fraction (the check has teeth)");
    check(m_beam > 0.0, "the pile carries mass too");
}

// ============================================================================================
// (f) SINGULAR MASS. Once a structure joins the dynamic system M is generally SINGULAR: a plate's
// rotation DOFs carry only rho_I (zero unless the plate has weight), and a weightless plate -- the
// GUI default, w = 0 -- contributes no mass at all while still owning DOFs. Those DOFs are zero ROWS
// of M. Newmark itself is fine with that (it only multiplies by M; the effective stiffness
// K_eff = K + a0 M + a1 C is non-singular because K stiffens every structural DOF) -- EXCEPT for the
// initial acceleration, which by default is obtained by FACTORIZING M. Handing a singular matrix to
// a direct solver's SPD path is not safe: measured with PARDISO it was an access violation,
// which /EHsc catch(...) does NOT catch,
// so it takes the whole application down with no message. That is exactly the crash class this
// project has already been bitten by twice.
//
// The fix is solve_newmark's optional a0_init, which the dynamic branch supplies in closed form:
// from rest, a(0) = -r a_g(0) solves M a(0) = F(0) = -M r a_g(0) EXACTLY for ANY M, singular or not.
// This test pins both halves: the hazard is real (zero rows exist), and the closed form is not an
// approximation (on a NON-singular M it reproduces the factorized answer to round-off).
// ============================================================================================
void test_singular_mass_initial_acceleration() {
    std::printf("-- (f) singular M (massless structural DOFs) + closed-form initial acceleration --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 3.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);

    DofMap dofs(mesh.node_count, 2);
    Model m = build_model(dofs, mesh, false, false);   // with_mass = false -> weightless plate
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();
    const double rho_soil = 20.0 / 9.81;

    auto build_M = [&](const Structures& s) {
        SparseMatrixBuilder bM(neq);
        katai::core::assemble_mass(mesh, dofs, {rho_soil}, bM);
        katai::core::assemble_structural_mass(mesh, dofs, s, bM);
        return bM.build();
    };
    const CsrMatrix M0 = build_M(m.structures);
    const CsrMatrix K = assemble_full_K(mesh, dofs, m);

    // The hazard: a weightless plate leaves zero rows in M.
    const Eigen::MatrixXd Md = dense(M0);
    int zero_rows = 0;
    for (int i = 0; i < neq; ++i) if (Md.row(i).cwiseAbs().maxCoeff() == 0.0) ++zero_rows;
    std::printf("   weightless plate: %d of %d rows of M are exactly zero -> M is singular\n",
                zero_rows, neq);
    check(zero_rows > 0, "a weightless plate leaves massless DOFs (M singular) -- the hazard is real");

    // The closed form, checked on a NON-SINGULAR M so the factorized route is available to compare.
    DofMap dofs2(mesh.node_count, 2);
    const Model m2 = build_model(dofs2, mesh, true, false);   // with_mass = true -> plate has weight
    for (int n : mesh.bottom_nodes) { dofs2.fix_node_component(n, 0); dofs2.fix_node_component(n, 1); }
    dofs2.finalize();
    const int neq2 = dofs2.equation_count();
    const CsrMatrix M2 = [&] {
        SparseMatrixBuilder bM(neq2);
        katai::core::assemble_mass(mesh, dofs2, {rho_soil}, bM);
        katai::core::assemble_structural_mass(mesh, dofs2, m2.structures, bM);
        return bM.build();
    }();
    const CsrMatrix K2 = assemble_full_K(mesh, dofs2, m2);
    SparseMatrixBuilder bC2(neq2);
    for (int r = 0; r < M2.rows; ++r)
        for (int p = M2.row_ptr[r]; p < M2.row_ptr[r + 1]; ++p)
            bC2.add_entry(r, M2.col_indices[p], 0.3 * M2.values[p]);
    for (int r = 0; r < K2.rows; ++r)
        for (int p = K2.row_ptr[r]; p < K2.row_ptr[r + 1]; ++p)
            bC2.add_entry(r, K2.col_indices[p], 1e-3 * K2.values[p]);
    const CsrMatrix C2 = bC2.build();

    const Eigen::VectorXd r2 = katai::core::seismic_influence_x(mesh, dofs2, m2.structures);
    const Eigen::VectorXd Mr2 = M2 * r2;
    const double amp = 1.5, freq = 2.0;
    auto ag = [&](double t) { return amp * std::cos(2 * kPi * freq * t); };   // a_g(0) != 0
    auto force = [&](int step) { return Eigen::VectorXd(-ag(step * 0.005) * Mr2); };
    const Eigen::VectorXd z2 = Eigen::VectorXd::Zero(neq2);

    // Route A: factorize M for a(0)  (the historical path -- needs M non-singular).
    const auto RA = katai::core::solve_newmark(M2, C2, K2, force, 0.005, 20, z2, z2, 0.5, 0.25);
    // Route B: closed-form a(0) = -r a_g(0), M never factorized.
    const auto RB = katai::core::solve_newmark(M2, C2, K2, force, 0.005, 20, z2, z2, 0.5, 0.25, {}, {},
                                               Eigen::VectorXd(-ag(0.0) * r2));
    double du = 0.0, umax = 0.0;
    for (size_t s = 0; s < RA.u.size(); ++s) {
        du = std::fmax(du, (RA.u[s] - RB.u[s]).cwiseAbs().maxCoeff());
        umax = std::fmax(umax, RA.u[s].cwiseAbs().maxCoeff());
    }
    std::printf("   non-singular M: max|u(factorized a0) - u(closed-form a0)| = %.3e (max|u| = %.3e)\n",
                du, umax);
    check(umax > 1e-12, "the comparison run actually moved");
    check(du < 1e-10 * umax, "closed-form a(0) = -r a_g(0) == factorizing M a(0) = F(0)");

    // And on the SINGULAR system the closed form runs and stays finite (this used to be the crash).
    SparseMatrixBuilder bC0(neq);
    for (int r = 0; r < M0.rows; ++r)
        for (int p = M0.row_ptr[r]; p < M0.row_ptr[r + 1]; ++p)
            bC0.add_entry(r, M0.col_indices[p], 0.3 * M0.values[p]);
    for (int r = 0; r < K.rows; ++r)
        for (int p = K.row_ptr[r]; p < K.row_ptr[r + 1]; ++p)
            bC0.add_entry(r, K.col_indices[p], 1e-3 * K.values[p]);
    const Eigen::VectorXd r0 = katai::core::seismic_influence_x(mesh, dofs, m.structures);
    const Eigen::VectorXd Mr0 = M0 * r0;
    auto force0 = [&](int step) { return Eigen::VectorXd(-ag(step * 0.005) * Mr0); };
    const Eigen::VectorXd z0 = Eigen::VectorXd::Zero(neq);
    const auto RS = katai::core::solve_newmark(M0, bC0.build(), K, force0, 0.005, 20, z0, z0, 0.5, 0.25,
                                               {}, {}, Eigen::VectorXd(-ag(0.0) * r0));
    bool finite = true;
    for (const auto& u : RS.u) if (!u.allFinite()) finite = false;
    std::printf("   singular M: %d steps integrated, all states finite = %s\n",
                (int)RS.u.size() - 1, finite ? "yes" : "NO");
    check(finite && RS.u.size() == 21, "Newmark integrates a SINGULAR-M (SSI) system without factorizing M");
}

// ============================================================================================
// (e) The geogrid: where the linear dynamic branch KNOWINGLY departs from the static path.
//
// The static geogrid is TENSION-ONLY (elements/geogrid.hpp axial_return: N_tr <= 0 -> N = 0 AND
// tangent = 0 -- a slack fabric carries nothing and stiffens nothing). A LINEAR dynamic system cannot
// carry a state-dependent stiffness, so the dynamic branch uses the taut tangent EA throughout. That
// is the standard linearisation about the static state: in a real reinforced-soil wall the geogrid is
// already taut under gravity, so its INCREMENTAL stiffness during shaking is EA.
//
// This test pins BOTH sides of that statement so neither can drift unnoticed:
//   - taut (tension): the two routes are the same system, to machine precision;
//   - slack (compression): the static path drops the geogrid entirely (== a model with no geogrid at
//     all), while the linear dynamic path keeps EA. The divergence is real, bounded and DELIBERATE.
// Same story, milder, for the anchor (no yield) and the interface (no Coulomb slip) -- v1 dynamics is
// a linearised analysis. Recorded in docs/validation/seismic-verification.md.
// ============================================================================================
void test_geogrid_linearization() {
    std::printf("-- (e) geogrid: linear dynamic branch vs tension-only static branch --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 2.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 2);

    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int b) { return mesh.x[a] < mesh.x[b]; });

    MaterialModel soil;
    soil.type = MaterialType::LinearElastic;
    soil.youngs_modulus = kSoilE; soil.poisson_ratio = kSoilNu;
    const std::vector<MaterialModel> mats{soil};
    const std::vector<katai::core::LinearElastic> elas{{kSoilE, kSoilNu}};

    Structures st;
    katai::core::geogrid::GeogridProps gp;
    gp.EA = 2.0e4; gp.Np = -1.0;
    for (int e = 0; 2 * e + 2 < (int)top.size(); ++e)
        st.geogrids.push_back(GeogridElement{{top[2 * e], top[2 * e + 2], top[2 * e + 1]}, gp});
    const Structures none;   // same soil, no geogrid -> the "fully slack" reference

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    // Equal and opposite horizontal forces at the two top corners: outward (sign=+1) stretches the
    // strip, inward (sign=-1) compresses it.
    auto load = [&](double sign) {
        Eigen::VectorXd f = Eigen::VectorXd::Zero(neq);
        const int l = dofs.equation(dofs.global_dof(top.front(), 0));
        const int r = dofs.equation(dofs.global_dof(top.back(), 0));
        if (l >= 0) f[l] = -sign * 20.0;
        if (r >= 0) f[r] = sign * 20.0;
        return f;
    };
    auto solve_static = [&](const Structures& s, const Eigen::VectorXd& f) {
        katai::core::NewtonOptions opt{1, 40, 1e-12};
        return katai::core::solve_nonlinear(
            mesh, dofs, mats, f,
            [](const CsrMatrix& A, const Eigen::VectorXd& b) { return solve_dense(A, b); },
            opt, {}, {}, s).displacement;
    };
    auto solve_linear = [&](const Structures& s, const Eigen::VectorXd& f) {
        SparseMatrixBuilder bK(neq);
        katai::core::assemble_stiffness(mesh, dofs, elas, bK);
        katai::core::assemble_structural_stiffness(mesh, dofs, s, bK);
        return solve_dense(bK.build(), f);
    };
    // Compare a static (total_dofs, global-indexed) field against a linear (equation-indexed) one.
    auto diff = [&](const Eigen::VectorXd& u_static, const Eigen::VectorXd& u_lin, double& umax) {
        double d = 0.0; umax = 0.0;
        for (int g = 0; g < dofs.total_dofs(); ++g) {
            const int eq = dofs.equation(g);
            if (eq < 0) continue;
            d = std::fmax(d, std::fabs(u_lin[eq] - u_static[g]));
            umax = std::fmax(umax, std::fabs(u_static[g]));
        }
        return d;
    };

    // --- Taut: identical systems.
    {
        const Eigen::VectorXd f = load(+1.0);
        double umax = 0.0;
        const double d = diff(solve_static(st, f), solve_linear(st, f), umax);
        std::printf("   taut  (tension):     max|u_lin - u_static| = %.3e  (max|u| = %.3e, rel %.2e)\n",
                    d, umax, d / umax);
        check(d < 1e-9 * umax, "taut geogrid: linear dynamic branch == tension-only static branch");
    }
    // --- Slack: static drops the geogrid; the linear branch keeps EA. Pin both halves.
    {
        const Eigen::VectorXd f = load(-1.0);
        double umax = 0.0, um2 = 0.0;
        const Eigen::VectorXd u_static = solve_static(st, f);
        const double d_slack = diff(u_static, solve_linear(none, f), um2);   // vs NO-geogrid linear
        std::printf("   slack (compression): static-with-geogrid vs no-geogrid = %.3e (rel %.2e)\n",
                    d_slack, d_slack / um2);
        check(d_slack < 1e-9 * um2, "slack geogrid carries nothing in the static path (== no geogrid)");

        const double d_lin = diff(u_static, solve_linear(st, f), umax);
        std::printf("   slack (compression): static vs LINEAR-with-geogrid   = %.3e (rel %.2e)"
                    "  <- deliberate v1 linearisation\n", d_lin, d_lin / umax);
        check(d_lin > 1e-3 * umax,
              "slack geogrid: linear branch keeps EA (known, documented v1 divergence)");
    }
}

// ============================================================================================
// (d) Plate-only cantilever: K phi = w^2 M phi vs the Euler-Bernoulli closed form.
// ============================================================================================
void test_plate_cantilever_frequency() {
    std::printf("-- (d) plate cantilever eigenfrequencies vs Euler-Bernoulli --\n");
    // Slender beam (L/d = 100) so shear deformation and rotary inertia -- which Timoshenko includes
    // and Euler-Bernoulli does not -- shift the frequencies by well under a percent.
    constexpr int nel = 8, nn = 2 * nel + 1;
    const double L = 10.0, d = 0.1, E = 3.0e7, rho = 2.4;   // rho [Mg/m^3]
    katai::core::plate::PlateProps pp;
    pp.EA = E * d; pp.EI = E * d * d * d / 12.0; pp.nu = 0.2;
    pp.rho_A = rho * d;                     // mass per unit length [Mg/m]
    pp.rho_I = pp.rho_A * d * d / 12.0;     // rotary inertia [Mg m]

    // A bare node line -- no soil elements at all, so K and M come ONLY from the plate.
    Mesh mesh;
    mesh.node_count = nn; mesh.element_count = 0; mesh.nodes_per_element = 6;
    mesh.x.resize(nn); mesh.y.assign(nn, 0.0);
    for (int i = 0; i < nn; ++i) mesh.x[i] = L * i / (nn - 1);

    DofMap dofs(mesh.node_count, 2);
    std::vector<int> rot(nn);
    for (int i = 0; i < nn; ++i) rot[i] = dofs.add_extra_dof();
    Structures st;
    for (int e = 0; e < nel; ++e)
        st.plates.push_back(PlateElement{{2 * e, 2 * e + 2, 2 * e + 1},
                                         {rot[2 * e], rot[2 * e + 2], rot[2 * e + 1]}, pp});
    dofs.fix_node_component(0, 0); dofs.fix_node_component(0, 1); dofs.fix(rot[0]);  // clamped end
    dofs.finalize();
    const int neq = dofs.equation_count();

    SparseMatrixBuilder bK(neq), bM(neq);
    katai::core::assemble_structural_stiffness(mesh, dofs, st, bK);
    katai::core::assemble_structural_mass(mesh, dofs, st, bM);
    const Eigen::MatrixXd K = dense(bK.build()), M = dense(bM.build());
    check(M.diagonal().minCoeff() > 0.0, "plate mass matrix is positive on every free DOF");

    // Generalized symmetric eigenproblem K phi = w^2 M phi.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(K, M);
    check(es.info() == Eigen::Success, "generalized eigenproblem solved");
    if (es.info() != Eigen::Success) return;

    // The axial (bar) modes live in the same spectrum; keep the BENDING ones by their frequency
    // ordering -- the two lowest are the first two flexural modes (axial modes of a slender beam are
    // far stiffer: f_axial,1 = (1/4L) sqrt(E/rho) is orders above f_bend,1).
    std::vector<double> f;
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        const double lam = es.eigenvalues()(i);
        if (lam > 1e-8) f.push_back(std::sqrt(lam) / (2 * kPi));
    }
    std::sort(f.begin(), f.end());
    check(f.size() >= 2, "at least two vibration modes found");
    if (f.size() < 2) return;

    // Blevins Table 8-1 (cantilever): f_n = (beta_n L)^2/(2 pi L^2) * sqrt(EI/rho_A).
    const double coef = std::sqrt(pp.EI / pp.rho_A) / (2 * kPi * L * L);
    const double f1_th = 1.87510407 * 1.87510407 * coef;
    const double f2_th = 4.69409113 * 4.69409113 * coef;
    std::printf("   f_1 = %.5f Hz  (Euler-Bernoulli %.5f, err %+.2f%%)\n",
                f[0], f1_th, 100 * (f[0] - f1_th) / f1_th);
    std::printf("   f_2 = %.5f Hz  (Euler-Bernoulli %.5f, err %+.2f%%)\n",
                f[1], f2_th, 100 * (f[1] - f2_th) / f2_th);
    // A slender Timoshenko beam sits just below Euler-Bernoulli (shear + rotary inertia soften it);
    // 1% / 2% bands hold that difference to its physical size rather than hiding a mass-matrix error.
    check(std::fabs(f[0] - f1_th) < 0.01 * f1_th, "cantilever f_1 = 1.8751^2/(2 pi L^2) sqrt(EI/rho_A)");
    check(std::fabs(f[1] - f2_th) < 0.02 * f2_th, "cantilever f_2 = 4.6941^2/(2 pi L^2) sqrt(EI/rho_A)");

    // Doubling the mass per length must drop every frequency by exactly sqrt(2) (f ~ 1/sqrt(rho_A)) --
    // an internal scaling law that pins rho_A as a mass, independent of the absolute value above.
    Structures st2 = st;
    for (auto& pe : st2.plates) { pe.props.rho_A *= 2.0; pe.props.rho_I *= 2.0; }
    SparseMatrixBuilder bM2(neq);
    katai::core::assemble_structural_mass(mesh, dofs, st2, bM2);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es2(K, dense(bM2.build()));
    std::vector<double> f2v;
    for (int i = 0; i < es2.eigenvalues().size(); ++i)
        if (es2.eigenvalues()(i) > 1e-8) f2v.push_back(std::sqrt(es2.eigenvalues()(i)) / (2 * kPi));
    std::sort(f2v.begin(), f2v.end());
    std::printf("   f_1(2 rho_A) = %.5f Hz   ratio f_1/f_1(2rho) = %.6f  (sqrt(2) = %.6f)\n",
                f2v[0], f[0] / f2v[0], std::sqrt(2.0));
    check(std::fabs(f[0] / f2v[0] - std::sqrt(2.0)) < 1e-6, "doubling rho_A scales frequencies by 1/sqrt(2)");
}

// ============================================================================================
// (g) SEISMIC STRUCTURAL FORCES -- the quasi-static limit.
//
// A dynamic run reports the wall's N/Q/M envelope over the shaking; those are the numbers an engineer
// designs the section with, so they need an oracle that does not come from the dynamic code at all.
// Shake the system SLOWLY (f << f_1): inertia and damping forces vanish next to the elastic ones, so
// the response degenerates to a sequence of STATIC solutions under the instantaneous body force
// -M r a_g(t). At the instant a_g is largest, the wall moment must therefore equal the moment from a
// STATIC solve under -M r a_g,peak -- computed here through the long-validated solve_nonlinear +
// plate_force_diagram path, sharing nothing with the Newmark driver.
//
// This pins the whole chain: the dynamic displacements, the force recovery from them, and the
// envelope accumulation. If the envelope silently reported the peak-DISPLACEMENT instant instead of
// the true maximum, or dropped a station, this check moves.
// ============================================================================================
void test_seismic_forces_quasi_static() {
    std::printf("-- (g) seismic wall force envelope vs the quasi-static limit --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 3.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);

    DofMap dofs(mesh.node_count, 2);
    const Model m = build_model(dofs, mesh, true, false);   // wall with weight, no geogrid
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    const double rho_soil = 20.0 / 9.81;
    SparseMatrixBuilder bM(neq);
    katai::core::assemble_mass(mesh, dofs, {rho_soil}, bM);
    katai::core::assemble_structural_mass(mesh, dofs, m.structures, bM);
    const CsrMatrix M = bM.build();
    const CsrMatrix K = assemble_full_K(mesh, dofs, m);
    const Eigen::VectorXd r = katai::core::seismic_influence_x(mesh, dofs, m.structures);
    const Eigen::VectorXd Mr = M * r;

    // Rayleigh damping tied to the two lowest modes -- generic, and it must not matter at f << f_1.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(K), dense(M));
    check(es.info() == Eigen::Success, "system eigenproblem solved (for f_1)");
    if (es.info() != Eigen::Success) return;
    const double f1 = std::sqrt(std::fmax(es.eigenvalues()(0), 1e-12)) / (2 * kPi);
    const auto ray = katai::core::rayleigh_from_modes(f1, 0.05, 5 * f1, 0.05);
    SparseMatrixBuilder bC(neq);
    for (int rr = 0; rr < M.rows; ++rr)
        for (int p = M.row_ptr[rr]; p < M.row_ptr[rr + 1]; ++p)
            bC.add_entry(rr, M.col_indices[p], ray.alpha * M.values[p]);
    for (int rr = 0; rr < K.rows; ++rr)
        for (int p = K.row_ptr[rr]; p < K.row_ptr[rr + 1]; ++p)
            bC.add_entry(rr, K.col_indices[p], ray.beta * K.values[p]);
    const CsrMatrix C = bC.build();

    // Quasi-static shaking: one full slow cycle at f = f_1/60, well resolved in time.
    const double amp = 2.0, freq = f1 / 60.0, dur = 1.0 / freq, nst = 600.0;
    const double dt = dur / nst;
    std::printf("   f_1 = %.3f Hz;  shaking at f = %.4f Hz (f_1/60) -> quasi-static\n", f1, freq);
    auto ag = [&](double t) { return amp * std::sin(2 * kPi * freq * t); };
    auto force = [&](int step) { return Eigen::VectorXd(-ag(step * dt) * Mr); };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);

    // Dynamic: accumulate the wall's |M| envelope exactly as build_problem's observer does.
    const std::vector<katai::core::PlateElement> chain(m.structures.plates.begin(),
                                                       m.structures.plates.end());
    double env_M = 0.0;
    auto observer = [&](int, double, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                        const Eigen::VectorXd&) {
        const Eigen::VectorXd uf = katai::core::expand_to_full(dofs, u);
        for (const auto& st : katai::core::plate_force_diagram(chain, mesh, dofs, uf))
            env_M = std::fmax(env_M, std::fabs(st.M));
    };
    katai::core::solve_newmark(M, C, K, force, dt, (int)nst, z, z, 0.5, 0.25, {}, observer,
                               Eigen::VectorXd(-ag(0.0) * r));

    // Static reference: the SAME body force at its peak, through the validated static solver.
    const Eigen::VectorXd f_static = -amp * Mr;
    katai::core::NewtonOptions opt{1, 30, 1e-12};
    const katai::core::NewtonResult nr = katai::core::solve_nonlinear(
        mesh, dofs, m.materials, f_static,
        [](const CsrMatrix& A, const Eigen::VectorXd& b) { return solve_dense(A, b); },
        opt, {}, {}, m.structures);
    check(nr.converged, "static reference solve converged");
    double ref_M = 0.0;
    for (const auto& st : katai::core::plate_force_diagram(chain, mesh, dofs, nr.displacement))
        ref_M = std::fmax(ref_M, std::fabs(st.M));

    std::printf("   wall max|M|:  dynamic envelope %.6g  vs  static under -M r a_peak %.6g  (err %+.2f%%)\n",
                env_M, ref_M, 100 * (env_M - ref_M) / ref_M);
    check(ref_M > 1e-9, "the static reference develops a real wall moment");
    // At f_1/60 the residual inertia scales as (f/f_1)^2 ~ 3e-4, so a 2% band is generous but still
    // far tighter than any wiring error (a dropped station or a snapshot-instead-of-envelope) survives.
    check(std::fabs(env_M - ref_M) < 0.02 * ref_M,
          "quasi-static seismic wall moment envelope = static moment under the same body force");

    // Linearity: the system is linear, so doubling the base amplitude must double the forces exactly.
    // Independent of the oracle above -- it checks the driver scales the way a linear system must.
    double env_M2 = 0.0;
    auto ag2 = [&](double t) { return 2.0 * amp * std::sin(2 * kPi * freq * t); };
    auto force2 = [&](int step) { return Eigen::VectorXd(-ag2(step * dt) * Mr); };
    auto observer2 = [&](int, double, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                         const Eigen::VectorXd&) {
        const Eigen::VectorXd uf = katai::core::expand_to_full(dofs, u);
        for (const auto& st : katai::core::plate_force_diagram(chain, mesh, dofs, uf))
            env_M2 = std::fmax(env_M2, std::fabs(st.M));
    };
    katai::core::solve_newmark(M, C, K, force2, dt, (int)nst, z, z, 0.5, 0.25, {}, observer2,
                               Eigen::VectorXd(-ag2(0.0) * r));
    std::printf("   amplitude x2 -> envelope ratio %.9f (linear system: exactly 2)\n", env_M2 / env_M);
    check(std::fabs(env_M2 / env_M - 2.0) < 1e-9, "doubling a_g doubles the force envelope (linearity)");

    // --- The wall<->soil Coulomb joint, same quasi-static oracle -------------------------------
    // The dynamic phase reports the joint's ELASTIC tau = ks du_s / sigma_n = kn du_n (that is what the
    // linear system solves; applying the Coulomb cap in post-processing would report a stress that is
    // NOT in equilibrium with the displacements the solver produced). This model's interface strength
    // is far above anything the test reaches, so the STATIC post-processor -- which DOES run the Coulomb
    // return -- stays on its elastic branch and the two must agree. That makes the validated static
    // interface_force_diagram a genuine oracle for the elastic branch AND pins the elastic flag itself:
    // if `elastic` accidentally applied sigma_n0 or a cap, this comparison moves.
    double env_tau = 0.0, env_sig = 0.0;
    auto obs_if = [&](int, double, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                      const Eigen::VectorXd&) {
        const Eigen::VectorXd uf = katai::core::expand_to_full(dofs, u);
        for (const auto& st : katai::core::interface_force_diagram(
                 m.structures.interfaces, 0, m.structures.interfaces.size(), mesh, dofs, uf, {}, true)) {
            env_tau = std::fmax(env_tau, std::fabs(st.tau));
            env_sig = std::fmax(env_sig, std::fabs(st.sigma_n));
        }
    };
    katai::core::solve_newmark(M, C, K, force, dt, (int)nst, z, z, 0.5, 0.25, {}, obs_if,
                               Eigen::VectorXd(-ag(0.0) * r));
    double ref_tau = 0.0, ref_sig = 0.0;
    for (const auto& st : katai::core::interface_force_diagram(
             m.structures.interfaces, 0, m.structures.interfaces.size(), mesh, dofs, nr.displacement, {}))
        { ref_tau = std::fmax(ref_tau, std::fabs(st.tau)); ref_sig = std::fmax(ref_sig, std::fabs(st.sigma_n)); }
    std::printf("   joint max|tau|:   dynamic envelope %.6g  vs  static %.6g  (err %+.2f%%)\n",
                env_tau, ref_tau, 100 * (env_tau - ref_tau) / ref_tau);
    std::printf("   joint max|sig_n|: dynamic envelope %.6g  vs  static %.6g  (err %+.2f%%)\n",
                env_sig, ref_sig, 100 * (env_sig - ref_sig) / ref_sig);
    check(ref_tau > 1e-9 && ref_sig > 1e-9, "the static reference develops real joint stresses");
    check(std::fabs(env_tau - ref_tau) < 0.02 * ref_tau,
          "quasi-static seismic joint tau envelope = static tau under the same body force");
    check(std::fabs(env_sig - ref_sig) < 0.02 * ref_sig,
          "quasi-static seismic joint sigma_n envelope = static sigma_n under the same body force");
}

// ============================================================================================
// (h) The interface's ELASTIC report is not the Coulomb report -- and must not be.
//
// The dynamic system solves the joint with k_n, k_s only. If the seismic output ran the Coulomb return
// (as the static post-processor does) it would cap tau at the strength and add the static sigma_n0 --
// reporting a stress the solver never produced, silently out of equilibrium with its own
// displacements. So `elastic` must differ from the default EXACTLY where the Coulomb branch bites, and
// nowhere else. Pin both directions on one displacement field.
// ============================================================================================
void test_interface_elastic_report() {
    std::printf("-- (h) interface elastic report vs the Coulomb report --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 3.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 3);
    DofMap dofs(mesh.node_count, 2);
    Model m = build_model(dofs, mesh, false, false);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    Eigen::VectorXd f = Eigen::VectorXd::Zero(neq);
    int corner = -1; double best = 1e300;
    for (int n : mesh.top_nodes) if (mesh.x[n] < best) { best = mesh.x[n]; corner = n; }
    f[dofs.equation(dofs.global_dof(corner, 0))] = 25.0;
    const Eigen::VectorXd u = solve_dense(assemble_full_K(mesh, dofs, m), f);
    const Eigen::VectorXd uf = katai::core::expand_to_full(dofs, u);
    const size_t nif = m.structures.interfaces.size();

    // High strength (the build_model default) -> the Coulomb branch never bites -> IDENTICAL reports.
    {
        const auto el = katai::core::interface_force_diagram(m.structures.interfaces, 0, nif, mesh, dofs, uf, {}, true);
        const auto co = katai::core::interface_force_diagram(m.structures.interfaces, 0, nif, mesh, dofs, uf, {}, false);
        double d = 0.0, mx = 0.0;
        for (size_t i = 0; i < el.size(); ++i) {
            d = std::fmax(d, std::fabs(el[i].tau - co[i].tau));
            mx = std::fmax(mx, std::fabs(co[i].tau));
        }
        std::printf("   strong joint (never slips): max|tau_elastic - tau_coulomb| = %.3e (max|tau| = %.3e)\n", d, mx);
        check(mx > 1e-9 && d < 1e-12 * mx,
              "on the elastic branch the two reports are the same number (the flag adds no bias)");
    }
    // Weak joint -> the Coulomb report caps tau; the elastic report must NOT (it mirrors the solver).
    {
        Structures weak = m.structures;
        for (auto& ie : weak.interfaces) { ie.props.c_i = 0.5; ie.props.phi_i = 0.0; ie.props.sigma_t = 0.0; }
        const auto el = katai::core::interface_force_diagram(weak.interfaces, 0, nif, mesh, dofs, uf, {}, true);
        const auto co = katai::core::interface_force_diagram(weak.interfaces, 0, nif, mesh, dofs, uf, {}, false);
        double max_el = 0.0, max_co = 0.0; bool any_capped = false;
        for (size_t i = 0; i < el.size(); ++i) {
            max_el = std::fmax(max_el, std::fabs(el[i].tau));
            max_co = std::fmax(max_co, std::fabs(co[i].tau));
            if (co[i].slipping) any_capped = true;
        }
        bool el_never_slips = true;
        for (const auto& st : el) if (st.slipping) el_never_slips = false;
        std::printf("   weak joint (c_i = 0.5): max|tau| elastic %.4g vs Coulomb-capped %.4g; capped = %s\n",
                    max_el, max_co, any_capped ? "yes" : "no");
        check(any_capped, "the weak joint does reach its Coulomb limit (the contrast is real)");
        check(max_co <= 0.5 + 1e-9, "the Coulomb report caps tau at c_i");
        check(max_el > max_co + 1e-9, "the elastic report does NOT cap -- it mirrors what the linear system solved");
        check(el_never_slips, "the elastic report never claims slip (the dynamic joint cannot slip)");
    }
}

}  // namespace

int main() {
    std::printf("Soil-structure interaction in the dynamic system (structural K + plate mass + r)\n\n");
    test_elastic_assembly_identity();
    std::printf("\n");
    test_rigid_mode_and_total_mass();
    std::printf("\n");
    test_plate_cantilever_frequency();
    std::printf("\n");
    test_geogrid_linearization();
    std::printf("\n");
    test_singular_mass_initial_acceleration();
    std::printf("\n");
    test_seismic_forces_quasi_static();
    std::printf("\n");
    test_interface_elastic_report();
    if (g_failures == 0)
        std::printf("\nOK: structural K = static tangent, K r = 0, r^T M r = total mass, "
                    "plate cantilever = Euler-Bernoulli, geogrid linearisation bounded, "
                    "singular-M start safe, seismic force envelope = quasi-static limit\n");
    else
        std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
