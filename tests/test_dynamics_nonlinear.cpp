// Nonlinear dynamic (Newmark + per-step Newton) solver -- analysis/dynamics_nonlinear.cpp. The linear
// solve_newmark treats f_int = K u; this generalizes it to f_int(u) from the constitutive return
// mapping, so the soil can yield (and, through the shared structural assembly, an interface can slip /
// a geogrid go slack) DURING shaking, exactly like PLAXIS Dynamics (Ref sec 11.10.4: a dynamic phase
// uses the same convergence criteria as a plastic phase).
//
// Per the project directive, self-consistency is not enough: each result is pinned by an INDEPENDENT
// route that shares no code with the driver.
//
//   (a) LINEAR LIMIT. With a linear-elastic material and no yielding structure, f_int(u) = K u, so the
//       nonlinear solver MUST reproduce the long-validated linear solve_newmark trajectory to round-off
//       (Newton converges in one step on a linear system). Independent oracle: solve_newmark itself.
//       This is what proves the new machinery reduces exactly to the verified linear path -- every
//       existing seismic verification therefore still covers the nonlinear solver in its linear regime.
//
//   (b) QUASI-STATIC LIMIT. Shaken slowly (f << f_1) inertia and damping vanish next to the elastic/
//       plastic internal forces, so the response degenerates to a sequence of STATIC equilibria under
//       the instantaneous body force -M r a_g(t). At peak a_g the dynamic field must equal the STATIC
//       solve under -M r a_peak -- computed through the validated solve_nonlinear, sharing nothing with
//       the Newmark driver. Run with a MOHR-COULOMB soil that genuinely yields, so this also pins the
//       nonlinear branch: the static oracle is itself elastoplastic, and a linear dynamic run cannot
//       match it (checked explicitly).
//
//   (d) PARENT STRUCTURAL STATE CARRY (Track 1a). Independent oracle: capacity ALGEBRA. A statically
//       preloaded anchor (N_s = 0.6 F_max) shaken with a quasi-static increment of +-0.6 F_max must
//       yield ONLY when the parent state is carried (total 1.2 F_max > F_max; the increment alone,
//       0.6 F_max, stays elastic) -- and must NOT yield even with carry when the increment is small
//       enough that the total stays below capacity (0.85 F_max). Plus the ZERO-FORCE IDENTITY: seeded
//       with the parent state and driven by nothing, the solver must stay EXACTLY at rest with the
//       committed structural state bit-for-bit equal to its seed (datum + seed + baseline mutually
//       consistent, no spurious start-up transient).
//
// Math: docs/references/dynamic-seismic-formulation.md. Linear SSI machinery: test_ssi_dynamics.
#include <katai/analysis/dynamics.hpp>
#include <katai/analysis/dynamics_nonlinear.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/structural_dynamics.hpp>   // seismic_influence_x
#include <katai/analysis/structural_forces.hpp>     // anchor_force (measures N like the solver)
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
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
// General (partial-pivot) LU: K_eff = K_T + a0 M + a1 C is NONsymmetric for non-associated plasticity,
// so an LDLT would be wrong for the MC test. Small systems -> pure-Eigen, MKL-free.
Eigen::VectorXd solve_lu(const CsrMatrix& A, const Eigen::VectorXd& b) {
    return dense(A).partialPivLu().solve(b);
}

constexpr double kE = 3.0e4, kNu = 0.3, kGamma = 18.0, kG = 9.81;

// ============================================================================================
// (a) Linear limit: nonlinear Newmark == linear solve_newmark, to round-off.
// ============================================================================================
void test_linear_limit() {
    std::printf("-- (a) linear limit: nonlinear Newmark reproduces linear solve_newmark --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 4.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 4);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    MaterialModel soil;
    soil.type = MaterialType::LinearElastic;
    soil.youngs_modulus = kE; soil.poisson_ratio = kNu;
    const std::vector<MaterialModel> mats{soil};
    const std::vector<katai::core::LinearElastic> elas{{kE, kNu}};
    const double rho = kGamma / kG;

    SparseMatrixBuilder bM(neq), bK(neq);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    katai::core::assemble_stiffness(mesh, dofs, elas, bK);
    const CsrMatrix M = bM.build(), K = bK.build();
    const auto ray = katai::core::rayleigh_from_modes(2.0, 0.05, 10.0, 0.05);
    SparseMatrixBuilder bC(neq);
    for (int r = 0; r < M.rows; ++r)
        for (int p = M.row_ptr[r]; p < M.row_ptr[r + 1]; ++p)
            bC.add_entry(r, M.col_indices[p], ray.alpha * M.values[p]);
    for (int r = 0; r < K.rows; ++r)
        for (int p = K.row_ptr[r]; p < K.row_ptr[r + 1]; ++p)
            bC.add_entry(r, K.col_indices[p], ray.beta * K.values[p]);
    const CsrMatrix C = bC.build();

    const Eigen::VectorXd r = katai::core::seismic_influence_x(mesh, dofs, {});
    const Eigen::VectorXd Mr = M * r;
    const double amp = 1.5, freq = 3.0;
    auto ag = [&](double t) { return amp * std::sin(2 * kPi * freq * t); };
    auto force = [&](int step) { return Eigen::VectorXd(-ag(step * 0.004) * Mr); };
    const double dt = 0.004; const int nst = 250;
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
    const Eigen::VectorXd a0 = -ag(0.0) * r;

    // Linear oracle: full-history solve_newmark, same a0 so the start is bit-identical.
    const auto RL = katai::core::solve_newmark(M, C, K, force, dt, nst, z, z, 0.5, 0.25, {}, {}, a0);

    // Nonlinear: collect the streamed trajectory.
    std::vector<Eigen::VectorXd> un;
    auto obs = [&](int, double, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                   const Eigen::VectorXd&) { un.push_back(u); };
    katai::core::NewmarkNonlinearOptions opt;
    opt.newton_tol = 1e-12;   // drive Newton to round-off so the comparison is with the exact linear step
    const auto RN = katai::core::solve_newmark_nonlinear(
        mesh, dofs, mats, M, C, force, dt, nst, z, z, a0, solve_lu, obs, {}, {}, {}, {}, {}, opt);

    check(RN.converged && RN.steps_completed == nst, "nonlinear solver integrated every step");
    check(un.size() == RL.u.size(), "trajectories have the same length");
    double du = 0.0, umax = 0.0;
    const size_t ncmp = std::min(un.size(), RL.u.size());
    for (size_t s = 0; s < ncmp; ++s) {
        du = std::fmax(du, (un[s] - RL.u[s]).cwiseAbs().maxCoeff());
        umax = std::fmax(umax, RL.u[s].cwiseAbs().maxCoeff());
    }
    std::printf("   max|u_nonlinear - u_linear| = %.3e   (max|u| = %.3e, relative %.2e)\n",
                du, umax, du / umax);
    check(umax > 1e-6, "the linear run actually moved (the check has teeth)");
    check(du < 1e-9 * umax, "nonlinear Newmark == linear solve_newmark on an elastic system (round-off)");
}

// ============================================================================================
// (b) Quasi-static limit with a YIELDING Mohr-Coulomb soil vs the static elastoplastic solve.
// ============================================================================================
void test_quasi_static_mc() {
    std::printf("-- (b) quasi-static limit: nonlinear seismic vs static elastoplastic solve (MC) --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 4.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 5, 5);

    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    // Tresca-type (phi = 0) soil: a horizontal body force builds shear; above c it yields near the base.
    MaterialModel soil;
    soil.type = MaterialType::MohrCoulomb;
    soil.youngs_modulus = kE; soil.poisson_ratio = kNu;
    soil.cohesion = 30.0; soil.friction_angle = 0.0; soil.dilatancy_angle = 0.0;
    const std::vector<MaterialModel> mats{soil};
    const std::vector<katai::core::LinearElastic> elas{{kE, kNu}};
    const double rho = kGamma / kG;

    SparseMatrixBuilder bM(neq), bK(neq);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    katai::core::assemble_stiffness(mesh, dofs, elas, bK);   // initial (elastic) K for Rayleigh + linear ref
    const CsrMatrix M = bM.build(), K0 = bK.build();
    const Eigen::VectorXd r = katai::core::seismic_influence_x(mesh, dofs, {});
    const Eigen::VectorXd Mr = M * r;

    // f_1 of the elastic system, to set a quasi-static shaking frequency.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(K0), dense(M));
    const double f1 = std::sqrt(std::fmax(es.eigenvalues()(0), 1e-12)) / (2 * kPi);
    const auto ray = katai::core::rayleigh_from_modes(f1, 0.05, 5 * f1, 0.05);
    SparseMatrixBuilder bC(neq);
    for (int rr = 0; rr < M.rows; ++rr)
        for (int p = M.row_ptr[rr]; p < M.row_ptr[rr + 1]; ++p)
            bC.add_entry(rr, M.col_indices[p], ray.alpha * M.values[p]);
    for (int rr = 0; rr < K0.rows; ++rr)
        for (int p = K0.row_ptr[rr]; p < K0.row_ptr[rr + 1]; ++p)
            bC.add_entry(rr, K0.col_indices[p], ray.beta * K0.values[p]);
    const CsrMatrix C = bC.build();

    // Amplitude tuned for PARTIAL yield: above the base-shear yield onset (a_g H rho ~ c) but well
    // below a collapse mechanism, so the static elastoplastic oracle converges at load_factor = 1 with
    // a real plastic zone (checked below via the linear/nonlinear gap).
    const double amp = 3.5, freq = f1 / 60.0, dur = 1.0 / freq;
    const int nst = 480; const double dt = dur / nst;
    std::printf("   f_1 = %.3f Hz;  shaking at f = %.4f Hz (f_1/60) -> quasi-static;  a_peak = %.1f m/s^2\n",
                f1, freq, amp);
    auto ag = [&](double t) { return amp * std::sin(2 * kPi * freq * t); };
    auto force = [&](int step) { return Eigen::VectorXd(-ag(step * dt) * Mr); };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);

    // Dynamic (nonlinear): keep BOTH the peak-|u| scalar (magnitude headline) and the deformed field at
    // the instant a_g is most POSITIVE -- there force = -a_g M r = -a_peak M r, matching the static solve
    // under -M r a_peak sign-for-sign, so the field comparison is not confounded by the response's overall
    // sign (the Tresca response is antisymmetric: u(-a) = -u(+a), so the wrong half-cycle would flip it).
    double peak = 0.0; double max_ag = 0.0; Eigen::VectorXd u_at_maxag;
    auto obs = [&](int, double t, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                   const Eigen::VectorXd&) {
        peak = std::fmax(peak, u.cwiseAbs().maxCoeff());
        const double g = ag(t);
        if (g > max_ag) { max_ag = g; u_at_maxag = u; }
    };
    const auto RN = katai::core::solve_newmark_nonlinear(
        mesh, dofs, mats, M, C, force, dt, nst, z, z, Eigen::VectorXd(-ag(0.0) * r),
        solve_lu, obs, {}, {}, {}, {});
    check(RN.converged, "nonlinear seismic solve converged over the whole record");

    // Static oracle: the same body force at its peak, through the validated elastoplastic solver.
    const Eigen::VectorXd f_static = -amp * Mr;
    katai::core::NewtonOptions nopt{8, 40, 1e-9};
    const katai::core::NewtonResult nr = katai::core::solve_nonlinear(
        mesh, dofs, mats, f_static, solve_lu, nopt);
    std::printf("   static solve load_factor = %.4f (1.0 = fully converged, <1 = collapse)\n", nr.load_factor);
    check(nr.converged && nr.load_factor > 0.999,
          "static elastoplastic reference converged (partial yield, not collapse)");

    // Compare the free-DOF fields (dynamic at min-a_g vs static under -M r a_peak). At f_1/60 the residual
    // inertia ~ (f/f_1)^2 ~ 3e-4, so a 2% band holds it well below any driver/recovery wiring error.
    double d = 0.0, us = 0.0;
    for (int g = 0; g < dofs.total_dofs(); ++g) {
        const int eq = dofs.equation(g);
        if (eq < 0) continue;
        d = std::fmax(d, std::fabs(u_at_maxag[eq] - nr.displacement[g]));
        us = std::fmax(us, std::fabs(nr.displacement[g]));
    }
    std::printf("   max|u|: dynamic peak %.6g  vs  static %.6g;  field diff at max-a_g rel %.2e\n",
                peak, us, d / us);
    check(us > 1e-6, "the static reference developed a real displacement");
    check(d < 0.02 * us, "quasi-static nonlinear seismic field == static elastoplastic field");

    // The nonlinearity actually engaged: a LINEAR (elastic) solve under the same peak force is
    // materially different from the yielding static/dynamic result. If the soil had stayed elastic this
    // gap would be ~0, and (b) would be a hollow linear check.
    const Eigen::VectorXd u_lin = solve_lu(K0, f_static);
    double dlin = 0.0;
    for (int g = 0; g < dofs.total_dofs(); ++g) {
        const int eq = dofs.equation(g);
        if (eq < 0) continue;
        dlin = std::fmax(dlin, std::fabs(u_lin[eq] - nr.displacement[g]));
    }
    std::printf("   nonlinearity witness: |u_elastic - u_plastic| / |u_plastic| = %.2e (must be >> 0)\n",
                dlin / us);
    check(dlin > 0.05 * us, "the MC soil genuinely yielded (linear elastic result differs materially)");
}

// ============================================================================================
// (c) CYCLIC HYSTERESIS -- a Mohr-Coulomb material dissipates ONLY when it touches the yield surface.
//
// PLAXIS Material Manual sec 3.5: Mohr-Coulomb is linear-elastic INSIDE the failure surface, so a cyclic
// stress path that stays within it stores and returns energy with NO hysteretic damping; only once the
// path TOUCHES the surface does plastic flow open the stress-strain loop and dissipate energy. The
// nonlinear dynamic solver must reproduce exactly this on/off behaviour.
//
// Independent oracle: energy conservation. Set the damping matrix to ZERO, so plasticity is the ONLY
// dissipation mechanism, and drive one full quasi-static cycle. The net work done by the external force,
// W = closed-integral F . du, then EQUALS the energy dissipated: for a below-yield (elastic) cycle the
// system returns to its start along a reversible path so W = 0; for a yielding cycle W = the plastic
// dissipation > 0. This is a line integral of the streamed force/displacement -- it shares no code with
// the solver.
// ============================================================================================
void test_cyclic_hysteresis() {
    std::printf("-- (c) cyclic hysteresis: MC dissipates ONLY on yielding (PLAXIS Mat.Man 3.5) --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 4.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 5, 5);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    MaterialModel soil;
    soil.type = MaterialType::MohrCoulomb;
    soil.youngs_modulus = kE; soil.poisson_ratio = kNu;
    soil.cohesion = 30.0; soil.friction_angle = 0.0; soil.dilatancy_angle = 0.0;   // Tresca (phi = 0)
    const std::vector<MaterialModel> mats{soil};
    const double rho = kGamma / kG;

    SparseMatrixBuilder bM(neq), bK(neq);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    katai::core::assemble_stiffness(mesh, dofs, {{kE, kNu}}, bK);
    const CsrMatrix M = bM.build(), K0 = bK.build();
    const CsrMatrix Czero = SparseMatrixBuilder(neq).build();   // ZERO damping: plasticity is the only sink
    const Eigen::VectorXd r = katai::core::seismic_influence_x(mesh, dofs, {});
    const Eigen::VectorXd Mr = M * r;

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(K0), dense(M));
    const double f1 = std::sqrt(std::fmax(es.eigenvalues()(0), 1e-12)) / (2 * kPi);
    const double freq = f1 / 60.0, dur = 1.0 / freq;   // exactly ONE full quasi-static cycle
    const int nst = 480; const double dt = dur / nst;
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);

    // Net work by the external force over the cycle, W = integral F . du (trapezoidal on the streamed field).
    auto run_cycle = [&](double amp, bool& ok) -> double {
        auto ag = [&](double t) { return amp * std::sin(2 * kPi * freq * t); };
        auto force = [&](int step) { return Eigen::VectorXd(-ag(step * dt) * Mr); };
        double W = 0.0; Eigen::VectorXd u_prev = z; double ag_prev = ag(0.0);
        auto obs = [&](int step, double t, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                       const Eigen::VectorXd&) {
            const double agt = ag(t);
            if (step > 0) {
                const Eigen::VectorXd F = -agt * Mr, Fp = -ag_prev * Mr;
                W += 0.5 * (F + Fp).dot(u - u_prev);   // 0.5 (F_k + F_{k-1}) . (u_k - u_{k-1})
            }
            u_prev = u; ag_prev = agt;
        };
        const auto R = katai::core::solve_newmark_nonlinear(
            mesh, dofs, mats, M, Czero, force, dt, nst, z, z, Eigen::VectorXd(-ag(0.0) * r),
            solve_lu, obs, {}, {}, {}, {});
        ok = R.converged;
        return W;
    };

    bool ok_e = false, ok_p = false;
    const double W_elastic = run_cycle(1.0, ok_e);   // below yield -> reversible, elastic
    const double W_plastic = run_cycle(3.5, ok_p);   // touches the yield surface -> plastic dissipation
    check(ok_e && ok_p, "both cyclic runs converged with zero damping");
    std::printf("   elastic cycle (amp 1.0): W = %.4e J/m   plastic cycle (amp 3.5): W = %.4e J/m\n",
                W_elastic, W_plastic);
    std::printf("   |W_elastic| / W_plastic = %.2e  (should be ~0: MC is elastic inside the surface)\n",
                std::fabs(W_elastic) / std::fmax(W_plastic, 1e-30));
    check(W_plastic > 1e-6, "the yielding cycle DISSIPATES energy (a real hysteresis loop forms)");
    check(std::fabs(W_elastic) < 0.02 * W_plastic,
          "the below-yield cycle dissipates ~nothing (no damping when the path stays inside the surface)");
}

// ============================================================================================
// (d) PARENT STRUCTURAL STATE CARRY (Track 1a) -- capacity-algebra oracle + zero-force identity.
//
// A fixed-end anchor grounds a surface node. A static pull develops N_s; F_max is then DEFINED as
// N_s / 0.6, so the preload sits at exactly 60% of capacity. A quasi-static harmonic pull of the
// same amplitude produces a dynamic increment of ~0.6 F_max (elastic transfer is linear). Then:
//   - WITHOUT the parent state, the increment (0.6 F_max < F_max) must NOT yield the anchor;
//   - WITH it, the total (1.2 F_max on the tension half-cycle) MUST yield -- the asymmetric onset
//     only the carried preload can produce;
//   - WITH it but a smaller shake (total 0.85 F_max), it must NOT yield: the onset is governed by
//     the TOTAL crossing capacity, not by the carry machinery itself.
// The plastic state is checked BITWISE (the elastic return branch hands the committed value back
// unchanged), so "did not yield" has no tolerance to hide behind.
// ============================================================================================
void test_parent_state_carry() {
    std::printf("-- (d) parent structural state carry: preloaded anchor, asymmetric yield onset --\n");
    const RectangularDomain domain{0.0, 0.0, 4.0, 4.0, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri6(domain, 4, 4);
    DofMap dofs(mesh.node_count, 2);
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    dofs.finalize();
    const int neq = dofs.equation_count();

    MaterialModel soil;
    soil.type = MaterialType::LinearElastic;   // keep the soil elastic: the anchor is the only cap
    soil.youngs_modulus = kE; soil.poisson_ratio = kNu;
    const std::vector<MaterialModel> mats{soil};
    const double rho = kGamma / kG;

    // Anchor: top-right surface node, fixed far end in +x -> a grounded axial spring (EA/L).
    int node_a = -1;
    for (int n = 0; n < mesh.node_count; ++n)
        if (std::fabs(mesh.x[n] - 4.0) < 1e-9 && std::fabs(mesh.y[n] - 4.0) < 1e-9) node_a = n;
    check(node_a >= 0, "found the anchor node at (4,4)");
    katai::core::AnchorElement an;
    an.node_a = node_a; an.node_b = -1;
    an.fixed_point = Eigen::Vector2d(12.0, 4.0);   // dir = +x, L_geom = 8
    an.EA = 2.4e5;                                 // kk = EA/8 = 3e4 kN/m per m
    katai::core::Structures S;
    S.anchors.push_back(an);

    // Static preload: pull node_a in -x (away from the fixed end) -> anchor tension N_s.
    const int eqx = dofs.equation(dofs.global_dof(node_a, 0));
    check(eqx >= 0, "the anchor node is free in x");
    const double P_s = 40.0;
    Eigen::VectorXd f_static = Eigen::VectorXd::Zero(neq);
    f_static(eqx) = -P_s;
    katai::core::NewtonOptions nopt{4, 40, 1e-11};
    // Pass 1 (unbounded anchor) only measures N_s so F_max can be DEFINED from it.
    katai::core::NewtonResult nr0 = katai::core::solve_nonlinear(mesh, dofs, mats, f_static, solve_lu,
                                                                 nopt, {}, {}, S);
    check(nr0.converged, "static preload solve converged (unbounded anchor)");
    const double N_s = katai::core::anchor_force(S.anchors[0], mesh, dofs, nr0.displacement, 0.0).N;
    check(N_s > 1.0, "the static pull really loads the anchor in tension");
    const double Fmax = N_s / 0.6;
    S.anchors[0].Fmax_tens = Fmax;   // preload = 60% of capacity, by construction
    // Pass 2, with the real capacity: identical (still elastic) -- but ITS committed state is the
    // honest parent state (same rule build_problem uses: the parent phase ran with these structures).
    const katai::core::NewtonResult nr = katai::core::solve_nonlinear(mesh, dofs, mats, f_static,
                                                                      solve_lu, nopt, {}, {}, S);
    check(nr.converged && nr.anchor_plastic.size() == 1 && nr.anchor_plastic[0] == 0.0,
          "the preloaded anchor is elastic (N_s = 0.6 F_max, no static yield)");
    std::printf("   N_s = %.4f kN/m  ->  F_max := N_s/0.6 = %.4f kN/m\n", N_s, Fmax);

    // Dynamic machinery: quasi-static harmonic pull at the same node, zero damping.
    SparseMatrixBuilder bM(neq), bK(neq);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    katai::core::assemble_stiffness(mesh, dofs, {{kE, kNu}}, bK);
    const CsrMatrix M = bM.build(), K0 = bK.build();
    const CsrMatrix Czero = SparseMatrixBuilder(neq).build();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(K0), dense(M));
    const double f1 = std::sqrt(std::fmax(es.eigenvalues()(0), 1e-12)) / (2 * kPi);
    const double freq = f1 / 60.0, dur = 1.0 / freq;
    const int nst = 240; const double dt = dur / nst;
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);

    // Parent state, exactly as build_problem carries it: equation-space displacement datum + the
    // committed plastic state of the SAME solve.
    katai::core::StructuralInit sinit;
    sinit.u_datum = Eigen::VectorXd::Zero(neq);
    for (int gd = 0; gd < dofs.total_dofs(); ++gd) {
        const int eq = dofs.equation(gd);
        if (eq >= 0) sinit.u_datum(eq) = nr.displacement(gd);
    }
    sinit.anchor_plastic = nr.anchor_plastic;

    auto run_dyn = [&](double P_d, const katai::core::StructuralInit& init, bool& ok, double& u_pk) {
        auto force = [&](int step) {
            Eigen::VectorXd f = Eigen::VectorXd::Zero(neq);
            f(eqx) = -P_d * std::sin(2 * kPi * freq * step * dt);
            return f;
        };
        double pk = 0.0;
        auto obs = [&](int, double, const Eigen::VectorXd& u, const Eigen::VectorXd&,
                       const Eigen::VectorXd&) { pk = std::fmax(pk, u.cwiseAbs().maxCoeff()); };
        const auto R = katai::core::solve_newmark_nonlinear(
            mesh, dofs, mats, M, Czero, force, dt, nst, z, z, z, solve_lu, obs,
            nr.gauss_states, {}, S, {}, init, {});
        ok = R.converged;
        u_pk = pk;
        return R;
    };

    bool ok1 = false, ok2 = false, ok3 = false, ok4 = false;
    double pk1 = 0.0, pk2 = 0.0, pk3 = 0.0, pk4 = 0.0;
    // (i) increment alone (no parent state): 0.6 F_max < F_max -> must stay elastic, bitwise.
    const auto Rn = run_dyn(P_s, {}, ok1, pk1);
    // (ii) carried parent state, same shake: total 1.2 F_max on the tension half-cycle -> must yield.
    const auto Rc = run_dyn(P_s, sinit, ok2, pk2);
    // (iii) carried parent state, small shake: total 0.85 F_max -> must stay elastic, bitwise.
    const auto Rs = run_dyn(P_s * 0.25 / 0.6, sinit, ok3, pk3);
    check(ok1 && ok2 && ok3, "all three quasi-static runs converged");
    check(!Rn.struct_state_carried && Rc.struct_state_carried && Rs.struct_state_carried,
          "struct_state_carried reports what the solver actually did");
    std::printf("   plastic U_p:  no-carry %.3e | carry+big %.3e | carry+small %.3e\n",
                Rn.anchor_plastic[0], Rc.anchor_plastic[0], Rs.anchor_plastic[0]);
    check(Rn.anchor_plastic[0] == 0.0,
          "WITHOUT the parent state the 0.6 F_max increment does NOT yield (bitwise zero U_p)");
    check(Rc.anchor_plastic[0] > 1e-8,
          "WITH the parent state the SAME shake yields (total 1.2 F_max -- asymmetric onset)");
    check(Rs.anchor_plastic[0] == nr.anchor_plastic[0],
          "carry alone does not fabricate yield: total 0.85 F_max stays elastic (U_p == seed, bitwise)");

    // (iv) ZERO-FORCE IDENTITY: seeded with the parent state and driven by nothing, the solver must
    // stay EXACTLY at rest (baseline consistency) and commit the seed back bit-for-bit.
    const auto Ri = run_dyn(0.0, sinit, ok4, pk4);
    check(ok4, "zero-force run converged");
    std::printf("   zero-force identity: peak |u| = %.3e (must be exactly 0)\n", pk4);
    check(pk4 == 0.0, "seeded + undriven -> the solver stays EXACTLY at rest (no start-up transient)");
    check(Ri.anchor_plastic[0] == nr.anchor_plastic[0],
          "seeded + undriven -> the committed anchor state IS the seed (bitwise)");
    check(Ri.steps_completed == nst, "the identity run integrated every step");
}

}  // namespace

int main() {
    std::printf("Nonlinear dynamic (Newmark + per-step Newton) solver\n\n");
    test_linear_limit();
    std::printf("\n");
    test_quasi_static_mc();
    std::printf("\n");
    test_cyclic_hysteresis();
    std::printf("\n");
    test_parent_state_carry();
    if (g_failures == 0)
        std::printf("\nOK: nonlinear Newmark reduces to linear solve_newmark, and its quasi-static "
                    "limit matches the static elastoplastic solve while genuinely yielding\n");
    else
        std::printf("\n%d CHECK(S) FAILED\n", g_failures);
    return g_failures == 0 ? 0 : 1;
}
