#include <katai/analysis/dynamics_nonlinear.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/internal_forces.hpp>   // SHARED internal-force/tangent assembly (COMMON with static)
#include <katai/fem/elements/element_traits.hpp>

namespace katai::core {
namespace {

// Add a CONSTANT CSR matrix like M/C scaled by s into the builder (for assembling K_eff = K_T + a0·M + a1·C).
inline void add_scaled(const math::CsrMatrix& A, double s, math::SparseMatrixBuilder& b) {
    for (int r = 0; r < A.rows; ++r)
        for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
            b.add_entry(r, A.col_indices[p], s * A.values[p]);
}

// Nonlinear Newmark solver templated over element (tri6/tri15) × kinematics
// (plane-strain/axisym). Internal force/tangent from the SHARED detail::InternalForceAssembler
// (bit-for-bit the same assembly as the static solve_nonlinear). The outer time integration +
// in-step Newton are kinematics-independent.
template <class E, class Kin>
NewmarkNonlinearResult solve_newmark_nonlinear_impl(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const math::CsrMatrix& M, const math::CsrMatrix& C,
    const std::function<Eigen::VectorXd(int)>& force, double dt, int nsteps,
    const Eigen::VectorXd& u0, const Eigen::VectorXd& v0, const Eigen::VectorXd& a0_init,
    const LinearSolve& linear_solve, const NewmarkStepObserver& on_step,
    const std::vector<GaussState>& initial_state, const std::vector<char>& active_element,
    const Structures& structures, const std::vector<MaterialProfile>& profile,
    const StructuralInit& init_struct, const NewmarkNonlinearOptions& options) {
    using Asm = detail::InternalForceAssembler<E, Kin>;
    const int neq = dofs.equation_count();
    constexpr int n_gp = E::kGaussCount;

    // Newmark integration constants (Chopra Table 15.2.2 / Bathe) — SAME as the linear solve_newmark.
    const double beta = options.beta, gamma = options.gamma;
    const double a0 = 1.0 / (beta * dt * dt), a1 = gamma / (beta * dt), a2 = 1.0 / (beta * dt);
    const double a3 = 1.0 / (2.0 * beta) - 1.0, a4 = gamma / beta - 1.0;
    const double a5 = dt * (gamma / (2.0 * beta) - 1.0);
    const double a6 = dt * (1.0 - gamma), a7 = dt * gamma;
    (void)a4; (void)a5;  // a4/a5 are used in the linear effective load; not needed in the nonlinear residual branch.

    NewmarkNonlinearResult result;
    using Clock = std::chrono::steady_clock;
    const auto t_start = Clock::now();
    auto elapsed = [](Clock::time_point t0) {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };

    // --- Committed (start of step) + trial (step attempt) states. The soil σ is CARRIED
    // from the parent phase (initial_state) → stress-dependent constitutive law + seismic
    // stress recovery. The structural PLASTIC state (anchor U_p, geogrid ε_p, interface
    // slip, embedded slip) is SEEDED from the parent (init_struct, TRACK 1a) or starts from
    // zero. A wrong size does NOT fall through silently: hard error — silently starting from
    // zero would present a result that ignores the static preload (unsafe-sided) as normal.
    const size_t n_states = static_cast<size_t>(mesh.element_count) * n_gp;
    if (initial_state.size() == n_states) result.gauss_states = initial_state;
    else                                  result.gauss_states.assign(n_states, GaussState{});
    // Seeding by the SHARED rule (nonlinear_solver.cpp seed_structural_state; same as the
    // static solver and structural_internal_force): empty=0, full size=copy, wrong size=hard error.
    const bool carried = seed_structural_state(
        structures, init_struct, result.anchor_plastic, result.geogrid_plastic,
        result.interface_slip, result.interface5_slip, result.embedded_skin_slip,
        result.embedded_foot_slip, result.plate_plastic, result.plate5_plastic);
    const bool has_datum = init_struct.u_datum.size() > 0;
    if (has_datum && init_struct.u_datum.size() != neq)
        throw std::invalid_argument("solve_newmark_nonlinear: init_struct.u_datum size is not "
                                    "equation_count (" + std::to_string(init_struct.u_datum.size()) +
                                    " != " + std::to_string(neq) + ")");
    result.struct_state_carried = carried || has_datum;

    std::vector<GaussState>& committed = result.gauss_states;
    std::vector<double>& anchor_c = result.anchor_plastic;
    std::vector<double>& geogrid_c = result.geogrid_plastic;
    std::vector<double>& iface_c = result.interface_slip;
    std::vector<double>& iface5_c = result.interface5_slip;
    std::vector<double>& eskin_c = result.embedded_skin_slip;
    std::vector<double>& efoot_c = result.embedded_foot_slip;
    std::vector<double>& plate_c = result.plate_plastic;
    std::vector<double>& plate5_c = result.plate5_plastic;
    std::vector<GaussState> trial = committed;
    std::vector<double> anchor_t = anchor_c, geogrid_t = geogrid_c, iface_t = iface_c,
                        iface5_t = iface5_c, eskin_t = eskin_c, efoot_t = efoot_c,
                        plate_t = plate_c, plate5_t = plate5_c;

    Asm fasm(mesh, dofs, materials, active_element, structures, profile);
    // The dynamic step duration is in seconds; time-dependent constitutive models
    // (SoftSoilCreep) read on the day scale. Over an earthquake duration creep is
    // practically negligible, but the time bookkeeping is kept CORRECT (the real share is
    // passed through instead of silently zeroing).
    fasm.dt_day = dt / 86400.0;
    typename Asm::State astate;
    astate.committed = &committed;   astate.trial = &trial;
    astate.anchor_c = &anchor_c;     astate.anchor_t = &anchor_t;
    astate.geogrid_c = &geogrid_c;   astate.geogrid_t = &geogrid_t;
    astate.iface_c = &iface_c;       astate.iface_t = &iface_t;
    astate.iface5_c = &iface5_c;     astate.iface5_t = &iface5_t;
    astate.eskin_c = &eskin_c;       astate.eskin_t = &eskin_t;
    astate.efoot_c = &efoot_c;       astate.efoot_t = &efoot_t;
    astate.plate_c = &plate_c;       astate.plate_t = &plate_t;
    astate.plate5_c = &plate5_c;     astate.plate5_t = &plate5_t;
    const typename Asm::Ramp no_ramp;   // dynamic path: NO prescribed-displacement ramp (relative frame)

    // Material tangent mode (same choice as consolidation_plastic): HS → continuum (cheap +
    // robust), otherwise consistent (MC closed-form / LE exact).
    bool any_hs = false;
    for (const auto& mm : materials) if (mm.type == MaterialType::HardeningSoil) any_hs = true;
    const TangentMode tmode = any_hs ? TangentMode::kContinuum : TangentMode::kConsistent;

    // State carried through the time integration: u (committed displacement), v, a.
    Eigen::VectorXd u = u0.size() == neq ? u0 : Eigen::VectorXd::Zero(neq);
    Eigen::VectorXd v = v0.size() == neq ? v0 : Eigen::VectorXd::Zero(neq);
    // TOTAL displacement seen by the structural elements = parent datum + dynamic increment
    // (TRACK 1a). Since the soil assembly never reads u_free (only du; stress from committed
    // σ), the datum shifts only the structural elements. Without a datum u_struct == u
    // (v1 bit-for-bit).
    Eigen::VectorXd u_struct = has_datum ? Eigen::VectorXd(init_struct.u_datum + u) : u;

    // BASELINE f_int₀ = ∫Bᵀσ_committed (+ initial structural internal force) — the internal
    // force at initial equilibrium; subtracted from the residual (header: rest-start + D7
    // consistency). At u=u0, du=0, committed state. With a datum the structural term is the
    // parent's REAL static force (with the seeded plastic state + datum) → the residual is
    // still exactly zero at t=0 and the dynamic response evolves ON TOP of the parent state.
    const Eigen::VectorXd baseline =
        fasm.assemble(u_struct, Eigen::VectorXd::Zero(neq), /*build_tangent=*/false, tmode, astate,
                      no_ramp, nullptr, &result.timings);

    // Initial acceleration: if a0_init is given, use IT (M is never factored → singular M in
    // SSI is safe). Otherwise solve M·a = F(0) − C·v0 − (f_int₀ − baseline) = F(0) − C·v0
    // (the difference is zero at u=u0). If M is singular, linear_solve may not handle it →
    // the seismic path should pass a CLOSED-FORM a0_init (header).
    Eigen::VectorXd a = a0_init.size() == neq ? a0_init
                                              : linear_solve(M, force(0) - C * v);

    result.converged = true;
    // Initial state (stream): first the committed structural state (seed), then kinematics —
    // so the observer can evaluate the t=0 diagram with the parent's real state.
    if (options.on_commit) options.on_commit(0, result);
    if (on_step) on_step(0, 0.0, u, v, a);

    // Tangent builder for Newton — cleared + rebuilt every time step/iteration.
    math::SparseMatrixBuilder builder(neq);

    // Evaluate the residual for a given increment du (+ build the tangent into the builder
    // if requested):
    //   a_{n+1} = a0·du − a2·v − a3·a,   v_{n+1} = v + a6·a + a7·a_{n+1}
    //   r = F(t_{n+1}) − M·a_{n+1} − C·v_{n+1} − (f_int(u+du) − baseline)
    // a1n/v1n (the n+1 values) are returned as outputs (for the commit).
    auto residual = [&](const Eigen::VectorXd& du, const Eigen::VectorXd& F1, bool build_tangent,
                        math::SparseMatrixBuilder* b, Eigen::VectorXd& a1n, Eigen::VectorXd& v1n) {
        const Eigen::VectorXd f_int =
            fasm.assemble(u_struct, du, build_tangent, tmode, astate, no_ramp, b, &result.timings);
        a1n = a0 * du - a2 * v - a3 * a;
        v1n = v + a6 * a + a7 * a1n;
        return Eigen::VectorXd(F1 - M * a1n - C * v1n - (f_int - baseline));
    };

    for (int step = 0; step < nsteps; ++step) {
        const Eigen::VectorXd F1 = force(step + 1);
        Eigen::VectorXd du = Eigen::VectorXd::Zero(neq);
        Eigen::VectorXd a1n = a, v1n = v;   // n+1 values (committed on convergence)
        bool step_ok = false;
        int stall = 0;
        for (int it = 0; it < options.max_newton; ++it) {
            builder.clear();
            const Eigen::VectorXd r = residual(du, F1, /*build_tangent=*/true, &builder, a1n, v1n);
            const double rnorm = r.norm();
            // Convergence scale: external + inertial force magnitude (absolute floor if both are zero).
            const double ref = std::max({F1.norm(), (M * a1n).norm(), 1e-30});
            ++result.total_newton_iterations;
            if (rnorm <= options.newton_tol * ref + 1e-9 * (baseline.norm() + 1.0)) { step_ok = true; break; }

            // K_eff = K_T + a0·M + a1·C (the tangent is already in the builder; add the CONSTANT M,C contributions).
            add_scaled(M, a0, builder);
            add_scaled(C, a1, builder);
            const auto t_csr = Clock::now();
            const math::CsrMatrix Keff = builder.build();
            result.timings.csr_build += elapsed(t_csr);
            const auto t_lin = Clock::now();
            const Eigen::VectorXd delta = linear_solve(Keff, r);
            result.timings.linear_solve += elapsed(t_lin);
            ++result.timings.n_solve;

            // Backtracking line search (Armijo): the full Newton step can overshoot in
            // plasticity (interface slip / soil yield → near-singular tangent); alpha is
            // halved until ||r|| drops. Same pattern as solve_nonlinear — on the nonlinear
            // dynamic residual.
            double alpha = 1.0;
            bool improved = false;
            Eigen::VectorXd a1t, v1t;
            for (int ls = 0; ls < 12; ++ls) {
                const Eigen::VectorXd rt = residual(du + alpha * delta, F1, false, nullptr, a1t, v1t);
                if (rt.norm() < (1.0 - 1.0e-4 * alpha) * rnorm) { improved = true; break; }
                alpha *= 0.5;
            }
            stall = improved ? 0 : stall + 1;
            du += alpha * delta;
            if (stall >= 4) break;   // no sustained descent → the step is not converging (Δt fixed; no sub-stepping)
        }
        if (!step_ok) { result.converged = false; break; }

        // The step converged → COMMIT the state (Gauss + structural + Newmark u,v,a).
        u += du;
        u_struct += du;
        v = v1n;
        a = a1n;
        committed = trial;
        anchor_c = anchor_t; geogrid_c = geogrid_t; iface_c = iface_t; iface5_c = iface5_t;
        eskin_c = eskin_t; efoot_c = efoot_t; plate_c = plate_t; plate5_c = plate5_t;
        ++result.steps_completed;
        if (options.on_commit) options.on_commit(step + 1, result);
        if (on_step) on_step(step + 1, (step + 1) * dt, u, v, a);
    }

    result.timings.total = elapsed(t_start);
    return result;
}

}  // namespace

NewmarkNonlinearResult solve_newmark_nonlinear(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const math::CsrMatrix& M, const math::CsrMatrix& C,
    const std::function<Eigen::VectorXd(int)>& force, double dt, int nsteps,
    const Eigen::VectorXd& u0, const Eigen::VectorXd& v0, const Eigen::VectorXd& a0_init,
    const LinearSolve& linear_solve, const NewmarkStepObserver& on_step,
    const std::vector<GaussState>& initial_state, const std::vector<char>& active_element,
    const Structures& structures, const std::vector<MaterialProfile>& profile,
    const StructuralInit& init_struct, const NewmarkNonlinearOptions& options) {
    const bool axi = options.kinematics == Kinematics::Axisymmetric;
    if (mesh.nodes_per_element == Tri15Element::kNodeCount) {
        if (axi)
            return solve_newmark_nonlinear_impl<Tri15Element, detail::AxisymKin>(
                mesh, dofs, materials, M, C, force, dt, nsteps, u0, v0, a0_init, linear_solve,
                on_step, initial_state, active_element, structures, profile, init_struct, options);
        return solve_newmark_nonlinear_impl<Tri15Element, detail::PlaneStrainKin>(
            mesh, dofs, materials, M, C, force, dt, nsteps, u0, v0, a0_init, linear_solve,
            on_step, initial_state, active_element, structures, profile, init_struct, options);
    }
    if (axi)
        return solve_newmark_nonlinear_impl<Tri6Element, detail::AxisymKin>(
            mesh, dofs, materials, M, C, force, dt, nsteps, u0, v0, a0_init, linear_solve,
            on_step, initial_state, active_element, structures, profile, init_struct, options);
    return solve_newmark_nonlinear_impl<Tri6Element, detail::PlaneStrainKin>(
        mesh, dofs, materials, M, C, force, dt, nsteps, u0, v0, a0_init, linear_solve,
        on_step, initial_state, active_element, structures, profile, init_struct, options);
}

}  // namespace katai::core
