#pragma once
// NONLINEAR dynamic (seismic) solver — Newmark-β time integration + per-step Newton-Raphson.
//
//     M u'' + C u' + f_int(u) = F(t),      C = αM + βK₀ (Rayleigh, CONSTANT with initial stiffness)
//
// The nonlinear generalization of the linear `solve_newmark` (dynamics.hpp): the internal
// force is no longer K·u but f_int(u) from the constitutive return mapping → the interface
// CAN slip, the geogrid can slacken, the anchor can yield, the soil can plastify (PLAXIS
// Dynamics is FULLY NONLINEAR too — Ref §11.10.4: the SAME convergence criteria as plastic
// analysis). At each time step the residual
//     r = F(t_{n+1}) − M·a_{n+1} − C·v_{n+1} − (f_int(u_{n+1}) − f_int₀)
// is driven to zero by Newton; a_{n+1},v_{n+1} are written in terms of u via the Newmark
// relations →
//     K_eff = K_T(u) + a0·M + a1·C   (REFACTORED every iteration; K_eff generally NONSYMMETRIC).
//
// BASELINE (f_int₀) — CONSISTENCY with the at-rest start + D7 superposition:
// The phase starts from REST (u0=v0=0, relative frame) and CONTINUES from the parent
// phase's stress state (initial_state = committed σ). But F(t) = −M·r·a_g carries only the
// DYNAMIC body force; it does NOT include the parent's static loads (gravity, K0). Those
// static loads already balance the committed σ → the initial internal force
// f_int₀ = ∫Bᵀσ_committed must be SUBTRACTED from the residual, otherwise there is a
// spurious imbalance at t=0 (the committed σ would set off motion). With the subtraction,
// a(0) = −r·a_g(0) holds EXACTLY at t=0 and the dynamic response superposes onto the static
// state as in D7. (Same pattern as the consolidation_plastic_impl Bbase.)
//
// v1 LIMITS (honest): (a) C = αM + βK₀ CONSTANT — the βK term of Rayleigh uses the INITIAL
// (elastic) stiffness (the common choice; to keep damping unchanged as the material softens
// and to keep the tangent simple). (b) The time step Δt is CONSTANT — if a step fails to
// converge, the solution stops there (converged=false), NO time sub-stepping. (c) NO
// modified Newton (freezing the tangent within a step) — every iteration refactors.
//
// Assembly via `detail::InternalForceAssembler` (analysis/internal_forces.hpp) → BIT-FOR-BIT
// the same constitutive/structural rules as the static solve_nonlinear (single source, no
// drift). Verification: test_dynamics_nonlinear (linear limit → solve_newmark exact;
// quasi-static limit → solve_nonlinear exact). Formulation:
// docs/references/dynamic-seismic-formulation.md.

#include <functional>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/dynamics.hpp>           // NewmarkStepObserver
#include <katai/analysis/nonlinear_solver.hpp>   // Structures, LinearSolve, Kinematics, NewtonResult::Timings
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>    // GaussState, MaterialModel, MaterialProfile
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// PARENT STRUCTURAL STATE: `StructuralInit` is now defined in analysis/nonlinear_solver.hpp
// (with the static generalization of Track 1a BOTH solvers use the same carry-over
// contract; single source). Its meaning in dynamics: structural elements are evaluated on
// the sum (u_datum + u_dynamic), from the parent's committed plastic state → the Coulomb
// cap / anchor capacity / geogrid slack is checked on the TOTAL (static+dynamic) effect;
// slip/yield onset is asymmetric and correct.

struct NewmarkNonlinearResult;   // for the on_commit signature (defined below)

struct NewmarkNonlinearOptions {
    double gamma = 0.5;         // Newmark γ (½ = no algorithmic damping)
    double beta = 0.25;         // Newmark β (¼ = average acceleration, unconditionally stable)
    int max_newton = 30;        // maximum Newton iterations per step
    double newton_tol = 1e-6;   // relative convergence threshold on ||r|| / ref
    Kinematics kinematics = Kinematics::PlaneStrain;
    // Streams the committed structural plastic state after each step is COMMITTED (first
    // call step=0 = seed/initial state); invoked IMMEDIATELY BEFORE on_step. So the
    // envelope post-processor can evaluate with the solver's REAL constitutive state
    // (committed slip / U_p / ε_p) (D6b rule: the post-processor's constitutive law must
    // MATCH the solver's). The passed reference is live for the duration of the solve.
    std::function<void(int step, const NewmarkNonlinearResult& committed)> on_commit;
};

struct NewmarkNonlinearResult {
    // Last committed state (for seismic stress recovery + phase chaining — unlike the
    // linear branch, the nonlinear branch CARRIES the plastic evolution). Sizes as in
    // solve_nonlinear.
    std::vector<GaussState> gauss_states;
    std::vector<double> anchor_plastic;
    std::vector<double> geogrid_plastic;
    std::vector<double> interface_slip;
    std::vector<double> interface5_slip;
    std::vector<double> embedded_skin_slip;
    std::vector<double> embedded_foot_slip;
    std::vector<double> plate_plastic;    // plate M-N hinge state [ε_p,κ_p]×Gauss (×6 / element)
    std::vector<double> plate5_plastic;   // 5-node plate (×10 / element)
    bool converged = true;            // did all time steps converge
    int steps_completed = 0;          // number of converged time steps
    int total_newton_iterations = 0;  // sum over all steps
    NewtonResult::Timings timings;    // assembly/solve time breakdown (same as solve_nonlinear)
    // Whether the parent structural datum/state was seeded (at least one StructuralInit
    // member given). Reporting decisions (superposition off, utilisation from the total)
    // must look at THIS — at what the solver actually did, not at the caller's intent.
    bool struct_state_carried = false;
};

// Newmark-β + per-step Newton. M, C CONSTANT (assembled; structural contributions included;
// C = αM+βK₀). `force(k)` = external force at t=k·Δt (for seismic, −M·r·a_g(k·Δt) + the
// free-field driver if any). u0/v0/a0_init: initial displacement/velocity/acceleration
// (seismic rest-start: u0=v0=0, a0_init=−r·a_g(0); if a0_init is not given it is solved
// from M·a=F(0)−C·v0−f_int₀ — this requires M to be NON-SINGULAR, which it usually is not
// in SSI → pass the closed-form a0_init). `linear_solve`: solves K_eff·δ=r (refactor every
// iteration; NONSYMMETRIC tangent → PARDISO mtype=11). `on_step`: each step (step,t,u,v,a)
// is STREAMED (constant memory). `initial_state`: committed EFFECTIVE Gauss state (parent
// phase; the f_int₀ baseline comes from it). active_element, structures, profile: as in
// solve_nonlinear. `init_struct`: parent STRUCTURAL datum + plastic state (TRACK 1a; empty
// = v1 bit-for-bit). Element type from mesh.nodes_per_element.
NewmarkNonlinearResult solve_newmark_nonlinear(
    const mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& materials,
    const math::CsrMatrix& M, const math::CsrMatrix& C,
    const std::function<Eigen::VectorXd(int)>& force, double dt, int nsteps,
    const Eigen::VectorXd& u0, const Eigen::VectorXd& v0, const Eigen::VectorXd& a0_init,
    const LinearSolve& linear_solve, const NewmarkStepObserver& on_step,
    const std::vector<GaussState>& initial_state = {},
    const std::vector<char>& active_element = {},
    const Structures& structures = {},
    const std::vector<MaterialProfile>& profile = {},
    const StructuralInit& init_struct = {},
    const NewmarkNonlinearOptions& options = {});

}  // namespace katai::core
