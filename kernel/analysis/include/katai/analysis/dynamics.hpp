#pragma once
// Structural/soil dynamics — the semi-discrete equation of motion  M u'' + C u' + K u = F(t)
// with Newmark-β time integration + Rayleigh damping. DIMENSION-AGNOSTIC: operates on the
// assembled M, C, K matrices; the same integrator drives both the 1D free-field column
// (verified here) and (a later phase) the 2D FE system. Formulation + verification:
// docs/references/dynamic-seismic-formulation.md (Newmark 1959; Chopra §15; Kramer §7).
// factor-once-solve-many; MKL-agnostic (the ConsolidationSolveFactory pattern).

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <katai/math/sparse_matrix.hpp>

namespace katai::core {

// Rayleigh damping C = α·M + β·K. Modal ratio ξ_n = α/(2ω_n) + β·ω_n/2. (α, β) is solved
// from two target (frequency [Hz], damping ratio) pairs (Chopra §11.4). With equal targets:
// ξ≈target within the band between them, rising outside.
struct RayleighCoefficients { double alpha = 0.0, beta = 0.0; };

inline RayleighCoefficients rayleigh_from_modes(double f1, double xi1, double f2, double xi2) {
    constexpr double kPi = 3.14159265358979323846;
    const double w1 = 2.0 * kPi * f1, w2 = 2.0 * kPi * f2;
    RayleighCoefficients r;
    const double den = w2 * w2 - w1 * w1;
    if (w1 <= 0.0 || w2 <= 0.0 || std::fabs(den) < 1e-30) return r;
    r.alpha = 2.0 * w1 * w2 * (xi1 * w2 - xi2 * w1) / den;
    r.beta = 2.0 * (xi2 * w2 - xi1 * w1) / den;
    return r;
}

// Linear solve factory for constant-Δt Newmark: given the effective stiffness K_eff,
// returns a closure that solves K_eff x = b (factor ONCE, back-substitute every step). If
// empty, dense Eigen LU (the MKL-free reference path, used in tests). Same contract as
// ConsolidationSolveFactory.
using DynamicsSolveFactory =
    std::function<std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(const math::CsrMatrix&)>;

struct NewmarkResult {
    std::vector<double> t;                 // time of each recorded state (size nsteps+1)
    std::vector<Eigen::VectorXd> u, v, a;  // displacement / velocity / acceleration per state
};

namespace detail {
inline Eigen::MatrixXd csr_to_dense(const math::CsrMatrix& A) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows, A.cols);
    for (int r = 0; r < A.rows; ++r)
        for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
            D(r, A.col_indices[p]) = A.values[p];
    return D;
}
}  // namespace detail

// Newmark-β linear time integration (§2). M u''+C u'+K u = F(t); γ=½,β=¼ (average
// acceleration, unconditionally stable, no algorithmic damping). `force(k)` = external
// force at t=k·Δt (for seismic, −M r a_g(k·Δt)). M,C,K share the same n×n size (for
// undamped pass an n×n zero matrix as C). Initial u0,v0; acceleration a0 from
// M a0 = F(0)−C v0−K u0. Returns: u/v/a of every state 0..nsteps (relative — in seismic the
// total acceleration is a_total = a + a_g). If `factory` is empty, dense LU (test path).
// Per-step observer (optional). If given, solve_newmark STREAMS each state (step, time, u, v, a) to it
// and does NOT accumulate the O(nsteps x ndof) history in NewmarkResult -- a large 2D mesh x many time
// steps would exhaust memory (bad_alloc -> crash in the GUI worker). If null (default), the full history
// is returned as before (existing tests / callers unchanged).
using NewmarkStepObserver =
    std::function<void(int step, double t, const Eigen::VectorXd& u,
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a)>;

// a0_init (optional): the INITIAL ACCELERATION. If not given, M a = F(0)−C v0−K u0 is
// solved — this requires M to be NON-SINGULAR. In a soil-structure (SSI) system M is
// USUALLY SINGULAR: massless structural DOFs (the plate rotational DOF; all DOFs of a w=0
// weightless plate) leave ZERO ROWS in M → factoring a singular M breaks the solver (in
// PARDISO SPD an access violation = NOT a C++ exception, catch(...) does NOT catch it).
// The time stepping uses M ONLY as a matvec and K_eff = K+a0M+a1C is non-singular (K gives
// every structural DOF stiffness) → a singular M is NOT a problem; only this initial solve
// is. For seismic base excitation starting from rest (u0=v0=0, F(0)=−M·r·a_g(0)) the right
// answer is known in CLOSED FORM:
//     a(0) = −r·a_g(0)     →  M a(0) = −M r a_g(0) = F(0), EXACT for EVERY M (even singular).
// The physics: at t=0 the system is not yet deformed; in the relative frame everything
// accelerates rigidly at −a_g(0) (total acceleration = a_rel + a_g = 0). When given, M is
// never factored.
inline NewmarkResult solve_newmark(
    const math::CsrMatrix& M, const math::CsrMatrix& C, const math::CsrMatrix& K,
    const std::function<Eigen::VectorXd(int)>& force, double dt, int nsteps,
    const Eigen::VectorXd& u0, const Eigen::VectorXd& v0, double gamma = 0.5, double beta = 0.25,
    const DynamicsSolveFactory& factory = {}, const NewmarkStepObserver& on_step = {},
    const Eigen::VectorXd& a0_init = {}) {
    const int n = M.rows;
    // Newmark integration constants (Chopra Table 15.2.2 / Bathe).
    const double a0 = 1.0 / (beta * dt * dt), a1 = gamma / (beta * dt), a2 = 1.0 / (beta * dt);
    const double a3 = 1.0 / (2.0 * beta) - 1.0, a4 = gamma / beta - 1.0;
    const double a5 = dt * (gamma / (2.0 * beta) - 1.0);
    const double a6 = dt * (1.0 - gamma), a7 = dt * gamma;

    // Factored solver for a CsrMatrix: the factory if given, else factor a dense LU once.
    auto make_solver = [&](const math::CsrMatrix& A)
        -> std::function<Eigen::VectorXd(const Eigen::VectorXd&)> {
        if (factory) return factory(A);
        auto lu = std::make_shared<Eigen::FullPivLU<Eigen::MatrixXd>>(detail::csr_to_dense(A));
        return [lu](const Eigen::VectorXd& b) { return lu->solve(b); };
    };

    // Effective stiffness K_eff = K + a0·M + a1·C (sum the three patterns → factor ONCE).
    math::SparseMatrixBuilder builder(n, n);
    auto add_scaled = [&](const math::CsrMatrix& A, double s) {
        for (int r = 0; r < A.rows; ++r)
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
                builder.add_entry(r, A.col_indices[p], s * A.values[p]);
    };
    add_scaled(K, 1.0);
    add_scaled(M, a0);
    add_scaled(C, a1);
    const auto solve_keff = make_solver(builder.build());

    // Initial state + acceleration: use a0_init if given (M is never factored → singular M
    // safe; see the note above), else solve M a = F(0) − C v0 − K u0 (old path, M must be
    // non-singular).
    Eigen::VectorXd u = u0, v = v0;
    Eigen::VectorXd a = a0_init.size() == n ? a0_init
                                            : make_solver(M)(force(0) - C * v - K * u);

    NewmarkResult R;
    auto record = [&](int step, double t) {
        if (on_step) { on_step(step, t, u, v, a); return; }   // stream (constant memory)
        R.t.push_back(t); R.u.push_back(u); R.v.push_back(v); R.a.push_back(a);
    };
    record(0, 0.0);
    for (int step = 0; step < nsteps; ++step) {
        // Effective load F_eff = F_{n+1} + M(a0 u+a2 v+a3 a) + C(a1 u+a4 v+a5 a).
        const Eigen::VectorXd Feff =
            force(step + 1) + M * (a0 * u + a2 * v + a3 * a) + C * (a1 * u + a4 * v + a5 * a);
        const Eigen::VectorXd u_new = solve_keff(Feff);
        const Eigen::VectorXd a_new = a0 * (u_new - u) - a2 * v - a3 * a;
        const Eigen::VectorXd v_new = v + a6 * a + a7 * a_new;
        u = u_new; v = v_new; a = a_new;
        record(step + 1, (step + 1) * dt);
    }
    return R;
}

}  // namespace katai::core
