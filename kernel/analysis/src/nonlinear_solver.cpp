#include <katai/analysis/nonlinear_solver.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include <katai/analysis/internal_forces.hpp>  // SHARED internal-force/tangent assembly (static+dynamic)
#include <katai/fem/assembly/assembler.hpp>  // expand_to_full
#include <katai/fem/elements/axisymmetric.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/math/solve_error.hpp>
#include <katai/math/thread_pool.hpp>

namespace katai::core {
namespace {

// gather_element, PlaneStrainKin, AxisymKin and the internal-force/tangent assembly now live
// in the SHARED analysis/internal_forces.hpp (katai::core::detail): the static Newton (this
// file) and the nonlinear dynamic Newmark+Newton (dynamics_nonlinear.cpp) use the SAME
// assembly → no drift.

// Turn one iterate's raw measurement into the convergence family (Scientific Manual §9.1).
// Called at every ACCEPTED iterate, so what it leaves behind is what the last such iterate
// satisfied -- including, when the solve fails, the last one that was reached.
void measure_convergence(const detail::LocalErrorProbe& p, const Eigen::VectorXd& f_int,
                         const Eigen::VectorXd& residual, double rnorm, double cf_norm,
                         double ref, double rtol, NewtonResult::Convergence& cv) {
    cv.measured = true;
    cv.tolerated = rtol;

    // CSP, Eq. 7-22. Guarded twice, for two different reasons: an increment that has not moved
    // yet has no energy to take a ratio of (report the elastic value, 1), and an UNLOADING
    // increment can produce a negative or greater-than-one ratio, which is not a stiffness
    // measure -- the parameter is defined on the interval and is clamped to it.
    cv.csp = p.energy_elastic > 0.0
                 ? std::min(1.0, std::max(0.0, p.energy_total / p.energy_elastic))
                 : 1.0;

    // Eq. 9-1. The denominator is the force scale the model itself is carrying, plus the
    // standing load from the phases before this one, weighted by how much stiffness is left.
    // The floor is 1e-6 of the phase's own load scale: it bites only when nothing is resisting
    // yet, and then a huge relative error is the honest reading.
    const double denom = std::max(f_int.norm() + cv.csp * cf_norm, 1e-6 * ref);
    cv.force_error = rnorm / denom;

    // Eq. 9-3/9-4. Only structural elements with rotational freedom have one.
    cv.has_moment = !p.moment_eq.empty() && p.m_ref > 0.0;
    if (cv.has_moment) {
        double worst = 0.0;
        for (int eq : p.moment_eq) worst = std::max(worst, std::fabs(residual(eq)));
        cv.moment_error = worst / p.m_ref;
    } else {
        cv.moment_error = 0.0;
    }

    cv.plastic_points = p.plastic;
    cv.plastic_inaccurate = p.plastic_inaccurate;
    cv.elastic_points = p.elastic_total;
    cv.nl_elastic_points = p.nl_elastic;
    cv.nl_elastic_inaccurate = p.nl_elastic_inaccurate;
    cv.worst_plastic_error = p.worst_plastic;
    cv.worst_nl_elastic_error = p.worst_nl_elastic;

    cv.iface_points = p.iface_plastic;
    cv.iface_inaccurate = p.iface_plastic_inaccurate;
    cv.worst_iface_error = p.worst_iface;

    // Eq. 9-9. The denominator is the source's, floor and all: the constitutive foot forces, or
    // one per cent of the declared capacities, or 1.0 -- whichever is largest. The middle term is
    // what stops a pile that carries almost nothing at its toe from reporting a huge relative
    // error while the capacity it was given says the force is negligible either way.
    cv.feet = p.feet;
    cv.foot_force_error =
        p.feet == 0 ? 0.0
                    : p.foot_num / std::max({p.foot_den_c, 0.01 * p.foot_den_max, 1.0});
}

// Newton solver templated over element (tri6/tri15) and kinematics (plane strain/axisym);
// uniform code, zero vtables. Internal-force/tangent assembly is delegated to the shared
// detail::InternalForceAssembler. The solve_nonlinear above selects.
template <class E, class Kin>
NewtonResult solve_nonlinear_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const std::vector<MaterialModel>& materials,
                                  const Eigen::VectorXd& f_ext,
                                  const LinearSolve& linear_solve,
                                  const NewtonOptions& options,
                                  const std::vector<GaussState>& initial_state,
                                  const std::vector<char>& active_element,
                                  const Structures& structures,
                                  const Eigen::VectorXd& presc,
                                  const Eigen::VectorXd& constant_force,
                                  const std::vector<MaterialProfile>& profile,
                                  const StructuralInit& init_struct) {
    const int neq = dofs.equation_count();
    constexpr int n_gp = E::kGaussCount;

    NewtonResult result;
    result.displacement = Eigen::VectorXd::Zero(dofs.total_dofs());
    const size_t n_states = static_cast<size_t>(mesh.element_count) * n_gp;
    // Committed Gauss state: seed with the supplied prestress (K0 procedure / previous
    // stage) when provided, else start from zero stress.
    if (initial_state.size() == n_states)
        result.gauss_states = initial_state;
    else
        result.gauss_states.assign(n_states, GaussState{});
    result.converged = true;

    // Converged total displacement (free) and the step trial state. trial is initialized
    // from committed: the states of passive elements (which never enter the assembly) are
    // preserved through commits too (for stress carry-over in staged construction).
    Eigen::VectorXd u_free = Eigen::VectorXd::Zero(neq);
    // TOTAL displacement seen by the structural elements = parent datum + this phase's u
    // (Track 1a). Since the soil assembly never reads u_free (only du), the datum shifts
    // only the structural elements; without a datum u_struct == u_free (old behaviour
    // bit-for-bit).
    const bool has_datum = init_struct.u_datum.size() > 0;
    if (has_datum && init_struct.u_datum.size() != neq)
        throw std::invalid_argument("solve_nonlinear: init_struct.u_datum size is not "
                                    "equation_count (" + std::to_string(init_struct.u_datum.size()) +
                                    " != " + std::to_string(neq) + ")");
    Eigen::VectorXd u_struct = has_datum ? Eigen::VectorXd(init_struct.u_datum)
                                         : Eigen::VectorXd::Zero(neq);
    std::vector<GaussState>& committed = result.gauss_states;
    std::vector<GaussState> trial = committed;

    // Structural path-dependent state (elastoplastic): anchor plastic elongation U_p,
    // geogrid plastic axial ε_p (×2 Gauss). committed/trial like gauss_states; SEEDED from
    // the parent (init_struct, the static generalization of Track 1a) or starts from zero
    // (empty = old behaviour).
    seed_structural_state(structures, init_struct, result.anchor_plastic, result.geogrid_plastic,
                          result.interface_slip, result.interface5_slip,
                          result.embedded_skin_slip, result.embedded_foot_slip,
                          result.plate_plastic, result.plate5_plastic);
    std::vector<double>& anchor_committed = result.anchor_plastic;
    std::vector<double>& geogrid_committed = result.geogrid_plastic;
    std::vector<double>& interface_committed = result.interface_slip;
    std::vector<double>& interface5_committed = result.interface5_slip;
    std::vector<double>& eskin_committed = result.embedded_skin_slip;
    std::vector<double>& efoot_committed = result.embedded_foot_slip;
    std::vector<double>& plate_committed = result.plate_plastic;
    std::vector<double>& plate5_committed = result.plate5_plastic;
    std::vector<double> anchor_trial = anchor_committed;
    std::vector<double> geogrid_trial = geogrid_committed;
    std::vector<double> interface_trial = interface_committed;
    std::vector<double> interface5_trial = interface5_committed;
    std::vector<double> eskin_trial = eskin_committed;
    std::vector<double> efoot_trial = efoot_committed;
    std::vector<double> plate_trial = plate_committed;
    std::vector<double> plate5_trial = plate5_committed;

    const double f_ext_norm = f_ext.norm();
    const double rtol = options.tolerance;

    // Nonzero Dirichlet (prescribed displacement ū): ramped 0→ū together with the external
    // load. cur_lambda (committed) and cur_target (step target) are read from inside
    // assemble; the fixed DOF's increment this step is (cur_target−cur_lambda)·ū, its total
    // value cur_target·ū.
    const bool has_presc = presc.size() == dofs.total_dofs();
    const bool has_cf = constant_force.size() == neq;
    double cur_lambda = 0.0, cur_target = 0.0;

    // HYBRID HS TANGENT POLICY: increments start with the fast CONTINUUM tangent; if an
    // increment fails to converge, before halving dlam the SAME increment is retried once
    // with the CONSISTENT (numerical, FD) tangent, and the rest of the solve stays
    // consistent (the hard regime has been entered — footing edge/low confining pressure;
    // measured: continuum stalls there at lf=0.84, consistent carries to 1.0, but is ~2.5×
    // more expensive in the easy regime). Active only when an HS material is present →
    // the LE/MC path is the old behaviour EXACTLY (TangentMode does not affect them, the
    // retry branch never opens).
    bool has_hs = false;
    for (const auto& mm : materials)
        if (mm.type == MaterialType::HardeningSoil) { has_hs = true; break; }
    bool hs_consistent_mode = false;

    using Clock = std::chrono::steady_clock;
    const auto t_start = Clock::now();
    auto elapsed = [](Clock::time_point t0) {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };

    // Internal force / consistent tangent are assembled by the SHARED assembler
    // (analysis/internal_forces.hpp): soil constitutive return mapping + every embedded structural
    // element, identical rules to the nonlinear-dynamic (Newmark+Newton) solver -> single source, no
    // drift. The assembler owns the per-solve element topology + fe/ke buffers (built once in its
    // constructor), so the old inline gather/buffers move into it.
    detail::InternalForceAssembler<E, Kin> fasm(mesh, dofs, materials, active_element,
                                                structures, profile);
    typename detail::InternalForceAssembler<E, Kin>::State astate;
    astate.committed = &committed;           astate.trial = &trial;
    astate.anchor_c = &anchor_committed;     astate.anchor_t = &anchor_trial;
    astate.geogrid_c = &geogrid_committed;   astate.geogrid_t = &geogrid_trial;
    astate.iface_c = &interface_committed;   astate.iface_t = &interface_trial;
    astate.iface5_c = &interface5_committed; astate.iface5_t = &interface5_trial;
    astate.eskin_c = &eskin_committed;       astate.eskin_t = &eskin_trial;
    astate.efoot_c = &efoot_committed;       astate.efoot_t = &efoot_trial;
    astate.plate_c = &plate_committed;       astate.plate_t = &plate_trial;
    astate.plate5_c = &plate5_committed;     astate.plate5_t = &plate5_trial;

    // Thin adapter so the Newton / line-search call sites below are unchanged: f_int (and, when
    // build_tangent, K_T into *builder) at increment du_free, evaluated at u_free + du_free. The
    // prescribed-displacement ramp (static only) is (cur_target - cur_lambda) * presc; hs_consistent_
    // mode selects the HS tangent (continuum vs consistent). Bit-for-bit the old inline assemble.
    auto assemble = [&](const Eigen::VectorXd& du_free, bool build_tangent,
                        math::SparseMatrixBuilder* builder,
                        detail::LocalErrorProbe* probe = nullptr) {
        astate.local = probe;
        const TangentMode tmode =
            hs_consistent_mode ? TangentMode::kConsistent : TangentMode::kContinuum;
        typename detail::InternalForceAssembler<E, Kin>::Ramp ramp;
        // factor drives the SOIL (an increment, against committed Gauss states); total drives the
        // STRUCTURAL elements, which are total-displacement formulations and need the whole share
        // of u_bar standing at the end of this step. See Ramp in internal_forces.hpp.
        if (has_presc) {
            ramp.presc = &presc;
            ramp.factor = cur_target - cur_lambda;
            ramp.total = cur_target;
        }
        // The time share is proportional to this increment's Δλ (SoftSoilCreep; time_interval=0 → 0, old path).
        fasm.dt_day = options.time_interval * std::max(0.0, cur_target - cur_lambda);
        fasm.substep_tol = options.substep_tolerance;
        return fasm.assemble(u_struct, du_free, build_tangent, tmode, astate, ramp, builder,
                             &result.timings);
    };

    const double cf_norm = has_cf ? constant_force.norm() : 0.0;
    const double ref = std::max({f_ext_norm, cf_norm, 1.0});
    const bool debug = std::getenv("KATAI_NL_DEBUG") != nullptr;
    // An increment must satisfy the LOCAL criteria as well as the global force residual before
    // it is called converged -- which is what the source of these criteria does. ON by default
    // since 2026-08-24; see NewtonOptions::enforce_local_criteria for the measurement that
    // decided it. The two environment variables are the whole-run overrides a study uses, and
    // the option is what a test and the phase's numerical controls drive.
    const bool enforce_local =
        std::getenv("KATAI_CONV_LOCAL") != nullptr
            ? true
            : (std::getenv("KATAI_CONV_NOLOCAL") != nullptr ? false
                                                            : options.enforce_local_criteria);

    // Adaptive (automatic) load incrementation. The external load is advanced from
    // 0 to f_ext by increments d_lambda; an increment that fails to converge is
    // retried with a halved size (sub-stepping) instead of aborting the analysis.
    // The collapse signal -- non-convergence even at the minimum increment -- is
    // thus decoupled from the (arbitrary) initial step count, making the limit/FoS
    // result robust across element orders and problem stiffness (cf. Crisfield;
    // PLAXIS/ABAQUS automatic stepping). When every increment converges at the
    // initial size (no cutback) the scheme reduces exactly to fixed N-step loading,
    // so the tri6 benchmarks are unchanged.
    result.iteration_limit = options.max_iterations;
    const double init_dlam = 1.0 / options.load_steps;
    const double min_dlam = init_dlam / 8.0;  // smaller -> declare collapse
    double lambda = 0.0;        // committed load fraction
    double dlam = init_dlam;    // current increment size
    result.converged = true;

    // The global tangent assembly buffers live ACROSS ITERATIONS: since the COO entry order
    // is deterministic (the assembly loops run over a fixed topology), the CSR pattern is
    // built once and later iterations only accumulate values (build_cached; falls back to a
    // full build by itself if the signature mismatches).
    math::SparseMatrixBuilder builder(neq);
    math::CsrPatternCache kt_cache;

    // The convergence family's memory (Scientific Manual §9.1). prev_sc carries sigma_c,j-1
    // and prev_deps carries Delta-eps_{j-1}; both are restarted at every increment ATTEMPT,
    // including a retried one, because sigma_0 is the increment's own committed state and a
    // cut-back increment starts its iteration afresh.
    detail::LocalErrorProbe probe;
    std::vector<GaussState> prev_sc;
    std::vector<double> prev_deps;
    std::vector<detail::PrevPoint> prev_iface, prev_iface5, prev_skin, prev_foot;
    probe.prev_sigma_c = &prev_sc;
    probe.prev_deps = &prev_deps;
    probe.prev_iface = &prev_iface;
    probe.prev_iface5 = &prev_iface5;
    probe.prev_skin = &prev_skin;
    probe.prev_foot = &prev_foot;
    probe.tolerated = rtol;
    result.convergence.tolerated = rtol;

    while (lambda < 1.0 - 1e-12) {
        if (dlam > 1.0 - lambda) dlam = 1.0 - lambda;  // do not overshoot
        const double target_lambda = lambda + dlam;
        Eigen::VectorXd target = target_lambda * f_ext;
        if (has_cf) target += constant_force;  // constant (non-ramped) load (geostatic gravity)
        cur_lambda = lambda;          // for the prescribed-displacement ramp (read by assemble)
        cur_target = target_lambda;

        Eigen::VectorXd du_free = Eigen::VectorXd::Zero(neq);
        prev_sc = committed;  // sigma_c,0 = sigma_0, where Fig. 9-1 starts
        std::fill(prev_deps.begin(), prev_deps.end(), 0.0);
        // A traction has no committed value stored anywhere -- the structural elements are
        // total-displacement formulations, so the state at the start of an increment is
        // whatever the first assembly of that increment computes. It is recorded there.
        bool first_iterate = true;
        // The last W residual norms of THIS increment, most recent last: the line search's
        // yardstick when the window is open. It starts empty at every increment (including a
        // retried one), so the memory never crosses a load step.
        std::vector<double> recent;
        const int ls_window = std::max(1, options.line_search_window);
        bool step_converged = false;
        // How this increment ended, if it did not converge. Set at each exit from the
        // iteration loop so the abandonment is named where it happens rather than
        // guessed afterwards from the iteration count.
        NewtonResult::Abandonment ended = NewtonResult::Abandonment::IterationBudget;
        int stall = 0;  // consecutive iterations that produced no descent
        for (int iter = 0; iter < options.max_iterations; ++iter) {
            builder.clear();
            probe.reset_counts();
            probe.first_iterate = first_iterate;
            first_iterate = false;
            // The family is gathered on every accepted iterate and NOT on the line search's
            // trial evaluations, which is both correct (only an accepted iterate's local error
            // means anything) and why it is free: measured against a build that skipped it
            // entirely, the same five cases run in 1.1/2.4/4.3/6.8/0.8 s with it and
            // 1.2/2.9/4.6/7.5/0.7 s without -- inside the noise, with identical iteration
            // counts and bit-identical answers.
            const Eigen::VectorXd f_int = assemble(du_free, true, &builder, &probe);
            const Eigen::VectorXd residual = target - f_int;
            const double rnorm = residual.norm();
            measure_convergence(probe, f_int, residual, rnorm, cf_norm, ref, rtol,
                                result.convergence);
            recent.push_back(rnorm);
            if ((int)recent.size() > ls_window) recent.erase(recent.begin());
            // What the trial step has to beat. With a window of 1 this is rnorm and the test
            // below is the strict Armijo condition, unchanged.
            const double gate = *std::max_element(recent.begin(), recent.end());
            if (debug) {
                const auto& c = result.convergence;
                std::fprintf(stderr,
                             "  lambda %.4f iter %d  rnorm=%.4e  rel=%.4e  | CSP=%.5f "
                             "force=%.3e  plastic %d/%d worst=%.3e  nl-el %d/%d worst=%.3e\n",
                             target_lambda, iter, rnorm, rnorm / ref, c.csp, c.force_error,
                             c.plastic_inaccurate, c.plastic_points, c.worst_plastic_error,
                             c.nl_elastic_inaccurate, c.nl_elastic_points,
                             c.worst_nl_elastic_error);
            }
            if (rnorm <= rtol * ref && (!enforce_local || result.convergence.local_ok())) {
                step_converged = true;
                ++result.total_iterations;
                break;
            }

            const auto t_csr = Clock::now();
            const math::CsrMatrix& kt = builder.build_cached(kt_cache);
            result.timings.csr_build += elapsed(t_csr);
            const auto t_lin = Clock::now();
            Eigen::VectorXd delta;
            bool solved = true;
            try {
                delta = linear_solve(kt, residual);
            } catch (const math::SingularSystem& refusal) {
                // The tangent is singular at this iterate and the solver refused to
                // return a vector that does not satisfy it. That is the normal state
                // at a limit load: once a collapse mechanism forms -- a pile whose
                // skin friction has fully plastified, a soil body at its bearing
                // capacity -- the stiffness is rank-deficient along the mechanism and
                // the out-of-balance force has a component along it, so no increment
                // can restore equilibrium. Measured on test_pile_capacity: rank 74 of
                // 75, the null vector a uniform axial translation of the pile, i.e.
                // the plunging mode itself.
                //
                // So this is a property of THIS increment, not a fatal error. Abandon
                // it and let the outer loop halve the increment, exactly as it does
                // for a persistently non-descending direction; cutting back to the
                // minimum increment is what turns this into the reported collapse
                // load. Before the solver verified its answers this path still worked,
                // but only by accident: the solver returned a finite meaningless
                // vector and the line search happened to reject it.
                //
                // Only SingularSystem is caught. A malformed request or a broken
                // backend raises SolveError instead and propagates, because turning
                // one of those into a cut-back would publish a collapse load that is
                // really a bug.
                solved = false;
                ++result.refused_solves;
                ended = NewtonResult::Abandonment::SolveRefused;
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f iter %d  SOLVE REFUSED: %s\n",
                                 target_lambda, iter, refusal.what());
            }
            result.timings.linear_solve += elapsed(t_lin);
            ++result.timings.n_solve;
            if (!solved) break;   // abandon the increment -> outer loop cuts back

            // Backtracking line search: the consistent tangent gives fast local
            // convergence but a full Newton step can overshoot far from the
            // solution (the perfectly-plastic tangent is singular/indefinite at
            // yielding points). Shrink alpha until the residual decreases
            // (Armijo). This globalizes Newton without sacrificing its rate.
            double alpha = 1.0;
            bool improved = false;
            for (int ls = 0; ls < 12; ++ls) {
                const Eigen::VectorXd fi = assemble(du_free + alpha * delta,
                                                    false, nullptr);
                if ((target - fi).norm() < (1.0 - 1.0e-4 * alpha) * gate) {
                    improved = true;
                    break;
                }
                alpha *= 0.5;
            }
            ++result.total_iterations;
            if (debug)
                std::fprintf(stderr, "      alpha=%.6f %s |delta|=%.4e |du|=%.4e\n", alpha,
                             improved ? "descent" : "NO-DESCENT",
                             delta.lpNorm<Eigen::Infinity>(),
                             du_free.lpNorm<Eigen::Infinity>());
            // Persistent no-descent means the Newton direction is useless for this
            // increment (too large -> singular/indefinite tangent); abandon it so
            // the outer loop cuts the increment rather than grinding to
            // max_iterations. But a SINGLE no-descent can be transient (e.g. at the
            // onset of contained plastic flow), so only abort after several
            // consecutive failures -- this keeps cutback cheap without spuriously
            // collapsing recoverable steps.
            stall = improved ? 0 : stall + 1;
            if (stall >= 4) {
                ended = NewtonResult::Abandonment::NoDescent;
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f iter %d  ABANDONED: four consecutive "
                                         "iterations without descent\n", target_lambda, iter);
                break;
            }
            du_free += alpha * delta;
        }

        if (step_converged) {
            hs_consistent_mode = false;  // increment closed → back to the cheap continuum path
            u_free += du_free;
            u_struct += du_free;
            committed = trial;  // commit the increment
            anchor_committed = anchor_trial;    // commit the structural elastoplastic state too
            geogrid_committed = geogrid_trial;
            interface_committed = interface_trial;
            interface5_committed = interface5_trial;
            eskin_committed = eskin_trial;
            efoot_committed = efoot_trial;
            plate_committed = plate_trial;
            plate5_committed = plate5_trial;
            lambda = target_lambda;
            result.load_factor = lambda;
            dlam = std::min(init_dlam, dlam * 1.5);  // allow recovery
        } else if (has_hs && !hs_consistent_mode) {
            // Hybrid tangent: continuum failed to converge this increment → retry the SAME
            // increment with the consistent (FD) tangent without touching dlam; the
            // remaining increments stay consistent too.
            hs_consistent_mode = true;
        } else {
            // The increment is abandoned: record WHY before the size is halved, so the
            // phase can report the reason it actually stopped for.
            if (ended == NewtonResult::Abandonment::NoDescent) ++result.no_descent;
            else if (ended == NewtonResult::Abandonment::IterationBudget) ++result.budget_exhausted;
            result.last_abandonment = ended;
            dlam *= 0.5;  // sub-stepping: shrink the increment and retry
            if (debug)
                std::fprintf(stderr, "  cutback at lambda %.4f: dlam %.6f (minimum %.6f), "
                                     "%d refused solve(s) so far\n",
                             lambda, dlam, min_dlam, result.refused_solves);
            if (dlam < min_dlam) {
                result.converged = false;
                break;  // true collapse: no equilibrium even at the smallest increment
            }
        }
    }

    result.displacement = expand_to_full(dofs, u_free);
    if (has_presc)  // prescribed displacement ū at fixed DOFs (full value; at lambda=1)
        for (int d = 0; d < dofs.total_dofs(); ++d)
            if (dofs.is_fixed(d)) result.displacement[d] += result.load_factor * presc(d);
    result.timings.total = elapsed(t_start);
    return result;
}

// STRUCTURAL f_s0 — with the shared assembler, using the ZERO committed Gauss + du=0 trick
// (the contract in the header): the soil loop gives Δε=0 → σ_trial = committed = 0 →
// identically zero contribution; the remaining pure structural force is ARITHMETICALLY
// IDENTICAL to what solve_nonlinear would produce with the same datum/seed (same code path).
template <class E, class Kin>
Eigen::VectorXd structural_f0_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                                   const std::vector<MaterialModel>& materials,
                                   const Structures& structures,
                                   const StructuralInit& init_struct) {
    const int neq = dofs.equation_count();
    const bool has_datum = init_struct.u_datum.size() > 0;
    if (has_datum && init_struct.u_datum.size() != neq)
        throw std::invalid_argument("structural_internal_force: u_datum size is not equation_count");
    std::vector<GaussState> committed(static_cast<size_t>(mesh.element_count) * E::kGaussCount);
    std::vector<GaussState> trial = committed;
    std::vector<double> anchor_c, geogrid_c, iface_c, iface5_c, eskin_c, efoot_c, plate_c, plate5_c;
    seed_structural_state(structures, init_struct, anchor_c, geogrid_c, iface_c, iface5_c,
                          eskin_c, efoot_c, plate_c, plate5_c);
    std::vector<double> anchor_t = anchor_c, geogrid_t = geogrid_c, iface_t = iface_c,
                        iface5_t = iface5_c, eskin_t = eskin_c, efoot_t = efoot_c,
                        plate_t = plate_c, plate5_t = plate5_c;
    const std::vector<char> no_mask;
    const std::vector<MaterialProfile> no_profile;
    detail::InternalForceAssembler<E, Kin> fasm(mesh, dofs, materials, no_mask, structures,
                                                no_profile);
    typename detail::InternalForceAssembler<E, Kin>::State st;
    st.committed = &committed; st.trial = &trial;
    st.anchor_c = &anchor_c;   st.anchor_t = &anchor_t;
    st.geogrid_c = &geogrid_c; st.geogrid_t = &geogrid_t;
    st.iface_c = &iface_c;     st.iface_t = &iface_t;
    st.iface5_c = &iface5_c;   st.iface5_t = &iface5_t;
    st.eskin_c = &eskin_c;     st.eskin_t = &eskin_t;
    st.efoot_c = &efoot_c;     st.efoot_t = &efoot_t;
    st.plate_c = &plate_c;     st.plate_t = &plate_t;
    st.plate5_c = &plate5_c;   st.plate5_t = &plate5_t;
    const typename detail::InternalForceAssembler<E, Kin>::Ramp no_ramp;
    const Eigen::VectorXd u = has_datum ? Eigen::VectorXd(init_struct.u_datum)
                                        : Eigen::VectorXd::Zero(neq);
    return fasm.assemble(u, Eigen::VectorXd::Zero(neq), /*build_tangent=*/false,
                         TangentMode::kNone, st, no_ramp, nullptr, nullptr);
}

} // namespace

bool seed_structural_state(const Structures& structures, const StructuralInit& init_struct,
                           std::vector<double>& anchor, std::vector<double>& geogrid,
                           std::vector<double>& iface_v, std::vector<double>& iface5_v,
                           std::vector<double>& eskin, std::vector<double>& efoot,
                           std::vector<double>& plate_p, std::vector<double>& plate5_p) {
    bool carried = false;
    auto seed = [&carried](std::vector<double>& dst, const std::vector<double>& src, size_t n,
                           const char* what) {
        if (src.empty()) { dst.assign(n, 0.0); return; }
        if (src.size() != n)
            throw std::invalid_argument(std::string("StructuralInit.") + what +
                                        " boyutu yapisal elemanlarla uyusmuyor (" +
                                        std::to_string(src.size()) + " != " + std::to_string(n) + ")");
        dst = src;
        carried = true;
    };
    size_t total_skin = 0;
    for (const auto& eb : structures.embedded_beams) total_skin += eb.skin.size();
    seed(anchor, init_struct.anchor_plastic, structures.anchors.size(), "anchor_plastic");
    seed(geogrid, init_struct.geogrid_plastic,
         structures.geogrids.size() * static_cast<size_t>(geogrid::kGaussCount), "geogrid_plastic");
    seed(iface_v, init_struct.interface_slip,
         structures.interfaces.size() * static_cast<size_t>(iface::kPointCount), "interface_slip");
    seed(iface5_v, init_struct.interface5_slip,
         structures.interfaces5.size() * static_cast<size_t>(iface::kPointCount5), "interface5_slip");
    seed(eskin, init_struct.embedded_skin_slip, total_skin, "embedded_skin_slip");
    seed(efoot, init_struct.embedded_foot_slip, structures.embedded_beams.size(),
         "embedded_foot_slip");
    seed(plate_p, init_struct.plate_plastic,
         structures.plates.size() * static_cast<size_t>(plate::kPlasticStateSize), "plate_plastic");
    seed(plate5_p, init_struct.plate5_plastic,
         structures.plates5.size() * static_cast<size_t>(plate::kPlasticStateSize5),
         "plate5_plastic");
    return carried;
}

Eigen::VectorXd structural_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                                          const std::vector<MaterialModel>& materials,
                                          const Structures& structures,
                                          const StructuralInit& init_struct,
                                          Kinematics kinematics) {
    const bool axi = kinematics == Kinematics::Axisymmetric;
    if (mesh.nodes_per_element == Tri15Element::kNodeCount) {
        if (axi)
            return structural_f0_impl<Tri15Element, detail::AxisymKin>(mesh, dofs, materials,
                                                                       structures, init_struct);
        return structural_f0_impl<Tri15Element, detail::PlaneStrainKin>(mesh, dofs, materials,
                                                                        structures, init_struct);
    }
    if (axi)
        return structural_f0_impl<Tri6Element, detail::AxisymKin>(mesh, dofs, materials,
                                                                  structures, init_struct);
    return structural_f0_impl<Tri6Element, detail::PlaneStrainKin>(mesh, dofs, materials,
                                                                   structures, init_struct);
}

NewtonResult solve_nonlinear(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<MaterialModel>& materials,
                             const Eigen::VectorXd& f_ext,
                             const LinearSolve& linear_solve,
                             const NewtonOptions& options,
                             const std::vector<GaussState>& initial_state,
                             const std::vector<char>& active_element,
                             const Structures& structures,
                             const Eigen::VectorXd& prescribed_displacement,
                             const Eigen::VectorXd& constant_force,
                             const std::vector<MaterialProfile>& profile,
                             const StructuralInit& init_struct) {
    // Element type (mesh) x kinematics (options) -> a monomorphic compile-time inner loop.
    const bool axi = options.kinematics == Kinematics::Axisymmetric;
    const Eigen::VectorXd& pd = prescribed_displacement;
    const Eigen::VectorXd& cf = constant_force;
    const std::vector<MaterialProfile>& pf = profile;
    if (mesh.nodes_per_element == Tri15Element::kNodeCount) {
        if (axi)
            return solve_nonlinear_impl<Tri15Element, detail::AxisymKin>(
                mesh, dofs, materials, f_ext, linear_solve, options, initial_state,
                active_element, structures, pd, cf, pf, init_struct);
        return solve_nonlinear_impl<Tri15Element, detail::PlaneStrainKin>(
            mesh, dofs, materials, f_ext, linear_solve, options, initial_state,
            active_element, structures, pd, cf, pf, init_struct);
    }
    if (axi)
        return solve_nonlinear_impl<Tri6Element, detail::AxisymKin>(
            mesh, dofs, materials, f_ext, linear_solve, options, initial_state,
            active_element, structures, pd, cf, pf, init_struct);
    return solve_nonlinear_impl<Tri6Element, detail::PlaneStrainKin>(
        mesh, dofs, materials, f_ext, linear_solve, options, initial_state,
        active_element, structures, pd, cf, pf, init_struct);
}

} // namespace katai::core

