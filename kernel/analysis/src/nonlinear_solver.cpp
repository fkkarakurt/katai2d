#include <katai/analysis/nonlinear_solver.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

#include <katai/analysis/excess_pore_force.hpp>
#include <katai/analysis/internal_forces.hpp>  // SHARED internal-force/tangent assembly (static+dynamic)
#include <katai/analysis/post/consolidation_recovery.hpp>  // store_excess_pore_pressure
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

// Turn one iterate's raw measurement into the convergence family (NewtonResult::Convergence).
// Called at every ACCEPTED iterate, so what it leaves behind is what the last such iterate
// satisfied -- including, when the solve fails, the last one that was reached.
void measure_convergence(const detail::LocalErrorProbe& p, const Eigen::VectorXd& f_int,
                         const Eigen::VectorXd& residual, double rnorm, double cf_norm,
                         double ref, double rtol, NewtonResult::Convergence& cv) {
    cv.measured = true;
    cv.tolerated = rtol;

    // CSP, the current stiffness parameter: the work of the increment over the work the same
    // strain would have done elastically. Guarded twice, for two different reasons: an increment
    // that has not moved yet has no energy to take a ratio of (report the elastic value, 1), and
    // an UNLOADING increment can produce a negative or greater-than-one ratio, which is not a
    // stiffness measure -- the parameter is defined on the interval and is clamped to it.
    cv.csp = p.energy_elastic > 0.0
                 ? std::min(1.0, std::max(0.0, p.energy_total / p.energy_elastic))
                 : 1.0;

    // The stiffness-weighted global force error, ||r|| / (||f_int|| + CSP*||f_const||). The
    // denominator is the force scale the model itself is carrying, plus the standing load from
    // the phases before this one, weighted by how much stiffness is left.
    // The floor is 1e-6 of the phase's own load scale: it bites only when nothing is resisting
    // yet, and then a huge relative error is the honest reading.
    const double denom = std::max(f_int.norm() + cv.csp * cf_norm, 1e-6 * ref);
    cv.force_error = rnorm / denom;
    // The default gate's own ratio, so that both are on the record at the same iterate.
    cv.global_error = rnorm / ref;

    // The moment error: the largest out-of-balance moment on a rotational equation over the sum
    // of absolute nodal moment contributions. Only structural elements with rotational freedom
    // have one.
    //
    // THE FLOOR IS AN ADDITION TO THAT RATIO, and it is here deliberately. Every other criterion
    // in this family has one -- the local error at a plastic stress point divides by
    // max(tau_max, c, 1 kPa), the foot force error by max(|F_c|, 1% of |F_max|, 1 kN) -- and the
    // reason for them is that a point carrying almost nothing must not report an enormous
    // relative error on a difference that is numerically nothing. A structure carrying almost no
    // moment is the same situation in different units, and the tree had no floor there because
    // the criterion was never consulted: the moment reference is a sum of ABSOLUTE nodal moments,
    // which cannot cancel, but it CAN be empty.
    //
    // Measured, on the day the criterion first gated (2026-08-25): a plate standing along the
    // whole of a line that is pushed down is driven UNIFORMLY -- every one of its nodes takes the
    // same settlement, so it translates without curving -- and its reference is 1.11e-12 kNm/m,
    // i.e. round-off. (This comment first called the plate "undriven", after K2D-A003; that
    // warning was stale and is retired -- the plate moves, it just does not bend.) The
    // ratio of one round-off to another came out at 2.06e-1 and refused the phase outright, at
    // load factor 0. The same plate under a strip load has a reference of 4.71e+02 kNm/m and an
    // error of 4.5e-14. Fourteen orders of magnitude of denominator separate a criterion that
    // measures something from one that measures nothing, and 1 kNm/m sits between them by ten.
    cv.moment_ref = p.m_ref;
    // has_moment is now about the MODEL, not about the loading: if something in it carries a
    // rotational freedom the criterion applies, and an unloaded structure reports a small error
    // rather than disappearing from the report. (Before the floor it had to mean "and the
    // reference is non-zero", which quietly excused exactly the case that needed saying.)
    cv.has_moment = !p.moment_eq.empty();
    if (cv.has_moment) {
        double worst = 0.0;
        for (int eq : p.moment_eq) worst = std::max(worst, std::fabs(residual(eq)));
        cv.moment_error = worst / std::max(p.m_ref, 1.0);
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

    // The embedded-beam foot force error, sum |F_eq - F_c| over the denominator below, floor and
    // all: the constitutive foot forces, or one per cent of the declared capacities, or 1.0 --
    // whichever is largest. The middle term is what stops a pile that carries almost nothing at
    // its toe from reporting a huge relative error while the capacity it was given says the force
    // is negligible either way.
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
    // A prescribed displacement that actually moves something: the predictor below is for it.
    const bool presc_active = has_presc && presc.lpNorm<Eigen::Infinity>() > 0.0;
    double cur_lambda = 0.0, cur_target = 0.0;

    // HYBRID HS TANGENT POLICY: increments start with the fast CONTINUUM tangent; if an
    // increment fails to converge, before halving dlam the SAME increment is retried once
    // with the CONSISTENT (numerical, FD) tangent, and the rest of the solve stays
    // consistent (the hard regime has been entered — footing edge/low confining pressure;
    // measured: continuum stalls there at lf=0.84, consistent carries to 1.0, but is ~2.5×
    // more expensive in the easy regime). Active only when an HS material is present →
    // the LE/MC path is the old behaviour EXACTLY (TangentMode does not affect them, the
    // retry branch never opens).
    //
    // WHICH MODELS QUALIFY is the question this flag actually asks: not "is Hardening Soil
    // here" but "is there a material whose integrator can be asked for a STRONGER tangent than
    // the one it just gave". It was written when Hardening Soil was the only such model and it
    // kept that name and that test afterwards, so a model added later got the retry only if
    // somebody remembered this line. Hoek-Brown did not: its integrator grew a finite-difference
    // consistent tangent and the retry never opened, so the tangent existed and was never
    // reached -- measured on a tunnel unloaded into a rock mass, where the phase stalled at
    // exactly the same load factor with the tangent as without it, and at 40 steps as at 150.
    //
    // NOT ENROLLED, and deliberately: Soft Soil and Soft Soil Creep also build a consistent
    // tangent by finite difference. Enrolling them is not a free improvement -- the retry
    // LATCHES for the rest of the phase, so it changes the iteration path of runs that would
    // have recovered anyway, and at a fixed tolerance a changed path can move a converged
    // answer. There is no failing soft-soil case in front of this line to measure that against,
    // and a speculative change to two verified models is not worth making blind.
    bool has_fd_tangent = false;
    for (const auto& mm : materials)
        if (mm.type == MaterialType::HardeningSoil || mm.type == MaterialType::HoekBrown) {
            has_fd_tangent = true;
            break;
        }
    bool hs_consistent_mode = false;
    // Set once the linear solver has refused a tangent: from then on the phase assembles with the
    // vertex floor (InternalForceAssembler::vertex_floor). A refusal used to go straight to cutting
    // the increment back, and three refusals in a row were reported as a collapse load. Measured on
    // a rigid wall translated away from a cohesionless backfill (c = 0, phi = psi = 41, tri6): the
    // phase stopped at 1% of a 4 mm movement and called it a collapse; with the floor it carries
    // the whole movement to 3.2% of Rankine's active thrust (0.16% on 15-noded elements, where the
    // old code happened to pass). Nothing had failed: points at the stress-free apex carry a zero
    // tangent, and a node surrounded by them has no stiffness at all. The answer does not depend
    // on the floor: 141.29 kN/m at 1e-3 and at 1e-1. A mechanism is not removed by it -- the floor adds
    // stiffness only where the exact tangent is identically zero, and the residual is unchanged --
    // so a genuine limit load still ends in a refusal or a stall, now on its own merits.
    bool vertex_floor_mode = false;
    // THE INITIAL-STIFFNESS FALLBACK. An increment whose Newton direction stopped descending --
    // every retry above spent -- is re-entered once with the ELASTIC tangent at its plastic
    // points, full steps and no demand that the residual fall at every step. That iteration is a
    // contraction, not a descent method: it converges where Newton's direction no longer descends,
    // which with non-associated flow it can fail to do on a non-smooth surface, but its residual
    // may rise for an iterate on the way. A monotone line search rejected exactly those iterates.
    // Measured on initial phases on sloping ground (Mohr-Coulomb and Soft Soil, psi = 0, 15-noded
    // elements): "did not converge" at 30-70% of the load with Newton, on any tolerance. It is
    // only ever entered for an increment Newton has already given up on, so every run Newton
    // closes is untouched; it latches off again when the increment closes.
    bool elastic_mode = false;
    // THE SMALL-STRAIN STIFFNESS OF AN INCREMENT, READ FROM A FIXED INCREMENT. In Hardening Soil
    // with small-strain stiffness the stiffness of a step is set by the direction of the step's
    // own strain increment against the strain history: a point that turns back answers with G0,
    // one that carries on with what is left of its curve. Read from the iterate, that makes the
    // internal force a function with slopes of order 1 / |d eps| at every point whose increment
    // is small beside its history, and with steps where the increment's principal directions
    // are not unique -- measured on a tunnel lining installed after grouting beside a piled
    // building (a phase that reverses the ground everywhere): Newton with the exact derivative,
    // and the elastic-operator pass after it, both wandered at 1 to 2% of the load scale with
    // corrections of micrometres, the same sixty stress points changing sides, and the phase
    // was given up at 1% of its load.
    // An increment that cannot be closed that way is solved with the stiffness read from a fixed
    // displacement increment -- the best iterate it reached -- which is a smooth problem that
    // Newton closes; the fixed increment is then replaced by the answer and the increment
    // solved again from it, kSmallRefresh times, so that the stiffness the answer stands on is
    // that of an increment one such pass away from itself. The phase keeps the mode, starting
    // each later increment from the last one's displacement scaled to its size.
    // A run that closes without it never enters it and is untouched.
    bool has_small = false;
    if (!options.srf_reduce)   // a Safety search runs on Mohr-Coulomb copies of the materials
        for (const auto& mm : materials)
            if (mm.type == MaterialType::HardeningSoil && mm.hs.G0_ref > 0.0) has_small = true;
    constexpr int kSmallRefresh = 2;
    bool small_freeze = false;
    int small_refresh = 0;          // refresh passes taken on the increment at hand
    bool small_warm = false;        // the next attempt starts from du_small
    Eigen::VectorXd du_small;
    // What the convergence family read on the Newton iterate the fallback took over from. The
    // fallback's own iterates are elastic-stiffness iterates, not the state the increment failed
    // in: when it fails too, the run reports the Newton one (its CSP is what names a collapse).
    auto newton_convergence = result.convergence;
    bool newton_convergence_kept = false;
    // THE STALL FLOOR: the same device as the vertex floor, raised while the iteration is creeping.
    // Near a limit state many stress points sit on the corners of their surfaces or toggle
    // between elastic and plastic from one iterate to the next; the tangent is exact on the side
    // each point is on, the full Newton step lands on the other side, and the line search can
    // only take a hair of it -- measured on a 30-degree cohesionless slope on 15-noded elements,
    // 200 iterations at steps of 1e-3 before the increment was accepted as stagnated, the global
    // residual long under its tolerance and the local criteria never. After three such steps in a
    // row the plastic points' tangents get a share of their elastic stiffness -- 0.1, doubling
    // while it still creeps, up to the elastic operator itself -- which turns the direction
    // towards the elastic-stiffness iteration that crosses those kinks without looking at them;
    // each good step (alpha >= 1/2) halves it again. Lifting it at once on the first good step put
    // the creep straight back, and on an axisymmetric footing the cycle wandered into a state
    // with no descent at all. Like the vertex floor it changes the path the
    // iteration takes, never the residual it has to reach.
    double stall_floor = 0.0;
    // Once an increment has been abandoned for a stalled line search, the phase judges every
    // attempt from then on against the worst of the last few residuals instead of the latest.
    bool ls_open = false;

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
    // Strength continuation (NewtonOptions::srf_reduce): the assembler reads WORKING copies whose
    // strength is rewritten for each factor the path visits; everything else reads the originals.
    const bool srf_mode = static_cast<bool>(options.srf_reduce);
    std::vector<MaterialModel> srf_materials;
    Structures srf_structures;
    std::vector<MaterialProfile> srf_profile;
    if (srf_mode) {
        srf_materials = materials;
        srf_structures = structures;
        srf_profile = profile;
    }
    detail::InternalForceAssembler<E, Kin> fasm(mesh, dofs, srf_mode ? srf_materials : materials,
                                                active_element,
                                                srf_mode ? srf_structures : structures,
                                                srf_mode ? srf_profile : profile);
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
        fasm.vertex_floor = std::max(vertex_floor_mode ? kVertexFloor : 0.0, stall_floor);
        fasm.elastic_tangent = elastic_mode;
        fasm.du_small = (small_freeze && du_small.size() == du_free.size()) ? &du_small : nullptr;
        return fasm.assemble(u_struct, du_free, build_tangent, tmode, astate, ramp, builder,
                             &result.timings);
    };

    const double cf_norm = has_cf ? constant_force.norm() : 0.0;
    const double ref = std::max({f_ext_norm, cf_norm, 1.0});
    const bool debug = std::getenv("KATAI_NL_DEBUG") != nullptr;
    // An increment must satisfy the LOCAL criteria as well as the global force residual before
    // it is called converged -- local accuracy is part of convergence, not a report on it. ON by
    // default since 2026-08-24; see NewtonOptions::enforce_local_criteria for the measurement
    // that decided it. The two environment variables are the whole-run overrides a study uses, and
    // the option is what a test and the phase's numerical controls drive.
    const bool enforce_local =
        std::getenv("KATAI_CONV_LOCAL") != nullptr
            ? true
            : (std::getenv("KATAI_CONV_NOLOCAL") != nullptr ? false
                                                            : options.enforce_local_criteria);
    // MEASUREMENT SEAM, not a setting: which GLOBAL criterion gates the step. The default is the
    // one this tree has always used -- ||r|| against a FIXED scale, max(||f_ext||, ||f_const||, 1)
    // -- and KATAI_CONV_CSPGATE swaps in the stiffness-weighted force error,
    // ||r|| / (||f_int|| + CSP*||f_const||), which is the form the criteria family already REPORTS.
    //
    // WHAT THE SEAM MEASURED (2026-08-25), because a seam with no reading is just an option. On
    // the footing walked to collapse the two ratios agree within ~1.5x wherever the run converges,
    // and the stiffness-weighted one is 1.4x (q=600) to 2.2x (q=900) STRICTER at collapse -- the
    // property the CSP normalisation exists for. But swapping the gate outright FAILS 5 of the
    // 142 fast tests, and the clearest of them equilibrates 0% of its load: a phase that
    // prestresses an anchor against ground that has not responded has almost no internal force,
    // so the weighted denominator is small and the ratio never reaches the tolerance. Its only
    // floor is 1e-6 of the phase's load scale -- six decades down. That is the SAME defect the
    // moment criterion had (see measure_convergence): a ratio in this family whose denominator
    // has no floor tied to what is being applied. The stiffness-weighted error cannot gate until
    // that floor is designed, so this stays a seam.
    const bool csp_gate = std::getenv("KATAI_CONV_CSPGATE") != nullptr;

    // Adaptive (automatic) load incrementation. The external load is advanced from
    // 0 to f_ext by increments d_lambda; an increment that fails to converge is
    // retried with a halved size (sub-stepping) instead of aborting the analysis.
    // The collapse signal -- non-convergence even at the minimum increment -- is
    // thus decoupled from the (arbitrary) initial step count, making the limit/FoS
    // result robust across element orders and problem stiffness (cf. Crisfield,
    // automatic load stepping). When every increment converges at the
    // initial size (no cutback) the scheme reduces exactly to fixed N-step loading,
    // so the tri6 benchmarks are unchanged.
    result.iteration_limit = options.max_iterations;
    const double init_dlam = 1.0 / options.load_steps;
    const double min_dlam = init_dlam * options.min_step_fraction;  // smaller -> collapse
    double lambda = 0.0;        // committed load fraction
    double dlam = init_dlam;    // current increment size
    result.converged = true;

    // The global tangent assembly buffers live ACROSS ITERATIONS: since the COO entry order
    // is deterministic (the assembly loops run over a fixed topology), the CSR pattern is
    // built once and later iterations only accumulate values (build_cached; falls back to a
    // full build by itself if the signature mismatches).
    math::SparseMatrixBuilder builder(neq);
    math::CsrPatternCache kt_cache;

    // The convergence family's memory (LocalErrorProbe). prev_sc carries sigma_c,j-1
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
    // THE FULL STEP'S ASSEMBLY, KEPT. The line search's first trial is the full Newton step, and
    // it is taken: 343 of the 371 iterations of a published staged excavation (Hardening Soil).
    // Its residual was assembled without the tangent, and the next iteration then assembled the
    // same displacement again, tangent and all -- every constitutive integration done twice. The
    // first trial now assembles the tangent too, into its own builder and with its own copy of
    // the probe's memory (the probe overwrites the previous iterate's stress with the current
    // one, which a trial that is turned down must not do), and an iteration that starts where an
    // accepted full step ended takes it as its own assembly. The same displacement under the same
    // controls is the same assembly, so the iterates are bit for bit those of before.
    math::SparseMatrixBuilder builder_t(neq);
    detail::LocalErrorProbe probe_t;
    std::vector<GaussState> prev_sc_t;
    std::vector<double> prev_deps_t;
    std::vector<detail::PrevPoint> prev_iface_t, prev_iface5_t, prev_skin_t, prev_foot_t;
    const auto point_probe = [&](detail::LocalErrorProbe& pb, bool trial) {
        pb.prev_sigma_c = trial ? &prev_sc_t : &prev_sc;
        pb.prev_deps = trial ? &prev_deps_t : &prev_deps;
        pb.prev_iface = trial ? &prev_iface_t : &prev_iface;
        pb.prev_iface5 = trial ? &prev_iface5_t : &prev_iface5;
        pb.prev_skin = trial ? &prev_skin_t : &prev_skin;
        pb.prev_foot = trial ? &prev_foot_t : &prev_foot;
    };

    // ===================================================================================
    // STRENGTH CONTINUATION BY ARC LENGTH (NewtonOptions::srf_reduce).
    //
    // A strength-reduction search that PRESCRIBES the factor s at every step solves, near the
    // factor of safety, for an equilibrium that exists but lies far away along a developing
    // mechanism: measured on a c-phi slope, the last 0.5% of s below the factor carried the
    // displacement from 9.5 to 78 mm, and a Newton iteration asked for it at a prescribed s
    // crawled for a thousand iterations or was given up, so the factor came out short of the
    // path's own limit. Here s is an UNKNOWN. Each step prescribes the length of the
    // displacement increment along the predicted direction instead (Riks' normal plane,
    // n . du = dl) and solves the bordered system for (du, s):
    //     r(du, s) = f_held - f_int(du; s) = 0,     K a = r,  K b = df_int/ds,
    //     du <- du + a - ds b,   ds = (n . (du + a) - dl) / (n . b),
    // df_int/ds by a one-sided difference of the internal force (one extra residual). Near the
    // limit the path flattens -- s barely rises while the displacement grows -- and each step
    // stays a well-posed small one. The factor of safety is where the path stops rising: s
    // falls (a limit point passed), or it has risen less than the resolution over each of the
    // last two doublings of the displacement (the plateau of a perfectly plastic mechanism,
    // where the rise over a doubling is what is left to rise).
    // Steps are accepted on exactly the criteria of the load-controlled solve (global residual
    // and the local ones), and the step length follows the iteration count of the last step.
    // ===================================================================================
    if (srf_mode) {
        const auto set_srf = [&](double sf) {
            srf_materials = materials;
            srf_structures = structures;
            srf_profile = profile;
            options.srf_reduce(srf_materials, srf_structures, srf_profile, sf);
        };
        const Eigen::VectorXd held =
            has_cf ? Eigen::VectorXd(constant_force) : Eigen::VectorXd::Zero(neq);
        const int max_it = std::max(options.max_iterations, 20);
        double s_cur = options.srf_start;
        double dl = -1.0;                 // the arc length (displacement increment norm)
        int flat = 0;
        int cut = 0;
        double s_ctrl_step = options.srf_first_step;   // for steps with no plastic response
        result.srf = s_cur;
        result.srf_path.push_back({s_cur, 0.0});
        set_srf(s_cur);
        // The last accepted arc step (displacement and factor), for the secant predictor below.
        Eigen::VectorXd sec_du;
        double sec_ds = 0.0;
        bool sec_ok = false;
        for (int stepn = 0; stepn < 2000; ++stepn) {
            if (s_cur >= options.srf_max) break;
            cur_lambda = 0.0;
            cur_target = 1.0;
            Eigen::VectorXd du = Eigen::VectorXd::Zero(neq);
            prev_sc = committed;
            std::fill(prev_deps.begin(), prev_deps.end(), 0.0);
            Eigen::VectorXd n;
            double s_try;
            bool arc = false;
            // ---- THE SECANT PREDICTOR. Once the path is on its arc, the step before is the best
            // account there is of where it goes: the same points go on yielding and the factor
            // rises by little for a lot of displacement. The tangent at the committed state does
            // not know that -- a point on its yield surface with no strain increment yet reads as
            // elastic, the predicted displacement per unit of factor is several times too small,
            // and the arc length is reached by a factor far beyond the limit. Measured on an
            // anchored excavation in Hardening Soil near its factor of 1.62: a predicted factor of
            // 1.86 for a step that converged at 1.6193, the first six of its 21 iterations spent
            // walking the factor back. The secant along the last accepted step at the current
            // arc length is the predictor; the tangent one starts the path, follows a step
            // controlled by the factor and retries a step that failed. With it the factor of that
            // excavation (1.630 either way) took 42 steps and 20 s instead of 189 and 119 s, and
            // the factor of a slope on Hardening Soil 51 s instead of 122 (1.832 against 1.831,
            // within its resolution).
            if (sec_ok) {
                const double scale = dl / sec_du.norm();
                du = scale * sec_du;
                n = du / du.norm();
                s_try = s_cur + scale * sec_ds;
                arc = true;
            } else {
                // ---- predictor at the committed state
                const double eps0 = 1.0e-6 * s_cur;
                set_srf(s_cur + eps0);
                const Eigen::VectorXd fe0 = assemble(du, false, nullptr);
                set_srf(s_cur);
                builder.clear();
                const Eigen::VectorXd f00 = assemble(du, true, &builder);
                const Eigen::VectorXd fs0 = (fe0 - f00) / eps0;
                Eigen::VectorXd v;
                bool solved = true;
                const auto t_csr0 = Clock::now();
                const math::CsrMatrix& k0 = builder.build_cached(kt_cache);
                result.timings.csr_build += elapsed(t_csr0);
                const auto t_lin0 = Clock::now();
                try {
                    v = linear_solve(k0, fs0);
                    ++result.timings.n_solve;
                } catch (const math::SingularSystem&) {
                    solved = false;
                }
                result.timings.linear_solve += elapsed(t_lin0);
                if (!solved || !v.allFinite()) break;
                const double vnorm = v.norm();
                // With no plastic response the path is vertical in s: such a step is controlled by s.
                arc = vnorm * options.srf_first_step > 1e-12 * std::max(1.0, u_free.norm()) &&
                      vnorm > 0.0;
                if (arc) {
                    if (dl < 0.0) dl = options.srf_first_step * vnorm;
                    const double ds_pred = dl / vnorm;
                    du = -ds_pred * v;
                    n = du / du.norm();
                    s_try = s_cur + ds_pred;
                } else {
                    s_try = s_cur + s_ctrl_step;
                }
            }
            // ---- corrector
            // The iterate's own assembly comes first and the factor's derivative only after the
            // convergence test, which (a) saves the derivative's residual on the iterate that
            // converges and (b) leaves the iterate's assembly the last one made when it does --
            // the trial states committed below are the last assembled. An iterate that starts
            // where an accepted full step ended takes that step's assembly as its own (the
            // load-controlled path below explains the kept assembly): the same displacement at
            // the same factor, the same bits.
            bool conv = false;
            int its = 0;
            bool kept = false;
            Eigen::VectorXd kept_fi;
            for (int it = 0; it < max_it; ++it) {
                its = it + 1;
                const double eps = 1.0e-6 * s_try;
                set_srf(s_try);
                const bool reuse = kept;
                kept = false;
                if (reuse) {
                    std::swap(builder, builder_t);
                    prev_sc.swap(prev_sc_t);
                    prev_deps.swap(prev_deps_t);
                    prev_iface.swap(prev_iface_t);
                    prev_iface5.swap(prev_iface5_t);
                    prev_skin.swap(prev_skin_t);
                    prev_foot.swap(prev_foot_t);
                    probe = probe_t;
                    point_probe(probe, false);
                } else {
                    builder.clear();
                    probe.reset_counts();
                    probe.first_iterate = (it == 0);
                }
                const Eigen::VectorXd fi = reuse ? std::move(kept_fi) : assemble(du, true, &builder, &probe);
                const Eigen::VectorXd r = held - fi;
                const double rnorm = r.norm();
                measure_convergence(probe, fi, r, rnorm, cf_norm, ref, rtol, result.convergence);
                if (debug)
                    std::fprintf(stderr, "  srf %.6f it %d rnorm=%.4e rel=%.4e plastic %d/%d %s\n",
                                 s_try, it, rnorm, rnorm / ref,
                                 result.convergence.plastic_inaccurate,
                                 result.convergence.plastic_points, arc ? "arc" : "s-ctrl");
                const bool global_ok =
                    csp_gate ? result.convergence.force_ok() : (rnorm <= rtol * ref);
                if (it > 0 && global_ok && (!enforce_local || result.convergence.enforced_ok())) {
                    conv = true;
                    break;
                }
                Eigen::VectorXd fe;
                if (arc) {
                    set_srf(s_try + eps);
                    fe = assemble(du, false, nullptr);
                    set_srf(s_try);
                }
                const auto t_csr = Clock::now();
                const math::CsrMatrix& kt = builder.build_cached(kt_cache);
                result.timings.csr_build += elapsed(t_csr);
                Eigen::VectorXd a, b;
                const auto t_lin = Clock::now();
                try {
                    a = linear_solve(kt, r);
                    ++result.timings.n_solve;
                    if (arc) {
                        b = linear_solve(kt, (fe - fi) / eps);
                        ++result.timings.n_solve;
                    }
                } catch (const math::SingularSystem&) {
                    result.timings.linear_solve += elapsed(t_lin);
                    break;
                }
                result.timings.linear_solve += elapsed(t_lin);
                if (!a.allFinite() || (arc && !b.allFinite())) break;
                double dsi = 0.0;
                if (arc) {
                    const double nb = n.dot(b);
                    if (std::fabs(nb) < 1e-300) break;
                    dsi = (n.dot(du) + n.dot(a) - dl) / nb;
                }
                const Eigen::VectorXd ddu = arc ? Eigen::VectorXd(a - dsi * b) : a;
                // A damped step when the full one throws the residual out: the residual at the
                // trial is compared with the current one, halving up to six times.
                double alpha = 1.0;
                for (int ls = 0; ls < 6; ++ls) {
                    set_srf(s_try + alpha * dsi);
                    const bool full = ls == 0;
                    if (full) {
                        builder_t.clear();
                        prev_sc_t = prev_sc;
                        prev_deps_t = prev_deps;
                        prev_iface_t = prev_iface;
                        prev_iface5_t = prev_iface5;
                        prev_skin_t = prev_skin;
                        prev_foot_t = prev_foot;
                        probe_t = probe;
                        point_probe(probe_t, true);
                        probe_t.reset_counts();
                        probe_t.first_iterate = false;
                    }
                    Eigen::VectorXd ft = full ? assemble(du + alpha * ddu, true, &builder_t, &probe_t)
                                              : assemble(du + alpha * ddu, false, nullptr);
                    const double rt = (held - ft).norm();
                    if (std::isfinite(rt) && rt < std::max(rnorm, 10.0 * rtol * ref)) {
                        if (full) { kept = true; kept_fi = std::move(ft); }
                        break;
                    }
                    alpha *= 0.5;
                }
                du += alpha * ddu;
                s_try += alpha * dsi;
            }
            result.total_iterations += its;
            if (conv) {
                // the trial states of the accepted iterate are the last assembled at s_try
                committed = trial;
                anchor_committed = anchor_trial;
                geogrid_committed = geogrid_trial;
                interface_committed = interface_trial;
                interface5_committed = interface5_trial;
                eskin_committed = eskin_trial;
                efoot_committed = efoot_trial;
                plate_committed = plate_trial;
                plate5_committed = plate5_trial;
                u_free += du;
                u_struct += du;
                const double s_prev = s_cur;
                s_cur = s_try;
                sec_ok = arc && du.norm() > 0.0;
                if (sec_ok) { sec_du = du; sec_ds = s_cur - s_prev; }
                result.srf = std::max(result.srf, s_cur);
                result.srf_path.push_back({s_cur, u_free.cwiseAbs().maxCoeff()});
                ++result.timings.n_steps;
                if (debug) {
                    // e = d ln(s - s0) / d ln u over the last doubling of the displacement
                    const double un = result.srf_path.back()[1];
                    double e_dbl = -1.0;
                    for (size_t j = result.srf_path.size(); j-- > 0;) {
                        const double uj = result.srf_path[j][1], sj = result.srf_path[j][0];
                        if (uj <= 0.5 * un && uj > 0.0 && sj > options.srf_start && s_cur > options.srf_start) {
                            e_dbl = std::log((s_cur - options.srf_start) / (sj - options.srf_start)) /
                                    std::log(un / uj);
                            break;
                        }
                    }
                    std::fprintf(stderr, "  SRF STEP %d: s %.6f -> %.6f in %d it, |du| %.3e maxu %.4e e %.5f%s\n",
                                 stepn, s_prev, s_cur, its, du.norm(), un, e_dbl, arc ? "" : " (s-ctrl)");
                }
                if (arc) {
                    if (s_cur < s_prev) { result.srf_limit = true; break; }
                    // THE PLATEAU. On a perfectly plastic mechanism the path approaches its
                    // limit F like s = F - a/u, so the rise over the last doubling of the
                    // displacement equals what is left to rise. When that rise is below the
                    // resolution the factor is known to the resolution, and the path stops.
                    // A path can flatten for a while and rise again, so two successive doublings
                    // must each rise less than the resolution.
                    const auto rise_over_doubling = [&](double u_end, double s_end) {
                        for (size_t j = result.srf_path.size(); j-- > 0;)
                            if (result.srf_path[j][1] <= 0.5 * u_end)
                                return result.srf_path[j][1] > 0.0
                                           ? std::make_pair(true, s_end - result.srf_path[j][0])
                                           : std::make_pair(false, 0.0);
                        return std::make_pair(false, 0.0);
                    };
                    const double u_now = result.srf_path.back()[1];
                    const auto r1 = rise_over_doubling(u_now, s_cur);
                    double s_half = s_cur;
                    for (size_t j = result.srf_path.size(); j-- > 0;)
                        if (result.srf_path[j][1] <= 0.5 * u_now) { s_half = result.srf_path[j][0]; break; }
                    const auto r2 = rise_over_doubling(0.5 * u_now, s_half);
                    flat = (r1.first && r2.first && r1.second < options.srf_resolution * s_cur &&
                            r2.second < options.srf_resolution * s_cur) ? 1 : 0;
                    if (flat) { result.srf_limit = true; break; }
                    // The next step aims at 20 iterations: measured on the slope above, 12
                    // took 117 steps and 37 s for the same factor (to 4e-5) that 20 reached in
                    // 35 steps and 15 s. A controller that grew the step by half whenever it
                    // converged in 25 was faster nowhere and moved the factor by 0.3% between
                    // tolerances.
                    dl *= std::clamp(std::sqrt(20.0 / std::max(its, 1)), 0.5, 2.0);
                } else {
                    s_ctrl_step *= 1.5;
                }
                cut = 0;
            } else {
                trial = committed;
                if (arc) dl *= 0.5; else s_ctrl_step *= 0.5;
                sec_ok = false;   // the retry predicts from the tangent
                set_srf(s_cur);
                if (++cut > 12 || (!arc && s_ctrl_step < options.srf_resolution * s_cur)) {
                    result.srf_limit = true;   // no equilibrium beyond: the limit
                    break;
                }
            }
        }
        set_srf(s_cur);
        lambda = 1.0;                     // the load loop below has nothing to do
        result.load_factor = 1.0;
        result.converged = true;
    }

    // The displacement of the last accepted increment and its size (see THE SECANT START).
    Eigen::VectorXd du_last;
    double dlam_last = 0.0;
    while (lambda < 1.0 - 1e-12) {
        if (dlam > 1.0 - lambda) dlam = 1.0 - lambda;  // do not overshoot
        const double target_lambda = lambda + dlam;
        Eigen::VectorXd target = target_lambda * f_ext;
        if (has_cf) target += constant_force;  // constant (non-ramped) load (geostatic gravity)
        cur_lambda = lambda;          // for the prescribed-displacement ramp (read by assemble)
        cur_target = target_lambda;

        Eigen::VectorXd du_free = Eigen::VectorXd::Zero(neq);
        // A phase in the fixed small-strain mode reads the stiffness of a new increment from the
        // last increment's displacement scaled to this one's size, until its own answer replaces
        // it; an attempt after a refresh starts from the answer it refreshes.
        if (small_freeze && small_refresh == 0 && !small_warm) {
            if (dlam_last > 0.0 && du_last.size() == neq) du_small = (dlam / dlam_last) * du_last;
            else if (du_small.size() != neq) du_small = Eigen::VectorXd::Zero(neq);
        }
        const bool warm_start = small_freeze && small_warm && du_small.size() == neq;
        if (warm_start) du_free = du_small;
        small_warm = false;
        prev_sc = committed;  // sigma_c,0 = sigma_0, where the iteration's stress walk starts
        std::fill(prev_deps.begin(), prev_deps.end(), 0.0);
        // Set when the predictor below moved the start. Such a start is not accepted as it stands:
        // its out-of-balance force can already be under the tolerance -- on a linear problem it
        // always is -- while carrying the rounding of the directional derivative, measured at
        // 1e-8 of a plate moment that one Newton correction returns to 1e-15. So at least one
        // correction is taken from it.
        bool predicted = false;
        stall_floor = 0.0;   // a new increment starts with the plain tangent
        if (presc_active) {
            // THE DISPLACEMENT PREDICTOR. A prescribed displacement enters an increment through
            // the fixed DOFs alone, so the first iterate is the committed field with only those
            // nodes moved: the whole increment is taken up by the one row of elements touching
            // them. Measured on a rigid footing pushed into sand (axisymmetric, 15-noded, c = 1,
            // phi = 30, psi = 0): that iterate carried an out-of-balance force 22 times the
            // reference, hundreds of stress points were returned from strains the converged
            // field never has, and Newton spent about 25 iterations walking out of that state
            // before converging quadratically in five -- so at 20 load steps the phase ran out
            // of iterations at 20% of the settlement, and at 100 steps at 37%.
            //
            // The predictor moves the free DOFs by the linear response to the increment first:
            // K_ff(committed) du = r(committed) - dlam * K_fp u_bar. K_fp u_bar is not assembled
            // anywhere -- the structural elements scatter free rows only -- so it is measured as
            // the directional derivative of the internal force along the prescribed motion, a
            // one-sided difference at the committed state. The start is KEPT only if its
            // out-of-balance force is below the plain start's, so no increment starts worse than
            // it did before, and what Newton converges to is decided by the residual, not by
            // where it started.
            const double tgt = cur_target;
            const double step = tgt - lambda;
            const double eps = 1.0e-5 * step;
            cur_target = lambda;  // the committed state: no prescribed increment yet
            builder.clear();
            const Eigen::VectorXd f0 = assemble(du_free, true, &builder);
            const math::CsrMatrix& k0 = builder.build_cached(kt_cache);
            cur_target = lambda + eps;
            const Eigen::VectorXd f_eps = assemble(du_free, false, nullptr);
            const Eigen::VectorXd r_lin = (target - f0) - (step / eps) * (f_eps - f0);
            cur_target = tgt;
            const double r_plain = (target - assemble(du_free, false, nullptr)).norm();
            try {
                const Eigen::VectorXd du0 = linear_solve(k0, r_lin);
                const double r_pred = (target - assemble(du0, false, nullptr)).norm();
                if (du0.allFinite() && r_pred < r_plain) {
                    du_free = du0;
                    predicted = true;
                }
            } catch (const math::SingularSystem&) {
                // No predictor at a singular committed tangent; the plain start stands.
            }
            ++result.timings.n_solve;
        }
        // THE SECANT START. The increment before this one is the best account there is of how
        // the ground is moving: under a load that keeps growing the same points go on yielding
        // and the field keeps its shape. The tangent at the committed state does not know that
        // -- a point lying on its yield surface with no strain increment yet reads as elastic --
        // so a start made from it puts the increment into ground far stiffer than it is.
        // Measured on a rigid footing pushed into sand (axisymmetric, 15-noded, 0.4 m): the
        // stiffness parameter at the first iterate was 0.73 against 0.09 at the converged one,
        // and ten of an increment's fifteen iterations went on finding the plastic zone again.
        // The last increment's displacement, scaled to this one's size, is tried as the start
        // and kept only where its out-of-balance force is the smallest of the starts on offer.
        // Under a load it does the same work: a strutted excavation in Hardening Soil took 624
        // iterations with it against 920, an embankment built undrained 1539 against 1897.
        if (!srf_mode && !warm_start && dlam_last > 0.0 &&
            du_last.size() == du_free.size()) {
            const Eigen::VectorXd du_sec = (dlam / dlam_last) * du_last;
            const double r_now = (target - assemble(du_free, false, nullptr)).norm();
            const double r_sec = (target - assemble(du_sec, false, nullptr)).norm();
            const bool keep = du_sec.allFinite() && std::isfinite(r_sec) && r_sec < r_now;
            if (keep) {
                du_free = du_sec;
                predicted = true;
            }
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: secant start %.4e against %.4e (%s)\n",
                             target_lambda, r_sec, r_now, keep ? "kept" : "dropped");
        }
        // A traction has no committed value stored anywhere -- the structural elements are
        // total-displacement formulations, so the state at the start of an increment is
        // whatever the first assembly of that increment computes. It is recorded there.
        bool first_iterate = true;
        // The last W residual norms of THIS increment, most recent last: the line search's
        // yardstick when the window is open. It starts empty at every increment (including a
        // retried one), so the memory never crosses a load step.
        std::vector<double> recent;
        // Not const: an increment that the monotone rule cannot close opens this once before it
        // is abandoned. See kStallEscalationWindow and the stall handler below.
        int ls_window = ls_open ? std::max(kStallEscalationWindow, options.line_search_window)
                                : std::max(1, options.line_search_window);
        bool step_converged = false;
        // The best iterate this attempt reached, for the stagnation acceptance below.
        double best_rnorm = std::numeric_limits<double>::infinity();
        double best_csp = 1.0;   // the stiffness parameter measured AT the best iterate
        Eigen::VectorXd best_du = du_free;
        // How this increment ended, if it did not converge. Set at each exit from the
        // iteration loop so the abandonment is named where it happens rather than
        // guessed afterwards from the iteration count.
        NewtonResult::Abandonment ended = NewtonResult::Abandonment::IterationBudget;
        int stall = 0;  // consecutive iterations that produced no descent
        int local_wait = 0;   // consecutive iterates in balance whose stress points are not settled
        bool early_cut = false;   // given up early, to be made smaller rather than retried
        bool early_stagnation = false;   // stagnation shown inside the acceptance band
        int creep = 0;  // consecutive iterations whose step was cut to under 1/4
        std::vector<double> progress;   // this attempt's residual norms, iterate by iterate
        std::vector<double> best_history;   // ... and the best of them so far, iterate by iterate
        stall_floor = 0.0;
        bool kept = false;              // the last accepted full step's assembly is in builder_t
        Eigen::VectorXd kept_fint;
        double kept_floor = 0.0;
        bool kept_consistent = false, kept_freeze = false;
        for (int iter = 0; iter < options.max_iterations; ++iter) {
            // Taken only under the controls it was made under: the stall floor moves after the line
            // search, and the tangent mode and the small-strain freeze between attempts.
            const bool reuse = kept && !elastic_mode &&
                               kept_floor == std::max(vertex_floor_mode ? kVertexFloor : 0.0, stall_floor) &&
                               kept_consistent == hs_consistent_mode && kept_freeze == small_freeze;
            kept = false;
            if (reuse) {
                std::swap(builder, builder_t);
                prev_sc.swap(prev_sc_t);
                prev_deps.swap(prev_deps_t);
                prev_iface.swap(prev_iface_t);
                prev_iface5.swap(prev_iface5_t);
                prev_skin.swap(prev_skin_t);
                prev_foot.swap(prev_foot_t);
                probe = probe_t;
                point_probe(probe, false);
            } else {
                builder.clear();
                probe.reset_counts();
                probe.first_iterate = first_iterate;
            }
            first_iterate = false;
            // The family is gathered on every accepted iterate and NOT on the line search's
            // trial evaluations, which is both correct (only an accepted iterate's local error
            // means anything) and why it is free: measured against a build that skipped it
            // entirely, the same five cases run in 1.1/2.4/4.3/6.8/0.8 s with it and
            // 1.2/2.9/4.6/7.5/0.7 s without -- inside the noise, with identical iteration
            // counts and bit-identical answers.
            const Eigen::VectorXd f_int =
                reuse ? std::move(kept_fint) : assemble(du_free, true, &builder, &probe);
            const Eigen::VectorXd residual = target - f_int;
            const double rnorm = residual.norm();
            if (rnorm < best_rnorm) {
                best_rnorm = rnorm;
                best_du = du_free;
            }
            measure_convergence(probe, f_int, residual, rnorm, cf_norm, ref, rtol,
                                result.convergence);
            if (rnorm <= best_rnorm) best_csp = result.convergence.csp;
            // The history is kept to the LARGEST window this increment could ever use, and only
            // the tail of it is consulted. Keeping just ls_window entries would make the stall
            // escalation below arrive with no memory at all -- it would be monotone for another
            // five iterations, which are exactly the five that were failing, and the escalation
            // measured as useless for that reason before this line was written.
            recent.push_back(rnorm);
            if ((int)recent.size() > kStallEscalationWindow) recent.erase(recent.begin());
            // What the trial step has to beat. With a window of 1 this is rnorm and the test
            // below is the strict Armijo condition, unchanged.
            const auto gate_first =
                recent.end() - std::min<std::ptrdiff_t>((std::ptrdiff_t)recent.size(), ls_window);
            const double gate = *std::max_element(gate_first, recent.end());
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
            const bool global_ok =
                csp_gate ? result.convergence.force_ok() : (rnorm <= rtol * ref);
            // AN INCREMENT IS SOLVED BEFORE IT IS ACCEPTED. At its first iterate the displacement
            // increment is zero and the out-of-balance force is the load just added. Where that
            // load is small against the scale the tolerance is a fraction of -- one fortieth of
            // an excavation stage against the weight of the whole model, with the 1e-2 tolerance
            // of the hardening models -- the first iterate already "converged", and the increment
            // was committed without the ground having moved: the load stayed in the residual
            // until enough of it had gathered to exceed the tolerance, and up to that much of the
            // stage was never applied at the end. The displacement of such a phase came out too
            // small, by more the smaller its load was against the model's weight. Only a first
            // iterate that is in balance to rounding -- a phase that changes nothing -- is kept.
            const bool solved_once = iter > 0 || rnorm <= 1e-10 * ref;
            // IN BALANCE, BUT THE STRESS POINTS HAVE NOT SETTLED. The local criteria read how far
            // the stress of a plastic point moved between two iterates, and a Newton step on a
            // ground that has softened to a few percent of its elastic stiffness is a large one:
            // the out-of-balance force is inside the tolerance and stays there, while the points
            // keep moving from one branch of their yield surface to another and the count of
            // unsettled ones never falls. Measured on an embankment built undrained on soft clay
            // (Hardening Soil): 0.84-0.95% out of balance against 1% tolerated, 96-100 points
            // unsettled against 58 allowed, for forty iterations, then cut back to nothing and
            // reported as stopped at a third of the load. What settles the points is small steps:
            // from here the increment is iterated with the elastic stiffness, whose corrections
            // are as small as the force still out of balance is.
            if (global_ok && solved_once && enforce_local && !elastic_mode && !presc_active &&
                !result.convergence.enforced_ok()) {
                if (++local_wait >= 3) {
                    elastic_mode = true;
                    if (debug)
                        std::fprintf(stderr, "  lambda %.4f iter %d: in balance but the stress "
                                             "points are not settled; continuing with the "
                                             "elastic stiffness\n", target_lambda, iter);
                }
            } else {
                local_wait = 0;
            }
            if (global_ok && solved_once &&
                (!enforce_local || result.convergence.enforced_ok()) &&
                !(predicted && iter == 0)) {
                // This iterate is the one about to be committed, so its integration is the one
                // the answer is walked along: record the guard here and nowhere else. Trial
                // iterates that were discarded are not part of the path and do not count.
                if (probe.integration_saturated > 0) {
                    ++result.convergence.saturated_increments;
                    result.convergence.saturated_points = std::max(
                        result.convergence.saturated_points, probe.integration_saturated);
                }
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
                if (debug) {
                    std::fprintf(stderr, "  lambda %.4f iter %d  SOLVE REFUSED: %s\n",
                                 target_lambda, iter, refusal.what());
                    // Which equations carry no stiffness at all: a row of zeros is a degree of
                    // freedom nothing holds, the commonest reason for a refusal that is not a
                    // mechanism.
                    std::vector<int> eq_gdof(neq, -1);
                    for (int g = 0; g < dofs.total_dofs(); ++g)
                        if (dofs.equation(g) >= 0) eq_gdof[dofs.equation(g)] = g;
                    int shown = 0;
                    for (int r = 0; r < kt.rows && shown < 12; ++r) {
                        double mx = 0.0;
                        for (int k = kt.row_ptr[r]; k < kt.row_ptr[r + 1]; ++k)
                            mx = std::max(mx, std::fabs(kt.values[k]));
                        if (mx > 0.0) continue;
                        const int g = eq_gdof[r];
                        const int nn = mesh.node_count;
                        if (g >= 0 && g < 2 * nn)
                            std::fprintf(stderr, "      zero row: eq %d = node %d comp %d at (%.4f, %.4f)\n",
                                         r, g / 2, g % 2, mesh.x[g / 2], mesh.y[g / 2]);
                        else
                            std::fprintf(stderr, "      zero row: eq %d = extra dof %d\n", r, g);
                        ++shown;
                    }
                }
            }
            result.timings.linear_solve += elapsed(t_lin);
            ++result.timings.n_solve;
            if (!solved) break;   // abandon the increment -> outer loop cuts back

            // Developer diagnostic (KATAI_NL_TANGENT_CHECK): is K_T the derivative of the
            // internal force along the direction just solved for? A Newton iteration whose full
            // step keeps being cut to alpha ~ 1e-3 near the solution is either walking a kinked
            // residual or carrying a tangent that is not its derivative; this tells which.
            if (debug && std::getenv("KATAI_NL_TANGENT_CHECK") != nullptr) {
                const Eigen::VectorXd kd = kt * delta;
                const double dn = delta.lpNorm<Eigen::Infinity>();
                for (double rel_h : {1e-4, 1e-6, 1e-8}) {
                    const double h = rel_h / std::max(dn, 1e-30);
                    const Eigen::VectorXd fp = assemble(du_free + h * delta, false, nullptr);
                    const Eigen::VectorXd fm = assemble(du_free - h * delta, false, nullptr);
                    const Eigen::VectorXd fd = (fp - fm) / (2.0 * h);
                    std::fprintf(stderr, "      tangent check h*|d|=%.0e: |fd-Kd|/|Kd| = %.3e\n",
                                 rel_h, (fd - kd).norm() / std::max(kd.norm(), 1e-300));
                }
            }

            // Backtracking line search: the consistent tangent gives fast local
            // convergence but a full Newton step can overshoot far from the
            // solution (the perfectly-plastic tangent is singular/indefinite at
            // yielding points). Shrink alpha until the residual decreases
            // (Armijo). This globalizes Newton without sacrificing its rate.
            //
            // The next trial is not simply half the last. The squared residual along the step,
            // phi(a) = |r(u + a delta)|^2 / 2, has phi(0) and -- delta being the Newton step --
            // phi'(0) = -|r|^2, and each rejected trial adds phi(a): the minimum of the quadratic
            // through the three is the next trial, safeguarded to [a/10, a/2] (Dennis & Schnabel
            // 1983, Alg. A6.3.1). The acceptance test is untouched, so it can only take fewer
            // evaluations to find the same kind of step, never accept a worse one; and its worst
            // case is the halving it replaces. Halving needed five residual evaluations to reach
            // the alpha = 1/16 a plastic mechanism typically asks for, and on a Hardening Soil
            // footing the line search's evaluations were three quarters of the whole run time.
            // The floor is the halving's own last trial, 2^-11, and a search that fails ends at
            // 2^-12 exactly as before.
            constexpr double kAlphaFloor = 1.0 / 2048.0;
            const double phi0 = 0.5 * rnorm * rnorm, slope0 = -rnorm * rnorm;
            double alpha = 1.0;
            bool improved = false;
            if (elastic_mode) {
                // The initial-stiffness iteration takes its full step; only a step that throws the
                // residual far out (a tenfold rise) is treated as a failure.
                const double rn = (target - assemble(du_free + delta, false, nullptr)).norm();
                if (std::isfinite(rn) && rn < 10.0 * rnorm) improved = true;
            }
            double alpha_rejected = 0.0;   // the smallest step the search turned down
            double rn_kept = 0.0;          // the residual at the step it took
            for (int ls = 0; ls < 12 && !improved; ++ls) {
                const bool full = ls == 0 && !elastic_mode;
                if (full) {
                    builder_t.clear();
                    prev_sc_t = prev_sc;
                    prev_deps_t = prev_deps;
                    prev_iface_t = prev_iface;
                    prev_iface5_t = prev_iface5;
                    prev_skin_t = prev_skin;
                    prev_foot_t = prev_foot;
                    probe_t = probe;
                    point_probe(probe_t, true);
                    probe_t.reset_counts();
                    probe_t.first_iterate = false;
                }
                const double floor_now = std::max(vertex_floor_mode ? kVertexFloor : 0.0, stall_floor);
                Eigen::VectorXd fi = full ? assemble(du_free + alpha * delta, true, &builder_t, &probe_t)
                                          : assemble(du_free + alpha * delta, false, nullptr);
                const double rn = (target - fi).norm();
                if (rn < (1.0 - 1.0e-4 * alpha) * gate) {
                    improved = true;
                    rn_kept = rn;
                    if (full) {
                        kept = true;
                        kept_fint = std::move(fi);
                        kept_floor = floor_now;
                        kept_consistent = hs_consistent_mode;
                        kept_freeze = small_freeze;
                    }
                    break;
                }
                if (alpha <= kAlphaFloor) break;
                alpha_rejected = alpha;
                const double phia = 0.5 * rn * rn;
                const double curv = phia - phi0 - slope0 * alpha;
                double next = curv > 0.0 ? -slope0 * alpha * alpha / (2.0 * curv) : 0.5 * alpha;
                if (!std::isfinite(next)) next = 0.5 * alpha;
                alpha = std::max(kAlphaFloor, std::clamp(next, 0.1 * alpha, 0.5 * alpha));
            }
            if (!improved) alpha = 0.5 * kAlphaFloor;
            // THE STEP IS LENGTHENED AGAIN WHERE THAT PAYS. Backtracking stops at the first step
            // that passes, and the interpolation that chooses the trials is built for a smooth
            // residual: across the kink of a yield surface it reads the rejected full step as
            // strong curvature and falls straight to a tenth of it. Measured on a strip footing
            // pushed into sand (15-noded, 0.4 m, phi = 30, psi = 0): step after step taken at
            // exactly 0.1 with the residual falling by exactly a tenth -- linear all the way to
            // the step taken, the kink somewhere in the nine tenths never looked at -- and the
            // creep that followed raised the stall floor, which gave the gain back: a cycle of six
            // iterations that moved the residual by 0.3%. A residual evaluation costs about a
            // fifteenth of an iteration, so up to three are spent between the step taken and the
            // one turned down, and the best is kept. A step can only get longer by this, and only
            // where the residual is smaller there.
            if (improved && !elastic_mode && alpha_rejected > alpha) {
                double lo = alpha, hi = alpha_rejected;
                for (int k = 0; k < 3 && hi > 1.5 * lo; ++k) {
                    const double mid = std::sqrt(lo * hi);
                    const double rn =
                        (target - assemble(du_free + mid * delta, false, nullptr)).norm();
                    if (std::isfinite(rn) && rn < rn_kept) {
                        lo = mid;
                        rn_kept = rn;
                    } else {
                        hi = mid;
                    }
                }
                alpha = lo;
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
            creep = (improved && alpha >= 0.25) ? 0 : creep + 1;
            if (elastic_mode) {
                // Its residual is not monotone, so it is judged on progress: not down by a
                // tenth over forty iterates is a stall.
                creep = 0;
                if (improved && progress.size() >= 40 &&
                    rnorm > 0.9 * progress[progress.size() - 40])
                    stall = 4;
            } else if (!improved && stall_floor > 0.0) {
                // The floor is a remedy for CREEP (a direction that descends a hair at a time),
                // not for no descent at all: with non-associated flow the stiffened tangent's
                // direction need not descend either, and raising the floor further only kept the
                // increment failing. Measured on a cohesive slope's K0 equilibrium step on 15-noded
                // elements (phi = 30, psi = 0): with the floor still climbing the step was abandoned
                // at 61% of the load, while the exact tangent carried it to the end. A step that
                // finds no descent hands the next one back to Newton's own direction.
                stall_floor = 0.0;
                creep = 0;
            } else if (improved && alpha >= 0.5) {
                stall_floor *= 0.5;                 // let go gradually: dropping it at once
                if (stall_floor < 0.01) stall_floor = 0.0;   // put the creep straight back
            } else if (creep >= 3) {
                stall_floor = std::min(1.0, stall_floor > 0.0 ? 2.0 * stall_floor : 0.1);
            }
            progress.push_back(rnorm);
            // AN INCREMENT THAT IS TOO LARGE IS MADE SMALLER BEFORE IT IS RETRIED. The retries
            // below -- the consistent tangent, the open line search, the elastic stiffness --
            // each start the same increment again from its beginning, and only when all of them
            // have failed is it halved. Where the increment is simply too large a step for the
            // state of the ground, that order spends the work first and finds the remedy last:
            // measured on an embankment built undrained on soft clay (Hardening Soil, 6-noded,
            // 1 m), the increment in which the fill's core reaches its unconfined strength took
            // 162 iterations, then 264 with the line search open, then a third pass, before it
            // was halved -- and the half took 27. For the models solved to a tolerance of 1e-3
            // or looser (the hardening family, whose iteration is the expensive one), an
            // increment still open after kEarlyCutIterations and not yet near the smallest size
            // is halved at once; the retries are kept for the increment that cannot be made
            // smaller. A load held constant (strength reduction) and an imposed displacement
            // are left as they were.
            constexpr int kEarlyCutIterations = 40;
            if (rtol >= 1.0e-3 && !srf_mode && !presc_active && dlam > 2.0 * min_dlam &&
                iter + 1 >= kEarlyCutIterations) {
                ended = NewtonResult::Abandonment::IterationBudget;
                early_cut = true;
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f iter %d: still open; the increment is "
                                         "halved before it is retried\n", target_lambda, iter);
                break;
            }
            // STAGNATION SHOWN, NOT AWAITED. The acceptance below keeps an increment at its best
            // iterate once every retry has been spent on it -- which, on a phase whose tolerance
            // is out of reach, is every increment, each ground through hundreds of iterations
            // first: measured on an embankment built on clay (Mohr-Coulomb, 15-noded), 3478
            // iterations for 20 increments, each accepted at about 1e-4 of the load scale in the
            // end, 396 s. An increment whose best iterate is already inside the acceptance band
            // and whose residual has not fallen by a tenth over its last thirty iterates is in the
            // state the retries end in; it is accepted there, on the same terms (counted, K2D-A019).
            // An iteration that is still converging -- quadratically or not -- never meets it.
            if (rtol < kStagnationAccept && progress.size() > kStagnationWindow &&
                best_rnorm <= kStagnationAccept * ref &&
                rnorm > 0.9 * progress[progress.size() - 1 - kStagnationWindow]) {
                ended = NewtonResult::Abandonment::NoDescent;
                early_stagnation = true;
                break;
            }
            // THE SAME OUTSIDE THE BAND. An attempt whose best iterate has not come down by a tenth
            // in thirty iterations will not close in the hundred and seventy it has left, wherever
            // it stands: measured on a strip footing pushed into sand (15-noded, 0.4 m, phi = 30,
            // psi = 0), five attempts in a row ran their whole budget of 200 at 5.5e-3 of the
            // reference, each one a cycle of six iterations that moved the residual by 0.3%. It is
            // ended there and handed to the remedies below -- the open window, the smaller
            // increment -- thirty iterations after it stopped instead of at the budget.
            best_history.push_back(best_rnorm);
            if (!elastic_mode && best_history.size() > kStagnationWindow &&
                best_rnorm > 0.9 * best_history[best_history.size() - 1 - kStagnationWindow]) {
                ended = NewtonResult::Abandonment::NoDescent;
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f iter %d: not down by a tenth in %d "
                                         "iterations; ended\n", target_lambda, iter,
                                 (int)kStagnationWindow);
                break;
            }
            if (stall >= 4) {
                ended = NewtonResult::Abandonment::NoDescent;
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f iter %d  ABANDONED: four consecutive "
                                         "iterations without descent%s\n", target_lambda, iter,
                                 ls_open ? " (with the non-monotone window already open)" : "");
                // WHERE the out-of-balance force that would not fall sits (debug only): the ten
                // largest entries with their node, and the share of the whole they carry. A stall
                // that lives in a handful of nodes is a local defect -- a stress point returning
                // a state its neighbours cannot balance -- and one spread over the mesh is not.
                if (debug) {
                    std::vector<int> eq_gdof(static_cast<size_t>(residual.size()), -1);
                    for (int g = 0; g < dofs.total_dofs(); ++g)
                        if (dofs.equation(g) >= 0) eq_gdof[(size_t)dofs.equation(g)] = g;
                    std::vector<int> order(static_cast<size_t>(residual.size()));
                    for (size_t k = 0; k < order.size(); ++k) order[k] = (int)k;
                    const size_t top = std::min<size_t>(10, order.size());
                    std::partial_sort(order.begin(), order.begin() + (long)top, order.end(),
                                      [&residual](int a, int b) {
                                          return std::fabs(residual(a)) > std::fabs(residual(b));
                                      });
                    double share = 0.0;
                    for (size_t k = 0; k < top; ++k) share += residual(order[k]) * residual(order[k]);
                    std::fprintf(stderr, "    the ten largest residual entries carry %.1f%% of "
                                         "its square:\n", 100.0 * share / residual.squaredNorm());
                    for (size_t k = 0; k < top; ++k) {
                        const int g = eq_gdof[(size_t)order[k]];
                        const int node = g >= 0 && g < 2 * mesh.node_count ? g / 2 : -1;
                        if (node >= 0)
                            std::fprintf(stderr, "      r=%+.4e  node %d (%.3f, %.3f) %s\n",
                                         residual(order[k]), node, mesh.x[node], mesh.y[node],
                                         g % 2 == 0 ? "x" : "y");
                        else
                            std::fprintf(stderr, "      r=%+.4e  structural dof %d\n",
                                         residual(order[k]), g);
                    }
                }
                break;
            }
            du_free += alpha * delta;
        }

        // A budget that ran out on an iteration which had stopped making progress -- its residual
        // not down by a tenth over its last thirty iterates, a rate no budget closes -- ended in a
        // stall, and is recorded as one. The iteration itself is not cut short (a slow increment
        // that is still converging keeps its whole budget); only the reading of how it ended
        // changes. With the stall floor a collapsing footing finds a hair of descent at every
        // iterate, so it used to end on the budget, and a budget is (rightly) never published as
        // a capacity: measured on the Prandtl footing of the corpus, which then refused to name
        // its limit load.
        if (!step_converged && ended == NewtonResult::Abandonment::IterationBudget &&
            progress.size() > 30 &&
            progress.back() > 0.9 * progress[progress.size() - 31])
            ended = NewtonResult::Abandonment::NoDescent;

        // An increment that stopped descending within options.stagnation_accept of the force
        // scale is kept before anything is retried (see the option).
        if (!step_converged && options.stagnation_accept > rtol &&
            ended == NewtonResult::Abandonment::NoDescent &&
            best_rnorm <= options.stagnation_accept * ref) {
            du_free = best_du;
            assemble(du_free, false, nullptr);  // the trial states of the iterate being kept
            result.convergence.csp = best_csp;
            ++result.stagnation_accepted;
            result.worst_accepted_residual =
                std::max(result.worst_accepted_residual, best_rnorm / ref);
            step_converged = true;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: stopped descending at %.3e of the reference; kept\n",
                             target_lambda, best_rnorm / ref);
        }

        // STAGNATION ACCEPTANCE. Only an increment that is about to be CUT BACK is looked at --
        // every retry above still gets its chance first, so a run that converged before this rule
        // existed converges on exactly the same path. See kStagnationAccept for the measurement.
        const bool would_cut_back =
            early_stagnation ||
            (!step_converged &&
             !(ended == NewtonResult::Abandonment::SolveRefused && !vertex_floor_mode) &&
             !(has_fd_tangent && !hs_consistent_mode) &&
             !(ended == NewtonResult::Abandonment::NoDescent && !ls_open) &&
             !(ended == NewtonResult::Abandonment::NoDescent && !elastic_mode && !presc_active));
        if (would_cut_back && ended != NewtonResult::Abandonment::SolveRefused &&
            rtol < kStagnationAccept && best_rnorm <= kStagnationAccept * ref) {
            du_free = best_du;
            assemble(du_free, false, nullptr);  // the trial states of the iterate being kept
            // ... and its stiffness parameter: the measurement standing in result.convergence is
            // the LAST iterate's, which on a stagnated increment can be one thrown far out -- read
            // as the committed state it reported a mechanism (0.004) where the kept iterate had 0.9.
            result.convergence.csp = best_csp;
            ++result.stagnation_accepted;
            result.worst_accepted_residual =
                std::max(result.worst_accepted_residual, best_rnorm / ref);
            step_converged = true;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: stagnated at %.3e of the reference; accepted\n",
                             target_lambda, best_rnorm / ref);
        }

        // AN IMPOSED DISPLACEMENT ALWAYS HAS AN ANSWER. Under a load a phase can end because the
        // ground cannot carry it; under a prescribed movement there is no such limit -- the ground
        // takes the movement and the reaction is whatever it turns out to be. An increment of it
        // that can be cut no further and still will not close is therefore never a collapse, only
        // an iteration that has stopped: with non-associated flow the out-of-balance force settles
        // on the kinks of the yield surfaces and Newton's direction, exact as it is, crosses them
        // at a hundredth of its length. Measured on a strip footing pushed into sand (15-noded,
        // 0.4 m, phi = 30, psi = 0): 199 iterations at 5.5e-3 of the reference with the tangent
        // matching its finite difference to 1e-7, and the phase given up at 59% of the settlement.
        // Such an increment is kept at its best iterate if that is within kLastResortAccept, and
        // reported with the force left out of balance (K2D-A019) like every other one that is.
        if (!step_converged && presc_active && !srf_mode && would_cut_back &&
            ended != NewtonResult::Abandonment::SolveRefused && 0.5 * dlam < min_dlam &&
            best_rnorm <= kLastResortAccept * ref) {
            du_free = best_du;
            assemble(du_free, false, nullptr);  // the trial states of the iterate being kept
            result.convergence.csp = best_csp;
            ++result.stagnation_accepted;
            result.worst_accepted_residual =
                std::max(result.worst_accepted_residual, best_rnorm / ref);
            step_converged = true;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: cannot be cut further; kept at %.3e of the "
                                     "reference\n", target_lambda, best_rnorm / ref);
        }

        // The answer of a fixed-stiffness pass becomes the increment the stiffness is read from,
        // and the increment is solved again from it (see small_freeze).
        if (step_converged && small_freeze && small_refresh < kSmallRefresh) {
            const double moved = (du_free - du_small).norm();
            if (moved > 1.0e-6 * std::max(du_free.norm(), 1.0e-300)) {
                if (debug)
                    std::fprintf(stderr, "  lambda %.4f: small-strain stiffness re-read from the "
                                         "answer (pass %d, moved %.3e of |du|)\n",
                                 target_lambda, small_refresh + 1,
                                 moved / std::max(du_free.norm(), 1.0e-300));
                du_small = du_free;
                small_warm = true;
                ++small_refresh;
                ++result.small_strain_passes;
                elastic_mode = false;
                continue;
            }
        }

        if (step_converged) {
            if (small_freeze && du_small.size() == du_free.size()) {
                ++result.small_strain_increments;
                result.small_strain_gap =
                    std::max(result.small_strain_gap,
                             (du_free - du_small).norm() / std::max(du_free.norm(), 1.0e-300));
            }
            small_refresh = 0;
            hs_consistent_mode = false;  // increment closed → back to the cheap continuum path
            elastic_mode = false;        // and back to Newton
            du_last = du_free;
            dlam_last = dlam;
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
            ++result.timings.n_steps;
            result.load_factor = lambda;
            dlam = std::min(init_dlam, dlam * 1.5);  // allow recovery
        } else if (has_small && !small_freeze &&
                   ended != NewtonResult::Abandonment::SolveRefused) {
            // The same increment again, with the small-strain stiffness read from the best
            // iterate this attempt reached (see small_freeze). Before anything is cut: the
            // increment was not shown to be too large, only not smooth.
            small_freeze = true;
            small_refresh = 0;
            du_small = best_du.allFinite() ? best_du : Eigen::VectorXd::Zero(neq);
            small_warm = true;
            elastic_mode = false;
            ++result.small_strain_freezes;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: retrying the increment with the small-strain "
                                     "stiffness read from its best iterate (dlam untouched)\n",
                             lambda);
        } else if (early_cut) {
            // Halved without a retry and without counting as an abandonment: the increment was
            // not shown to be unsolvable, only too large (see kEarlyCutIterations).
            elastic_mode = false;
            dlam *= 0.5;
            if (small_freeze) { du_small *= 0.5; small_refresh = 0; }
            if (debug)
                std::fprintf(stderr, "  early cut at lambda %.4f: dlam %.6f (minimum %.6f)\n",
                             lambda, dlam, min_dlam);
        } else if (ended == NewtonResult::Abandonment::SolveRefused && !vertex_floor_mode) {
            // A refused solve is first asked whether it was a vertex and not a mechanism: the SAME
            // increment is retried with the floor on and dlam untouched, and the phase keeps the
            // floor, the same latch as the two retries below and for the same reason.
            vertex_floor_mode = true;
            ++result.vertex_floor_retries;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: tangent refused; retrying the increment with "
                                     "the vertex floor (dlam untouched)\n", lambda);
        } else if (has_fd_tangent && !hs_consistent_mode) {
            // Hybrid tangent: continuum failed to converge this increment → retry the SAME
            // increment with the consistent (FD) tangent without touching dlam; the
            // remaining increments stay consistent too.
            hs_consistent_mode = true;
        } else if (ended == NewtonResult::Abandonment::NoDescent && !ls_open) {
            // A STALL IS WHAT THE NON-MONOTONE RULE EXISTS FOR, so the increment is re-entered
            // with the window open before its size is touched -- the same shape as the hybrid
            // tangent above, and for the same reason: when an increment cannot be closed the
            // cheap way, try the stronger tool on THAT increment rather than paying for it
            // everywhere. Paying for it everywhere was measured and costs accuracy at a fixed
            // tolerance: binding the local criteria at the shipped tolerance lands 0.147% from the
            // four-decades-tighter answer with the window always open, against 0.05% without it.
            // The remaining increments keep the window, exactly as they keep the consistent
            // tangent, because a phase that needed it once is a phase that is in that regime.
            ls_open = true;
            ++result.line_search_escalations;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: retrying the increment with a line-search "
                                     "window of %d (dlam untouched)", lambda,
                             kStallEscalationWindow);
        } else if (ended == NewtonResult::Abandonment::NoDescent && !elastic_mode && !presc_active) {
            // NOT FOR A PRESCRIBED DISPLACEMENT. The elastic-stiffness iteration was brought in for
            // load-driven phases (initial phases on sloping ground, where Newton's direction stops
            // descending) and was measured there. Under an imposed displacement it is the wrong
            // remedy: the increment that Newton cannot close is closed by making it smaller, and the
            // slowly converging elastic iteration put in front of the cut-back sent the phase down a
            // different path -- an axisymmetric footing pushed into sand (15-noded, 0.25 m, c = 1,
            // phi = 30, psi = 0) stopped at 13% of its settlement where without it, as before the
            // retry existed, it reaches the full settlement at 604.5 kN.
            elastic_mode = true;
            newton_convergence = result.convergence;
            newton_convergence_kept = true;
            ++result.elastic_retries;
            if (debug)
                std::fprintf(stderr, "  lambda %.4f: retrying the increment with the elastic "
                                     "tangent (dlam untouched)\n", lambda);
        } else {
            // The increment is abandoned: record WHY before the size is halved, so the
            // phase can report the reason it actually stopped for.
            if (ended == NewtonResult::Abandonment::NoDescent) ++result.no_descent;
            else if (ended == NewtonResult::Abandonment::IterationBudget) ++result.budget_exhausted;
            result.last_abandonment = ended;
            if (elastic_mode && newton_convergence_kept) result.convergence = newton_convergence;
            dlam *= 0.5;  // sub-stepping: shrink the increment and retry
            if (small_freeze) { du_small *= 0.5; small_refresh = 0; }
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
    if (std::getenv("KATAI_NL_PROFILE") != nullptr) {
        const auto& t = result.timings;
        std::fprintf(stderr,
                     "[nl-profile] total %.3f s | tangent %.3f s (%d) | residual %.3f s (%d) | "
                     "csr %.3f s | linear %.3f s (%d) | other %.3f s | neq %d | %s\n",
                     t.total, t.assemble_tangent, t.n_tangent, t.assemble_residual, t.n_residual,
                     t.csr_build, t.linear_solve, t.n_solve,
                     t.total - t.assemble_tangent - t.assemble_residual - t.csr_build -
                         t.linear_solve,
                     neq, result.converged ? "converged" : "not converged");
    }
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


// The excess pore pressure's force (excess_pore_force.hpp). Written beside the internal force it
// has to agree with: the same Gauss rule and weight (gauss weight x Kin weight), the same B, the same
// pore direction m, and Kw/n from the stress point's profiled material -- the lines the assembler
// adds to sigma', kept separate because the baseline must not re-integrate the constitutive law.
template <class E, class Kin>
void excess_pore_force_impl(const mesh::Mesh& mesh, const DofMap& dofs,
                            const std::vector<MaterialModel>& materials,
                            const std::vector<MaterialProfile>& profile,
                            const std::vector<GaussState>& states,
                            const std::vector<char>& active_element, Eigen::VectorXd& rhs,
                            bool every_dof = false, std::vector<double>* at_gauss = nullptr) {
    const auto gauss = E::gauss_points();
    const int n_gp = E::kGaussCount;
    const typename Kin::Strain mvec = Kin::pore_vector();
    typename E::NodeCoords coords;
    std::array<int, E::kDofCount> edofs;
    for (int e = 0; e < mesh.element_count; ++e) {
        if (!active_element.empty() && !active_element[e]) continue;
        const int e_mat = mesh.element_material[e];
        const MaterialModel& mat = materials[e_mat];
        const MaterialProfile prof = e_mat < (int)profile.size() ? profile[e_mat] : MaterialProfile{};
        bool gathered = false;
        for (int g = 0; g < n_gp; ++g) {
            const GaussState& s = states[(size_t)e * n_gp + g];
            if (s.pw_carried == 0.0 && !(mat.undrained && s.eps_vol_und != 0.0)) continue;
            if (!gathered) { detail::gather_element<E>(mesh, dofs, e, coords, edofs); gathered = true; }
            double pw = 0.0;
            if (mat.undrained && s.eps_vol_und != 0.0) {
                double kwn = mat.kw_over_n(mat.undrained_poisson);
                if (!prof.uniform()) {
                    const typename E::ShapeValues sh = E::shape_functions(gauss[g].xi, gauss[g].eta);
                    double y = 0.0;
                    for (int i = 0; i < E::kNodeCount; ++i) y += sh(i) * coords(i, 1);
                    MaterialModel mg = mat;
                    mg.youngs_modulus = profile_at(mat.youngs_modulus, prof.E_inc, prof.y_ref, y);
                    kwn = mg.kw_over_n(mg.undrained_poisson);
                }
                pw = kwn * s.eps_vol_und;
            }
            if (s.pw_carried != 0.0) pw += s.pw_carried;
            if (at_gauss) {   // the pressure itself, point by point; no force wanted
                (*at_gauss)[(size_t)e * n_gp + g] = pw;
                continue;
            }
            const auto grad = Kin::template gradients<E>(coords, gauss[g].xi, gauss[g].eta);
            const typename Kin::Strain sig = pw * mvec;
            const auto fe = ((gauss[g].weight * grad.weight) * grad.B.transpose() * sig).eval();
            for (int a = 0; a < E::kDofCount; ++a) {
                if (every_dof) {   // by global DOF, supported ones included (the reactions)
                    if (edofs[a] >= 0 && edofs[a] < rhs.size()) rhs(edofs[a]) += fe(a);
                    continue;
                }
                const int eq = dofs.equation(edofs[a]);
                if (eq >= 0) rhs(eq) += fe(a);
            }
        }
    }
}

template <class E, class Kin>
void excess_pore_force_every_dof(const mesh::Mesh& mesh, const DofMap& dofs,
                                 const std::vector<MaterialModel>& materials,
                                 const std::vector<MaterialProfile>& profile,
                                 const std::vector<GaussState>& states,
                                 const std::vector<char>& active_element, Eigen::VectorXd& full) {
    excess_pore_force_impl<E, Kin>(mesh, dofs, materials, profile, states, active_element, full, true);
}

} // namespace

std::vector<double> excess_pore_at_gauss(const mesh::Mesh& mesh, const DofMap& dofs,
                                         const std::vector<MaterialModel>& materials,
                                         const std::vector<MaterialProfile>& profile,
                                         const std::vector<GaussState>& states,
                                         const std::vector<char>& active_element) {
    const bool t15 = mesh.nodes_per_element == Tri15Element::kNodeCount;
    std::vector<double> pw(states.size(), 0.0);
    if (states.size() != (size_t)mesh.element_count *
                             (t15 ? Tri15Element::kGaussCount : Tri6Element::kGaussCount))
        return pw;
    Eigen::VectorXd unused;
    // The pressure does not depend on the kinematics (only its force does).
    if (t15)
        excess_pore_force_impl<Tri15Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile, states, active_element, unused, false, &pw);
    else
        excess_pore_force_impl<Tri6Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile, states, active_element, unused, false, &pw);
    return pw;
}

void add_excess_pore_nodal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                                 const std::vector<MaterialModel>& materials,
                                 const std::vector<MaterialProfile>& profile,
                                 const std::vector<GaussState>& states,
                                 const std::vector<char>& active_element, bool axisymmetric,
                                 Eigen::VectorXd& full) {
    const bool t15 = mesh.nodes_per_element == Tri15Element::kNodeCount;
    if (states.size() != (size_t)mesh.element_count *
                             (t15 ? Tri15Element::kGaussCount : Tri6Element::kGaussCount))
        return;
    if (t15 && axisymmetric)
        excess_pore_force_every_dof<Tri15Element, detail::AxisymKin>(mesh, dofs, materials, profile, states, active_element, full);
    else if (t15)
        excess_pore_force_every_dof<Tri15Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile, states, active_element, full);
    else if (axisymmetric)
        excess_pore_force_every_dof<Tri6Element, detail::AxisymKin>(mesh, dofs, materials, profile, states, active_element, full);
    else
        excess_pore_force_every_dof<Tri6Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile, states, active_element, full);
}

void add_excess_pore_force(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<MaterialModel>& materials,
                           const std::vector<MaterialProfile>& profile,
                           const std::vector<GaussState>& states,
                           const std::vector<char>& active_element, bool axisymmetric,
                           Eigen::VectorXd& rhs) {
    if (states.size() != (size_t)mesh.element_count *
                             (mesh.nodes_per_element == Tri15Element::kNodeCount
                                  ? Tri15Element::kGaussCount : Tri6Element::kGaussCount))
        return;   // no committed state of this mesh: nothing is carried
    if (mesh.nodes_per_element == Tri15Element::kNodeCount) {
        if (axisymmetric)
            excess_pore_force_impl<Tri15Element, detail::AxisymKin>(mesh, dofs, materials, profile,
                                                                    states, active_element, rhs);
        else
            excess_pore_force_impl<Tri15Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile,
                                                                         states, active_element, rhs);
        return;
    }
    if (axisymmetric)
        excess_pore_force_impl<Tri6Element, detail::AxisymKin>(mesh, dofs, materials, profile, states,
                                                               active_element, rhs);
    else
        excess_pore_force_impl<Tri6Element, detail::PlaneStrainKin>(mesh, dofs, materials, profile,
                                                                    states, active_element, rhs);
}

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

// ---- StructureLaws: the structures alone, by the static solver's assembler -------------------
struct StructureLaws::Impl {
    virtual ~Impl() = default;
    virtual void add(const Eigen::VectorXd& v, Eigen::VectorXd& f, math::SparseMatrixBuilder* K) = 0;
    virtual void set_soil_pore(const std::vector<double>& nodal_excess,
                               const std::function<double(int, double)>& share) = 0;
    virtual void commit() = 0;
    virtual StructuralInit committed() const = 0;
    virtual void reset(const StructuralInit& s) = 0;
    virtual void set_offset(const Eigen::VectorXd& v0) = 0;
};

namespace {
template <class E, class Kin>
struct StructureLawsImpl final : StructureLaws::Impl {
    using Asm = detail::InternalForceAssembler<E, Kin>;
    const mesh::Mesh& mesh;
    std::vector<char> active;               // the soil's: a joint reads the ground beside it
    std::vector<MaterialProfile> no_profile;
    Asm fasm;
    // The soil's stress points, only for the excess pore pressure a joint reads beside it: the
    // assembler is told to leave the soil alone (structures_only), so nothing else is read.
    std::vector<GaussState> soil_c, soil_t;
    std::vector<double> anchor_c, geogrid_c, iface_c, iface5_c, eskin_c, efoot_c, plate_c, plate5_c;
    std::vector<double> anchor_t, geogrid_t, iface_t, iface5_t, eskin_t, efoot_t, plate_t, plate5_t;
    Eigen::VectorXd u_datum, offset;
    typename Asm::State st;

    StructureLawsImpl(const mesh::Mesh& m, const DofMap& dofs,
                      const std::vector<MaterialModel>& materials, const Structures& structures,
                      const StructuralInit& init, const std::vector<char>& act)
        : mesh(m), active(act), fasm(m, dofs, materials, active, structures, no_profile) {
        fasm.structures_only = true;
        const int neq = dofs.equation_count();
        u_datum = init.u_datum.size() == neq ? init.u_datum : Eigen::VectorXd::Zero(neq);
        offset = Eigen::VectorXd::Zero(neq);
        soil_c.assign(static_cast<size_t>(m.element_count) * E::kGaussCount, GaussState{});
        soil_t = soil_c;
        seed_structural_state(structures, init, anchor_c, geogrid_c, iface_c, iface5_c, eskin_c,
                              efoot_c, plate_c, plate5_c);
        anchor_t = anchor_c; geogrid_t = geogrid_c; iface_t = iface_c; iface5_t = iface5_c;
        eskin_t = eskin_c; efoot_t = efoot_c; plate_t = plate_c; plate5_t = plate5_c;
        st.committed = &soil_c;    st.trial = &soil_t;
        st.anchor_c = &anchor_c;   st.anchor_t = &anchor_t;
        st.geogrid_c = &geogrid_c; st.geogrid_t = &geogrid_t;
        st.iface_c = &iface_c;     st.iface_t = &iface_t;
        st.iface5_c = &iface5_c;   st.iface5_t = &iface5_t;
        st.eskin_c = &eskin_c;     st.eskin_t = &eskin_t;
        st.efoot_c = &efoot_c;     st.efoot_t = &efoot_t;
        st.plate_c = &plate_c;     st.plate_t = &plate_t;
        st.plate5_c = &plate5_c;   st.plate5_t = &plate5_t;
    }
    void add(const Eigen::VectorXd& v, Eigen::VectorXd& f, math::SparseMatrixBuilder* K) override {
        const typename Asm::Ramp no_ramp;
        const Eigen::VectorXd u = u_datum + offset + v;
        f += fasm.assemble(u, Eigen::VectorXd::Zero(u.size()), /*build_tangent=*/K != nullptr,
                           TangentMode::kContinuum, st, no_ramp, K, nullptr);
    }
    void set_soil_pore(const std::vector<double>& nodal_excess,
                       const std::function<double(int, double)>& share) override {
        store_excess_pore_pressure(mesh, nodal_excess, active, share, soil_c);
    }
    void commit() override {
        anchor_c = anchor_t; geogrid_c = geogrid_t; iface_c = iface_t; iface5_c = iface5_t;
        eskin_c = eskin_t; efoot_c = efoot_t; plate_c = plate_t; plate5_c = plate5_t;
    }
    StructuralInit committed() const override {
        StructuralInit s;
        s.u_datum = u_datum;
        s.anchor_plastic = anchor_c; s.geogrid_plastic = geogrid_c;
        s.interface_slip = iface_c; s.interface5_slip = iface5_c;
        s.embedded_skin_slip = eskin_c; s.embedded_foot_slip = efoot_c;
        s.plate_plastic = plate_c; s.plate5_plastic = plate5_c;
        return s;
    }
    void reset(const StructuralInit& s) override {
        anchor_c = s.anchor_plastic; geogrid_c = s.geogrid_plastic;
        iface_c = s.interface_slip; iface5_c = s.interface5_slip;
        eskin_c = s.embedded_skin_slip; efoot_c = s.embedded_foot_slip;
        plate_c = s.plate_plastic; plate5_c = s.plate5_plastic;
        anchor_t = anchor_c; geogrid_t = geogrid_c; iface_t = iface_c; iface5_t = iface5_c;
        eskin_t = eskin_c; efoot_t = efoot_c; plate_t = plate_c; plate5_t = plate5_c;
    }
    void set_offset(const Eigen::VectorXd& v0) override {
        if (v0.size() == offset.size()) offset = v0;
    }
};
}  // namespace

StructureLaws::StructureLaws(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<MaterialModel>& materials,
                             const Structures& structures, const StructuralInit& init,
                             const std::vector<char>& active, Kinematics kinematics) {
    const bool axi = kinematics == Kinematics::Axisymmetric;
    if (mesh.nodes_per_element == Tri15Element::kNodeCount) {
        if (axi) impl_ = std::make_unique<StructureLawsImpl<Tri15Element, detail::AxisymKin>>(mesh, dofs, materials, structures, init, active);
        else impl_ = std::make_unique<StructureLawsImpl<Tri15Element, detail::PlaneStrainKin>>(mesh, dofs, materials, structures, init, active);
    } else {
        if (axi) impl_ = std::make_unique<StructureLawsImpl<Tri6Element, detail::AxisymKin>>(mesh, dofs, materials, structures, init, active);
        else impl_ = std::make_unique<StructureLawsImpl<Tri6Element, detail::PlaneStrainKin>>(mesh, dofs, materials, structures, init, active);
    }
}
StructureLaws::~StructureLaws() = default;
void StructureLaws::add(const Eigen::VectorXd& v, Eigen::VectorXd& f, math::SparseMatrixBuilder* K) {
    impl_->add(v, f, K);
}
void StructureLaws::set_soil_pore(const std::vector<double>& nodal_excess,
                                  const std::function<double(int, double)>& share) {
    impl_->set_soil_pore(nodal_excess, share);
}
void StructureLaws::commit() { impl_->commit(); }
StructuralInit StructureLaws::committed() const { return impl_->committed(); }
void StructureLaws::reset(const StructuralInit& state) { impl_->reset(state); }
void StructureLaws::set_offset(const Eigen::VectorXd& v0) { impl_->set_offset(v0); }

} // namespace katai::core

