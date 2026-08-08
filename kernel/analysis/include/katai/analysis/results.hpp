#pragma once
// Result types of a staged analysis -- what a solved phase hands to whoever
// asked for it (Stage B1: extracted from the application driver, because a
// result is engine vocabulary; a front end only displays it).
//
// Everything here is engine-typed on purpose. The one field that used to leak
// the project schema into the result -- the design approach -- carries the
// engine's own enum (design_code.hpp); the schema value is mapped once, at the
// driver seam, by to_core_design_approach.

#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/design_code.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/analysis/structural_forces.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// One thing the run did that the file does not literally say: an object clipped to the
// soil it could attach to, a parameter the selected model does not read, a fallback taken.
// The rule this vocabulary exists to enforce: an input may be treated differently from the
// way it was written, but never in silence. An object that would contribute NOTHING is not
// a diagnostic at all -- it is a refusal (SolveResult::ok = false with the message), because
// a load or a wall that does nothing is a modelling error, and a plausible answer to the
// wrong model is worse than no answer.
//
// `code` is a stable machine tag: front ends may reword nothing, tests match on it, and a
// user can grep a log for it. Codes are never reused once retired.
enum class DiagnosticSeverity {
    Note,      // worth stating; the run is exactly what the file asks for
    Warning,   // the run is defensible, but it is not literally what the file asks for
    Refusal    // the run stopped here; SolveResult::message carries the same sentence
};

struct Diagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::Warning;
    std::string code;      // e.g. "K2D-W101" -- stable, greppable, never reworded
    std::string subject;   // the object's name as the user wrote it (may be empty)
    std::string message;   // engineer-readable: what was found, and what was done instead
};

// One structural element's internal-force OUTPUT (PLAXIS Output -> Structures -> M/Q/N):
// force stations along the element + the |max| envelope, computed from the converged solution
// by the validated post-processors in structural_forces.hpp.
// kind: 0 = plate / embedded wall (N, Q, M), 1 = anchor (single N), 2 = geogrid (N, tension-only).
struct StructForce {
    std::string name;
    int kind = 0;
    std::vector<ForceStation> stations;
    bool yielded = false;   // anchor at F_max / plate M-N hinge formed (committed-state event)
    double max_N = 0.0, max_Q = 0.0, max_M = 0.0;
    // Dynamic (seismic) phase only: the stations hold the ENVELOPE over the shaking -- each station's
    // N/Q/M is max|.| over every time step, NOT one instant. The extremes at different stations
    // generally occur at DIFFERENT times, so this is a design envelope, not a state of equilibrium.
    // ux/uy still come from the peak-response instant (they only drive the deformed overlay).
    // NEW FIELDS GO AT THE END: results_io writes these positionally.
    bool envelope = false;
    // Dynamic phase: true when the parent phase's STATIC action was superposed, so the stations hold
    // the TOTAL design action (static + dynamic) -- what the section is designed for. False = the
    // dynamic action alone (no usable parent state); the user must superpose by hand.
    bool superposed = false;
};

// One interface's results (PLAXIS Output -> Interfaces): shear tau, normal effective stress sigma_n,
// relative slip and opening along the joint, from the converged solution via the post-processor in
// structural_forces.hpp. any_slip = a Coulomb limit (plastic slip) was reached somewhere.
struct InterfaceResult {
    std::string name;
    std::vector<InterfaceStation> stations;
    double max_abs_tau = 0.0, max_abs_sigma_n = 0.0, max_abs_slip = 0.0;
    bool any_slip = false;
    // Dynamic (seismic) phase only, like StructForce::envelope: the stations hold max|.| over the
    // shaking, and they are the ELASTIC branch (no Coulomb cap, no sigma_n0) because that is what the
    // linear dynamic system solves. any_slip is therefore always false here -- the dynamic interface
    // does not slip by construction. NEW FIELDS GO AT THE END (results_io writes positionally).
    bool envelope = false;
    bool superposed = false;        // the parent's static state was added -> stations are TOTAL action
    // Coulomb demand/capacity |tau| / tau_max along the joint (only when superposed -- tau_max needs
    // the TOTAL normal stress). The elastic dynamic joint cannot slip, so this is the only way the run
    // can reveal that the real joint WOULD have: > 1 means the linearisation is violated there.
    // Reported as a MEASURE, not a verdict: near the ground surface sigma_n -> 0, so tau_max -> c_i and
    // a short zone at the top exceeds capacity under almost any real shaking. That is genuine local
    // slip, not grounds for condemning the whole analysis -- so give the engineer both the worst value
    // and HOW MUCH of the joint is over, and let them judge. Capped at kUtilCap where tau_max = 0 (the
    // joint is in tension = separated, so it carries no shear at all and the ratio is unbounded).
    double max_utilisation = 0.0;
    double over_fraction = 0.0;     // fraction of the joint's stations with utilisation > 1
    // v5: this envelope came from the NONLINEAR (Coulomb) dynamic branch -- the joint's slip WAS
    // checked instant by instant, so any_slip is a real finding either way ([SLIPPING] or [bonded]).
    // False on a linear envelope: no slip check exists there at all, and the UI must say
    // "[elastic, no slip check]" instead of claiming "bonded" (a claim the analysis never tested).
    bool slip_checked = false;
};

// Raw converged structural state of a STATIC-family phase, for the phase chain (Track 1a): a
// NONLINEAR Dynamic child seeds its structural elements from this (displacement datum incl. the
// extra structural DOFs + the committed plastic state), so its Coulomb / yield / slack caps act on
// the TOTAL (static + dynamic) action instead of a zero-preload increment -- slip and yield ONSET
// are then correct where a static preload exists. In-memory chain carrier only: results_io does NOT
// persist it (a result reloaded from disk simply cannot be carried from; the dynamic phase then
// falls back to the increment-from-zero path and says so in its message).
struct StructCarryState {
    Eigen::VectorXd full_disp;   // total_dofs, incl. extra structural DOFs (NewtonResult.displacement)
    std::vector<double> anchor_plastic, geogrid_plastic, interface_slip, interface5_slip,
                        embedded_skin_slip, embedded_foot_slip,   // committed plastic state (solver order)
                        plate_plastic, plate5_plastic;            // plate M-N hinge state ([eps_p,kap_p]xGauss)
};

struct SolveResult {
    katai::mesh::Mesh mesh;                     // the mesh actually solved (split for embedded walls)
    Eigen::VectorXd disp;                      // full displacement (2 * node_count)
    NodalStressField stress;                   // nodal EFFECTIVE [sxx, syy, sxy] (pore on rhs)
    std::vector<double> pore;                  // nodal pore pressure u >= 0 (hydrostatic), size node_count
    bool ok = false;
    std::string message;
    double max_disp = 0.0;
    // Fraction of the applied load that reached equilibrium. 1.0 on full convergence; on a failed
    // (collapse) solve it is the highest equilibrated level = the incremental limit load (PLAXIS
    // SumMstage). Lets the GUI report "equilibrated X% of the load" and enables limit-load checks.
    double load_factor = 1.0;
    // Factor of safety from a Safety (phi-c reduction) analysis; < 0 when not a Safety run.
    double fos = -1.0;
    // True when no failure mechanism developed up to the srf_max cap: `fos` is then only a LOWER
    // BOUND (e.g. "FoS > 3.0"), not a definitive factor of safety. The GUI must say ">" not "=".
    bool fos_lower_bound = false;
    // Design approach (EC7 / TBDY 2018) echoed from the phase, for the report + verdict, as the
    // ENGINE's enum -- the schema's value is mapped once at the driver seam. When a material-factored
    // approach (DA1-C2, DA3) is set, `fos` on a Safety run is the OVER-DESIGN FACTOR (ODF >= 1.0 =>
    // the design satisfies the ULS), not the characteristic factor of safety.
    DesignApproach design_approach = DesignApproach::None;
    // True when the K0 phase had to resolve a genuine geostatic imbalance (non-level ground
    // surface / layers / water table -- PLAXIS "plastic nil-step"). The displacement field is
    // then the equilibrium redistribution, not zero.
    bool nil_step = false;
    // Structural force diagrams (one per drawn structural line; empty when no structures).
    std::vector<StructForce> struct_forces;
    // Interface results (one per interface line / wall interface; empty when no interfaces). tau /
    // sigma_n / slip along each Coulomb joint -- shows what the soil-structure interface is doing.
    std::vector<InterfaceResult> interface_forces;
    // Support reactions [kN/m]: the nodal internal force at the FIXED dofs (zero at free
    // dofs), from the committed Gauss stresses -- what each support or imposed displacement
    // exerts. Layout node*2 + component, like `disp`. SOIL contribution (a structural end
    // force landing on a fixed node is not included in v1). Static-family phases only;
    // empty otherwise. In-memory result (not in the .res file yet).
    Eigen::VectorXd reaction;
    // Staged construction: element activity of this phase (empty = everything active). The GUI
    // dims/hides the excavated region; passive elements are excluded from the stress recovery.
    std::vector<char> active;
    // Calculation time breakdown + Newton iteration count of the main solve (PLAXIS-style
    // calculation report). Zero for paths that aggregate many solves (Safety bisection).
    NewtonResult::Timings timings;
    int iterations = 0;
    // Consolidation (time-dependent) phase: settlement-time curve for the PLAXIS U-t plot. Empty
    // for non-consolidation results. consol_time[k] is the time [day] and consol_settlement[k] the
    // maximum vertical settlement |u_y| [m] of step k (consol_excess_pore[k] the max excess pore).
    std::vector<double> consol_time;
    std::vector<double> consol_settlement;
    std::vector<double> consol_excess_pore;
    // Final nodal degree of saturation S (size node_count) for transient/fully-coupled unsaturated
    // flow; empty otherwise. 1.0 in the saturated zone, < 1 where suction develops (van Genuchten).
    std::vector<double> saturation;
    // Dynamic (seismic) phase results (empty otherwise): the surface (monitor) TOTAL horizontal
    // acceleration time-history a_x(t) [m/s^2], its 5%-damped elastic response spectrum, and the TBDY
    // 2018 design spectrum -- both over dyn_period [s], in m/s^2 -- for the code-comparison overlay.
    std::vector<double> dyn_time, dyn_surface_ax;
    std::vector<double> dyn_period, dyn_response_sa, dyn_design_sa;
    // EC8 (EN 1998-1) elastic design spectrum [m/s^2] over dyn_period, when the phase enabled the
    // overlay (empty otherwise). Recommended EN parameters; national annexes may differ (stated).
    std::vector<double> dyn_design_sa_ec8;
    double dyn_peak_surface_a = 0.0;   // peak |surface total a_x| (= response-spectrum PGA) [m/s^2]
    // Model fundamental frequency f1 [Hz] (small-strain elastic K, fixed base; inverse-power
    // iteration on (K, M)) -- computed for every Dynamic phase so the engineer never has to
    // ESTIMATE it by the quarter-wave travel-time hand rule, which mis-places f1 by 20%+ on a
    // strong impedance-contrast profile (measured: 2.34 vs the true 3.005 Hz -- see
    // docs/validation/site-response-benchmark.md). Drives the Rayleigh-band note in the message.
    // 0 = could not be computed (the run proceeds; the note is simply absent).
    double dyn_model_f1 = 0.0;
    // Raw structural state for the phase chain (Track 1a; see StructCarryState). Set by the
    // static-family tail; empty for Dynamic / Safety / consolidation results (no carriable state).
    StructCarryState struct_state;
    // Everything the run did differently from what the file literally says (see Diagnostic):
    // clipped geometry, an ignored parameter, a fallback taken. Empty on a clean run. Front
    // ends print this list; they do not interpret it. NEW FIELDS GO AT THE END (results_io
    // writes positionally) -- diagnostics belong to the run, so they are not persisted.
    std::vector<Diagnostic> diagnostics;
};

} // namespace katai::core
