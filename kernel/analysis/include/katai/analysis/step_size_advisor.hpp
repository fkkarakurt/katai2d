#pragma once
// Time-step adequacy advisors (Stage B7: extracted from the application
// driver). Both return a complete user-facing warning string, or "" when the
// step is adequate -- warn in the editor, never silently mis-integrate. The
// engine takes neutral inputs; whether a phase is Dynamic or Consolidation,
// and how a schema material resolves to an oedometer modulus, is decided once
// at the caller's seam.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/constants.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Neutral description of a dynamic phase's time discretisation.
struct DynamicStepInput {
    double duration = 1.0;     // [s]; <= 0 falls back to 1 s, as the driver always did
    int time_steps = 1;        // clamped to [1, 20000]
    bool record = false;       // accelerogram input (the record's own dt is the criterion)
    int record_samples = 0;    // loaded samples when record
    double record_dt = 0.0;    // record sampling [s]; <= 0 falls back to 0.02
    bool ricker = false;       // Ricker pulse: energy reaches ~2.5x the central frequency
    double frequency = 0.0;    // central/harmonic frequency [Hz]; <= 0 falls back to 2
    bool nonlinear = false;    // nonlinear dynamic: the step also controls plasticity + Newton
};

// Time-step adequacy for a Dynamic (seismic) phase. The base motion a_g(t) is SAMPLED pointwise at
// t = k*dt with no anti-aliasing: too few steps per cycle and the solver is driven by a signal that is
// not the motion the user asked for. Newmark (gamma=1/2, beta=1/4) is unconditionally stable, so
// nothing diverges -- the answer is simply wrong and looks normal. Standard practice (and the phase
// editor's own guidance) is >= ~20 steps per period of the highest frequency of interest; below 2
// per period the signal is aliased outright (Nyquist). Returns "" when the step is adequate.
inline std::string dynamic_step_warning(const DynamicStepInput& in) {
    const double dur = in.duration > 0.0 ? in.duration : 1.0;
    const int nst = std::clamp(in.time_steps, 1, 20000);
    // Accelerogram input: the meaningful step criterion is the RECORD's own sampling -- a solver
    // step coarser than record_dt skips samples (peaks are simply never seen by the interpolant).
    if (in.record) {
        if (in.record_samples < 2)
            return "Accelerogram input is selected but no record is loaded -- import one "
                   "(two-column t/a or one acceleration per line), or pick another waveform.";
        const double rdt = in.record_dt > 0.0 ? in.record_dt : 0.02;
        const double sdt = dur / nst;
        if (sdt > rdt * 1.0001) {
            char rb[300];
            std::snprintf(rb, sizeof(rb),
                          "Time step (%.4g s) is coarser than the record's sampling (%.4g s): the "
                          "record is UNDERSAMPLED and its peaks are skipped. Use at least %d time "
                          "steps (dt <= record dt).%s",
                          sdt, rdt, (int)std::ceil(dur / rdt),
                          in.nonlinear
                              ? " Nonlinear dynamic: the step also controls plasticity accuracy and "
                                "Newton convergence."
                              : "");
            return rb;
        }
        return "";
    }
    const double f = in.frequency > 0.0 ? in.frequency : 2.0;
    const double f_max = in.ricker ? 2.5 * f : f;
    const double per_cycle = (double)nst / (dur * f_max);   // steps per period of f_max
    if (per_cycle >= 20.0) return "";
    // Nonlinear path: the step ALSO controls plasticity accuracy and per-step Newton convergence,
    // not just input-signal sampling -- say so (no invented second threshold: a measured nonlinear
    // step criterion needs a refinement study, tracked; a non-converged step already stops the run
    // and is reported honestly in the phase message).
    const char* nl_note = in.nonlinear
        ? " Nonlinear dynamic: the step also controls plasticity accuracy and Newton convergence --"
          " prefer finer steps than this minimum."
        : "";
    char buf[320];
    if (per_cycle < 2.0)
        std::snprintf(buf, sizeof(buf),
                      "Time step too coarse: %.1f steps per cycle at %.3g Hz -- BELOW the Nyquist limit "
                      "(2). The base motion is aliased: the solver is driven by a different signal than "
                      "the one requested. Use at least %d time steps (>= 20 per cycle).",
                      per_cycle, f_max, (int)std::ceil(20.0 * dur * f_max));
    else
        std::snprintf(buf, sizeof(buf),
                      "Time step coarse: %.1f steps per cycle at %.3g Hz. Use >= 20 per cycle so the "
                      "base motion and the response are resolved -- suggest at least %d time steps.",
                      per_cycle, f_max, (int)std::ceil(20.0 * dur * f_max));
    return std::string(buf) + nl_note;
}

// One material's contribution to the Vermeer-Verruijt critical step, already resolved to physical
// quantities: the caller's seam decides how a constitutive model yields its oedometer modulus
// (constant-E closed form, HS Eoed_ref, or the Soft Soil ln-law tangent at a reference stress).
struct ConsolidationStepMaterial {
    double Eoed = 0.0;      // oedometer modulus [kPa]
    double k_y = 0.0;       // vertical permeability [m/day]
    double porosity = 0.0;  // n = e / (1 + e)
};

// Vermeer & Verruijt (1981) first-time-step guidance for a time-dependent phase: the equal step
// dt = duration/steps must not fall far below dt_crit = h^2*gamma_w/(eta*k)*(1/Eoed + n/Kw)
// (docs/references/consolidation-formulation.md sec. 4). Below it the backward-Euler diffusion dt*H
// cannot damp the element-local pore mode, so the equal-order tri6 (P2-P2) pore field can show a
// checkerboard at the near-undrained instant (docs/validation/lbb-undrained-checkerboard.md).
// Settlement / consolidation curves stay reliable regardless -- this warns only about the EARLY-TIME
// PORE field. h = mean element size; eta = 40 (tri6) / 80 (tri15); the most-restrictive dt_crit over
// the materials actually present in the mesh binds. `materials` is indexed by material id, exactly
// like the mesh's element_material. Returns "" when the step is adequate.
inline std::string consolidation_step_warning(const std::vector<ConsolidationStepMaterial>& materials,
                                              const katai::mesh::Mesh& mesh,
                                              double duration_day, int time_steps) {
    if (mesh.element_count == 0 || time_steps < 1 || duration_day <= 0.0) return {};

    // Mean element size h (area-based, like the flow solver's free-surface transition width).
    double area_sum = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        area_sum += 0.5 * std::fabs((mesh.x[b] - mesh.x[a]) * (mesh.y[c] - mesh.y[a]) -
                                    (mesh.x[c] - mesh.x[a]) * (mesh.y[b] - mesh.y[a]));
    }
    const double h = std::sqrt(2.0 * area_sum / std::max(1, mesh.element_count));
    const bool tri15 = mesh.nodes_per_element == 15;
    const double eta = tri15 ? 80.0 : 40.0;
    constexpr double kWaterBulk = 2.0e6;   // water bulk modulus [kPa] (matches the consolidation solve)

    std::vector<char> used(materials.size(), 0);
    for (int e = 0; e < mesh.element_count; ++e) {
        const int m = mesh.element_material[e];
        if (m >= 0 && m < (int)materials.size()) used[m] = 1;
    }
    double dt_crit = 0.0;
    for (size_t mi = 0; mi < materials.size(); ++mi) {
        if (!used[mi]) continue;
        const ConsolidationStepMaterial& M = materials[mi];
        dt_crit = std::max(dt_crit, consolidation_critical_dt(
                                        h, eta, M.Eoed, M.k_y, M.porosity, kWaterBulk, kGammaWater));
    }
    const double dt = duration_day / time_steps;
    if (dt_crit <= 0.0 || dt >= dt_crit) return {};

    char buf[440];
    std::snprintf(buf, sizeof(buf),
        "First time step dt = %.3g day is below the Vermeer-Verruijt critical step dt_crit = %.3g day "
        "(%s elements). Settlements stay reliable, but the early-time pore-pressure field can oscillate "
        "(checkerboard). Use fewer / larger steps to raise dt above dt_crit%s.",
        dt, dt_crit, tri15 ? "15-noded" : "6-noded",
        tri15 ? "" : ", or switch to 15-noded elements (Mesh tab)");
    return buf;
}

} // namespace katai::core
