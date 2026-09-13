#pragma once
// Consolidation phase strategy (Stage B9). Time-dependent Biot consolidation,
// the "Consolidation" phase: the configuration's load increment dF = f - B (the
// staged-construction imbalance) is applied at t = 0+, generating an undrained excess
// pore pressure that dissipates over the phase's time interval -- the classic
// settlement-time (Terzaghi U-t) development. v1: soil-only.
//
// Same contract as the transient-flow pilot: neutral inputs resolved at the
// caller's seam, refusal messages engine-owned and byte-identical, material
// names carried as diagnostic labels. Two differences worth stating:
//   - the linear solver enters as factories built by the composition root
//     (one factorization, back-substitution every step), never named here;
//   - on success this strategy does NOT set R.ok or R.message -- unlike the
//     flow-only pilot, a consolidation phase falls through to the driver's
//     common result tail, exactly as it always has.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/consolidation.hpp>
#include <katai/analysis/constants.hpp>
#include <katai/analysis/hydraulic_boundary.hpp>
#include <katai/analysis/post/consolidation_recovery.hpp>
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/analysis/results.hpp>
#include <katai/analysis/phase_solver/coupled_structures.hpp>  // the shared structural report
#include <katai/analysis/structural_dynamics.hpp>   // assemble_structural_stiffness (the shared assembler)
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// One material's consolidation description, resolved at the caller's seam.
// Constitutive behaviour (incl. the Soft Soil Creep detection) comes from the
// engine's own `models` table, not from here.
struct ConsolidationPhaseMaterial {
    std::string name;              // diagnostic label for honest refusals
    double kx = 1.0, ky = 1.0;     // permeability [m/day]
    double porosity = 0.3;         // n = e / (1 + e)
    bool nonporous = false;        // refused: a non-porous region would silently behave water-filled
    // Oedometer modulus E_oed [kPa] of the material, resolved at the caller's seam (how a schema
    // model yields one is schema knowledge: HS uses Eoed_ref, Soft Soil the NC tangent, the rest
    // the closed form). Only the AUTOMATIC first time step reads it -- it is the same dt_crit the
    // editor already warns below. 0 = unknown, which makes the automatic step refuse rather than
    // guess.
    double eoed = 0.0;
    // Undrained (C) is a TOTAL stress material: it has no pore pressures to consolidate, so a
    // consolidation calculation has nothing to act on in it. This build refuses the phase rather
    // than solving part of the mesh in total stress and the rest in effective stress with one
    // Kw/n between them.
    bool total_stress = false;
};

// The phase's neutral configuration.
struct ConsolidationPhase {
    std::vector<ConsolidationPhaseMaterial> materials;  // by material id
    std::vector<FlowEdge> flow_edges;                   // B4 vocabulary, Closed edges included
    bool have_flow_bcs = false;                         // unfiltered declaration scan (B4 rule)
    std::vector<char> active;                           // element activity; empty = everything active
    // Nodes of the DRAINS active in this phase. In a consolidation phase a drain sets the excess
    // pore pressure to zero and its specified head is not used (docs/k2d-format.md, `hydros`) --
    // which is exactly what this solver's drained-node mask means, so a drain enters here rather
    // than as a new kind of boundary condition.
    std::vector<int> drain_nodes;
    // Structural elements that sit in this phase's coupled solve. Plates, anchors and geogrids
    // enter it; interfaces and embedded beams are refused below, and the refusal says why.
    bool has_embedded_beams = false;     // embedded beam rows (skin + foot springs)
    // An interface SPLITS the mesh, so the joint has two node sets at one place and what the water
    // does between them is an INPUT (docs/k2d-format.md, `flow_barrier`): fully permeable = one
    // pressure (the two share a pore equation, `pore_tie` below), impermeable = two. The third,
    // semi-permeable, is a CONDUCTANCE dh/R between them and is not expressible as either -- it is
    // refused rather than rounded to whichever neighbour looks closer.
    bool has_semi_permeable_interface = false;
    const std::vector<int>* pore_tie = nullptr;              // node -> the node it shares a pore eq with
    const std::vector<IfaceDiag>* iface_diagrams = nullptr;  // joints to report tau / sigma_n / slip for
    // The structural elements themselves, the lines to report forces for, and the parent phase's
    // total displacement. The last one matters because structural elements are TOTAL-displacement
    // formulated: a wall installed in an earlier phase carries its force at the parent's datum,
    // and a consolidation phase solves for the INCREMENT from there. Reporting the increment as
    // the wall's moment would understate it by exactly the parent's share.
    const Structures* structures = nullptr;
    const std::vector<DiagSpec>* diagrams = nullptr;
    const Eigen::VectorXd* carry_full = nullptr;
    double duration_day = 1.0;                          // <= 0 falls back to 1 day
    int time_steps = 25;                                // clamped to [1, 2000]
    double yscale = 1.0;                                // model height, for the top-drain tolerance
    // How the phase ENDS (katai/analysis/results.hpp ConsolidationStop; docs/k2d-format.md,
    // `cstop`). TimeInterval uses duration_day / time_steps above and is the historical path,
    // untouched. The other two ignore both -- a time interval does not apply to a phase that ends
    // on a state -- and march until the state is reached, reporting the time it took.
    ConsolidationStop stop = ConsolidationStop::TimeInterval;
    double stop_min_pore = 1.0;   // [kPa] MinExcessPore threshold on |p|max (default 1 kPa)
    double stop_degree = 0.9;     // [-]  DegreeOfConsolidation target, as a PRESSURE ratio (default 90%)
    double first_dt = 0.0;        // [day] first time step; <= 0 = automatic (Vermeer-Verruijt dt_crit)
    int max_steps = 1000;         // safety cap on the march; reaching it REFUSES, it does not report
};

// Solve the phase. Fills the settlement-time series, the final displacement and
// the recovered nodal effective stresses in R, and returns the final committed
// Gauss state in `committed_out` (the phase chain carries it). Returns false on
// an honest refusal or non-convergence, with R.message set; on success the
// caller's common result tail completes R.
inline bool solve_consolidation_phase(
    const katai::mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<MaterialModel>& models, const std::vector<MaterialProfile>& profiles,
    const std::vector<GaussState>& init, const Eigen::VectorXd& dF, bool nonlinear_soil,
    const ConsolidationPhase& in,
    const ConsolidationSolveFactory& le_factory, const ConsolidationSolveFactory& plastic_factory,
    SolveResult& R, std::vector<GaussState>& committed_out) {
    // WHAT A STRUCTURE IS IN A COUPLED PHASE. Plates, anchors and geogrids add stiffness to
    // mechanical degrees of freedom the soil already owns. An interface adds that too, and one
    // thing more: it splits the mesh, so the joint carries two pore pressures where the ground
    // carried one, and which of them the water sees is the interface's own cross-permeability
    // input -- fully permeable (the default: the two share a pore equation, so continuity and the
    // flux balance across the joint hold by construction) or impermeable (two separate pressures,
    // which is what the bare split gives). Both are values of that input.
    //
    // Semi-permeable is the third, and it is not a blend of the two: it is a conductance, the
    // joint passing q_n = dh / R per unit area, which needs a term in H that this solver does not
    // have. Refusing it costs the model that asks for it; rounding it to permeable or impermeable
    // would cost every model that asks for it, silently, and in whichever direction the rounding
    // happened to go.
    if (in.has_semi_permeable_interface) {
        R.message = "A semi-permeable interface (cross permeability with a hydraulic resistance) "
                    "is in this consolidation phase. That is a conductance across the joint -- it "
                    "passes q = dh / R per unit area -- and this solver carries only the two ends "
                    "of that scale, fully permeable and impermeable. Set the interface to one of "
                    "those for a consolidation phase, or run the stage as a Plastic phase with a "
                    "groundwater-flow field, which does read the resistance.";
        return false;
    }
    if (in.has_embedded_beams) {
        R.message = "This consolidation phase contains an embedded beam (pile row). Its skin and "
                    "foot resistance follow the EFFECTIVE stress around the pile, which changes "
                    "throughout consolidation, so the elastic springs this phase would use are "
                    "not the pile's behaviour here. Use a Plastic phase for the pile, or remove "
                    "it from the consolidating stage.";
        return false;
    }
    const size_t nmat = in.materials.size();
    // Soft Soil Creep x consolidation: honest refusal. End-to-end verification WAS attempted and
    // failed (2026-07-21, test_soft_soil_gui): on a fast-draining column the settlement DEPENDS on
    // the time-step size and the creep tail came out 37% short of the exact mu* H ln(100) oracle
    // that the verified Plastic-phase path holds to +0.1%. Until the Biot x creep interaction is
    // resolved: apply loads in a timed Plastic phase, model holding as a timed no-change Plastic
    // phase (both pinned against the isotach closed forms).
    for (size_t mi = 0; mi < nmat; ++mi) {
        bool used_m = false;
        for (int e = 0; e < mesh.element_count && !used_m; ++e)
            if ((in.active.empty() || in.active[e]) && mesh.element_material[e] == (int)mi) used_m = true;
        if (used_m && models[mi].type == MaterialType::SoftSoilCreep) {
            R.message = "Material '" + in.materials[mi].name + "' is Soft Soil Creep: the "
                        "consolidation (Biot) x creep interaction is not verified yet -- "
                        "measured settlements DEPEND on the time-step size on this path, so "
                        "it is refused rather than shipped silently wrong. Apply loads in a "
                        "Plastic phase WITH its real time interval, and model holding "
                        "periods as a chained no-change Plastic phase with the holding time "
                        "(both verified against the isotach closed forms).";
            return false;
        }
    }
    // Non-porous x consolidation: honest refusal. The Biot solver gives every active element a
    // pore-pressure DOF -- a non-porous region would silently carry pore pressure and pass water
    // (concrete gaining pore stiffness, drainage paths running through it). The fluid stiffness is
    // per material now, but a stiffness is not an exclusion: region-based pore-DOF exclusion is
    // the separate work item this refusal waits for.
    for (size_t mi = 0; mi < nmat; ++mi) {
        bool used_np = false;
        for (int e = 0; e < mesh.element_count && !used_np; ++e)
            if ((in.active.empty() || in.active[e]) && mesh.element_material[e] == (int)mi) used_np = true;
        if (used_np && in.materials[mi].nonporous) {
            R.message = "Material '" + in.materials[mi].name + "' is Non-porous: consolidation "
                        "/ fully-coupled phases give every element a pore-pressure DOF, so a "
                        "non-porous region would silently carry pore pressure and pass water. "
                        "Model the concrete with Drained + high stiffness in "
                        "consolidation phases, or keep Non-porous to static/dynamic phases.";
            return false;
        }
        if (used_np && in.materials[mi].total_stress) {
            R.message = "Material '" + in.materials[mi].name + "' is Undrained (C), a total "
                        "stress analysis: it carries no pore pressure, so there is nothing in it "
                        "to consolidate, and a consolidation calculation has nothing to act on in "
                        "it. Solving the phase would put part of the mesh "
                        "in total stress and the rest in effective stress. Give the material "
                        "effective parameters with Drained or Undrained (A)/(B) for the "
                        "consolidating phases.";
            return false;
        }
    }
    // Permeability per material [m/day]; every used material needs k > 0.
    std::vector<Permeability> cperm(nmat, {1.0, 1.0});
    std::vector<char> used(nmat, 0);
    for (int e = 0; e < mesh.element_count; ++e)
        if (in.active.empty() || in.active[e]) used[mesh.element_material[e]] = 1;
    for (size_t mi = 0; mi < nmat; ++mi) {
        cperm[mi] = {in.materials[mi].kx, in.materials[mi].ky};
        if (used[mi] && (cperm[mi].kx <= 0.0 || cperm[mi].ky <= 0.0)) {
            R.message = "Material '" + in.materials[mi].name + "' has no permeability -- set "
                        "kx, ky (> 0) on its Groundwater tab before a consolidation phase.";
            return false;
        }
    }
    // Pore-fluid stiffness Kw/n from the real water bulk modulus (Verruijt), one value per
    // material from that material's own porosity (PoreFluidStiffness, consolidation.hpp): a
    // layered model is answered the same whatever order its materials are listed in.
    constexpr double kWaterBulk = 2.0e6;   // bulk modulus of water [kPa]
    std::vector<double> kw_by_material(nmat, kWaterBulk / 0.3);
    for (size_t mi = 0; mi < nmat; ++mi)
        kw_by_material[mi] = kWaterBulk / std::max(0.05, in.materials[mi].porosity);
    const PoreFluidStiffness kw_over_n(std::move(kw_by_material), kWaterBulk / 0.3);

    // Drainage boundary (engine service, B4): prescribed-head / seepage edges drain; with no
    // declared flow BCs the model top drains; inactive-only nodes carry no pore DOF.
    // (not const: an active drain adds its own nodes to the same mask -- the drain and the
    // draining boundary say the same thing to this solver, that the excess pore pressure there
    // is zero.)
    std::vector<char> drained =
        flow_drained_nodes(in.flow_edges, in.have_flow_bcs, mesh, in.active, in.yscale);
    for (int n : in.drain_nodes)
        if (n >= 0 && n < mesh.node_count) drained[n] = 1;
    // Two nodes that share a pore pressure share its boundary condition too: if either side of a
    // permeable joint drains, the joint drains. Without this the mask could say "drained" on one
    // side and "free" on the other of the SAME pressure, and the numbering would have to pick one
    // -- silently, and differently depending on which side the splitter numbered first.
    if (in.pore_tie)
        for (int n = 0; n < mesh.node_count; ++n) {
            const int rep = (*in.pore_tie)[n];
            if (rep >= 0 && rep < mesh.node_count && (drained[n] || drained[rep]))
                drained[n] = drained[rep] = 1;
        }

    // The structural elements' elastic stiffness, assembled ONCE with the shared assembler -- the
    // same element matrices and DOF mapping solve_nonlinear uses, which is what makes "the wall in
    // this phase is the wall in the Plastic phase" a fact rather than a hope. It enters the coupled
    // system as a matrix (katai/analysis/consolidation.hpp), so the core stays structure-agnostic.
    math::CsrMatrix struct_k_storage;
    const math::CsrMatrix* struct_k = nullptr;
    if (in.structures && (!in.structures->plates.empty() || !in.structures->plates5.empty() ||
                          !in.structures->anchors.empty() || !in.structures->geogrids.empty() ||
                          !in.structures->interfaces.empty() || !in.structures->interfaces5.empty())) {
        math::SparseMatrixBuilder kb(dofs.equation_count());
        assemble_structural_stiffness(mesh, dofs, *in.structures, kb);
        struct_k_storage = kb.build();
        struct_k = &struct_k_storage;
    }

    // Time stepping. The phase's own interval (the historical path) OR a state criterion below.
    const double duration = in.duration_day > 0.0 ? in.duration_day : 1.0;
    const int nsteps = std::clamp(in.time_steps, 1, 2000);
    const double dt = duration / nsteps;

    // One run of the core from the CURRENT state: n steps of dt, with the phase's load increment
    // applied in the first of them when `with_load`. Both cores are restartable BY CONSTRUCTION --
    // the linear-elastic skeleton is time-invariant, and the elastoplastic one rebuilds its
    // baseline internal force from the committed Gauss state it is handed -- so a march with a
    // CHANGING dt is a SEQUENCE of fixed-dt runs, and the verified fixed-dt core is not touched at
    // all. That equality is asserted, not assumed: test_consolidation_stop pins one run of 120
    // steps against four runs of 30.
    struct Chunk {
        ConsolidationResult series;            // times RELATIVE to the chunk start; [0] = its start
        std::vector<GaussState> committed;     // elastoplastic path: state at the chunk END only
        bool ok = false;
    };
    std::vector<double> p_state(mesh.node_count, 0.0);   // excess pore at the chunk start (0 at t=0)
    std::vector<GaussState> state = init;                // committed effective Gauss state
    auto run_chunk = [&](double step_dt, int n, bool with_load) {
        Chunk c;
        const Eigen::VectorXd* load = with_load ? &dF : nullptr;
        if (nonlinear_soil) {
            // Elastoplastic skeleton: effective stress from the constitutive return mapping; each
            // time step is a monolithic coupled Newton. K_T is nonsymmetric for non-associated
            // flow -> the caller's plastic factory selects the matching solver. The previous
            // phase's committed effective stresses (`init`) seed the state; dF drives it.
            ConsolidationPlasticResult r = solve_consolidation_plastic(
                mesh, dofs, models, cperm, kGammaWater, kw_over_n, drained, state, p_state,
                step_dt, n, in.active, load, plastic_factory, 40, 1e-6, profiles, struct_k,
                in.pore_tie);
            c.ok = r.converged;
            c.series = std::move(r.series);
            c.committed = std::move(r.committed);
        } else {
            c.series = solve_consolidation(mesh, dofs, models, cperm, kGammaWater, kw_over_n,
                                           drained, p_state, step_dt, n, in.active, load,
                                           le_factory, profiles, struct_k, in.pore_tie);
            c.ok = !c.series.displacement.empty();
        }
        return c;
    };
    const char* const kFailedMessage =
        nonlinear_soil
            ? "Elastoplastic consolidation did not converge in a time step (the load increment may "
              "exceed the soil capacity, or the time step is too large -- use a smaller load "
              "increment or more, smaller steps)."
            : "Consolidation produced no result.";

    // The phase's own settlement-time curve, in ABSOLUTE phase time: max vertical settlement |u_y|
    // and max |excess pore| per recorded step, plus the running PEAK excess pore pressure -- the
    // reference the degree of consolidation is a ratio of.
    double pore_peak = 0.0;
    Eigen::VectorXd v_total = Eigen::VectorXd::Zero(dofs.total_dofs());
    double t_now = 0.0;
    int steps_done = 0;
    bool load_pending = true;
    R.consol_time.push_back(0.0);          // the phase's own start: no displacement, no excess pore
    R.consol_settlement.push_back(0.0);
    R.consol_excess_pore.push_back(0.0);

    auto step_pmax = [&](const Chunk& c, int s) {
        double pmax = 0.0;
        for (int n = 0; n < mesh.node_count; ++n)
            pmax = std::fmax(pmax, std::fabs(c.series.pore[s](n)));
        return pmax;
    };
    // Append EVERY step of a chunk and adopt its end state. Only ever called with a whole chunk:
    // the committed Gauss state exists at its end and nowhere else, so "stop half way through" is
    // spelled as a shorter RE-RUN (identical arithmetic, same state), never as a partial commit.
    auto commit = [&](const Chunk& c) {
        const int n = (int)c.series.times.size() - 1;
        for (int s = 1; s <= n; ++s) {
            double smax = 0.0;
            for (int nd = 0; nd < mesh.node_count; ++nd) {
                const int d = dofs.global_dof(nd, 1);
                smax = std::fmax(smax, std::fabs(v_total(d) + c.series.displacement[s](d)));
            }
            const double pmax = step_pmax(c, s);
            R.consol_time.push_back(t_now + c.series.times[s]);
            R.consol_settlement.push_back(smax);
            R.consol_excess_pore.push_back(pmax);
            pore_peak = std::fmax(pore_peak, pmax);
        }
        v_total += c.series.displacement[n];
        for (int nd = 0; nd < mesh.node_count; ++nd) p_state[nd] = c.series.pore[n](nd);
        R.excess_pore = p_state;   // the field itself, not only its maximum (results.hpp)
        if (nonlinear_soil) state = c.committed;
        t_now += c.series.times[n];
        steps_done += n;
        load_pending = false;
    };

    if (in.stop == ConsolidationStop::TimeInterval) {
        Chunk c = run_chunk(dt, nsteps, true);
        if (!c.ok) { R.message = kFailedMessage; return false; }
        commit(c);
        R.consol_stop_met = true;   // a time interval always ends: it is a duration, not a target
    } else {
        // --- Ending on a STATE ---------------------------------------------------------------
        // No duration is given, so the step size cannot come from one. The march starts from the
        // Vermeer-Verruijt critical step (docs/references/consolidation-formulation.md sec. 4) and
        // then DOUBLES every kStepsPerLevel steps, which holds dt/t ~ 1/kStepsPerLevel for the
        // whole run: a log-time march, the shape dissipation actually has and the shape its own
        // curve is read in. Both constants below were MEASURED on the 1-D column whose answer is
        // known in closed form (KV-CON-003); neither was chosen for looking reasonable.
        //
        //   kStepsPerLevel -- the reported stop time is first-order in 1/N, and it is the ANSWER
        //   of this phase, not a detail of its curve. Measured (stop at 1 kPa of the 10 kPa
        //   generated; closed form 14.566 day, a fine equal-step grid gives 14.60):
        //       N =   8 -> 15.367 (+5.5%)     N =  32 -> 14.778 (+1.5%)    N = 128 -> 14.624 (+0.4%)
        //       N =  16 -> 15.014 (+3.1%)     N =  64 -> 14.675 (+0.8%)
        //   The halving is clean. 32 buys 1.5% for 203 steps where 8 costs 5.5% for 67 -- and the
        //   steps INSIDE a level are back-substitutions on one factorisation (dt is constant
        //   there), so N costs far less than its step count suggests: the factorisations are the
        //   levels, and there are ~log2(t_end/dt_0) of those whatever N is.
        constexpr int kStepsPerLevel = 32;    // steps per doubling of dt
        constexpr int kRefineSubsteps = 16;   // pieces the crossing step is re-run in
        const int max_steps = std::clamp(in.max_steps, 1, 20000);

        // Nothing can dissipate without somewhere to dissipate TO. A criterion that can never be
        // met would otherwise be answered by exhausting the step budget on a model whose defect
        // is structural, and the report would blame the budget.
        bool any_drained = false;
        for (char d : drained) if (d) { any_drained = true; break; }
        if (!any_drained) {
            R.message = "This consolidation phase is asked to run until the excess pore pressure "
                        "reaches a target, but the model has no drainage boundary: every boundary "
                        "is closed, so the excess pore pressure redistributes and never leaves. "
                        "Give the model a Prescribed head / Seepage face edge (or a drain), or end "
                        "the phase on a time interval.";
            return false;
        }

        // The automatic first time step. dt_crit is a STABILITY bound on the pore field, not an
        // accuracy bound on the pressure the ratio is measured against, and the difference shows:
        // at dt_crit exactly, the first step of the column overshoots the undrained pressure it
        // generates by 12.8% (11.275 kPa for a 10 kPa surcharge) -- the drainage-boundary layer is
        // thinner the smaller the step, and an equal-order pore field cannot resolve it. That peak
        // IS the denominator of the degree of consolidation, so the overshoot would land in the
        // answer. Measured on the same column against the exact 9.9983 kPa: 1x dt_crit -> 11.275,
        // 2x -> 10.034, 4x -> 9.9983, 64x -> 9.9980. Four is where it is clean, and dt_crit is a
        // lower bound, so multiplying it is allowed.
        constexpr double kFirstStepOverCritical = 4.0;
        double dt_level = in.first_dt;
        if (dt_level <= 0.0) {
            const double h = mean_element_size(mesh);
            const double eta = consolidation_eta(mesh);
            for (size_t mi = 0; mi < nmat; ++mi)
                if (used[mi])
                    dt_level = std::fmax(dt_level,
                                         consolidation_critical_dt(h, eta, in.materials[mi].eoed,
                                                                   in.materials[mi].ky,
                                                                   in.materials[mi].porosity,
                                                                   kWaterBulk, kGammaWater));
            dt_level *= kFirstStepOverCritical;
            if (!(dt_level > 0.0)) {
                R.message = "This consolidation phase ends on a target rather than on a time "
                            "interval, so the solver must choose the first time step itself -- and "
                            "the Vermeer-Verruijt critical step cannot be formed here (a material "
                            "has no oedometer modulus: nu = 0.5, or a model whose stiffness is not "
                            "resolved). Enter a first time step for the phase, or end it on a time "
                            "interval.";
                return false;
            }
        }

        // Is the criterion met at this state? MinExcessPore is an ABSOLUTE threshold on |p| (it
        // applies to suction as much as to pressure, which is why the magnitude is taken).
        // DegreeOfConsolidation is defined as a PRESSURE ratio (see ConsolidationStop):
        // |p|max(t) <= (1 - U) * |p|max,initial.
        // A ratio needs a denominator that is a pressure. Where the staged change generated
        // nothing (a phase that activates no load, or one whose change is already in the ground),
        // the peak is round-off, and 10% of round-off is a target the march would chase until its
        // step budget ran out and then REFUSE -- a budget message for a modelling fact. Below the
        // floor the criterion is simply true, and the run says so (K2D-A013 at the tail).
        constexpr double kPoreFloor = 1e-9;   // [kPa] -- a micropascal is not an excess pore pressure
        auto criterion_met = [&](double pmax, double peak) {
            return in.stop == ConsolidationStop::MinExcessPore
                       ? pmax <= in.stop_min_pore
                       : pmax <= std::fmax((1.0 - in.stop_degree) * peak, kPoreFloor);
        };
        // First step of a chunk that meets it (1-based), or -1. The peak is carried forward: the
        // reference of the ratio is the largest excess pore pressure the phase has reached, which
        // for a load applied at t = 0+ is the undrained pressure that load generated.
        auto scan = [&](const Chunk& c) {
            double peak = pore_peak;
            const int n = (int)c.series.times.size() - 1;
            for (int s = 1; s <= n; ++s) {
                const double pmax = step_pmax(c, s);
                peak = std::fmax(peak, pmax);
                if (criterion_met(pmax, peak)) return s;
            }
            return -1;
        };

        while (!R.consol_stop_met && steps_done < max_steps) {
            const int n = std::min(kStepsPerLevel, max_steps - steps_done);
            Chunk c = run_chunk(dt_level, n, load_pending);
            if (!c.ok) { R.message = kFailedMessage; return false; }
            const int hit = scan(c);
            if (hit < 0) { commit(c); dt_level *= 2.0; continue; }
            if (hit > 1) {   // land exactly on the step BEFORE the crossing (same dt, same state)
                Chunk c2 = run_chunk(dt_level, hit - 1, load_pending);
                if (!c2.ok) { R.message = kFailedMessage; return false; }
                commit(c2);
            }
            // Localise the crossing INSIDE that one step: re-run it in kRefineSubsteps pieces and
            // end at the first that meets the criterion. Without this the answer would be the step
            // AFTER the crossing, which this late in a log-time march is ~12% of t away -- and the
            // answer of this phase IS a time.
            const double dt_fine = dt_level / kRefineSubsteps;
            Chunk cr = run_chunk(dt_fine, kRefineSubsteps, load_pending);
            if (!cr.ok) { R.message = kFailedMessage; return false; }
            const int m = scan(cr);
            if (m < 0) { commit(cr); continue; }   // the finer path has not crossed yet: march on
            if (m < kRefineSubsteps) {
                Chunk cf = run_chunk(dt_fine, m, load_pending);
                if (!cf.ok) { R.message = kFailedMessage; return false; }
                commit(cf);
            } else {
                commit(cr);
            }
            R.consol_stop_met = true;
        }

        if (!R.consol_stop_met) {
            const double left = R.consol_excess_pore.back();
            char buf[600];
            if (in.stop == ConsolidationStop::MinExcessPore)
                std::snprintf(buf, sizeof(buf),
                    "The consolidation phase did not reach its target: after %d time steps and "
                    "t = %.6g day the maximum excess pore pressure is still %.4g kPa, above the "
                    "%.4g kPa asked. Raise the phase's maximum number of steps, or end the phase "
                    "on a time interval -- the run is refused rather than reported as if the "
                    "target had been met.",
                    steps_done, t_now, left, in.stop_min_pore);
            else
                std::snprintf(buf, sizeof(buf),
                    "The consolidation phase did not reach its target: after %d time steps and "
                    "t = %.6g day the maximum excess pore pressure is %.4g kPa of the %.4g kPa "
                    "generated -- a degree of consolidation of %.1f%% (pressure ratio) against the "
                    "%.1f%% asked. Raise the phase's maximum number of steps, or end the phase on "
                    "a time interval -- the run is refused rather than reported as if the target "
                    "had been met.",
                    steps_done, t_now, left, pore_peak,
                    100.0 * (pore_peak > 0.0 ? 1.0 - left / pore_peak : 1.0),
                    100.0 * in.stop_degree);
            R.message = buf;
            return false;
        }
    }

    // What the phase reached, in the definition it was asked in (results.hpp ConsolidationStop).
    // pore_peak = 0 means the staged change generated no excess pore pressure at all: the ratio has
    // no denominator, the criterion is trivially true, and the caller warns rather than printing a
    // "90% consolidated" that describes nothing.
    R.consol_stop = in.stop;
    R.consol_pore_reference = pore_peak;
    R.consol_degree_reached =
        pore_peak > 0.0 ? 1.0 - R.consol_excess_pore.back() / pore_peak : 1.0;
    // A target the phase was ALREADY inside when it started is met at the first time step, and the
    // time that comes out is the size of that step -- not a consolidation time. The number is not
    // wrong; the reading of it would be, so the run says so instead of letting the report present
    // "reached the target in 1.6e-05 day" as a settlement rate.
    if (in.stop != ConsolidationStop::TimeInterval) {
        const bool never_rose = in.stop == ConsolidationStop::MinExcessPore
                                    ? pore_peak <= in.stop_min_pore
                                    : pore_peak <= 1e-9;   // the kPoreFloor of the march above
        if (never_rose)
            add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A013", "consolidation",
                           "this phase was asked to run until the excess pore pressure fell to a "
                           "target, but the staged change never generated excess pore pressure "
                           "above it (peak " + std::to_string(pore_peak) +
                           " kPa), so the target was already met at the first time step. The time "
                           "reported is that step, not a consolidation time. Check that the load / "
                           "excavation this phase is meant to apply is actually switched on in it.");
    }

    // Final effective Gauss state: from the return mapping (elastoplastic) or recovered from the
    // total displacement of the linear-elastic skeleton (sigma' = init + D B v_final).
    std::vector<GaussState> committed;
    if (nonlinear_soil) {
        committed = std::move(state);
    } else {
        committed = init;
        recover_consolidation_stress(mesh, dofs, models, v_total, in.active, committed);
    }
    R.disp = v_total.head(mesh.node_count * 2);
    R.stress = recover_nodal_stresses_from_gauss(mesh, committed, in.active);
    R.load_factor = 1.0;
    R.iterations = steps_done;

    // The structures, reported by the one function that reports them for both coupled phases
    // (phase_solver/coupled_structures.hpp): the diagrams from the TOTAL displacement -- this
    // phase's increment plus the parent's datum, because structural elements are
    // total-displacement formulated -- and the three findings the elastic branch owes the reader.
    if (in.structures) {
        Eigen::VectorXd disp_total = v_total;
        if (in.carry_full && in.carry_full->size() == v_total.size()) disp_total += *in.carry_full;
        report_coupled_structures(*in.structures, in.diagrams, in.iface_diagrams, mesh, dofs,
                                  disp_total, R);
    }

    committed_out = std::move(committed);
    return true;
}

} // namespace katai::core
