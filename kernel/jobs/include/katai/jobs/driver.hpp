#pragma once
// The analysis driver (layer 2, katai/jobs): resolve a katai::model project into a
// finite-element problem, run its phases through the engine's phase strategies, and
// return displacement / stress / structural results. Every front end -- GUI, CLI,
// scripting -- runs a project through this seam, so a .k2d means the same run
// everywhere (Stage D; the physics itself left this file for the engine in Stage B).
//
// Historical name build_problem.hpp; relocated with its schema-resolution
// content intact once phase solving lived engine-side. Public spellings remain
// katai::app pending the deferred namespace decision (see ARCHITECTURE.md).

// Section 5.2, batch 2: the include set below is the DECLARATION closure only. The
// twenty-odd engine headers the solve bodies consume (phase strategies, dynamics,
// consolidation, interfaces, assembly, solvers, ...) moved to src/driver.cpp with the
// bodies -- every consumer of this header used to parse all of them on every build.
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/constants.hpp>            // gamma_w, g, utilisation cap (Stage B2: engine-owned)
#include <katai/analysis/results.hpp>              // SolveResult family (Stage B1: engine-owned)
#include <katai/analysis/step_size_advisor.hpp>    // dt adequacy warnings (Stage B7)
#include <katai/materials/material_model.hpp>      // GaussState (PhaseIO), MaterialParams
#include <katai/materials/profile_builder.hpp>     // depth-gradient profiles (Stage B6)
#include <katai/materials/registry.hpp>            // constitutive catalogue (name -> entry)
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

namespace katai::app {

// The result types live in the engine now (katai/analysis/results.hpp, Stage B1):
// a result is engine vocabulary, and a front end only displays it. Aliased here so
// every existing katai::app spelling -- the application, results_io and the tests --
// keeps meaning the same types while the driver itself is still being dissolved.
// One deliberate change rode along: SolveResult::design_approach now carries the
// ENGINE's DesignApproach (same enumerators, same file-stable values); the schema's
// enum is mapped once, where the phase configuration enters (to_core_design_approach).
using katai::core::StructForce;
using katai::core::InterfaceResult;
using katai::core::StructCarryState;
using katai::core::SolveResult;

// The analysis constants live in the engine now (katai/analysis/constants.hpp,
// Stage B2); aliased for the application's and tests' existing spellings.
using katai::core::kGammaWater;
using katai::core::kGravity;
using katai::core::kUtilCap;

// Time-step advisors live in the engine now (Stage B7, katai/analysis/
// step_size_advisor.hpp) over neutral inputs; the wrappers here keep the old
// signatures, apply the phase-type gates (schema knowledge), and copy fields
// across. The warning texts are the engine's, byte for byte -- the editor and
// the tests pin them.
inline std::string dynamic_step_warning(const model::Phase& ph) {
    if (ph.type != model::PhaseType::Dynamic) return "";
    katai::core::DynamicStepInput in;
    in.duration = ph.duration;
    in.time_steps = ph.time_steps;
    in.record = ph.seismic_wave == model::SeismicWave::Record;
    in.record_samples = (int)ph.accel_record.size();
    in.record_dt = ph.record_dt;
    in.ricker = ph.seismic_wave == model::SeismicWave::Ricker;
    in.frequency = ph.seismic_freq;
    in.nonlinear = ph.dynamic_nonlinear;
    return katai::core::dynamic_step_warning(in);
}

// Depth-varying profiles (PLAXIS E'_inc / c'_inc) are built by the engine now (Stage B6,
// katai/materials/profile_builder.hpp), gated by the catalogue's profile_E / profile_c
// flags instead of an inline model-type test here -- which inputs a model reads is
// catalogue knowledge, stated once. The driver wrapper lives after the constitutive
// seam below (it resolves entries by name).

// EVERY deformation phase now carries the profile into its own assembly and evaluates E'(y) / c'(y) at
// the stress point: Plastic/K0 + Safety via solve_nonlinear, Dynamic via assemble_stiffness,
// Consolidation and Fully-coupled via their own coupled assemblies. (Transient flow solves head only --
// no stiffness, so a profile is irrelevant there rather than ignored.) Nothing to refuse any more.
//
// One subtlety worth stating: phi-c reduction (Safety) divides c'_inc by srf along with c'_ref, because
// c(y) = c'_ref + c'_inc (y_ref - y) is ONE strength -- factoring only the reference would leave the deep
// soil at full gradient strength and report an FoS that is too high. E'_inc is a stiffness, so phi-c
// reduction leaves it alone (strength_reduction.cpp).

// Vermeer-Verruijt first-step guidance lives in the engine now (Stage B7); this seam
// resolves each schema material to the physical quantities the advisor needs. How a
// model yields its oedometer modulus is per-model knowledge kept here with the other
// schema mappings: Soft Soil (Creep) uses the NC oedometric tangent d(sigma'_v)/d(eps1)
// = sigma'_v/lambda*, evaluated at the 100 kPa reference (dt_crit is only a WARNING
// threshold; the stress-dependent true value is the phase solve's own business); HS
// uses Eoed_ref; everything else the constant-E closed form.
inline std::string consolidation_step_warning(const model::Project& pr, const katai::mesh::Mesh& mesh,
                                              const model::Phase& ph) {
    using PT = model::PhaseType;
    if (ph.type != PT::Consolidation && ph.type != PT::FullyCoupled) return {};
    std::vector<katai::core::ConsolidationStepMaterial> mats(pr.materials.size());
    for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
        const auto& M = pr.materials[mi];
        const bool hs = M.model == model::SoilModel::HardeningSoil || M.model == model::SoilModel::HSsmall;
        mats[mi].Eoed = (M.model == model::SoilModel::SoftSoil ||
                         M.model == model::SoilModel::SoftSoilCreep)
                            ? 100.0 / std::max(M.lam_star, 1e-6)
                            : (hs ? M.Eoedref : katai::core::oedometer_modulus(M.E, M.nu));
        mats[mi].k_y = M.ky;
        mats[mi].porosity = M.e_init / (1.0 + M.e_init);
    }
    return katai::core::consolidation_step_warning(mats, mesh, ph.duration, ph.time_steps);
}

// Initial-stress phase (PLAXIS): K0 procedure (geostatic equilibrium, ~zero self-weight
// displacement) or gravity loading (self-weight produces settlement from a stress-free start).
// Consolidation is a chained, time-dependent (Biot) phase: the configuration's load increment is
// applied at t=0 and the resulting excess pore pressure dissipates over the phase's time interval.
enum class InitialPhase { K0Procedure, GravityLoading, Safety, Consolidation, TransientFlow, FullyCoupled, Dynamic };

// The .k2d `initial_procedure` resolved to the driver's InitialPhase. Every front end
// resolves through here, so a checked-in file means the same run everywhere.
inline InitialPhase initial_phase_from(model::InitialProcedure p) {
    switch (p) {
        case model::InitialProcedure::GravityLoading: return InitialPhase::GravityLoading;
        case model::InitialProcedure::Safety:         return InitialPhase::Safety;
        default:                                      return InitialPhase::K0Procedure;
    }
}

// Staged-construction plumbing for one phase of a multi-phase run (solve_phases). `config`
// selects the active objects (null = everything active = classic single-phase). For a CHAINED
// phase (after the initial one) `init_states` carries the previous phase's committed Gauss
// stresses and the solve ramps the configuration imbalance d = f(active) - f_int(committed,
// active) -- the PLAXIS SumMstage: excavation unloading, fill weight, structure installation
// all emerge from that single rule. `out_states` returns this phase's committed stresses.
// The numerical controls the driver normally derives from the material class: the tolerated
// force residual and the number of load increments (see driver.cpp -- Hardening Soil runs at a
// PLAXIS-realistic 1%, Mohr-Coulomb at 1e-6, a linear problem at 1e-10). Deriving them is a
// convenience, not a contract, and it has a cost a reviewer is entitled to ask about: is a
// published number a physics result, or an artefact of the tolerance it happened to be computed
// at? That question can only be answered by re-running the same problem at other tolerances,
// which is what this override is for (tests/test_tolerance_independence.cpp).
//
// Zero means "keep the derived default", so an untouched PhaseIO reproduces today's runs
// bit-for-bit. Exposing these per phase in the .k2d contract is a separate step (the hardening
// plan's WP-3); until then they are a jobs-layer seam, not a file-format promise.
struct NumericalControls {
    double tolerance = 0.0;   // tolerated relative force residual; 0 = by material class
    int steps = 0;            // load increments; 0 = by material class
};

struct PhaseIO {
    const model::Phase* config = nullptr;
    NumericalControls numeric;   // all-zero: the driver chooses, exactly as before
    const std::vector<katai::core::GaussState>* init_states = nullptr;
    std::vector<katai::core::GaussState>* out_states = nullptr;
    bool chained = false;
    // The PARENT phase's converged result. A Dynamic phase needs it: its own solve is the LINEAR
    // dynamic INCREMENT about the static state (soil unloads/reloads elastically at small strain), so
    // the action a section is actually designed for is static + dynamic. Reporting the increment alone
    // leaves the superposition -- and every strength check -- to the user. PLAXIS instead continues the
    // parent phase's state and reports the TOTAL (Sci sec 6.4; Ref sec 7.9.2.3, 9.4.5).
    // Null (or a parent that is itself Dynamic) = no static state to add; the phase then reports the
    // dynamic action alone and says so.
    const SolveResult* prev = nullptr;
};

// The construction seam between the project schema and the constitutive
// catalogue: the schema enum resolves to a canonical registry name, and the
// schema's material fields are copied into the registry's neutral parameter
// block. Everything constitutive that used to be inlined here -- parameter
// wiring, K0^NC memory, cap calibration, the Undrained (B) Tresca override and
// the honest drainage refusals -- now lives with the models in
// katai/materials/registry.hpp; this function only ferries fields across.
inline const char* constitutive_name(model::SoilModel sm) {
    switch (sm) {
        case model::SoilModel::LinearElastic: return "LinearElastic";
        case model::SoilModel::MohrCoulomb:   return "MohrCoulomb";
        case model::SoilModel::HardeningSoil: return "HardeningSoil";
        case model::SoilModel::HSsmall:       return "HSsmall";
        case model::SoilModel::SoftSoil:      return "SoftSoil";
        case model::SoilModel::SoftSoilCreep: return "SoftSoilCreep";
    }
    return "";
}

// The field-by-field copy itself (to_material_params) moved to src/driver.cpp with the
// bodies -- nothing outside the driver constructs registry parameters from the schema.
katai::core::MaterialParams to_material_params(const model::Material& m);

// Depth-varying profile table, indexed by material id exactly like `mats` / `models`
// (engine: profile_builder.hpp; the catalogue's flags gate which gradients apply).
inline std::vector<katai::core::MaterialProfile> build_profiles(const model::Project& pr) {
    std::vector<katai::core::MaterialProfile> prof(pr.materials.size());
    for (size_t i = 0; i < pr.materials.size(); ++i) {
        const auto& m = pr.materials[i];
        if (const katai::core::ModelEntry* entry = katai::core::find_model(constitutive_name(m.model)))
            prof[i] = katai::core::build_profile(*entry, m.E_inc, m.c_inc, m.y_ref);
    }
    return prof;
}

// The single-phase (or one-phase-of-a-chain) solve: resolve the schema into the neutral
// problem, run the phase family through the engine strategies, and return fields,
// diagrams and messages. The BODY is compiled once in kernel/jobs/src/driver.cpp
// (section 5.2): as a header-inline function every driver test paid ~220 s recompiling
// it; as a library symbol the suite pays for it exactly once.
SolveResult solve_gravity_le(const model::Project& pr, const katai::mesh::Mesh& mesh_in,
                             InitialPhase phase = InitialPhase::K0Procedure,
                             const Eigen::VectorXd* flow_head = nullptr,
                             const PhaseIO& io = {});

// Multi-phase staged construction (PLAXIS Phases): runs the INITIAL phase (project.initial
// activation, `init_phase` = K0 / gravity) and then every project phase in order, carrying the
// committed Gauss stresses forward (SumMstage chaining). Phase displacements are INCREMENTAL
// (each phase starts from zero displacement, PLAXIS-style); the caller may accumulate them.
// Stops at the first failed phase (its honest message is in the last result).
// `on_phase`, when set, is called BEFORE each phase solve with (current 0-based index, total phase
// count, phase display name). It lets a GUI report "Phase k of n: <name>" live from a worker thread
// (purely informational; must be thread-safe on the caller's side).
// `cancelled`, when set, is polled at every phase boundary (before and right after `on_phase`, so
// the progress receiver itself may cancel the announced phase). Returning true stops the run and
// the COMPLETED prefix is returned -- cooperative cancellation: a running phase is never
// interrupted. The Job object (katai/jobs/job.hpp) is the intended consumer.
using PhaseProgress = std::function<void(int, int, const std::string&)>;
std::vector<SolveResult> solve_phases(const model::Project& pr,
                                      const katai::mesh::Mesh& mesh_in,
                                      InitialPhase init_phase = InitialPhase::K0Procedure,
                                      const PhaseProgress& on_phase = {},
                                      const std::function<bool()>& cancelled = {});

} // namespace katai::app
