#pragma once
// The published facade (Stage D, layer 3): the ONE header a front end -- the CLI, the
// GUI, the Python binding, any future consumer -- programs against. It re-exports, by
// name, exactly the vocabulary the section 7.3 boundary sanctions: TOTAL over the
// input contract, CLOSED over the implementation.
//
//   build      Project and every schema component type and enum -- construct in code
//   files      load_project / save_project / project_from_json / project_to_json
//              (reader notes surface through the optional Issue vector)
//   validate   validate_project -> ValidationReport (Issue, Severity)
//   run        Job / JobState / PhaseProgress / PhaseTiming -- submit, observe,
//              cancel. Job::run() executes the input contract itself, so a consumer
//              of this surface cannot solve an invalid project.
//   results    MeshResult / SolveResult (displacements, stresses, pore pressures,
//              structural diagrams, factor of safety, time series)
//
// What is deliberately ABSENT is the boundary: no stiffness or element state, no
// Gauss points, no solver iterations, no detail:: namespace, no engine header. A
// consumer that needs those is not a front end and does not sit on this surface.
// test_api_surface pins sufficiency: the whole workflow -- construct, save, load,
// validate, run, read, refuse -- compiles and runs against this header alone.
//
// New API is born in the convention namespace (katai::api); the exported names keep
// their historical spellings underneath until the deferred rename (ARCHITECTURE.md).

#include <katai/api/version.hpp>              // kVersion/kVersionDate/kAppName: the ONE identity
#include <katai/io/project_io.hpp>
#include <katai/io/results_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/materials/material_model.hpp>   // the undrained-stiffness conversions (MMM 2.4)
#include <katai/jobs/job.hpp>
#include <katai/linsolve/direct_solver.hpp>   // backend_name(): result provenance

namespace katai::api {

// -- build: the input contract's types ------------------------------------------------
using model::Project;
using model::Material;
using model::PlateMaterial;
using model::AnchorMaterial;
using model::GeogridMaterial;
using model::EmbeddedBeamMaterial;
using model::SoilPolygon;
using model::StructElement;
using model::Load;
using model::PrescribedDisp;
using model::Phase;
using model::BoundaryConditions;
using model::MeshSettings;
using model::SoilModel;
using model::Drainage;
using model::PhaseType;
using model::SeismicWave;
using model::DesignApproach;
using model::StructKind;
using model::LoadKind;
using model::BCType;
using model::FlowBCType;
using model::InitialProcedure;

// -- files ----------------------------------------------------------------------------
using model::kProjectFileVersion;
using model::load_project;
using model::save_project;
using model::project_from_json;
using model::project_to_json;
// The .res results file (versioned binary; staleness-guarded by the model hash --
// canonically fnv1a64(project_to_json(project)), the same in every front end).
using app::kResultsFileVersion;
using app::save_results;
using app::load_results;
using app::fnv1a64;

// -- validate -------------------------------------------------------------------------
using io::Issue;
using io::Severity;
using io::ValidationReport;
using io::validate_project;

// -- run ------------------------------------------------------------------------------
using jobs::Job;
using jobs::JobState;
using jobs::PhaseProgress;
using jobs::PhaseTiming;

// -- results --------------------------------------------------------------------------
using app::MeshResult;
using app::SolveResult;
// Everything a run did that its file does not literally say -- clipped geometry, a fallback
// taken, an object refused. A front end that hides this list hides the one report that is
// made against the MESH, which the input contract (validated before meshing) cannot see.
using app::Diagnostic;
using app::DiagnosticSeverity;

// -- derived input ---------------------------------------------------------------------
// What a material's undrained stiffness inputs MEAN, computed the way the solve computes
// them (PLAXIS MMM section 2.4). This is not solver state: it is the reading of an input,
// and a front end that lets an engineer type nu_u or Skempton's B has to be able to show
// the other two numbers, because they are the same statement in the units the engineer
// happens to think in. Whichever of the two was entered, all four fields come back filled.
struct UndrainedStiffness {
    double k_eff = 0.0;    // K' of the effective elastic pair the material's model reads
    double nu_u = 0.0;     // equivalent undrained Poisson ratio (Eq. 2-55 when B was entered)
    double kw_over_n = 0.0;  // pore-fluid bulk stiffness (Eq. 2-50) [kN/m2]
    double skempton_B = 0.0; // Skempton's B (Eq. 2-57)
};
inline UndrainedStiffness undrained_stiffness(const Material& m) {
    // The Hardening Soil family's elasticity is the unload/reload pair; E and nu are boxes
    // it never reads. The same choice the constitutive catalogue makes when it builds Kw/n.
    const bool hs = m.model == SoilModel::HardeningSoil || m.model == SoilModel::HSsmall;
    const double E = hs ? m.Eurref : m.E;
    const double nu = hs ? m.nu_ur : m.nu;
    UndrainedStiffness r;
    if (!(E > 0.0) || !(nu > -1.0 && nu < 0.5)) return r;
    r.k_eff = E / (3.0 * (1.0 - 2.0 * nu));
    r.nu_u = m.und_mode == 1 ? core::undrained_poisson_from_skempton(m.skempton_B, nu) : m.nu_u;
    if (!(r.nu_u < 0.5)) return r;
    r.kw_over_n = 3.0 * (r.nu_u - nu) / ((1.0 - 2.0 * r.nu_u) * (1.0 + nu)) * r.k_eff;
    r.skempton_B = core::skempton_from_kw_over_n(r.kw_over_n, r.k_eff);
    return r;
}

// -- provenance -----------------------------------------------------------------------
// The linear-solver backend actually linked into this program (never what the build
// merely wanted): a published result is traceable to the solver that produced it.
using linsolve::backend_name;

}  // namespace katai::api
