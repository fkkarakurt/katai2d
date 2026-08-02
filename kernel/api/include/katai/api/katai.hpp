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

// -- provenance -----------------------------------------------------------------------
// The linear-solver backend actually linked into this program (never what the build
// merely wanted): a published result is traceable to the solver that produced it.
using linsolve::backend_name;

}  // namespace katai::api
