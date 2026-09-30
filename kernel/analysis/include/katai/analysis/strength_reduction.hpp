#pragma once
// Shear strength reduction method (phi-c reduction) for slope stability -- P1.7.
//
// The factor of safety is the strength reduction factor SRF at which a slope,
// under its full self-weight, can no longer reach static equilibrium. The
// strength parameters are factored as
//     c_f = c / SRF,   phi_f = atan(tan(phi) / SRF)
// (Griffiths & Lane 1999), and for each trial SRF a nonlinear Mohr-Coulomb
// analysis under gravity is attempted. Failure to converge is the standard
// indicator of collapse; the critical SRF is bracketed by bisection.

#include <functional>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

struct StrengthReductionOptions {
    double srf_min = 0.5;        // assumed stable lower bracket
    double srf_max = 2.5;        // assumed unstable upper bracket
    int bisection_iterations = 12;
    NewtonOptions newton{};      // per-trial nonlinear solve controls
    // INCREMENTAL strength reduction (safety_analysis): the ground is brought to equilibrium
    // under its own weight ONCE, at full strength, and the strength is then reduced step by step
    // from each converged state -- the step growing while it converges and halving when it does
    // not -- until a step of fos_resolution x SRF fails. The factor of safety is the middle of
    // that last bracket. Off (or KATAI_SRM_BISECTION set, or anything the incremental path does
    // not carry) falls back to the bisection, which solves every trial from the unstressed state.
    bool incremental = true;
    double msf_step = 0.1;           // the first strength-reduction increment
    double fos_resolution = 1e-3;    // relative width of the final bracket
};

// Bisect the strength reduction factor to the slope's factor of safety. The base
// material is a single Mohr-Coulomb model; gravity_load is the free-DOF body-force
// vector (e.g. from assemble_gravity). linear_solve should back a non-symmetric
// system (non-associated flow gives an unsymmetric tangent).
double factor_of_safety(const mesh::Mesh& mesh, const DofMap& dofs,
                        const Eigen::VectorXd& gravity_load,
                        const MaterialModel& base, const LinearSolve& linear_solve,
                        const StrengthReductionOptions& options = {});

// Multi-material phi-c reduction + the FAILURE MECHANISM (for the GUI). Each Mohr-Coulomb / Hardening
// Soil material's strength is factored by the trial SRF (c_f=c/SRF, phi_f=atan(tan phi/SRF)); the slope
// is re-solved under gravity from the unstressed state each trial. Returns the factor of safety and the
// NewtonResult at the highest STILL-STABLE SRF -- its displacement localizes along the slip surface, so
// the GUI can show WHERE the slope is failing (the Safety phase's incremental displacement).
struct SafetyResult {
    double fos = 0.0;          // factor of safety = critical strength reduction factor
    NewtonResult mechanism;    // solve at the highest stable SRF (displacement = slip surface)
    bool ok = false;           // a valid result was obtained
    // True when a collapse was actually observed during the bisection, so the factor of safety is
    // bracketed (finite). False means the slope stayed stable all the way to srf_max -- NO failure
    // mechanism developed, and `fos` is only a LOWER BOUND (the cap), not a real factor of safety.
    // (A laterally-confined block under self-weight, an over-restrained boundary, or a model with no
    // destabilizing load never collapses; reporting the cap as a definitive FoS would be misleading.)
    bool bracketed = false;
    // Started from a parent phase, and the ground did not carry this phase's own configuration
    // change at full strength (the SRF = 1 step): no factor of safety is defined.
    bool start_failed = false;
};
// profile: optional depth gradient (materials/material_model.hpp MaterialProfile), passed
// through to solve_nonlinear as-is. If c'_inc is GIVEN, the φ-c reduction divides c'_ref
// and the gradient carries c'(y) → strength growing with depth enters the FoS. Empty ⇒ old
// behaviour exactly.
//
// structures: the phase's active structural elements, solved in every trial together with the
// soil. gravity_load must then already hold their self-weight -- a structure present without its
// weight is a different model. What the reduction touches is decided by what a strength IS:
//  - an INTERFACE is a soil strength (c_i = R_inter c, tan φ_i = R_inter tan φ), so it is reduced
//    exactly as the soil is -- c_i / SRF, tan φ_i / SRF, and its tension cut-off / SRF;
//  - a plate's M_p / N_p, an anchor's and a geogrid's capacity and an embedded beam's skin and base
//    resistance are the STRUCTURE's own declared capacities, so they stay at their input values.
// With an empty initial_state the search starts from the unstressed state, and every interface
// starts unstressed with it: its wished-in-place normal stress σ_n0 belongs to a geostatic seed
// this search does not have, and kept alone it would push the two faces apart with nothing in
// the soil to balance it. A prestressed anchor has the same problem with no such fix -- its
// lock-off force presupposes a ground that has already moved -- and is refused by the caller.
//
// baseline + init_struct: START FROM A PARENT PHASE. Given initial_state (the parent's committed
// stresses) together with baseline (the parent's internal force in equation space -- soil plus
// structures at init_struct, exactly the constant force a chained phase holds), the search does
// not re-solve the self-weight: the SRF = 1 state is the parent's equilibrium plus this phase's
// own configuration change (gravity_load - baseline, ramped as a chained phase ramps it; a nil
// step when nothing changed), and the strength is reduced incrementally from it with the structures carried from
// init_struct (datum + plastic state). Interfaces keep their sigma_n0, a lock-off force is part
// of the carried anchor state, and undrained soil keeps the pore pressures the phases built.
SafetyResult safety_analysis(const mesh::Mesh& mesh, const DofMap& dofs,
                             const Eigen::VectorXd& gravity_load,
                             const std::vector<MaterialModel>& materials,
                             const LinearSolve& linear_solve,
                             const StrengthReductionOptions& options = {},
                             const std::vector<GaussState>& initial_state = {},
                             const std::vector<char>& active_element = {},
                             const std::vector<MaterialProfile>& profile = {},
                             const Structures& structures = {},
                             const Eigen::VectorXd& baseline = {},
                             const StructuralInit& init_struct = {});

} // namespace katai::core
