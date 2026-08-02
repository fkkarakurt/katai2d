#pragma once
// Shear strength reduction method (phi-c reduction) for slope stability -- P1.7.
//
// The factor of safety is the strength reduction factor SRF at which a slope,
// under its full self-weight, can no longer reach static equilibrium. The
// strength parameters are factored as
//     c_f = c / SRF,   phi_f = atan(tan(phi) / SRF)
// (Rocscience/Griffiths & Lane), and for each trial SRF a nonlinear Mohr-Coulomb
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
// the GUI can show WHERE the slope is failing (PLAXIS "Safety" phase, incremental displacement).
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
};
// profile: optional depth gradient (materials/material_model.hpp MaterialProfile), passed
// through to solve_nonlinear as-is. If c'_inc is GIVEN, the φ-c reduction divides c'_ref
// and the gradient carries c'(y) → strength growing with depth enters the FoS. Empty ⇒ old
// behaviour exactly.
SafetyResult safety_analysis(const mesh::Mesh& mesh, const DofMap& dofs,
                             const Eigen::VectorXd& gravity_load,
                             const std::vector<MaterialModel>& materials,
                             const LinearSolve& linear_solve,
                             const StrengthReductionOptions& options = {},
                             const std::vector<GaussState>& initial_state = {},
                             const std::vector<char>& active_element = {},
                             const std::vector<MaterialProfile>& profile = {});

} // namespace katai::core
