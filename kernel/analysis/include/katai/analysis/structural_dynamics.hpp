#pragma once
// Soil-structure interaction (SSI) in the DYNAMIC (seismic) analysis -- the structural
// elements' ELASTIC stiffness and CONSISTENT mass into the dynamic system's K and M.
//
// The dynamic system (analysis/dynamics.hpp solve_newmark):
//     M u'' + C u' + K u = -M.r.a_g(t),      C = alpha.M + beta.K  (Rayleigh)
// The soil contribution comes from assemble_stiffness / assemble_mass (assembly/assembler.hpp).
// This header adds the STRUCTURE's contribution to the SAME K and M -> soil and structure solve
// in ONE system = soil-structure interaction (the counterpart of PLAXIS 2D Dynamics'
// plate+interface seismic analysis).
//
// STIFFNESS: assemble_structural_stiffness is the ELASTIC branch, at u=0, of solve_nonlinear's
// structural tangent (analysis/nonlinear_solver.cpp) -- the same elements, the same DOF mapping,
// the same element matrices. That IDENTITY is not an assumption: test_ssi_dynamics
// cross-verifies it independently -- a static linear solve with the K built here must reproduce
// the verified solve_nonlinear path to machine precision.
//
// v1 CONSTITUTIVE LIMIT (honest): the dynamic branch is built on a LINEAR-ELASTIC skeleton (the
// soil itself is LE too), so the structural contributions are taken elastic as well:
//   - interface: elastic (k_n, k_s) -- no Coulomb SLIP (seismic interface slip is a later
//     increment),
//   - anchor:    elastic EA/L -- no yield (F_max),
//   - geogrid:   elastic EA -- no tension-only cut.
// These limits are recorded in docs/validation/seismic-verification.md.
// EMBEDDED BEAM (pile row): IS assembled in the dynamic branch (below; the driver's Dynamic
// branch passes it too) -- the old "out of v1 scope / early-refuse" note was STALE (an audit
// finding); the skin/foot coupling is taken elastic (the same class as the linear-skeleton
// limit above).
//
// MASS: assemble_structural_mass gives mass to PLATES only (plate::PlateProps rho_A, rho_I):
//     M_translation = int rho_A N'N ds   (into the u_x and u_y blocks separately, uncoupled)
//     M_rotation    = int rho_I N'N ds   (the phi DOF)
// rho_A = w/g [Mg/m], rho_I = rho_A.d^2/12 [Mg.m] (d = equivalent thickness).
// Anchor/geogrid/interface are MASSLESS (as in PLAXIS: those elements carry no weight). A
// massless DOF (e.g. phi with rho_I=0) is fine: K_eff = K + a0.M + a1.C takes a nonsingular
// contribution from K on those rows.
//
// Mathematics + verification: docs/references/dynamic-seismic-formulation.md sec 10 (SSI),
// docs/validation/seismic-verification.md. Test: test_ssi_dynamics.

#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <katai/analysis/nonlinear_solver.hpp>   // Structures + the element structs
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/geogrid.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {


// Adds the structural elements' ELASTIC stiffness to the builder (ON TOP of the soil K -- the
// same builder). Scope: plates, plates5, anchors, geogrids, interfaces, interfaces5 and the
// embedded beams (elastic; see the constitutive limit in the header intro).
// Definition in kernel/analysis/src/structural_dynamics.cpp (section 5.2).
void assemble_structural_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                                          const Structures& structures,
                                          math::SparseMatrixBuilder& builder);

// Adds the structural CONSISTENT mass to the builder (ON TOP of the soil M -- the same builder).
// Only plates carry mass: rho_A (translation, x and y blocks separately) + rho_I (rotation).
// rho_A=rho_I=0 (the default) -> no contribution. Anchor/geogrid/interface are massless.
// Definition in kernel/analysis/src/structural_dynamics.cpp (section 5.2).
void assemble_structural_mass(const mesh::Mesh& mesh, const DofMap& dofs,
                                     const Structures& structures,
                                     math::SparseMatrixBuilder& builder);

// Consistent nodal forces of the structural SELF-WEIGHT (into the STATIC external load vector)
// -- the static counterpart of the mass twin (assemble_structural_mass): line weight
// w = rho_A.g [kN/m/m] in the vertical (-y) direction,
//     f_y(i) -= int rho_A.g . N_i ds
// The SAME elements, SAME Gauss points, SAME DOF mapping (independent-wall trans_dof included;
// the embedded beam with its own node coordinates + own extra DOFs). NO contribution to the
// rotation DOF -- a uniform line load does no work on phi (translation/rotation interpolate
// independently in Timoshenko). On a straight uniform 3-node element this is exactly the
// Simpson distribution (L/6, L/6, 2L/3); on the 5-node one Boole (7,32,12,32,7)/90.L.
// rho_A = 0 (the default w=0) -> no contribution (bit-for-bit the old behaviour).
// Anchor/geogrid/interface are weightless in PLAXIS too. The PLAXIS w input contract: for a
// wall embedded in soil the user enters w with the overlapping soil deducted
// (w = gamma_plate.d - gamma_soil.d) -- the entered value is applied as given (PLAXIS does the
// same). Definition in kernel/analysis/src/structural_dynamics.cpp (section 5.2).
void assemble_structural_weight(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const Structures& structures, double g,
                                       Eigen::VectorXd& f);

// The influence vector r of HORIZONTAL base motion: the displacement every DOF takes when the
// base translates by a unit horizontally. In a rigid translation ALL HORIZONTAL TRANSLATION
// DOFs are 1, verticals and ROTATIONS 0 ->
//     r = 1 on {soil node u_x} UNION {INDEPENDENT structural translation u_x (wall trans_dof)}.
// The seismic body force F = -M.r.a_g uses it. Skipping an independent wall DOF means the wall
// takes NO inertial force from the base -> silently wrong SSI; test_ssi_dynamics catches it
// with the rigid-mode check. Definition in kernel/analysis/src/structural_dynamics.cpp.
Eigen::VectorXd seismic_influence_x(const mesh::Mesh& mesh, const DofMap& dofs,
                                           const Structures& structures);

}  // namespace katai::core
