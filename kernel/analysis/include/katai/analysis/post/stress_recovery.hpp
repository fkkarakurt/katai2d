#pragma once
// Stress recovery (post-processing) — P0.7.
//
// In FEM the stress is most accurate at the Gauss (integration) points:
//     σ = D · B · u_e   (Voigt: [sxx, syy, sxy])
// For contours/visualization and node-based output the Gauss values are extrapolated to
// the nodes (3 Gauss points → linear field → 6 nodes), then averaged across elements
// (nodal smoothing). Reproduces a linear stress field EXACTLY.

#include <array>
#include <vector>

#include <Eigen/Core>

#include <katai/fem/elements/tri6.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/materials/material_model.hpp>  // GaussState
#include <katai/mesh/mesh.hpp>

namespace katai::core {

using ElementDisplacement = Eigen::Matrix<double, 12, 1>;

// Stresses at an element's 3 Gauss points (Voigt [sxx, syy, sxy]).
std::array<Eigen::Vector3d, tri6::kGaussCount> gauss_point_stresses(
    const tri6::NodeCoords& coords, const Eigen::Matrix3d& constitutive,
    const ElementDisplacement& element_displacement);

// Extrapolates the values at the 3 Gauss points to the 6 nodes (linear field).
std::array<Eigen::Vector3d, tri6::kNodeCount> extrapolate_gauss_to_nodes(
    const std::array<Eigen::Vector3d, tri6::kGaussCount>& gauss_values);

// Node-averaged stress field.
struct NodalStressField {
    std::vector<Eigen::Vector3d> stress;  // node_count entries of [sxx, syy, sxy]
};

// Recovers the nodal stresses for the whole mesh. full_displacement layout is
// node*2 + component (same as the expand_to_full output).
NodalStressField recover_nodal_stresses(
    const mesh::Mesh& mesh, const std::vector<LinearElastic>& materials,
    const Eigen::VectorXd& full_displacement);

// Recovers the nodal stresses directly from the SOLVER's committed Gauss states (stress =
// effective σ'; K0 pre-stress + plastic history INCLUDED — unlike recover_nodal_stresses it
// does not re-derive from displacement alone). Each element's Gauss values are extrapolated
// to the nodes (tri6: 3 Gauss→6 nodes linear; tri15: 12 Gauss→15 nodes least-squares
// cubic), then averaged across elements. gauss_states order matches the solver
// (e*kGaussCount + g). The CORRECT path to the true stress field in staged construction,
// the K0 procedure and nonlinear (MC/HS) analyses.
// active_element: staged-construction mask (empty = all) — passive (excavated) elements do
// not join the nodal average; nodes touching only passive elements stay zero.
NodalStressField recover_nodal_stresses_from_gauss(
    const mesh::Mesh& mesh, const std::vector<GaussState>& gauss_states,
    const std::vector<char>& active_element = {});

// Full-DOF nodal internal force from the committed Gauss stresses, F = sum_e B^T sigma
// (axisymmetric: r-weighted per radian, hoop included). At a FIXED dof this is the
// SUPPORT REACTION the constraint exerts (soil contribution; loads applied directly on
// a fixed node are not subtracted). Layout node*2 + component; activity-masked like the
// recovery above. Plane strain and axisymmetric, tri6 and tri15.
Eigen::VectorXd nodal_internal_force_from_gauss(
    const mesh::Mesh& mesh, const std::vector<GaussState>& gauss_states,
    bool axisymmetric, const std::vector<char>& active_element = {});

} // namespace katai::core
