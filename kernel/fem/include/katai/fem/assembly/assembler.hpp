#pragma once
// Global assembly — accumulates element contributions into the reduced (free-DOF) system.
// Fixed DOFs are eliminated (see DofMap). The produced matrix is symmetric positive definite.

#include <functional>
#include <vector>

#include <Eigen/Core>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/materials/material_model.hpp>  // GaussState
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Global stiffness: scatters each element Ke into the free-DOF equations.
// materials is indexed by material id. builder produces a square system of size
// equation_count() (free-free contributions; fixed-DOF rows/columns are dropped →
// zero-Dirichlet, the right-hand side is unaffected).
// active_element: optional activity mask (empty = all active). Passive (excavated / not
// yet placed) elements do NOT enter the assembly — the SAME contract as assemble_gravity.
// Omitting it in a staged phase leaves excavated soil in the system = a silently wrong
// (unsafe-sided) result.
// profile: optional depth gradient (empty OR all uniform() ⇒ the old constant-E path
// BIT-FOR-BIT). If given, Ke = Σ_g w_g·J_g·Bᵀ·D(y_g)·B (E per Gauss point) — see
// MaterialProfile.
//
// gauss_elastic: optional PER-GAUSS-POINT elastic property (size element_count·kGaussCount,
// index e·ngp+g). If given it replaces both `materials` AND `profile`. Purpose: when the
// STRESS STATE, not material/depth, sets the stiffness (Hardening Soil / HSsmall:
// E_ur = E_ur,ref·((c·cosφ+σ3·sinφ)/(c·cosφ+p_ref·sinφ))^m), the caller can compute it
// from the committed stress and pass it directly. Empty ⇒ no effect.
void assemble_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                        const std::vector<LinearElastic>& materials,
                        math::SparseMatrixBuilder& builder,
                        const std::vector<char>& active_element = {},
                        const std::vector<MaterialProfile>& profile = {},
                        const std::vector<LinearElastic>& gauss_elastic = {});

// Global consistent mass matrix: scatters each element Mₑ = ∫ ρ NᵀN dA into the free-DOF
// equations (x and y blocks separate; no x-y coupling). density = the mass density ρ = γ/g
// indexed by material id. Dynamic analysis (M u'' + C u' + K u = F, dynamics.hpp).
// Total mass is conserved: 1ᵀ M 1 = ρ·area (Σ N_i = 1 partition of unity). SPD. Element
// type from mesh.nodes_per_element (tri6/tri15). Plane strain; axisym mass in a later phase.
// active_element: the SAME contract as assemble_stiffness (empty = all active). A passive
// element gives NO mass either — otherwise excavated soil generates seismic inertial force.
void assemble_mass(const mesh::Mesh& mesh, const DofMap& dofs,
                   const std::vector<double>& density,
                   math::SparseMatrixBuilder& builder,
                   const std::vector<char>& active_element = {});

// Water-table-aware consistent mass: ρ_sat if the Gauss point is BELOW the water table
// (=γ_sat/g, grains + pore water = TOTAL mass), ρ_unsat above. The mass counterpart of
// assemble_gravity_phreatic — since the TOTAL mass carries the seismic inertial force
// (−M·r·a_g), using γ_unsat below water understates the inertia (unsafe-sided).
// water_table_y(x) → the water-table elevation.
void assemble_mass_phreatic(const mesh::Mesh& mesh, const DofMap& dofs,
                            const std::vector<double>& rho_unsat,
                            const std::vector<double>& rho_sat,
                            const std::function<double(double)>& water_table_y,
                            math::SparseMatrixBuilder& builder,
                            const std::vector<char>& active_element = {});

// Gravity body force: f += ∫ Nᵀ (0, -γ) dA. unit_weight (γ = ρg, down +) is indexed by
// material id. rhs must have size equation_count().
// If active_element is given (empty = all active), only active elements contribute —
// active soil weight in staged construction (see staged_construction.hpp).
void assemble_gravity(const mesh::Mesh& mesh, const DofMap& dofs,
                      const std::vector<double>& unit_weight,
                      Eigen::VectorXd& rhs,
                      const std::vector<char>& active_element = {});

// Phreatic (water-table-aware) gravity body force: at each Gauss point
// γ = (y ≤ y_wt(x)) ? γ_sat : γ_unsat (saturated unit weight below the water table, moist
// above). Total-stress equilibrium: used together with the hydrostatic pore-pressure load
// (assemble_pore_pressure_load), the computed stress becomes EFFECTIVE (the buoyancy
// γ' = γ_sat − γ_w emerges naturally; see docs/references/effective-stress-formulation.md
// §3). gamma_unsat/gamma_sat are indexed by material id. water_table_y(x) returns the
// phreatic surface elevation at that x (no water → a very small value → moist everywhere).
// active_element: the staged-construction mask (empty = all active).
void assemble_gravity_phreatic(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<double>& gamma_unsat,
                               const std::vector<double>& gamma_sat,
                               const std::function<double(double)>& water_table_y,
                               Eigen::VectorXd& rhs,
                               const std::vector<char>& active_element = {});

// Lysmer-Kuhlemeyer viscous absorbing (dashpot) boundary: accumulates the damping matrix
// C_b = ∫_Γ Nᵀ D_c N dΓ along an edge chain, D_c = c_n (n⊗n) + c_t (t⊗t) (normal/tangent
// dashpot; c_n=ρ·Vp, c_t=ρ·Vs). Added to the global damping C, the outgoing wave undergoes
// ~reflection-free radiation (radiation damping). The SIGN of the normal is irrelevant
// (D_c is invariant under n→−n). ordered_boundary_nodes is an edge chain (tri6:
// corner,mid,corner,...; tri15: 5 per edge). Dynamics (D3, dynamics.hpp);
// docs/references/dynamic-seismic-formulation.md §8.
void assemble_boundary_dashpot(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<int>& ordered_boundary_nodes,
                               double c_n, double c_t,
                               math::SparseMatrixBuilder& builder);

// Uniform surface traction: f += ∫ Nᵀ (tx, ty) ds along ordered edge nodes.
// ordered_boundary_nodes is a boundary edge chain (e.g. mesh.top_nodes); triples
// (corner, mid, corner) form the tri6 edges. For a pressure p on top,
// (traction_x, traction_y) = (0, -p).
void assemble_surface_traction(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<int>& ordered_boundary_nodes,
                               double traction_x, double traction_y,
                               Eigen::VectorXd& rhs);

// LINEARLY VARYING surface/line traction (distributed load): the traction (tx_i, ty_i)
// given at each chain node is interpolated along the edge with the shape functions and the
// consistent nodal forces f += ∫ Nᵢ (Σⱼ Nⱼ tⱼ) ds accumulate. node_traction_x/y must have
// the same size as ordered_boundary_nodes. The general form of the uniform version
// (above); for the PLAXIS distributed load (q1→q2 linear between the ends). The nodes must
// form an edge chain (corner, mid, corner, ...) (existing on the surface, or a line
// conformed into the mesh).
void assemble_surface_traction_varying(const mesh::Mesh& mesh, const DofMap& dofs,
                                       const std::vector<int>& ordered_boundary_nodes,
                                       const std::vector<double>& node_traction_x,
                                       const std::vector<double>& node_traction_y,
                                       Eigen::VectorXd& rhs);

// --- Axisymmetric (r-z) variants (P1.6) ----------------------------------------
// Coordinates are interpreted as (r, z); integration is r-weighted (per radian).
// The element type is picked from mesh.nodes_per_element (tri6/tri15).

// Global axisymmetric stiffness (Ke = ∫ Bᵀ D B r dA, D 4x4 including the hoop).
void assemble_axisym_stiffness(const mesh::Mesh& mesh, const DofMap& dofs,
                               const std::vector<LinearElastic>& materials,
                               math::SparseMatrixBuilder& builder);

// Axisymmetric surface traction: f += ∫ Nᵀ (t_r, t_z) r ds (ordered edge nodes).
// For an internal pressure p on an inner surface, (t_r, t_z) = (p, 0).
void assemble_axisym_traction(const mesh::Mesh& mesh, const DofMap& dofs,
                              const std::vector<int>& ordered_boundary_nodes,
                              double traction_r, double traction_z,
                              Eigen::VectorXd& rhs);

// Axisymmetric body force (r-weighted self-weight): f_z += -γ ∫ N_i r dA.
void assemble_axisym_gravity(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<double>& unit_weight, Eigen::VectorXd& rhs,
                             const std::vector<char>& active_element = {});

// Axisymmetric consistent nodal internal force: F += ∫ Bᵀ σ r dA (σ including the hoop =
// GaussState::stress_zz). For the K0 procedure's geostatic baseline (constant_force). The
// r-weighted/4-component counterpart of the plane-strain assemble_internal_force.
void assemble_axisym_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                                    const std::vector<GaussState>& states, Eigen::VectorXd& rhs,
                                    const std::vector<char>& active_element = {});

// Axisymmetric varying (linear q1→q2) surface traction: f += ∫ N_i (Σ N_j t_j) r ds.
void assemble_axisym_traction_varying(const mesh::Mesh& mesh, const DofMap& dofs,
                                      const std::vector<int>& ordered_boundary_nodes,
                                      const std::vector<double>& tx, const std::vector<double>& ty,
                                      Eigen::VectorXd& rhs);

// Predefined (steady-state) pore-pressure load: f += ∫ Bᵀ u·m dV (plane strain
// m=[1,1,0]). Equilibrium is set up in total stress, the constitutive law produces
// effective stress σ'; with this term added to the right-hand side the computed stress is
// effective (see docs/references/effective-stress-formulation.md). pore(x, y) returns the
// pore pressure (≥ 0) at the Gauss point. Element type from mesh.nodes_per_element.
// active_element: staged-construction mask (empty = all active) — the pore load applies
// only to active soil (an excavated region has no pore pressure).
void assemble_pore_pressure_load(const mesh::Mesh& mesh, const DofMap& dofs,
                                 const std::function<double(double, double)>& pore,
                                 Eigen::VectorXd& rhs,
                                 const std::vector<char>& active_element = {});

// Internal force of the initial (committed) stress: f += ∫ Bᵀ σ0 dV (plane strain,
// free-DOF). σ0 = states[e*kGaussCount+g].stress (Voigt [sxx,syy,sxy]). If active_element
// is given, only active elements contribute. The K0 procedure's "baseline" equilibrium
// force (wished-in-place: constant_force = this + the interface σ_n0 forces → u=0 in the
// seeded state). Element type from mesh.nodes_per_element.
void assemble_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<GaussState>& states,
                             Eigen::VectorXd& rhs,
                             const std::vector<char>& active_element = {});

// Expands the free-DOF solution to the full displacement vector (total_dofs; fixed DOFs = 0).
Eigen::VectorXd expand_to_full(const DofMap& dofs,
                               const Eigen::VectorXd& free_solution);

} // namespace katai::core
