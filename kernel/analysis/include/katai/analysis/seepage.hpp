#pragma once
// Steady-state confined (saturated) groundwater flow — Darcy + continuity (Laplace).
// Governing equation ∇·(k ∇h) = 0; hydraulic head h = y + u/γ_w. Galerkin FE: the same
// tri6/tri15 shape functions, 1 DOF per node (head h). The resulting conductivity matrix
// H = ∫ Gᵀ k G dΩ is symmetric positive definite (PARDISO SPD).
// Inhomogeneous Dirichlet (prescribed head) is moved to the right-hand side by "lifting".
// Formulation and verification plan: docs/references/seepage-formulation.md.

#include <functional>
#include <vector>

#include <Eigen/Core>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/water_retention.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Anisotropic permeability — principal axes aligned with x,y. Unit: length/time.
// For isotropic use kx == ky.
struct Permeability {
    double kx = 1.0;
    double ky = 1.0;
};

// Accumulates the global conductivity matrix H = ∫ Gᵀ diag(kx,ky) G dΩ into the builder
// (free-free entries only) and applies the inhomogeneous Dirichlet lift to rhs:
//   H_ff h_f = Q_f − H_fp h_p   ⇒   rhs[eq_a] -= H_e(a,b)·head_prescribed[b].
// `dofs` is a 1-DOF-per-node (head h) map. `permeability` is indexed by the element
// material id. `head_prescribed` holds the prescribed head at fixed (Dirichlet) nodes
// (size = dofs.total_dofs() = node_count; free-node values are ignored). An impermeable
// boundary q_n=0 is natural (do nothing). The element type is picked from
// mesh.nodes_per_element.
void assemble_seepage(const mesh::Mesh& mesh, const DofMap& dofs,
                      const std::vector<Permeability>& permeability,
                      const std::vector<double>& head_prescribed,
                      math::SparseMatrixBuilder& builder, Eigen::VectorXd& rhs);

// Seepage → deformation coupling: builds the pore-pressure load from the steady-state
// nodal head field. At each Gauss point the pressure head ψ = h_gp − y_gp is interpolated
// (h and y with the same shape functions), u = γ_w·max(0, ψ) (suction cut-off), then
// f += ∫ Bᵀ u·m dV (plane strain m=[1,1,0]) into the MECHANICAL dof map.
// In the confined (fully saturated, ψ≥0) case u is exactly interpolated. The computed
// effective stress σ'=σ_total+u·m thus reflects the seepage pore field (Terzaghi). `head`
// is indexed by node (1-DOF flow solution); `dofs` is the mechanical (2-DOF) map. PLAXIS's
// flow→effective-stress hand-over. The element type is picked from mesh.nodes_per_element.
// active_element: staged/porosity mask (empty = all) — NonPorous elements receive NO
// pore-pressure load (PLAXIS: a non-porous material has neither initial nor excess pore
// pressure).
void assemble_pore_load_from_head(const mesh::Mesh& mesh, const DofMap& dofs,
                                  const Eigen::VectorXd& head, double gamma_w,
                                  Eigen::VectorXd& rhs,
                                  const std::vector<char>& active_element = {});

// Head-based saturation gravity: at each Gauss point γ = (ψ = h−y ≥ 0) ? γ_sat : γ_unsat,
// with h interpolated by the SAME shape functions as the pore-load kernel (saturation and
// pore pressure form a consistent pair — Terzaghi). The flow-solution version of
// assemble_gravity_phreatic: the phreatic surface is the ψ=0 contour of the flow head
// field, not a water-table polyline. `dofs` is the mechanical 2-DOF map.
void assemble_gravity_from_head(const mesh::Mesh& mesh, const DofMap& dofs,
                                const std::vector<double>& gamma_unsat,
                                const std::vector<double>& gamma_sat,
                                const Eigen::VectorXd& head, Eigen::VectorXd& rhs);

// Consistent nodal flux Q = K·h: each element's He·h_e is scattered to the global nodes.
// At source-free (free) nodes Q≈0; at prescribed-head (Dirichlet) nodes Q = the boundary
// discharge ("flux reaction"). The sum of Q over an inlet boundary = the seepage quantity;
// the sum over ALL nodes = 0 (mass conservation — the row sums of K are zero, constant
// head ⇒ zero flux). `head` is the full nodal field (size node_count). The element type is
// picked from mesh.nodes_per_element.
Eigen::VectorXd compute_nodal_flux(const mesh::Mesh& mesh,
                                   const std::vector<Permeability>& permeability,
                                   const Eigen::VectorXd& head);

// Prescribed normal flux (Neumann) boundary condition: Q_i += ∫ N_i q_n ds along an edge
// chain. q_n = specific discharge normal to the boundary (Darcy), inflow POSITIVE (unit:
// length/time; volumetric discharge per length at unit plane thickness). For wells/
// infiltration/prescribed-discharge boundaries. Impermeable q_n=0 is already natural (this
// function is not needed). ordered_boundary_nodes is an ordered edge chain (tri6:
// corner,mid,corner,...; tri15: 5 nodes per edge). `dofs` is the 1-DOF (head) map.
void assemble_seepage_flux(const mesh::Mesh& mesh, const DofMap& dofs,
                           const std::vector<int>& ordered_boundary_nodes,
                           double q_n, Eigen::VectorXd& rhs);

// --- Unconfined (free-surface) seepage ------------------------------------------
// SPD linear solve callback (PARDISO); keeps the kernel decoupled from the solver (see
// nonlinear_solver.hpp LinearSolve).
using SeepageLinearSolve =
    std::function<Eigen::VectorXd(const math::CsrMatrix&, const Eigen::VectorXd&)>;

// Unsaturated-zone relative permeability + Picard iteration settings.
struct UnconfinedOptions {
    double k_min = 1e-4;      // relative-permeability floor in the unsaturated zone (k·k_min)
    double transition = 0.1;  // pressure-head ψ transition width (length) — ~1-2× the
                              // near-surface element size; smaller = more accurate discharge.
    // Unsaturated k_rel model. nullptr → simple linear ramp (k_min..1 over `transition`
    // width); if provided (per material, indexed by element_material), van Genuchten/Mualem
    // k_rel(ψ) is used — CONSISTENT with transient/fully-coupled flow (water_retention.hpp).
    // Suction ψ=−(h−y)=y−h; saturated (h≥y) → k_rel=1; floor = WaterRetention.k_rel_min.
    // `transition` remains only for the free-surface active-set dead band (it no longer
    // shapes k_rel).
    const std::vector<WaterRetention>* retention = nullptr;
    int max_iter = 500;
    double tol = 1e-5;        // convergence: max |Δh|/scale of the relaxed head (saturated
                              // zone). ~1e-5 = the mesh-resolution limit of the free surface
                              // (physically settled); machine precision is not expected.
    double relax = 0.15;      // under-relaxation — strong (breaks the free-surface
                              // flip-flop; 0.3+ can get stuck in a limit cycle). For hard
                              // problems (seepage face) drop to 0.05 → the floor drops
                              // ∝ relax × raw amplitude.
};

struct UnconfinedResult {
    Eigen::VectorXd head;     // full nodal head field (node_count)
    int iterations = 0;
    double residual = 0.0;    // last relative |Δh|
    bool converged = false;
};

// Unconfined steady-state seepage. Fixed mesh, variable relative permeability k_rel(ψ)
// (pressure head ψ=h−y; ψ≥0 saturated → k_rel=1, ψ<0 unsaturated → drops to k_min): the
// phreatic surface forms by itself on the ψ=0 contour (Bathe & Khoshgoftaar 1979; the
// relative-permeability approach of the PLAXIS flow module). Picard iteration
// (under-relaxed) — k_rel is updated each step, the SPD system is solved via
// `linear_solve`. `head_prescribed`: Dirichlet nodal values (reservoir head; free nodes
// are ignored).
UnconfinedResult solve_unconfined_seepage(
    const mesh::Mesh& mesh, const DofMap& dofs,
    const std::vector<Permeability>& permeability,
    const std::vector<double>& head_prescribed,
    const SeepageLinearSolve& linear_solve, const UnconfinedOptions& options = {});

// Unconfined discharge recovery: the k_rel-weighted (effective permeability) version of
// compute_nodal_flux — flux ~0 in the unsaturated zone. Q = K_eff·h.
Eigen::VectorXd compute_nodal_flux(const mesh::Mesh& mesh,
                                   const std::vector<Permeability>& permeability,
                                   const Eigen::VectorXd& head,
                                   const UnconfinedOptions& options);

// Unconfined + SEEPAGE FACE. On a freely draining boundary (dam downstream face,
// excavation face; above tailwater) water exits at atmospheric pressure: below the exit
// point h=y (zero pressure, discharging), above it impermeable (unsaturated). The exit
// point is NOT known A PRIORI → active-set iteration: a face node with ψ=h−y>0 becomes
// discharging (h=y is pinned); a discharging node with inward flux (Q>0) is released.
// k_rel is nested inside the phreatic-surface iteration (the DofMap is rebuilt each step).
// `fixed_nodes`/`fixed_values`: reservoir Dirichlet (parallel arrays); `seepage_nodes`:
// candidate discharge nodes. The exit point emerges by itself.
// Verification: Charny q=k(h₁²−h₂²)/(2L) (exact INCLUDING the seepage face). (Casagrande
// exit point.)
UnconfinedResult solve_unconfined_seepage_face(
    const mesh::Mesh& mesh, const std::vector<int>& fixed_nodes,
    const std::vector<double>& fixed_values, const std::vector<int>& seepage_nodes,
    const std::vector<Permeability>& permeability,
    const SeepageLinearSolve& linear_solve, const UnconfinedOptions& options = {});

} // namespace katai::core
