#pragma once
// Transient (time-dependent) saturated groundwater flow -- storage + Darcy continuity.
// PLAXIS 2D 2025.1 Scientific Manual sec 3.1.1-3.3 (Eq 3-10, 3-31..3-37). HEAD form h:
//   S_s dh/dt = div(k grad h) + q,   S_s = n.gamma_w/K_w (specific storage, [1/length]).
// Semi-discrete (Galerkin, the same tri6/tri15 N):  S.(dh/dt) + H.h = q,
//   H = int G' diag(kx,ky) G dOmega  (= the seepage conductivity; SPD),
//   S = int N' S_s N dOmega           (storage/capacitance; consistent-mass-like, SPD).
// Backward Euler (alpha=1, the PLAXIS default, unconditionally stable):
//   (S + dt.H).h^{n+1} = S.h^n + dt.q^{n+1}.
// Fixed dt -> A = (S+dt.H) is factorized ONCE (factor-once-solve-many; A SPD -> PARDISO mtype=2).
// Time-varying Dirichlet head (rapid drawdown) via `head_bc(node,t)`; the boundary SET is fixed
// (only the value varies in time) -> factor-once survives; the lift moves to the RHS every step.
// The steady-state limit (dh/dt->0, fixed BCs) reduces bit-close to `assemble_seepage`
// (regression). The unsaturated (van Genuchten S_e + Mualem k_rel) generalisation = W2;
// fully-coupled (Bishop chi) = W3. Oracle (W1): 1D linear diffusion -- the sudden-head erfc step
// (Carslaw-Jaeger/Bear) + Terzaghi-matched U(Tv).
// Formulation: docs/references/transient-unsaturated-flow-formulation.md.

#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/seepage.hpp>          // Permeability
#include <katai/fem/elements/element_traits.hpp>   // Tri6Element / Tri15Element
#include <katai/materials/water_retention.hpp> // van Genuchten / Mualem (W2)
#include <katai/math/sparse_matrix.hpp>             // SparseMatrixBuilder / CsrMatrix
#include <katai/mesh/mesh.hpp>

namespace katai::core {

struct TransientFlowResult {
    std::vector<double> times;             // time of every step (t=0 included; size = nsteps+1)
    std::vector<Eigen::VectorXd> head;     // full nodal head h per step (size = node_count)
};

// Time-varying Dirichlet head BC: returns the prescribed head at `node` at time `t` (a fixed BC
// ignores t). Only nodes with `is_prescribed[node]` are constrained (the set is fixed, the value
// varies in time).
using HeadBoundary = std::function<double(int node, double t)>;

// Factory solving the SPD system (S+dt.H) factor-once-solve-many (CSR -> (rhs -> x)); e.g.
// PARDISO mtype=2. Empty = dense Eigen LU (the MKL-free reference path -- tests). Structurally
// compatible with consolidation.hpp.
using TransientFlowSolveFactory =
    std::function<std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(const math::CsrMatrix&)>;


// Transient saturated flow (head form). `perm`/`specific_storage` by material id.
// `is_prescribed` (size node_count; empty = all free) marks the Dirichlet head nodes;
// `head_bc(n,t)` is their (possibly time-varying) prescribed head. `initial_head` (size
// node_count) is the t=0 head field. `flux` (optional, size node_count) is a time-constant
// nodal source/flux (well Q>0 = inflow; for an edge flux the caller may accumulate int N'q_bar
// nodally). `solve_factory` empty = dense LU (small meshes/tests). Element type from
// mesh.nodes_per_element. The steady-state limit reduces bit-close to assemble_seepage.
// Definition in kernel/analysis/src/transient_flow.cpp (section 5.2).
TransientFlowResult solve_transient_flow(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& specific_storage, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head, double dt, int nsteps,
    const std::vector<double>* flux = nullptr,
    const TransientFlowSolveFactory& solve_factory = {});

// --- W2: TRANSIENT UNSATURATED (RICHARDS) FLOW -- mass-conservative modified Picard (Celia 1990)
// Richards (head form, mixed/mass-conservative):  dtheta/dt + S_s.S_e.dh/dt = div[k_rel(psi).k.grad h] + q,
//   theta = n.S(psi) (volumetric water content), psi = suction = y-h (h_p=h-y pressure head;
//   h_p<0 => unsaturated), k_rel, S, dS/dpsi from van Genuchten/Mualem (water_retention.hpp).
//   In the saturated region (psi<=0) theta is constant => only the S_s storage remains ->
//   reduces (consistently) to W1. Backward Euler + Celia modified Picard: each time step solves
//   [M_C/dt + H^m].delta = q - H^m.h^m - (1/dt)[int N(theta^m-theta^n) + M_ss(h^m-h^n)]
//   (delta->0 => the mass form / conservative on the mesh),
//   M_C = int N' C N (C = dtheta/dh + S_s.S_e), M_ss = int N' S_s.S_e N, H = int G' k_rel.k G.
//   ALL terms are interpolated at the Gauss point (theta, k_rel, capacitance) -> CONSISTENT
//   mass. The mixed (theta-based) formulation = global mass conservation
//   (1'int N(theta^{n+1}-theta^n) = int dtheta); NO lumping (avoids the tri6 corner-node
//   int N_corner = 0 degeneracy). Element-generic; mixed materials Gauss-based. Verification:
//   saturated Terzaghi + mass conservation + Philip sqrt(t).
struct UnsaturatedFlowResult {
    std::vector<double> times;                  // size nsteps+1
    std::vector<Eigen::VectorXd> head;          // nodal total head h (node_count)
    std::vector<Eigen::VectorXd> saturation;    // nodal degree of saturation S (volume-weighted)
    std::vector<double> water_volume;           // total int theta dV (Gauss-consistent) -- mass check
    std::vector<double> darcy_influx;           // net boundary Darcy inflow that step (pos = inflow)
    int max_picard_iters = 0;
    bool converged = true;
};


// Transient unsaturated (Richards) flow, head form, mass-conservative modified Picard
// (Celia 1990). By material id: `perm` saturated k_sat, `porosity` n, `ss_sat` saturated
// specific storage (= n.gamma_w/Kw), `retention` van Genuchten/Mualem. BCs:
// `is_prescribed` + `head_bc(node,t)`; `initial_head` (node_count). `flux` (optional) nodal
// source. `solve_factory` empty = dense LU (refactorized every Picard iteration -- H is
// nonlinear). In the saturated limit (psi>=0 everywhere) reduces to W1 (consistent storage).
// Element type from the mesh.
// Definition in kernel/analysis/src/transient_flow.cpp (section 5.2).
UnsaturatedFlowResult solve_transient_unsaturated_flow(
    const mesh::Mesh& mesh, const std::vector<Permeability>& perm,
    const std::vector<double>& porosity, const std::vector<double>& ss_sat,
    const std::vector<WaterRetention>& retention, const std::vector<char>& is_prescribed,
    const HeadBoundary& head_bc, const std::vector<double>& initial_head, double dt, int nsteps,
    const std::vector<double>* flux = nullptr, const TransientFlowSolveFactory& solve_factory = {},
    int max_picard = 60, double picard_tol = 1e-7);

}  // namespace katai::core
