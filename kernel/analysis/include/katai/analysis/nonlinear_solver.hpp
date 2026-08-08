#pragma once
// Nonlinear (Newton-Raphson) static solver — P1.1.
//
// Incremental-iterative scheme: the external load is applied in N steps; at each
// step the equilibrium
//   r(u) = f_ext - f_int(u) = 0
// is solved with the consistent tangent K_T. Internal force f_int = ∫ Bᵀσ dA,
// tangent K_T = ∫ Bᵀ D_T B dA (from the material-point integration at each Gauss
// point). Stress/history lives at the Gauss points (committed); once a step
// converges it is committed → ready for path-dependent plasticity.
//
// The linear solver is supplied FROM OUTSIDE via a callback (LinearSolve): K_T·δ = r.
// This keeps katai_core independent of MKL/PARDISO (modularity). For an LE material
// a single step reproduces the linear solution to round-off in ~2 iterations.

#include <functional>
#include <vector>

#include <Eigen/Core>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/embedded_beam.hpp>
#include <katai/fem/elements/geogrid.hpp>
#include <katai/fem/elements/interface.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Given K (CSR, equation_count×equation_count) and r, returns δ.
using LinearSolve =
    std::function<Eigen::VectorXd(const math::CsrMatrix&, const Eigen::VectorXd&)>;

// Analysis kinematics: plane strain (3 strain components) or axisymmetric (4,
// including hoop). In axisym the mesh is interpreted as (r, z).
enum class Kinematics { PlaneStrain, Axisymmetric };

struct NewtonOptions {
    int load_steps = 1;       // (initial) number of steps the external load is split into
    int max_iterations = 50;  // maximum iterations per step
    double tolerance = 1e-9;  // relative convergence threshold on ||r|| / ||f_ext||
    Kinematics kinematics = Kinematics::PlaneStrain;
    // Time interval of the phase [days] — for time-dependent constitutive models
    // (SoftSoilCreep); distributed over the increments in proportion to Δλ (PLAXIS: time
    // advances together with SumMstage). 0 = timeless phase (no creep accumulates; all old
    // callers bit-for-bit).
    double time_interval = 0.0;
};

// Plate (structural wall/beam) embedded in soil — 3-node Timoshenko beam (see
// elements/plate.hpp). Translational DOFs (u_x,u_y) are SHARED with the mesh nodes (assembly
// summation → soil-structure interaction); the rotational DOF (φ) is plate-specific
// (obtained via DofMap::add_extra_dof). nodes: [end A (ξ=−1), end B (ξ=+1), middle (ξ=0)]
// mesh node indices (tri6 edge order). rot_dof: global index of each node's rotational DOF.
// With props.plastic() off the element is ELASTIC (f_int = K·u_total, bit-for-bit the old
// path); with it on, the M-N hinge (Mp/Np diamond, §10).
struct PlateElement {
    std::array<int, 3> nodes;     // mesh nodes [A, B, middle] — for GEOMETRY (coordinates)
    std::array<int, 3> rot_dof;   // global indices of the rotational DOFs
    plate::PlateProps props;
    // Translational DOFs: by default SHARED with the mesh nodes (plate-in-soil). For an
    // embedded WALL (barrier), INDEPENDENT extra DOFs: if trans_dof[2k+c]≥0 that one is used
    // (instead of the mesh node) → the wall is detached from the mesh and tied to both sides
    // with interfaces. -1 = share the mesh node.
    std::array<int, 6> trans_dof = {-1, -1, -1, -1, -1, -1};  // [A_x,A_y, B_x,B_y, mid_x,mid_y]
};

// Anchor — one-directional AXIAL spring (normal force only, NO rotation). PLAXIS MMM Eq 18-1:
// N = (EA/L)·U. The two types share the SAME element technology (PLAXIS Reference Manual):
//  - node-to-node: between two mesh nodes (strut/internal support); L = equivalent length
//    (≤0 ⇒ geometric distance). The spring acts along the geometry direction.
//  - fixed-end: one mesh node + a fixed far end (node_b<0 ⇒ fixed_point); ground anchor (bond
//    zone outside the mesh, assumed fixed). Direction = node_a→fixed_point.
// Uses translational DOFs (shared with mesh nodes) → no DofMap change needed.
// ELASTOPLASTIC: N develops elastically with EA/L, admissible N ∈ [−Fmax_comp, +Fmax_tens]
// (both positive magnitudes; ≤0 ⇒ unbounded = purely elastic, old behaviour bit-for-bit). At
// yield the force plateaus and permanent elongation U_p accumulates (support capacity).
// Math: structural-plate-formulation.md §7a.
struct AnchorElement {
    int node_a = -1;                       // structure-end mesh node
    int node_b = -1;                        // node-to-node: the other mesh node; fixed-end: <0
    Eigen::Vector2d fixed_point{0.0, 0.0};  // fixed-end far-end coordinate (node_b<0)
    double EA = 0.0;                        // axial stiffness E·A
    double L = 0.0;                         // equivalent length (stiffness EA/L); ≤0 ⇒ geometric distance
    double Fmax_tens = -1.0;                // max tensile force; ≤0 ⇒ unbounded (elastic)
    double Fmax_comp = -1.0;                // max compressive force (positive); ≤0 ⇒ unbounded
    // Prestress (lock-off force) per metre of wall, tension-positive. An anchor or strut is
    // almost never installed slack: it is tensioned against the wall, and that force is what
    // holds the excavation before any further movement occurs. The anchor then behaves as an
    // elastic spring FROM that state — N = N0 + (EA/L)·(U − U_p) with U measured from the
    // installation datum the phase chain carries — which is the PLAXIS reading of a prestressed
    // anchor: the lock-off force is applied once, and afterwards the force follows the wall.
    // 0 ⇒ installed slack, i.e. what every KATAI anchor was before this field existed.
    double prestress = 0.0;
};

// Geogrid — 3-node axial membrane (tension-only + optional N_p; see elements/geogrid.hpp).
// Translational DOFs (u_x,u_y) are SHARED with the mesh nodes; NO bending/rotation. nodes:
// tri6 edge [A, B, middle]. Math: structural-plate-formulation.md §8.
struct GeogridElement {
    std::array<int, 3> nodes;     // mesh node indices [A, B, middle]
    geogrid::GeogridProps props;
};

// Interface — zero-thickness soil-structure Coulomb interface (see elements/interface.hpp).
// soil_nodes: mesh nodes (soil side, base DOFs) [A,B,middle]; struct_dof: global indices of
// the extra DOFs of the coincident structure-side nodes [A_ux,A_uy, B_ux,B_uy, mid_ux,mid_uy]
// (via DofMap::add_extra_dof, 2 per node). The structure (plate/load) attaches to struct_dof;
// the interface ties the two sides with k_n,k_s + Coulomb. Math: interface-formulation.md.
struct InterfaceElement {
    std::array<int, 3> soil_nodes;        // mesh nodes [A, B, middle]
    std::array<int, 6> struct_dof;        // structure-side extra DOFs [A_ux,A_uy,B_ux,B_uy,mid_ux,mid_uy]
    iface::InterfaceProps props;
    // Initial normal stress (per Newton-Cotes point, nc_points order q=0:A, 1:middle, 2:B).
    // K0 install / staged carry-over: carries the K0 horizontal stress on a split wall line
    // at Δu_n=0 → no spurious "installation" movement. Default 0 = old behaviour bit-for-bit.
    // (interface-formulation §6.)
    std::array<double, 3> sigma_n0 = {0.0, 0.0, 0.0};
};

// 5-NODE (quartic) plate — sits on a tri15 edge (counterpart of the PLAXIS 15-node soil
// element). Same semantics as the 3-node one: nodes=geometry (mesh nodes), rot_dof=rotational
// extra DOFs, trans_dof=translations (≥0 extra-DOF/embedded wall, −1 share the mesh node).
// 15 DOFs. Math: structural-plate-formulation.md.
struct PlateElement5 {
    std::array<int, 5> nodes;
    std::array<int, 5> rot_dof;
    plate::PlateProps props;
    std::array<int, 10> trans_dof = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};  // [n_i_x, n_i_y]×5
};

// 5-NODE interface (tri15 edge, 5-point Newton-Cotes). Same as the 3-node one; σ_n0 per NC point.
struct InterfaceElement5 {
    std::array<int, 5> soil_nodes;
    std::array<int, 10> struct_dof;   // structure-side extra DOFs [n_i_ux, n_i_uy]×5
    iface::InterfaceProps props;
    std::array<double, 5> sigma_n0 = {0.0, 0.0, 0.0, 0.0, 0.0};
};

// Bundle of structural elements embedded in the soil (plate + anchor + geogrid + interface;
// 3- and 5-node).
struct Structures {
    std::vector<PlateElement> plates;
    std::vector<AnchorElement> anchors;
    std::vector<GeogridElement> geogrids;
    std::vector<InterfaceElement> interfaces;
    std::vector<PlateElement5> plates5;          // tri15-edge structural elements
    std::vector<InterfaceElement5> interfaces5;
    std::vector<ebeam::EmbeddedBeam> embedded_beams;  // embedded beam (pile row, mesh-nonconforming skin)
};

// PARENT STRUCTURAL STATE (Track 1a; first in nonlinear dynamics, now in the static chain
// too). The soil marches each increment from committed σ and therefore inherits the parent's
// stress state via `initial_state`; structural elements, however, are TOTAL-displacement
// formulated (f = f(u_total, plastic state)) → the parent's structural state carries over
// only if (a) the converged displacement DATUM and (b) the committed plastic states are
// provided. Then the Coulomb cap / anchor capacity / geogrid slack is checked on the TOTAL
// effect and an UNCHANGED (nil) phase is a true no-op — without the carry-over, the SumMstage
// imbalance re-ramps the parent's structural tractions every phase (measured: wall M drifted
// 32% across a nil phase). Empty members = start from zero (old behaviour BIT-FOR-BIT).
// A WRONG size does not fall through silently: std::invalid_argument (silently starting from
// zero would present a result that ignores the static preload — unsafe-sided — as normal).
struct StructuralInit {
    // Parent's converged displacement, in EQUATION space (equation_count). Structural
    // elements are evaluated on the sum (u_datum + u_phase); the soil assembly NEVER reads
    // u_free (it only reads the increment du, stress is carried in committed σ) → the datum
    // does not affect the soil. Size 0 = zero.
    Eigen::VectorXd u_datum;
    // Committed structural plastic states (same size/order as NewtonResult). Empty = from zero.
    std::vector<double> anchor_plastic, geogrid_plastic, interface_slip, interface5_slip,
                        embedded_skin_slip, embedded_foot_slip;
    // Plate M-N hinge state [ε_p,κ_p]×Gauss (plate::kPlasticStateSize(5) per element; §10).
    std::vector<double> plate_plastic, plate5_plastic;
};

// Unpack the init_struct plastic seeds into the target vectors — the static Newton, the
// nonlinear dynamics and structural_internal_force use the SAME rule (single source): empty
// member → n zeros; full size → copy; wrong size → std::invalid_argument. Returns whether at
// least one member was seeded.
bool seed_structural_state(const Structures& structures, const StructuralInit& init_struct,
                           std::vector<double>& anchor, std::vector<double>& geogrid,
                           std::vector<double>& iface, std::vector<double>& iface5,
                           std::vector<double>& eskin, std::vector<double>& efoot,
                           std::vector<double>& plate_p, std::vector<double>& plate5_p);

// STRUCTURAL internal force f_s0 = the internal force of the structural elements under
// (u_datum + seeded plastic state), in EQUATION space. The baseline (constant_force) of a
// chained phase MUST INCLUDE this so that residual(0) = 0 holds by the parent's OWN
// equilibrium and the ramp = f − B carries only the true configuration change. Computed with
// the shared assembler (ARITHMETICALLY IDENTICAL to what solve_nonlinear will use; NOT a
// separate formula): with committed Gauss = ZERO and du = 0 the soil loop produces exactly
// zero (Δε=0 → σ_trial=committed=0), leaving the pure structural force.
Eigen::VectorXd structural_internal_force(const mesh::Mesh& mesh, const DofMap& dofs,
                                          const std::vector<MaterialModel>& materials,
                                          const Structures& structures,
                                          const StructuralInit& init_struct,
                                          Kinematics kinematics = Kinematics::PlaneStrain);

struct NewtonResult {
    Eigen::VectorXd displacement;            // total_dofs (fixed DOFs = 0)
    std::vector<GaussState> gauss_states;    // element_count * 3 (converged)
    // Structural path-dependent state (elastoplastic). anchor_plastic: committed plastic
    // elongation U_p per anchor. geogrid_plastic: committed plastic axial ε_p per geogrid ×
    // 2 Gauss. In multi-phase runs it can be fed back into the next phase (like gauss_states).
    std::vector<double> anchor_plastic;      // anchors.size()
    std::vector<double> geogrid_plastic;     // geogrids.size() * 2
    std::vector<double> interface_slip;      // interfaces.size() * 3 (Newton-Cotes node pairs)
    std::vector<double> interface5_slip;     // interfaces5.size() * 5 (5-node interface)
    std::vector<double> embedded_skin_slip;  // all embedded beam skin points (axial plastic slip)
    std::vector<double> embedded_foot_slip;  // embedded_beams.size() (foot axial plastic slip)
    std::vector<double> plate_plastic;       // plates.size() * 6  ([ε_p,κ_p]×3 Gauss; M-N hinge)
    std::vector<double> plate5_plastic;      // plates5.size() * 10 ([ε_p,κ_p]×5 Gauss)
    bool converged = false;
    int total_iterations = 0;                // sum over all steps
    // Highest fraction of f_ext that was successfully equilibrated. Equals 1.0
    // on full convergence; on failure it is the last converged load level, which
    // for perfect plasticity approaches the collapse (limit) load from below ->
    // the basis of incremental limit analysis (e.g. Prandtl bearing capacity).
    double load_factor = 0.0;
    // Increments abandoned because the linear solver refused to answer: the tangent
    // was singular along a collapse mechanism, so no Newton direction exists. Normal
    // at the limit load and recovered from by cutting the increment back, but counted
    // rather than swallowed -- a solve that could not be performed is information
    // about the model, and an analysis that reports a load factor after several of
    // them is reporting a collapse it should be able to name.
    int refused_solves = 0;
    // Wall-clock breakdown of the computation (seconds) + call counters. The counterpart of
    // PLAXIS's calculation-time report; the base measurement for performance studies.
    // Instrumentation is at iteration granularity (one chrono call per iteration) → the cost
    // is negligible, always on.
    struct Timings {
        double assemble_tangent = 0.0;   // K_T + f_int assembly (1 per iteration)
        double assemble_residual = 0.0;  // f_int only (line-search trials)
        double csr_build = 0.0;          // COO→CSR compile (SparseMatrixBuilder::build)
        double linear_solve = 0.0;       // LinearSolve callback (factorization + solve)
        double total = 0.0;              // solve_nonlinear total
        int n_tangent = 0;               // number of tangent assemblies
        int n_residual = 0;              // number of residual-only assemblies (line search)
        int n_solve = 0;                 // number of linear solves
    } timings;
};

// f_ext: free-DOF (equation_count) external load vector (e.g. the output of
// assemble_surface_traction). materials is indexed by the mesh.element_material ids.
//
// initial_state: optional committed INITIAL Gauss state (pre-stress — the K0 procedure or a
// previous phase). If empty, starts from zero stress. If given, its size must be
// element_count*kGaussCount; when the pre-stress balances the full external load (geostatic
// K0 + gravity), load_steps=1 yields zero displacement in the first phase (see
// initial-stress-k0.md).
// active_element: optional element activity mask (empty = all active). Passive elements do
// not enter the assembly (staged construction — excavation/fill); their committed stresses
// are preserved. Nodes orphaned by passive elements must be fixed in the DofMap (see
// staged_construction.hpp fix_inactive_nodes), otherwise the system is singular.
// structures: optional structural elements embedded in the soil (plate + anchor).
// Translational DOFs are shared with mesh nodes; plate rotational DOFs must have been added
// to the DofMap via add_extra_dof.
//
// prescribed_displacement: optional NONZERO Dirichlet (prescribed displacement) — size
// total_dofs, meaningful only on FIXED DOFs (free-DOF values are ignored). Ramped together
// with the external load from 0 to the full value (load_factor·ū); since the internal force
// is computed from the total displacement (including ū), the K_fp·ū term enters
// automatically. For wall/boundary translation, prescribed settlement, active/passive earth
// pressure (Rankine) experiments. result.displacement contains ū at the fixed DOFs.
// constant_force: optional CONSTANT (non-ramped) free-DOF external load (equation_count).
// Applied IN FULL at every load step (target = constant_force + λ·f_ext). To hold the
// initial equilibrium (geostatic gravity, balancing the pre-stress) constant and ramp only
// f_ext/the prescribed displacement: this way the prescribed displacement is applied
// incrementally in plasticity while gravity stays constant (no spurious imbalance on
// cutback). Empty = 0.
// profile: optional depth gradient (by material id; see materials/material_model.hpp
// MaterialProfile). Empty OR all uniform() ⇒ the old constant-E/c path BIT-FOR-BIT. If
// given, E' and c' are evaluated PER STRESS (Gauss) POINT (as PLAXIS does).
// init_struct: parent STRUCTURAL datum + plastic state (static generalization of Track 1a;
// empty = old behaviour bit-for-bit). If given, constant_force MUST INCLUDE
// structural_internal_force(init_struct) — otherwise residual(0) ≠ 0 and the first
// increment ramps a spurious imbalance.
NewtonResult solve_nonlinear(const mesh::Mesh& mesh, const DofMap& dofs,
                             const std::vector<MaterialModel>& materials,
                             const Eigen::VectorXd& f_ext,
                             const LinearSolve& linear_solve,
                             const NewtonOptions& options = {},
                             const std::vector<GaussState>& initial_state = {},
                             const std::vector<char>& active_element = {},
                             const Structures& structures = {},
                             const Eigen::VectorXd& prescribed_displacement = {},
                             const Eigen::VectorXd& constant_force = {},
                             const std::vector<MaterialProfile>& profile = {},
                             const StructuralInit& init_struct = {});

} // namespace katai::core
