#pragma once
// Biot coupled consolidation (time-dependent deformation + pore-water flow) -- PLAXIS 2D 2025.1
// Scientific Manual sec 4 (Biot 1956). Linear-elastic skeleton, fully implicit (alpha=1) time
// integration. Block system (Eq 4-18/19):  [K  L; L' -S*] [dv;dp] = [df; dt.H.p_n],  S* = dt.H + S.
//   K=int B'MB, L=int B'm N, H=int G'(k/gamma_w)G, S=int (n/Kw)N'N.  Element-generic (tri6/tri15);
//   the same N serves both u and p.
// Pore pressure p is the EXCESS (steady part removed), tension-positive. Consolidation (df=0): the
// initial excess pore pressure p0 dissipates and settlement develops through K dv = -L dp ->
// Terzaghi U(Tv). Formulation: docs/references/consolidation-formulation.md.

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <katai/fem/assembly/assembler.hpp>        // expand_to_full
#include <katai/fem/elements/element_traits.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/analysis/seepage.hpp>           // Permeability
#include <katai/materials/material_model.hpp>
#include <katai/math/sparse_matrix.hpp>             // SparseMatrixBuilder / CsrMatrix
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Oedometer (constrained, 1-D) modulus E_oed = E(1−ν)/((1+ν)(1−2ν)) — the confined stiffness that
// governs consolidation (c_v = k·E_oed/γ_w). Degenerate (ν→0.5) returns 0.
inline double oedometer_modulus(double E, double nu) {
    const double d = (1.0 + nu) * (1.0 - 2.0 * nu);
    return d > 1e-30 ? E * (1.0 - nu) / d : 0.0;
}

// Vermeer & Verruijt (1981) critical (minimum) time step for the coupled u-p consolidation solve
// (PLAXIS 2D Sci.Man §4.4; docs/references/consolidation-formulation.md §4):
//     Δt_crit = h²·γ_w / (η·k_y) · (1/E_oed + n/K_w)
// A step far below this leaves the backward-Euler diffusion Δt·H unable to damp the element-local
// pore mode, so the equal-order tri6 (P2-P2) pore field can show a checkerboard at the near-undrained
// instant (docs/validation/lbb-undrained-checkerboard.md). h = element size [m]; η = 40 (tri6) /
// 80 (tri15); k_y vertical permeability [length/time]; n porosity; K_w water bulk modulus. Returns 0
// on degenerate input (so callers treat it as "no constraint").
inline double consolidation_critical_dt(double h, double eta, double Eoed, double k_y, double n,
                                        double Kw, double gamma_w) {
    if (h <= 0.0 || eta <= 0.0 || k_y <= 0.0 || Eoed <= 0.0) return 0.0;
    const double storage = 1.0 / Eoed + (Kw > 0.0 ? n / Kw : 0.0);
    return h * h * gamma_w / (eta * k_y) * storage;
}

// Mean element size h [m] of the mesh (area-based, the same measure the free-surface transition
// width uses): h = sqrt(2*A_mean) from the corner triangle areas. It is the h that enters
// consolidation_critical_dt, and it has ONE definition because two definitions of the same length
// would let the editor warn about a step the solver then chooses differently.
inline double mean_element_size(const mesh::Mesh& mesh) {
    if (mesh.element_count == 0) return 0.0;
    double area_sum = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        area_sum += 0.5 * std::fabs((mesh.x[b] - mesh.x[a]) * (mesh.y[c] - mesh.y[a]) -
                                    (mesh.x[c] - mesh.x[a]) * (mesh.y[b] - mesh.y[a]));
    }
    return std::sqrt(2.0 * area_sum / std::max(1, mesh.element_count));
}

// Interpolation-order factor eta of the Vermeer-Verruijt criterion: 40 for the 6-noded, 80 for the
// 15-noded triangle (PLAXIS 2D Sci.Man sec. 4.4).
inline double consolidation_eta(const mesh::Mesh& mesh) {
    return mesh.nodes_per_element == 15 ? 80.0 : 40.0;
}

struct ConsolidationResult {
    std::vector<double> times;                 // time of every step
    std::vector<Eigen::VectorXd> displacement; // full DOF displacement per step (total_dofs)
    std::vector<Eigen::VectorXd> pore;         // nodal excess pore pressure per step (node_count)
};

// Factory for the fixed-dt coupled solve: the combined symmetric-indefinite system matrix.
// Given A = [K L; L' -(dt.H+S)] IN FULL -- a backend that needs one triangle extracts it
// itself -- returns a closure that solves A x = b. The caller factorizes ONCE and each
// time step only back-substitutes. Empty = dense Eigen LU (the MKL-free reference path --
// test_consolidation). The core stays MKL-independent (the callback is built at the driver).
using ConsolidationSolveFactory =
    std::function<std::function<Eigen::VectorXd(const Eigen::VectorXd&)>(const math::CsrMatrix&)>;


// --- ELASTOPLASTIC (MC/HS) Biot consolidation (PLAXIS Sci.Man sec 4.3) ------------------------
// The nonlinear generalisation of the LE core above: effective stress comes from the constitutive
// return mapping (integrate_point) and the tangent K_T is state-dependent; every time step
// iterates a MONOLITHIC coupled Newton
//   [K_T L; L' -S*][dv;dp] = [r_u; r_p],  r_u = df - df_int(dv) - L.dp,  r_p = dt.H.p_n - (L' dv - S* dp)
// until convergence (the hot solve_nonlinear loop is NOT modified -- the return mapping is called
// in isolation). Pore = an open DOF (NO undrained wrapper); water compressibility sits in
// S = int (n/Kw)N'N (near-incompressible at t=0+ -> the undrained plastic response emerges by
// itself). `initial_state`: committed EFFECTIVE Gauss states (K0/previous phase; size elem*ngp).
// A non-converged step returns `converged=false` (never a hang).
struct ConsolidationPlasticResult {
    ConsolidationResult series;                // times / displacement / pore (same as LE)
    std::vector<GaussState> committed;         // final effective Gauss states (for the phase chain)
    bool converged = true;                     // did every step converge
};


// TIED PORE DEGREES OF FREEDOM ACROSS A SPLIT SEAM. `pore_tie` (optional, node_count) maps a node
// to the node whose pore equation it SHARES, or -1 for a node that owns its own. It exists because
// an interface splits the mesh: the two sides of the joint are two node sets at the same place, and
// what the water does between them is an INPUT, not a consequence of the splitting.
//   fully permeable (the default, and what every model written before the field said) -- one
//     pressure at the joint: the two nodes share an equation, so continuity and the flux balance
//     across the seam both hold by construction rather than by a constraint equation;
//   impermeable -- two pressures: no tie, which is what the bare split already gives.
// PLAXIS states the same two in the same words (Reference Table 5-2, Scientific Manual sec. 3.4).
// Semi-permeable is a third thing -- a conductance dh/R between the two -- and is not this
// parameter; the phase strategy refuses it rather than rounding it to one of the two neighbours.
// null = every node owns its pore equation, and the numbering is BIT-FOR-BIT what it was.
//
// STRUCTURAL STIFFNESS IN A COUPLED SOLVE. `struct_k` (optional, equation_count square) is added
// to the MECHANICAL block of the coupled system: A = [K + K_s   L; Lt  -(dt H + S)]. It arrives as a
// MATRIX rather than as a structural system on purpose -- this core knows soil, water and time, and
// nothing about plates or anchors, exactly as it knows nothing about which linear backend solves it.
// The caller's seam builds it with the shared assembler (assemble_structural_stiffness), which is
// the same element matrices and the same DOF mapping solve_nonlinear uses. It is a CONSTANT
// contribution, so the structural response here is ELASTIC: no anchor yield, no geogrid tension
// cut-off, no plate hinge -- the same declared limit the dynamic branch carries, and the phase
// strategy checks the elastic forces against their capacities so the limit cannot be reached in
// silence. null = no structural elements, and then every arithmetic operation below is BIT-FOR-BIT
// what it was before this parameter existed.
//
// Linear-elastic Biot consolidation. dofs: translations (2/node), finalized (lateral/base BCs).
// materials: LinearElastic; perm: k by material id; gamma_w, kw_over_n=Kw/n (same as the undrained
// wrapper); drained_node[n]=1 -> p=0 at that node (drainage boundary); initial_pore[n] = initial
// excess pore pressure (size node_count); dt fixed time step, nsteps step count. `active`
// (optional, size = element_count; empty = all active) is the staged excavation/fill mask --
// passive elements enter none of the K/L/H/S assemblies. `load_increment` (optional, size =
// equation_count, in the mechanical DOF equation space), when given, is applied at the t=0+ step
// (staged surcharge/fill -> instant undrained excess pore pressure; later steps dissipate it).
// null = only the initial excess pore pressure dissipates (classic Terzaghi). `solve_factory`
// (optional): when given, the sparse PARDISO sym-indefinite path (factor-once-solve-many); empty =
// dense Eigen LU (the MKL-free reference, small meshes). Element type from mesh.nodes_per_element.
// profile (optional): depth gradient, E'(y) per stress point (materials/material_model.hpp
// MaterialProfile). Empty OR uniform() => the constant-E path BIT-FOR-BIT.
// Definition in kernel/analysis/src/consolidation.cpp (section 5.2).
ConsolidationResult solve_consolidation(const mesh::Mesh& mesh, const DofMap& dofs,
                                        const std::vector<MaterialModel>& materials,
                                        const std::vector<Permeability>& perm,
                                        double gamma_w, double kw_over_n,
                                        const std::vector<char>& drained_node,
                                        const std::vector<double>& initial_pore,
                                        double dt, int nsteps,
                                        const std::vector<char>& active = {},
                                        const Eigen::VectorXd* load_increment = nullptr,
                                        const ConsolidationSolveFactory& solve_factory = {},
                                        const std::vector<MaterialProfile>& profile = {},
                                        const math::CsrMatrix* struct_k = nullptr,
                                        const std::vector<int>* pore_tie = nullptr);

// Elastoplastic (MC/HS) Biot consolidation -- monolithic coupled Newton (described above).
// `initial_state`: committed EFFECTIVE Gauss states (K0/previous phase; size elem*ngp; empty =
// zero stress). `solve_factory`: the factory solving the coupled system -- with non-associated
// plasticity K_T is NONSYMMETRIC, so **RealNonsymmetric (mtype=11)** must be supplied
// (factor-once-solve-many; refactorized every iteration). With LE materials it reduces to the
// same answer as solve_consolidation (Terzaghi). Element type from the mesh.
// Definition in kernel/analysis/src/consolidation.cpp (section 5.2).
ConsolidationPlasticResult solve_consolidation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps, const std::vector<char>& active,
    const Eigen::VectorXd* load_increment, const ConsolidationSolveFactory& solve_factory,
    int max_newton = 40, double newton_tol = 1e-6,
    const std::vector<MaterialProfile>& profile = {},
    const math::CsrMatrix* struct_k = nullptr,
    const std::vector<int>* pore_tie = nullptr);

}  // namespace katai::core
