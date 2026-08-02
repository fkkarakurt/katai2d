#pragma once
// Biot coupled consolidation (time-dependent deformation + pore-water flow) -- PLAXIS 2D 2025.1
// Scientific Manual sec 4 (Biot 1956). Linear-elastic skeleton, fully implicit (alpha=1) time
// integration. Block system (Eq 4-18/19):  [K  L; L' -S*] [dv;dp] = [df; dt.H.p_n],  S* = dt.H + S.
//   K=int B'MB, L=int B'm N, H=int G'(k/gamma_w)G, S=int (n/Kw)N'N.  Element-generic (tri6/tri15);
//   the same N serves both u and p.
// Pore pressure p is the EXCESS (steady part removed), tension-positive. Consolidation (df=0): the
// initial excess pore pressure p0 dissipates and settlement develops through K dv = -L dp ->
// Terzaghi U(Tv). Formulation: docs/references/consolidation-formulation.md.

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
                                        const std::vector<MaterialProfile>& profile = {});

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
    const std::vector<MaterialProfile>& profile = {});

}  // namespace katai::core
