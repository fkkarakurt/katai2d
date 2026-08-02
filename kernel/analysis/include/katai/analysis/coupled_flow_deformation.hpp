#pragma once
// Fully-coupled flow-deformation (PLAXIS 2D's most general analysis; Ref sec 7.4.4, Sci.Man sec 3
// Eq 3-8 + sec 4). The UNSATURATED generalisation of Biot coupled consolidation
// (consolidation.hpp) -- W3. Linear-elastic skeleton, monolithic saddle point +
// time-varying hydraulic BCs + van Genuchten/Mualem unsaturated retention.
//
// Three generalisations (in the saturated limit ALL vanish -> reduces to consolidation
// BIT-FOR-BIT):
//   1. SATURATION COUPLING: continuity S_eff.m'de/dt (Eq 3-8) -> the L block scales by chi=S_eff.
//   2. BISHOP EFFECTIVE STRESS: sigma = sigma' + chi.p.m (chi=S_eff); saturated chi=1 => classic
//      Terzaghi. On an LE skeleton chi affects only the coupling L (the sigma'=D.eps skeleton
//      response is independent of p).
//   3. UNSATURATED STORAGE + PERMEABILITY: S = int N'(S.n/Kw + n.dS/dp_w)N (compressibility +
//      moisture capacity), H = int G'(k_rel.k/gamma_w)G (Mualem k_rel). Suction psi=-p/gamma_w
//      (p>=0 => psi<=0 saturated S_eff=1).
//
// Unsaturated -> coefficients depend on p => NONLINEAR: every time step iterates PICARD (lag the
// coefficients, re-assemble+solve). The saddle A=[K  L_chi; L_chi'  -(dt.H_kr+S_st)] is
// SYMMETRIC-INDEFINITE (PARDISO mtype=-2; empty = dense LU). In the saturated region the
// coefficients are p-independent (S_eff=1, k_rel=1, moisture=0) => Picard converges in 1
// iteration and A = consolidation's A => bit-for-bit regression. Formulation:
// docs/references/transient-unsaturated-flow-formulation.md (W3); consolidation-formulation.md
// (the Biot core), effective-stress-formulation.md (Bishop).

#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/consolidation.hpp>     // ConsolidationResult, ConsolidationSolveFactory
#include <katai/analysis/seepage.hpp>           // Permeability
#include <katai/fem/assembly/assembler.hpp>         // expand_to_full
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/element_traits.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/materials/water_retention.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

struct CoupledFlowResult {
    ConsolidationResult series;                 // times / displacement / pore (consolidation-compatible)
    std::vector<Eigen::VectorXd> saturation;    // nodal degree of saturation per step, S = S_res+(1-S_res).S_e
    int max_picard_iters = 0;
    bool converged = true;
};


// Fully-coupled flow-deformation (W3). LE skeleton; `retention`/`porosity` by material id
// (unsaturated van Genuchten/Mualem + Bishop chi=S_eff). The other parameters are IDENTICAL to
// solve_consolidation (easy regression): `drained_node[n]=1` -> p=0 boundary; `initial_pore`
// initial pore (tension/suction: psi=-p/gamma_w; p>=0 saturated); `load_increment` applied at
// t=0+; `active` staged mask; `solve_factory` sym-indefinite (mtype=-2) factor-once-solve-many,
// empty = dense LU (reference). In the saturated limit reduces to solve_consolidation
// BIT-FOR-BIT. Definition in kernel/analysis/src/coupled_flow_deformation.cpp (section 5.2).
CoupledFlowResult solve_coupled_flow_deformation(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<double>& initial_pore,
    double dt, int nsteps, const std::vector<char>& active = {},
    const Eigen::VectorXd* load_increment = nullptr,
    const ConsolidationSolveFactory& solve_factory = {},
    int max_picard = 50, double picard_tol = 1e-10,
    const std::vector<MaterialProfile>& profile = {});

// --- ELASTOPLASTIC (MC/HS) FULLY-COUPLED FLOW-DEFORMATION (W3 follow-up) ------------------------
// The UNION of W3 (LE coupled, above) and elastoplastic consolidation (consolidation_plastic) =
// PLAXIS 2D's MOST GENERAL analysis: unsaturated (van Genuchten/Mualem + Bishop chi=S_eff) Biot
// consolidation on an MC/HS plastic skeleton. The skeleton sigma' comes from the return mapping
// (integrate_point) and the tangent K_T is state-dependent; the flow coefficients (S_eff, k_rel,
// storage) depend on pore. Every time step iterates a MONOLITHIC Newton-Picard:
//   [K_T  L_chi;  L_chi'  -(dt.H_kr+S_st)] [dv; dp] = r,
//   r_u = df - (f_int(dv) - Bbase) - L_chi.dp,  r_p = dt.H_kr.p_n - (L_chi' dv - (dt.H_kr+S_st) dp).
// K_T comes from the plastic return mapping (Newton); L_chi/H_kr/S_st are evaluated at the latest
// pore estimate (Picard lag). Two limits verify it: (a) SATURATED (retention=sat) -> reduces to
// solve_consolidation_plastic BIT-FOR-BIT (S_eff=k_rel=1, dS/dp_w=0 => same L/H/S blocks + same
// Newton); (b) LE skeleton -> converges to the SAME fixed point as solve_coupled_flow_deformation
// (W3). K_T is NONSYMMETRIC under non-associated flow => solve_factory RealNonsymmetric
// (mtype=11). Formulation: transient-unsaturated-flow-formulation.md (W3);
// consolidation-formulation.md sec 4.3 (the elastoplastic Biot core),
// effective-stress-formulation.md.
struct CoupledFlowPlasticResult {
    ConsolidationResult series;                 // times / displacement / pore (consolidation-compatible)
    std::vector<Eigen::VectorXd> saturation;    // nodal degree of saturation per step, S = S_res+(1-S_res).S_e
    std::vector<GaussState> committed;          // final effective Gauss states (for the phase chain)
    int max_newton_iters = 0;
    bool converged = true;
};


// Elastoplastic (MC/HS) fully-coupled flow-deformation (W3 follow-up). Unsaturated van
// Genuchten/Mualem + Bishop chi=S_eff + plastic skeleton (integrate_point return mapping).
// Parameters are the union of solve_coupled_flow_deformation (W3) and
// solve_consolidation_plastic: `retention`/`porosity` by material id (unsaturated SWCC);
// `initial_state` committed EFFECTIVE Gauss states (K0/previous phase; empty = zero stress);
// `load_increment` applied at t=0+ (staged surcharge -> instant undrained response); `active`
// staged mask; `solve_factory` must be **RealNonsymmetric (mtype=11)** because K_T is
// NONSYMMETRIC under non-associated flow. With LE materials it reduces to
// solve_coupled_flow_deformation (W3), in the saturated limit to solve_consolidation_plastic.
// Definition in kernel/analysis/src/coupled_flow_deformation.cpp (section 5.2).
CoupledFlowPlasticResult solve_coupled_flow_deformation_plastic(
    const mesh::Mesh& mesh, const DofMap& dofs, const std::vector<MaterialModel>& materials,
    const std::vector<Permeability>& perm, const std::vector<WaterRetention>& retention,
    const std::vector<double>& porosity, double gamma_w, double kw_over_n,
    const std::vector<char>& drained_node, const std::vector<GaussState>& initial_state,
    const std::vector<double>& initial_pore, double dt, int nsteps,
    const std::vector<char>& active = {}, const Eigen::VectorXd* load_increment = nullptr,
    const ConsolidationSolveFactory& solve_factory = {},
    int max_newton = 40, double newton_tol = 1e-6,
    const std::vector<MaterialProfile>& profile = {});

}  // namespace katai::core
