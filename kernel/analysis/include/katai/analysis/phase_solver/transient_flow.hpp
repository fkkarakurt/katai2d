#pragma once
// Transient groundwater-flow phase strategy (Stage B9, pilot). The first phase
// of the staged-analysis driver to live in the engine as a strategy: a plain
// function over neutral inputs that fills the engine-owned SolveResult. This
// file fixes the contract the remaining phase strategies follow -- and it is
// the call target the headless execution layer (katai/jobs) will use, so its
// vocabulary is written for every front end, not for the GUI.
//
// Contract rules (B9 plan, roadmap section 4.6):
//   - inputs are physical quantities resolved once at the caller's seam;
//   - a material NAME may enter: it is the diagnostic label every front end
//     needs for an honest refusal, not schema structure;
//   - refusal messages are engine-owned and byte-identical to what the phase
//     editor has always shown;
//   - the strategy fills everything about the outcome except the mesh handoff,
//     which stays with the caller (it owns the working mesh's lifetime).

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/constants.hpp>
#include <katai/analysis/hydraulic_boundary.hpp>
#include <katai/analysis/results.hpp>
#include <katai/analysis/seepage.hpp>
#include <katai/analysis/transient_flow.hpp>
#include <katai/materials/water_retention.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// One material's flow description, resolved at the caller's seam.
struct TransientFlowMaterial {
    std::string name;              // diagnostic label for honest refusals
    double kx = 1.0, ky = 1.0;     // permeability [m/day]
    WaterRetention retention;      // van Genuchten / Mualem SWCC
    double porosity = 0.3;         // n = e / (1 + e)
    bool nonporous = false;        // impermeable barrier, holds no water (zero storage)
};

// The phase's neutral configuration. Defaults reproduce the driver's historical
// no-configuration behaviour; the >0 and clamp rules live here so every caller
// gets the same ones.
struct TransientFlowPhase {
    std::vector<TransientFlowMaterial> materials;  // by material id
    std::vector<FlowEdge> flow_edges;              // B4 vocabulary; Closed edges included
    std::vector<char> active;                      // element activity; empty = everything active
    std::vector<double> fallback_head;             // per node: initial head where no Head BC prescribes one
    double duration_day = 1.0;                     // <= 0 falls back to 1 day
    int time_steps = 25;                           // clamped to [1, 2000]
};

// Solve the phase and fill R (flow-only fields: the pore/head evolution, the
// final saturation, zero displacement). No deformation: the pore/head field
// evolves under storage + the flow BCs; unsaturated van Genuchten/Mualem
// reduces to saturated where psi <= 0. Returns R.ok.
inline bool solve_transient_flow_phase(const katai::mesh::Mesh& mesh,
                                       const TransientFlowPhase& in, SolveResult& R) {
    const size_t nmat = in.materials.size();
    std::vector<Permeability> cperm(nmat, {1.0, 1.0});
    std::vector<WaterRetention> ret(nmat);
    std::vector<double> poros(nmat, 0.3), ss(nmat, 0.0);
    std::vector<char> used(nmat, 0);
    for (int e = 0; e < mesh.element_count; ++e)
        if (in.active.empty() || in.active[e]) used[mesh.element_material[e]] = 1;
    constexpr double kWaterBulk = 2.0e6;   // water bulk modulus [kPa] (matches the consolidation solve)
    for (size_t mi = 0; mi < nmat; ++mi) {
        const TransientFlowMaterial& Mt = in.materials[mi];
        cperm[mi] = {Mt.kx, Mt.ky};
        ret[mi] = Mt.retention;
        poros[mi] = Mt.porosity;
        ss[mi] = poros[mi] * kGammaWater / kWaterBulk;   // specific storage S_s = n*gamma_w/Kw
        // NonPorous: impermeable barrier with no storage (holds no water) -- the same rule
        // as the steady flow solve.
        if (Mt.nonporous) {
            cperm[mi] = {1e-8, 1e-8};
            ss[mi] = 0.0;
            continue;
        }
        if (used[mi] && (cperm[mi].kx <= 0.0 || cperm[mi].ky <= 0.0)) {
            R.message = "Material '" + Mt.name + "' has no permeability -- set kx, ky (> 0) on "
                        "its Groundwater tab before a transient flow phase.";
            return false;
        }
    }
    std::vector<char> is_presc; std::vector<double> head_val;
    flow_head_nodes(in.flow_edges, mesh, is_presc, head_val);
    bool any = false; for (char c : is_presc) if (c) { any = true; break; }
    if (!any) {
        R.message = "Transient flow needs at least one prescribed-head boundary: right-click a "
                    "soil edge in Flow conditions and set a head before a transient flow phase.";
        return false;
    }
    std::vector<double> h0(mesh.node_count);
    for (int n = 0; n < mesh.node_count; ++n)
        h0[n] = is_presc[n] ? head_val[n] : in.fallback_head[n];
    const double duration = in.duration_day > 0.0 ? in.duration_day : 1.0;
    const int nsteps = std::clamp(in.time_steps, 1, 2000);
    const double dt = duration / nsteps;
    HeadBoundary hb = [&head_val](int n, double) { return head_val[n]; };
    const UnsaturatedFlowResult fr = solve_transient_unsaturated_flow(
        mesh, cperm, poros, ss, ret, is_presc, hb, h0, dt, nsteps);
    if (!fr.converged || fr.head.empty()) {
        R.message = "Transient flow did not converge (the SWCC may be very steep -- use more, "
                    "smaller time steps).";
        return false;
    }
    R.disp = Eigen::VectorXd::Zero(mesh.node_count * 2);
    R.stress.stress.assign(mesh.node_count, Eigen::Vector3d::Zero());
    for (size_t s = 0; s < fr.times.size(); ++s) {
        double pmax = 0.0;
        for (int n = 0; n < mesh.node_count; ++n)
            pmax = std::fmax(pmax, kGammaWater * std::fmax(0.0, fr.head[s][n] - mesh.y[n]));
        R.consol_time.push_back(fr.times[s]);
        R.consol_settlement.push_back(0.0);
        R.consol_excess_pore.push_back(pmax);
    }
    { const Eigen::VectorXd& Sf = fr.saturation.back();
      R.saturation.assign(Sf.data(), Sf.data() + Sf.size()); }
    R.pore.assign(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        R.pore[n] = kGammaWater * std::fmax(0.0, fr.head.back()[n] - mesh.y[n]);
    R.load_factor = 1.0; R.iterations = nsteps; R.max_disp = 0.0;
    R.active = in.active;
    R.ok = true;
    R.message = "Transient flow solved over " + std::to_string(duration) + " day(s).";
    return true;
}

} // namespace katai::core
