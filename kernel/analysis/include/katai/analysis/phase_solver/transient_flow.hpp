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
    double gamma_water = kGammaWater;   // the unit weight of water [kN/m^3]
    std::vector<TransientFlowMaterial> materials;  // by material id
    std::vector<FlowEdge> flow_edges;              // B4 vocabulary; Closed edges included
    std::vector<char> active;                      // element activity; empty = everything active
    std::vector<double> fallback_head;             // per node: initial head where no Head BC prescribes one
    // Prescribed boundary flux, node by node (the integral of q N_i along each Flux edge, inflow
    // positive, accumulated by the caller exactly as the steady flow solve does); empty = none.
    // A flux edge used to reach this solver as a Closed one, so a rainfall infiltration was
    // silently zero.
    std::vector<double> nodal_flux;
    // Seepage faces, drains and wells' h_min (HeadLimit), node by node on this mesh. The nodes of
    // the Head edges are not listed here: the strategy makes them limits itself where the water
    // can stand below them.
    std::vector<HeadLimit> limits;
    // Semi-permeable walls and interfaces (LeakySeam), on this solver's own node numbering.
    std::vector<LeakySeam> seams;
    // The head the phase starts from, node by node: the head the phase before it ended with.
    // Empty = the prescribed heads on their nodes and the fallback elsewhere (the water table).
    std::vector<double> initial_head;
    // The phase's own water level, node by node (the elevation of that water above each node);
    // non-empty = the Head edges take it as their head. The water outside the ground is what a
    // Head edge stands for, and a phase that states its water level states it: a drawdown is a
    // phase whose level is lower. Empty = the heads the edges carry in the project.
    std::vector<double> head_level;
    // The faces an excavation exposes, node by node: the head of the phase's water there (NaN =
    // not such a face). They are water boundaries like the Head edges -- held at that head where
    // the water stands above them, seepage faces where it does not -- which is how the steady
    // flow of a staged phase treats them; a face the drawing gave no flow condition would
    // otherwise be closed, and a pit pumped dry would take no water from the ground around it.
    std::vector<double> face_head;
    // How the Head edges reach the heads of this phase from the heads it starts from: at its first
    // instant (false -- an instantaneous change), or linearly over its duration (true -- a
    // reservoir lowered at a rate). Read only with an initial head.
    bool ramp_heads = false;
    // The sparse direct solver of each Picard iteration (the matrix H + M_C/dt is SPD), injected by
    // the caller as every phase strategy's is. Empty = a dense LU, which is for the small meshes of
    // the unit tests only: it was what every model got until 2026-10-23, and it grows with the cube
    // of the node count.
    TransientFlowSolveFactory solve_factory;
    double duration_day = 1.0;                     // <= 0 falls back to 1 day
    bool axisymmetric = false;                     // x is the radius: the flow per radian
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
        ss[mi] = poros[mi] * in.gamma_water / kWaterBulk;   // specific storage S_s = n*gamma_w/Kw
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
    // THE GROUND OF THIS PHASE. An excavated region holds no soil for the water to move through
    // or be stored in, so the flow is solved on the active elements alone -- as the steady flow of
    // a staged phase is. Until 2026-10-23 every element took part, and the soil an excavation had
    // removed went on conducting and storing water in a transient flow after it. A node no active
    // element touches belongs to no equation; it is held at the head it starts from.
    katai::mesh::Mesh amesh;
    const bool partial = !in.active.empty() &&
                         std::find(in.active.begin(), in.active.end(), 0) != in.active.end();
    if (partial) {
        amesh = mesh;
        amesh.connectivity.clear();
        amesh.element_material.clear();
        amesh.element_count = 0;
        const int npe = mesh.nodes_per_element;
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!in.active[e]) continue;
            for (int k = 0; k < npe; ++k) amesh.connectivity.push_back(mesh.node_of(e, k));
            amesh.element_material.push_back(mesh.element_material[e]);
            ++amesh.element_count;
        }
    }
    const katai::mesh::Mesh& fmesh = partial ? amesh : mesh;
    std::vector<char> touched(mesh.node_count, 0);
    for (int e = 0; e < fmesh.element_count; ++e)
        for (int k = 0; k < fmesh.nodes_per_element; ++k) touched[fmesh.node_of(e, k)] = 1;
    std::vector<char> is_presc; std::vector<double> head_val;
    flow_head_nodes(in.flow_edges, mesh, is_presc, head_val);
    for (int n = 0; n < mesh.node_count; ++n)
        if (!touched[n]) is_presc[n] = 0;   // set below, once the starting head is known
    if ((int)in.head_level.size() == mesh.node_count)
        for (int n = 0; n < mesh.node_count; ++n)
            if (is_presc[n]) head_val[n] = in.head_level[n];
    if ((int)in.face_head.size() == mesh.node_count)
        for (int n = 0; n < mesh.node_count; ++n)
            if (touched[n] && !is_presc[n] && std::isfinite(in.face_head[n])) {
                is_presc[n] = 1;
                head_val[n] = in.face_head[n];
            }
    bool any = false; for (char c : is_presc) if (c) { any = true; break; }
    if (!any && in.limits.empty()) {
        R.message = "Transient flow needs at least one prescribed-head boundary: right-click a "
                    "soil edge in Flow conditions and set a head before a transient flow phase.";
        return false;
    }
    const double duration = in.duration_day > 0.0 ? in.duration_day : 1.0;
    const int nsteps = std::clamp(in.time_steps, 1, 2000);
    const double dt = duration / nsteps;
    // The state the phase starts from: the head the phase before it ended with, or -- with none
    // -- the prescribed heads on their nodes and the water table elsewhere.
    const bool from_parent = (int)in.initial_head.size() == mesh.node_count;
    std::vector<double> h0(mesh.node_count);
    for (int n = 0; n < mesh.node_count; ++n)
        h0[n] = from_parent ? in.initial_head[n] : in.fallback_head[n];
    const bool ramp = from_parent && in.ramp_heads;
    HeadBoundary hb = [&head_val, &h0, ramp, duration](int n, double t) {
        if (!ramp) return head_val[n];
        const double f = std::clamp(t / duration, 0.0, 1.0);
        return h0[n] + f * (head_val[n] - h0[n]);
    };
    // THE WATER BOUNDARIES. A Head edge stands for the water outside the ground, so where that
    // water is below a node of the edge the node is not under it: it is a seepage face, open
    // while water leaves the ground there and closed otherwise. Held at the head instead -- as
    // the edge once was over its whole length -- it would draw the ground's water out by
    // suction, through a face that is open to the air. A node the water stands above for the
    // whole phase keeps its prescribed head; one it falls below (or rises above) during the phase
    // follows it.
    std::vector<HeadLimit> limits;
    for (const HeadLimit& L : in.limits)
        if (L.node >= 0 && L.node < mesh.node_count && touched[L.node]) limits.push_back(L);
    for (int n = 0; n < mesh.node_count; ++n) {
        if (!is_presc[n]) continue;
        const double lowest = std::min(hb(n, dt), hb(n, duration));
        if (lowest >= mesh.y[n]) continue;   // under the water for the whole phase
        is_presc[n] = 0;
        HeadLimit L; L.node = n; L.follows_bc = true;
        limits.push_back(L);
    }
    // With no head to start from, a node held at its head for the whole phase starts there, as
    // it always has; the nodes the water leaves start from the water table like the ground.
    if (!from_parent)
        for (int n = 0; n < mesh.node_count; ++n)
            if (is_presc[n]) h0[n] = head_val[n];
    for (int n = 0; n < mesh.node_count; ++n)
        if (!touched[n]) { is_presc[n] = 1; head_val[n] = h0[n]; }   // outside the ground: inert
    const bool has_flux = (int)in.nodal_flux.size() == mesh.node_count;
    // An inert node is held at its starting head whatever the ramp says.
    const HeadBoundary hb_ground = [&hb, &touched, &h0](int n, double t) {
        return touched[n] ? hb(n, t) : h0[n];
    };
    // A seam with a side the phase has excavated joins ground to no ground: left out.
    std::vector<LeakySeam> seams;
    for (const LeakySeam& sm : in.seams)
        if (touched[sm.a0] && touched[sm.a1] && touched[sm.b0] && touched[sm.b1]) seams.push_back(sm);
    const UnsaturatedFlowResult fr = solve_transient_unsaturated_flow(
        fmesh, cperm, poros, ss, ret, is_presc, hb_ground, h0, dt, nsteps, has_flux ? &in.nodal_flux : nullptr,
        in.solve_factory, 60, 1e-7, limits.empty() ? nullptr : &limits, seams.empty() ? nullptr : &seams,
        in.axisymmetric);
    if (!fr.converged || fr.head.empty()) {
        R.message = "Transient flow did not converge (the SWCC may be very steep -- use more, "
                    "smaller time steps).";
        return false;
    }
    if (!fr.limits_settled) {
        R.message = "The seepage faces, wells and drains of this transient flow phase did not "
                    "settle into a set in every time step: their one-sided conditions are cycling. "
                    "Use more, smaller time steps, and check that a well's h_min is below the head "
                    "the ground can supply it at.";
        return false;
    }
    R.disp = Eigen::VectorXd::Zero(mesh.node_count * 2);
    R.stress.stress.assign(mesh.node_count, Eigen::Vector3d::Zero());
    for (size_t s = 0; s < fr.times.size(); ++s) {
        double pmax = 0.0;
        for (int n = 0; n < mesh.node_count; ++n)
            pmax = std::fmax(pmax, in.gamma_water * std::fmax(0.0, fr.head[s][n] - mesh.y[n]));
        R.consol_time.push_back(fr.times[s]);
        R.consol_settlement.push_back(0.0);
        R.consol_excess_pore.push_back(pmax);
    }
    { const Eigen::VectorXd& Sf = fr.saturation.back();
      R.saturation.assign(Sf.data(), Sf.data() + Sf.size()); }
    R.pore.assign(mesh.node_count, 0.0);
    for (int n = 0; n < mesh.node_count; ++n)
        R.pore[n] = in.gamma_water * std::fmax(0.0, fr.head.back()[n] - mesh.y[n]);
    R.flow_head.assign(fr.head.back().data(), fr.head.back().data() + mesh.node_count);
    R.flow_head_transient = true;
    R.load_factor = 1.0; R.iterations = nsteps; R.max_disp = 0.0;
    R.active = in.active;
    R.ok = true;
    R.message = "Transient flow solved over " + std::to_string(duration) + " day(s).";
    return true;
}

} // namespace katai::core
