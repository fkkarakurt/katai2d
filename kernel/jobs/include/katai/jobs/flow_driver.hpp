#pragma once
// The groundwater-flow driver (layer 2, katai/jobs; historical name build_flow.hpp).
// Turn a project's groundwater-flow boundary conditions into a steady-state seepage solve
// on the FE mesh (a steady-state groundwater flow calculation). Confined and unconfined regimes are
// both handled by the variable-k_rel free-surface solver with seepage-face active sets (a fully
// saturated domain simply keeps k_rel = 1 everywhere). Returns the nodal head / pore fields for
// display and for coupling into the mechanical solve (assemble_pore_load_from_head +
// assemble_gravity_from_head), plus the boundary discharge and a mass-balance check.
//
// Core formulation + validation: docs/references/seepage-formulation.md (Darcy/Laplace, Charny,
// seepage face); the GUI-path numbers are pinned in tests/test_seepage_gui.cpp.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/seepage.hpp>   // Permeability
#include <katai/mesh/mesh.hpp>
#include <katai/model/project.hpp>

namespace katai::app {

struct FlowResult {
    Eigen::VectorXd head;        // nodal hydraulic head h [m] (size node_count)
    std::vector<double> pore;    // nodal pore pressure u = gamma_w * max(0, h - y) [kPa]
    double discharge = 0.0;      // TOTAL inflow into the domain [m3/day per m]: through
                                 // prescribed-head boundaries and prescribed-flux edges alike.
                                 // AXISYMMETRIC: [m3/day] over the whole ring (2 pi x per radian)
    double balance_err = 0.0;    // |sum of all nodal fluxes| / inflow (mass conservation, ~0)
    // What the hydraulic conditions took OUT of the ground [m3/day per m], and how many well
    // nodes are held at h_min. A well limited by h_min extracts what the ground can give at that
    // head, not what it was asked for -- a difference the run has to be able to state.
    double hydro_discharge = 0.0;
    int hydro_limited = 0;
    // The head on the FAR side of a flow barrier, at the same caller node (equal to `head` where
    // there is no barrier). An impermeable screen makes the head discontinuous across its line,
    // and one number per node cannot say that -- this is the other one.
    std::vector<double> head_far;
    int iterations = 0;
    bool ok = false;
    std::string message;
};

// Water unit weight [kN/m^3] (= build_problem.hpp kGammaWater; distinct name so the two headers
// can be included in either order / independently).
inline constexpr double kFlowGammaWater = 9.81;

// The hydraulic conditions drawn inside the model (`hydros`), resolved onto a mesh -- shared by the
// steady and the transient groundwater flow, so a well or a drain means the same thing in both.
//   level      the nodes that MAY be held at a head, and which way each holds: a drain's nodes
//              from above (a Vacuum drain `always`), an extraction well's nodes at its h_min
//              from below;
//   well_flux  each well's discharge spread along its line, every edge carrying a share in
//              proportion to (the conductivity across the line x the edge's length), so the
//              water comes out of the layers that can deliver it (negative = extraction).
// A line that does not lie on the mesh, or a well of zero length, is refused (ok = false).
struct HydroNode { int node; double level; bool from_above; bool always; std::size_t owner; };
struct HydroResolution {
    std::vector<HydroNode> level;
    std::vector<double> well_flux;
    bool any_well = false, any_drain = false;
    bool ok = true;
    std::string message;
};
HydroResolution resolve_hydros(const model::Project& pr, const katai::mesh::Mesh& mesh,
                               const std::vector<katai::core::Permeability>& perm,
                               const model::Phase* phase);

// Definition in kernel/jobs/src/flow_driver.cpp (section 5.2). `phase` selects which wells and
// drains are switched on (hydraulic conditions are activated per calculation phase); nullptr
// means the project as drawn, with every hydraulic condition active -- which is what the standalone
// "calculate groundwater flow" action asks for and what every caller written before wells existed
// gets, unchanged.
FlowResult solve_groundwater_flow(const model::Project& pr, const katai::mesh::Mesh& mesh,
                                  const model::Phase* phase = nullptr);

}  // namespace katai::app
