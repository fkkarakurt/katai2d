#pragma once
// Per-line structural force diagrams (Stage B9). One drawn structure line -- a
// plate, an embedded wall, an anchor, a geogrid, an embedded beam -- owns a
// contiguous slice of a Structures vector; a DiagSpec remembers that slice so
// a named N/Q/M diagram along the structure can be produced from any
// full global-DOF displacement field. IfaceDiag is the same bookkeeping for
// Coulomb interfaces (tau / sigma_n / slip per line).
//
// The evaluators are shared by the static result tail (converged solution +
// committed plastic state) and the Dynamic phase's per-step seismic envelope,
// so the two paths can never report different quantities from the same
// displacements. `elastic` selects the ELASTIC branch of every recovery: the
// linear dynamic system solves anchors uncapped, geogrids at full EA and
// plates elastic, so its report must too -- a capped / tension-cut diagram
// would silently under-state the demand the linear solve actually carried
// (the D6b rule).

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/results.hpp>            // StructForce, InterfaceResult
#include <katai/analysis/structural_forces.hpp>  // validated per-element recoveries
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// Which contiguous slice of a Structures vector one drawn line owns.
// kind: 0 plate / embedded wall (tri6) · 1 anchor · 2 geogrid ·
// 3 embedded beam · 5 embedded wall (tri15).
// anchor_spacing: the solver's anchor force is per-metre of wall (EA/Ls); the
// REPORT is per-anchor [kN] (the input convention -- the user enters Fmax per
// anchor, so a per-metre report would be a xLs trap in capacity comparisons;
// audit finding). 1.0 for every other kind.
struct DiagSpec { int kind; std::string name; size_t begin, end; double anchor_spacing = 1.0; };

// Same bookkeeping for interfaces: which contiguous slice of
// structures.interfaces (order 6) or structures.interfaces5 (order 15) a wall
// or a standalone line owns.
struct IfaceDiag { std::string name; int order; size_t begin, end; };

namespace detail {
// The field a drawn line's elements read. Structural elements installed in a later phase read the
// total displacement less the datum of their installation (Structures::install_datum); so must the
// report, or it would show the ground's deformation before the element existed as force in it.
// Applied only to a TOTAL field -- a phase increment or a dynamic relative field already starts at
// the installation -- and only to a line installed at a non-zero datum; otherwise `disp` itself.
inline const Eigen::VectorXd& installed_field(const Structures& structures, int install,
                                              const DofMap& dofs, const Eigen::VectorXd& disp,
                                              bool disp_is_total, Eigen::VectorXd& buf) {
    if (!disp_is_total || install < 0 || (size_t)install >= structures.install_datum.size())
        return disp;
    const Eigen::VectorXd& datum = structures.install_datum[(size_t)install];
    const Eigen::VectorXd* fixed = (size_t)install < structures.install_fixed.size()
                                       ? &structures.install_fixed[(size_t)install]
                                       : nullptr;
    buf = disp;
    for (int g = 0; g < dofs.total_dofs() && g < (int)buf.size(); ++g) {
        const int eq = dofs.equation(g);
        if (eq >= 0) buf[g] -= datum(eq);
        else if (fixed && g < fixed->size()) buf[g] -= (*fixed)(g);
    }
    return buf;
}
}  // namespace detail

// One line's force diagram from a FULL global-DOF displacement field.
// anchor_plastic / geogrid_plastic may be empty (elastic / dynamic).
// plate_plastic / plate5_plastic: the FULL committed M-N hinge state vectors
// (solver order); the plate branches slice their own chain out. Empty =>
// elastic recovery (the D6b rule above).
inline StructForce force_diagram(const DiagSpec& sp, const Structures& structures,
                                 const mesh::Mesh& mesh, const DofMap& dofs,
                                 const Eigen::VectorXd& disp_all,
                                 const std::vector<double>& anchor_plastic,
                                 const std::vector<double>& geogrid_plastic, bool elastic = false,
                                 const std::vector<double>& plate_plastic = {},
                                 const std::vector<double>& plate5_plastic = {},
                                 bool disp_is_total = true, bool from_state = false) {
    StructForce d; d.name = sp.name;
    // The line's installation cohort, read from its first element (a drawn line is installed whole).
    int install = -1;
    if (sp.kind == 0 && sp.begin < structures.plates.size()) install = structures.plates[sp.begin].install;
    else if (sp.kind == 5 && sp.begin < structures.plates5.size()) install = structures.plates5[sp.begin].install;
    else if (sp.kind == 1 && sp.begin < structures.anchors.size()) install = structures.anchors[sp.begin].install;
    else if (sp.kind == 3 && sp.begin < structures.embedded_beams.size()) install = structures.embedded_beams[sp.begin].install;
    else if (sp.kind == 2 && sp.begin < structures.geogrids.size()) install = structures.geogrids[sp.begin].install;
    Eigen::VectorXd buf;
    const Eigen::VectorXd& field =
        detail::installed_field(structures, install, dofs, disp_all, disp_is_total, buf);
    d.kind = (sp.kind == 3 || sp.kind == 5) ? 0 : sp.kind;   // beam / tri15 wall report N/Q/M like a plate
    if (sp.kind == 0) {            // plate / embedded wall (tri6): N, Q (Barlow), M along the chain
        const std::vector<PlateElement> chain(
            structures.plates.begin() + sp.begin, structures.plates.begin() + sp.end);
        std::vector<double> ps_slice;
        constexpr size_t S3 = (size_t)plate::kPlasticStateSize;
        if ((!elastic || from_state) && plate_plastic.size() >= sp.end * S3)
            ps_slice.assign(plate_plastic.begin() + sp.begin * S3,
                            plate_plastic.begin() + sp.end * S3);
        d.stations = plate_force_diagram(chain, mesh, dofs, field, ps_slice, 5, !elastic);
    } else if (sp.kind == 5) {     // tri15 embedded wall: 5-node plate N/Q/M (plate_force_diagram overload)
        const std::vector<PlateElement5> chain(
            structures.plates5.begin() + sp.begin, structures.plates5.begin() + sp.end);
        std::vector<double> ps_slice;
        constexpr size_t S5 = (size_t)plate::kPlasticStateSize5;
        if ((!elastic || from_state) && plate5_plastic.size() >= sp.end * S5)
            ps_slice.assign(plate5_plastic.begin() + sp.begin * S5,
                            plate5_plastic.begin() + sp.end * S5);
        d.stations = plate_force_diagram(chain, mesh, dofs, field, ps_slice, 9, !elastic);
    } else if (sp.kind == 1) {     // anchor: one axial force (mirrors the solver exactly)
        const auto& an = structures.anchors[sp.begin];
        const double Up = sp.begin < anchor_plastic.size() ? anchor_plastic[sp.begin] : 0.0;
        const auto af = anchor_force(an, mesh, dofs, field, Up, elastic);
        ForceStation st;
        // Report per-anchor: the solver's per-metre force produces the correct
        // yield flag; the displayed number, times Ls, is in the same unit as
        // the capacity the user entered.
        st.x = mesh.x[an.node_a]; st.y = mesh.y[an.node_a]; st.N = af.N * sp.anchor_spacing;
        st.ux = field[anchor_dof(an, dofs, 0)];   // for deformed-mesh overlay
        st.uy = field[anchor_dof(an, dofs, 1)];
        d.stations.push_back(st);
        d.yielded = af.yielded;
    } else if (sp.kind == 3) {     // embedded beam (pile row): N, Q, M along the pile
        d.stations = embedded_beam_force_diagram(structures.embedded_beams[sp.begin], dofs, field);
    } else {                       // geogrid: N at the Gauss stations, arc-length accumulated
        double s_off = 0.0;
        for (size_t gi = sp.begin; gi < sp.end; ++gi) {
            const auto& ge = structures.geogrids[gi];
            constexpr size_t S = (size_t)geogrid::kStateStride;
            std::array<double, S> ep{};
            if (geogrid_plastic.size() >= S * (gi + 1))
                for (size_t q = 0; q < S; ++q) ep[q] = geogrid_plastic[S * gi + q];
            auto st = ge.node_count() == 3 && !ge.props.axisym   // the historical element, to the bit
                          ? geogrid_force_diagram(ge, mesh, dofs, field, {ep[0], ep[1]}, elastic)
                          : geogrid_force_diagram_general(ge, mesh, dofs, field, ep, elastic);
            for (auto& s2 : st) { s2.s += s_off; d.stations.push_back(s2); }
            s_off += std::hypot(mesh.x[ge.nodes[1]] - mesh.x[ge.nodes[0]],
                                mesh.y[ge.nodes[1]] - mesh.y[ge.nodes[0]]);
        }
    }
    return d;
}

// One Coulomb joint's tau / sigma_n / relative-movement stations. `elastic`
// reports the ELASTIC branch (tau = ks du_s, sigma_n = kn du_n, no Coulomb
// cap, no sigma_n0) -- which is what the linear dynamic system actually
// solves, so its output stays in equilibrium with its own displacements.
inline InterfaceResult force_diagram(const IfaceDiag& is, const Structures& structures,
                                     const mesh::Mesh& mesh, const DofMap& dofs,
                                     const Eigen::VectorXd& disp,
                                     const std::vector<double>& slip3,
                                     const std::vector<double>& slip5, bool elastic,
                                     bool disp_is_total = true, bool from_state = false,
                                     const Eigen::VectorXd* parent_total = nullptr) {
    InterfaceResult ir; ir.name = is.name;
    int install = -1;
    if (is.order == 15 && is.begin < structures.interfaces5.size())
        install = structures.interfaces5[is.begin].install;
    else if (is.order != 15 && is.begin < structures.interfaces.size())
        install = structures.interfaces[is.begin].install;
    Eigen::VectorXd buf;
    const Eigen::VectorXd& field =
        detail::installed_field(structures, install, dofs, disp, disp_is_total, buf);
    // The displacement the phase started from, measured from the same installation datum.
    Eigen::VectorXd pbuf;
    const Eigen::VectorXd* pfield = nullptr;
    if (from_state && parent_total && parent_total->size() == disp.size())
        pfield = &detail::installed_field(structures, install, dofs, *parent_total, true, pbuf);
    ir.stations = (is.order == 15)
        ? interface_force_diagram(structures.interfaces5, is.begin, is.end, mesh,
                                  dofs, field, slip5, elastic, from_state, pfield)
        : interface_force_diagram(structures.interfaces, is.begin, is.end, mesh,
                                  dofs, field, slip3, elastic, from_state, pfield);
    return ir;
}

} // namespace katai::core
