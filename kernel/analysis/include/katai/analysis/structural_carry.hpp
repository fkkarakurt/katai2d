#pragma once
// Parent structural-state carry (Track 1a; engine-owned since Stage B9).
// Soil carries its committed stresses through `initial_state`, but structural
// elements are total-displacement formulated (f = f(u_total, plastic state)) --
// a chained phase inherits the parent's structural state only if it is given
// (a) the converged displacement DATUM and (b) the committed plastic states.
// Without them every phase restarts its structures at u = 0 with zero plastic
// memory, and the staged-construction imbalance re-ramps the parent's structural
// tractions -- measured: an unchanged nil phase drifted the wall moment by 32%.
//
// build_structural_init carries a state to an IDENTICAL structure set only;
// build_structural_carry matches it structure by structure, so that a phase which
// activates or removes a structure still continues every other one.

#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>  // Structures, StructuralInit
#include <katai/analysis/results.hpp>           // StructCarryState
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/geogrid.hpp>       // geogrid::kGaussCount
#include <katai/fem/elements/interface.hpp>     // iface::kPointCount(5)
#include <katai/fem/elements/plate.hpp>         // plate::kPlasticStateSize(5)

namespace katai::core {

// Fill `out` from the parent's carried structural state. Returns true only when
// every member matches the current structure set exactly; on false, `out` is
// untouched and the caller must not consume the carry.
inline bool build_structural_init(const StructCarryState& ps, const Structures& structures,
                                  const DofMap& dofs, StructuralInit& out) {
    size_t skin_total = 0;
    for (const auto& eb : structures.embedded_beams) skin_total += eb.skin.size();
    if (ps.full_disp.size() != (Eigen::Index)dofs.total_dofs() ||
        ps.anchor_plastic.size() != structures.anchors.size() ||
        ps.geogrid_plastic.size() !=
            structures.geogrids.size() * (size_t)geogrid::kGaussCount ||
        ps.interface_slip.size() !=
            structures.interfaces.size() * (size_t)iface::kPointCount ||
        ps.interface5_slip.size() !=
            structures.interfaces5.size() * (size_t)iface::kPointCount5 ||
        ps.embedded_skin_slip.size() != skin_total ||
        ps.embedded_foot_slip.size() != structures.embedded_beams.size() ||
        ps.plate_plastic.size() !=
            structures.plates.size() * (size_t)plate::kPlasticStateSize ||
        ps.plate5_plastic.size() !=
            structures.plates5.size() * (size_t)plate::kPlasticStateSize5)
        return false;
    // Datum in EQUATION space: the solvers' structural loops read free DOFs
    // only (a fixed DOF contributes zero there, exactly as it does in every
    // static solve).
    out.u_datum = Eigen::VectorXd::Zero(dofs.equation_count());
    for (int g = 0; g < dofs.total_dofs(); ++g) {
        const int eq = dofs.equation(g);
        if (eq >= 0) out.u_datum(eq) = ps.full_disp(g);
    }
    out.anchor_plastic = ps.anchor_plastic;
    out.geogrid_plastic = ps.geogrid_plastic;
    out.interface_slip = ps.interface_slip;
    out.interface5_slip = ps.interface5_slip;
    out.embedded_skin_slip = ps.embedded_skin_slip;
    out.embedded_foot_slip = ps.embedded_foot_slip;
    out.plate_plastic = ps.plate_plastic;
    out.plate5_plastic = ps.plate5_plastic;
    return true;
}

// What the carry across a phase boundary did, for the phase to act on and to say.
struct StructuralCarry {
    bool carried = false;         // the phase starts from the parent's structural state
    int kept = 0;                 // drawn structures continued from the parent
    int installed = 0;            // drawn structures installed in this phase, on the deformed ground
    Eigen::VectorXd full_datum;   // the datum in THIS phase's global-DOF numbering (total_dofs)
    std::vector<char> new_anchor; // per structures.anchors: installed in this phase
    // Per child record: 1 = installed in this phase (not found in the parent). Filled whether or not
    // anything was carried -- with no parent state every structure is installed.
    std::vector<char> record_installed;
};

// THE CARRY ACROSS A CHANGED STRUCTURE SET. build_structural_init above accepts the parent's state
// only when the structure set is identical, and until 2026-09 any other case dropped the whole
// carry: activating one anchor reset every other structure's force to what the new phase alone put
// into it (measured: a wall's moment 52% low in a phase whose only change was an anchor carrying
// nothing). Here the parent's state is matched by DRAWN STRUCTURE (StructCarryRecord), because the
// numbers of the extra DOFs and the positions in the plastic-state vectors both move when a
// structure before them appears or disappears:
//   - kept (in both phases, same element counts): its plastic state, its extra-DOF datum and its
//     installation cohort are carried over;
//   - installed (new in this phase): zero plastic state, and the datum of this phase as its
//     installation datum, so it starts stress-free on the ground as the parent left it;
//   - removed (in the parent only): nothing is built, and its force leaves the baseline -- released
//     as part of this phase's staged change, which is what removing a strut does.
// The mesh DOFs keep their numbers across phases (node_dofs). Nothing is carried when no structure is
// kept: the phase then starts its structures from zero exactly as before. `child_records` are this
// phase's own records; their `install` is set here, together with every element's.
inline StructuralCarry build_structural_carry(const StructCarryState& ps,
                                              std::vector<StructCarryRecord>& child_records,
                                              Structures& structures, const DofMap& dofs,
                                              int node_dofs, StructuralInit& out) {
    StructuralCarry plan;
    plan.record_installed.assign(child_records.size(), 1);
    if (ps.full_disp.size() == 0) return plan;
    // A parent state with no records has none to match by: an identical set is still carried as
    // before (build_structural_init checks it), anything else is not.
    if (ps.records.empty()) {
        if (!child_records.empty() || !build_structural_init(ps, structures, dofs, out)) return plan;
        plan.carried = true;
        plan.full_datum = ps.full_disp;
        plan.new_anchor.assign(structures.anchors.size(), 0);
        return plan;
    }
    if (ps.node_dofs != node_dofs || ps.full_disp.size() < node_dofs) return plan;

    const auto len = [](const std::array<size_t, 2>& r) { return r[1] - r[0]; };
    const auto same_shape = [&](const StructCarryRecord& a, const StructCarryRecord& b) {
        return len(a.plates) == len(b.plates) && len(a.plates5) == len(b.plates5) &&
               len(a.anchors) == len(b.anchors) && len(a.geogrids) == len(b.geogrids) &&
               len(a.interfaces) == len(b.interfaces) && len(a.interfaces5) == len(b.interfaces5) &&
               len(a.embedded) == len(b.embedded) && len(a.skin) == len(b.skin) &&
               a.extra_dof[1] - a.extra_dof[0] == b.extra_dof[1] - b.extra_dof[0];
    };
    std::vector<const StructCarryRecord*> parent_of(child_records.size(), nullptr);
    for (size_t c = 0; c < child_records.size(); ++c)
        for (const auto& p : ps.records)
            if (p.si == child_records[c].si && same_shape(p, child_records[c])) {
                parent_of[c] = &p;
                break;
            }
    for (size_t c = 0; c < parent_of.size(); ++c) {
        (parent_of[c] ? plan.kept : plan.installed) += 1;
        plan.record_installed[c] = parent_of[c] ? 0 : 1;
    }
    if (plan.kept == 0) { plan.installed = 0; return plan; }

    // This phase's global DOF -> the parent's, for every DOF that exists in both.
    const int ntot = dofs.total_dofs();
    std::vector<int> c2p(ntot, -1);
    for (int g = 0; g < node_dofs && g < ntot; ++g) c2p[g] = g;
    for (size_t c = 0; c < child_records.size(); ++c) {
        const auto* p = parent_of[c];
        if (!p) continue;
        for (int k = 0; k < child_records[c].extra_dof[1] - child_records[c].extra_dof[0]; ++k) {
            const int gc = child_records[c].extra_dof[0] + k, gp = p->extra_dof[0] + k;
            if (gc < ntot && gp < ps.full_disp.size()) c2p[gc] = gp;
        }
    }
    const auto map_full = [&](const Eigen::VectorXd& parent_full) {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(ntot);
        for (int g = 0; g < ntot; ++g)
            if (c2p[g] >= 0 && c2p[g] < parent_full.size()) v[g] = parent_full[c2p[g]];
        return v;
    };
    const auto to_equations = [&](const Eigen::VectorXd& full) {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(dofs.equation_count());
        for (int g = 0; g < ntot; ++g) {
            const int eq = dofs.equation(g);
            if (eq >= 0) v(eq) = full[g];
        }
        return v;
    };
    plan.full_datum = map_full(ps.full_disp);
    out = StructuralInit{};
    out.u_datum = to_equations(plan.full_datum);

    // Plastic state, sized for this phase's structures and filled structure by structure.
    size_t skin_total = 0;
    for (const auto& eb : structures.embedded_beams) skin_total += eb.skin.size();
    out.anchor_plastic.assign(structures.anchors.size(), 0.0);
    out.geogrid_plastic.assign(structures.geogrids.size() * (size_t)geogrid::kGaussCount, 0.0);
    out.interface_slip.assign(structures.interfaces.size() * (size_t)iface::kPointCount, 0.0);
    out.interface5_slip.assign(structures.interfaces5.size() * (size_t)iface::kPointCount5, 0.0);
    out.embedded_skin_slip.assign(skin_total, 0.0);
    out.embedded_foot_slip.assign(structures.embedded_beams.size(), 0.0);
    out.plate_plastic.assign(structures.plates.size() * (size_t)plate::kPlasticStateSize, 0.0);
    out.plate5_plastic.assign(structures.plates5.size() * (size_t)plate::kPlasticStateSize5, 0.0);
    const auto copy_slice = [](const std::vector<double>& from, const std::array<size_t, 2>& pr,
                               std::vector<double>& to, const std::array<size_t, 2>& cr,
                               size_t per) {
        for (size_t k = 0; k < (pr[1] - pr[0]) * per; ++k)
            if (pr[0] * per + k < from.size() && cr[0] * per + k < to.size())
                to[cr[0] * per + k] = from[pr[0] * per + k];
    };

    // Installation cohorts: the parent's that a kept structure still uses, then one for this phase.
    std::vector<int> cohort_of_parent(ps.install_datum.size(), -1);
    std::vector<Eigen::VectorXd> cohorts;
    int this_phase = -1;
    for (size_t c = 0; c < child_records.size(); ++c) {
        auto& cr = child_records[c];
        const auto* p = parent_of[c];
        int install = -1;
        if (p) {
            copy_slice(ps.anchor_plastic, p->anchors, out.anchor_plastic, cr.anchors, 1);
            copy_slice(ps.geogrid_plastic, p->geogrids, out.geogrid_plastic, cr.geogrids,
                       (size_t)geogrid::kGaussCount);
            copy_slice(ps.interface_slip, p->interfaces, out.interface_slip, cr.interfaces,
                       (size_t)iface::kPointCount);
            copy_slice(ps.interface5_slip, p->interfaces5, out.interface5_slip, cr.interfaces5,
                       (size_t)iface::kPointCount5);
            copy_slice(ps.embedded_skin_slip, p->skin, out.embedded_skin_slip, cr.skin, 1);
            copy_slice(ps.embedded_foot_slip, p->embedded, out.embedded_foot_slip, cr.embedded, 1);
            copy_slice(ps.plate_plastic, p->plates, out.plate_plastic, cr.plates,
                       (size_t)plate::kPlasticStateSize);
            copy_slice(ps.plate5_plastic, p->plates5, out.plate5_plastic, cr.plates5,
                       (size_t)plate::kPlasticStateSize5);
            if (p->install >= 0 && (size_t)p->install < ps.install_datum.size()) {
                int& mapped = cohort_of_parent[(size_t)p->install];
                if (mapped < 0) {
                    mapped = (int)cohorts.size();
                    cohorts.push_back(to_equations(map_full(ps.install_datum[(size_t)p->install])));
                }
                install = mapped;
            }
        } else {
            if (this_phase < 0) {
                this_phase = (int)cohorts.size();
                cohorts.push_back(out.u_datum);
            }
            install = this_phase;
        }
        cr.install = install;
        for (size_t i = cr.plates[0]; i < cr.plates[1]; ++i) structures.plates[i].install = install;
        for (size_t i = cr.plates5[0]; i < cr.plates5[1]; ++i) structures.plates5[i].install = install;
        for (size_t i = cr.anchors[0]; i < cr.anchors[1]; ++i) structures.anchors[i].install = install;
        for (size_t i = cr.geogrids[0]; i < cr.geogrids[1]; ++i) structures.geogrids[i].install = install;
        for (size_t i = cr.interfaces[0]; i < cr.interfaces[1]; ++i)
            structures.interfaces[i].install = install;
        for (size_t i = cr.interfaces5[0]; i < cr.interfaces5[1]; ++i)
            structures.interfaces5[i].install = install;
        for (size_t i = cr.embedded[0]; i < cr.embedded[1]; ++i)
            structures.embedded_beams[i].install = install;
    }
    structures.install_datum = std::move(cohorts);
    plan.new_anchor.assign(structures.anchors.size(), 0);
    for (size_t c = 0; c < child_records.size(); ++c)
        if (!parent_of[c])
            for (size_t i = child_records[c].anchors[0]; i < child_records[c].anchors[1]; ++i)
                plan.new_anchor[i] = 1;
    plan.carried = true;
    return plan;
}

// The installation datums of a phase's structures as FULL global-DOF vectors, for the carry state
// the next phase reads (StructCarryState::install_datum). Fixed DOFs hold 0: a datum only exists on
// an equation.
inline std::vector<Eigen::VectorXd> install_datum_full(const Structures& structures,
                                                       const DofMap& dofs) {
    std::vector<Eigen::VectorXd> full;
    for (const auto& d : structures.install_datum) {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(dofs.total_dofs());
        for (int g = 0; g < dofs.total_dofs(); ++g) {
            const int eq = dofs.equation(g);
            if (eq >= 0 && eq < d.size()) v[g] = d(eq);
        }
        full.push_back(std::move(v));
    }
    return full;
}

} // namespace katai::core
