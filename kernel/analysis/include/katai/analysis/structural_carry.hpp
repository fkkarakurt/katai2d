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

// The fixed-DOF part of a full global-DOF displacement (Structures::fixed_datum): its value on
// every DOF this phase holds, zero on the others. An all-zero part is returned EMPTY, so a chain
// whose fixed DOFs never moved carries nothing new.
inline Eigen::VectorXd fixed_part(const Eigen::VectorXd& full, const DofMap& dofs) {
    Eigen::VectorXd v = Eigen::VectorXd::Zero(dofs.total_dofs());
    bool any = false;
    for (int g = 0; g < dofs.total_dofs() && g < full.size(); ++g)
        if (dofs.equation(g) < 0 && full[g] != 0.0) { v[g] = full[g]; any = true; }
    return any ? v : Eigen::VectorXd();
}

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
        structures.fixed_datum = fixed_part(ps.full_disp, dofs);
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
    structures.fixed_datum = fixed_part(plan.full_datum, dofs);

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
    std::vector<Eigen::VectorXd> cohorts, cohorts_fixed;
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
                    const Eigen::VectorXd full = map_full(ps.install_datum[(size_t)p->install]);
                    cohorts.push_back(to_equations(full));
                    cohorts_fixed.push_back(fixed_part(full, dofs));
                }
                install = mapped;
            }
        } else {
            if (this_phase < 0) {
                this_phase = (int)cohorts.size();
                cohorts.push_back(out.u_datum);
                cohorts_fixed.push_back(structures.fixed_datum);
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
    structures.install_fixed = std::move(cohorts_fixed);
    plan.new_anchor.assign(structures.anchors.size(), 0);
    for (size_t c = 0; c < child_records.size(); ++c)
        if (!parent_of[c])
            for (size_t i = child_records[c].anchors[0]; i < child_records[c].anchors[1]; ++i)
                plan.new_anchor[i] = 1;
    plan.carried = true;
    return plan;
}

// The installation datums of a phase's structures as FULL global-DOF vectors, for the carry state
// the next phase reads (StructCarryState::install_datum): the equation-space datum on a free DOF,
// and the fixed-DOF datum (Structures::install_fixed) on a fixed one.
inline std::vector<Eigen::VectorXd> install_datum_full(const Structures& structures,
                                                       const DofMap& dofs) {
    std::vector<Eigen::VectorXd> full;
    for (size_t k = 0; k < structures.install_datum.size(); ++k) {
        const Eigen::VectorXd& d = structures.install_datum[k];
        const Eigen::VectorXd* fx =
            k < structures.install_fixed.size() ? &structures.install_fixed[k] : nullptr;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(dofs.total_dofs());
        for (int g = 0; g < dofs.total_dofs(); ++g) {
            const int eq = dofs.equation(g);
            if (eq >= 0 && eq < d.size()) v[g] = d(eq);
            else if (eq < 0 && fx && g < fx->size()) v[g] = (*fx)(g);
        }
        full.push_back(std::move(v));
    }
    return full;
}

// A PLATE WHOSE SECTION CHANGES KEEPS ITS FORCES. A plate is evaluated on its total
// deformation, N = EA (eps - eps_p) and M = EI (kap - kap_p), so giving it another stiffness
// would multiply the forces it already carries by the ratio of the two -- a lining twice as
// stiff as the shield it replaces would start with twice the ring force. The part of the
// deformation the new section does not answer for is therefore moved into the state the
// element keeps at its integration points: eps_p' = eps - (EA_old / EA_new)(eps - eps_p), and
// the same for the curvature. N and M are then what they were, to round-off, and every further
// deformation is taken with the new stiffness. The shear force is kept the same way, through
// the shear strain offset at the reduced shear points.
// `full_datum` is the parent's total displacement in this phase's DOF numbering.
inline void rescale_plate_state(const mesh::Mesh& mesh, const DofMap& dofs,
                                const Structures& structures, const Eigen::VectorXd& full_datum,
                                const StructCarryRecord& rec, double EA_old, double EI_old,
                                double GA_old, StructuralInit& init) {
    const auto rel = [&](int gdof, int install) -> double {
        const double u = gdof >= 0 && gdof < full_datum.size() ? full_datum[gdof] : 0.0;
        if (install < 0) return u;
        const int eq = dofs.equation(gdof);
        if (eq >= 0)
            return (size_t)install < structures.install_datum.size()
                       ? u - structures.install_datum[(size_t)install](eq) : u;
        const Eigen::VectorXd* fx = (size_t)install < structures.install_fixed.size()
                                        ? &structures.install_fixed[(size_t)install] : nullptr;
        return u - (fx && gdof < fx->size() ? (*fx)(gdof) : 0.0);
    };
    const auto moved = [](double strain, double strain_p, double k_old, double k_new) {
        return k_new > 0.0 ? strain - (k_old / k_new) * (strain - strain_p) : strain_p;
    };
    for (size_t i = rec.plates[0]; i < rec.plates[1] && i < structures.plates.size(); ++i) {
        const PlateElement& pe = structures.plates[i];
        plate::NodeCoords X;
        plate::Dof u;
        for (int k = 0; k < 3; ++k) {
            X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]];
            for (int c = 0; c < 2; ++c) {
                const int g = pe.trans_dof[2 * k + c] >= 0 ? pe.trans_dof[2 * k + c]
                                                           : dofs.global_dof(pe.nodes[k], c);
                u(3 * k + c) = rel(g, pe.install);
            }
            u(3 * k + 2) = rel(pe.rot_dof[k], pe.install);
        }
        const size_t off = i * (size_t)plate::kPlasticStateSize;
        if (off + (size_t)plate::kPlasticStateSize > init.plate_plastic.size()) continue;  // silent-drop-ok: no state to move
        for (int q = 0; q < plate::kBendGaussCount; ++q) {
            const auto e = plate::detail::edge_kin(X, plate::kBendGaussXi[q]);
            Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
            plate::detail::b_rows(e, Be, Bk, Bg);
            double& ep = init.plate_plastic[off + (size_t)plate::kPlasticPerPoint * q];
            double& kp = init.plate_plastic[off + (size_t)plate::kPlasticPerPoint * q + 1];
            ep = moved((Be * u)(0), ep, EA_old, pe.props.EA);
            kp = moved((Bk * u)(0), kp, EI_old, pe.props.EI);
        }
        constexpr double g2 = 0.5773502691896257;
        const std::array<double, 2> xi2{-g2, g2};
        for (int q = 0; q < plate::kShearGaussCount; ++q) {
            const auto e = plate::detail::edge_kin(X, xi2[(size_t)q]);
            Eigen::Matrix<double, 1, 9> Be, Bk, Bg;
            plate::detail::b_rows(e, Be, Bk, Bg);
            double& g0 = init.plate_plastic[off + (size_t)plate::kShearStateBegin + q];
            g0 = moved((Bg * u)(0), g0, GA_old, pe.props.shear_rigidity());
        }
    }
    for (size_t i = rec.plates5[0]; i < rec.plates5[1] && i < structures.plates5.size(); ++i) {
        const PlateElement5& pe = structures.plates5[i];
        plate::NodeCoords5 X;
        plate::Dof5 u;
        for (int k = 0; k < 5; ++k) {
            X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]];
            for (int c = 0; c < 2; ++c) {
                const int g = pe.trans_dof[2 * k + c] >= 0 ? pe.trans_dof[2 * k + c]
                                                           : dofs.global_dof(pe.nodes[k], c);
                u(3 * k + c) = rel(g, pe.install);
            }
            u(3 * k + 2) = rel(pe.rot_dof[k], pe.install);
        }
        const size_t off = i * (size_t)plate::kPlasticStateSize5;
        if (off + (size_t)plate::kPlasticStateSize5 > init.plate5_plastic.size()) continue;  // silent-drop-ok: no state to move
        for (int q = 0; q < plate::kBendGaussCount5; ++q) {
            const auto e = plate::detail::edge_kin5(X, plate::kBendGaussXi5[q]);
            Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
            plate::detail::b_rows5(e, Be, Bk, Bg);
            double& ep = init.plate5_plastic[off + (size_t)plate::kPlasticPerPoint * q];
            double& kp = init.plate5_plastic[off + (size_t)plate::kPlasticPerPoint * q + 1];
            ep = moved((Be * u)(0), ep, EA_old, pe.props.EA);
            kp = moved((Bk * u)(0), kp, EI_old, pe.props.EI);
        }
        constexpr double a4 = 0.3399810435848563, b4 = 0.8611363115940526;
        const std::array<double, 4> xi4{-b4, -a4, a4, b4};
        for (int q = 0; q < plate::kShearGaussCount5; ++q) {
            const auto e = plate::detail::edge_kin5(X, xi4[(size_t)q]);
            Eigen::Matrix<double, 1, 15> Be, Bk, Bg;
            plate::detail::b_rows5(e, Be, Bk, Bg);
            double& g0 = init.plate5_plastic[off + (size_t)plate::kShearStateBegin5 + q];
            g0 = moved((Bg * u)(0), g0, GA_old, pe.props.shear_rigidity());
        }
    }
}

} // namespace katai::core
