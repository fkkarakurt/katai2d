#pragma once
// WHAT A COUPLED PHASE OWES THE READER ABOUT ITS STRUCTURES -- written once, called by both the
// consolidation and the fully-coupled strategies, because it is one statement and not two.
//
// Both phases solve their structural elements ELASTICALLY: the coupled system takes their
// stiffness (katai/analysis/consolidation.hpp, `struct_k`) and no return mapping runs on them. For
// a plate or an anchor that is a limit on the CAPACITY -- they are elastic until they hinge or
// yield. For a geogrid it is not a capacity at all but the element's own behaviour, tension-only,
// which the elastic branch does not have. For an interface it is the Coulomb slip, which is what
// an interface is FOR. None of the three may be reached in silence, so each is measured against
// the input the engineer entered and reported when it is passed:
//
//   K2D-A014  a plate / anchor / geogrid past its capacity, with the utilisation in the units the
//             capacity was entered in (per anchor, not per metre of wall);
//   K2D-A015  a geogrid in compression, where the elastic branch made it push back on the soil
//             instead of going slack -- the unsafe direction;
//   K2D-A016  an interface past its Coulomb capacity, where a joint that cannot slip is stiffer
//             than the real one, so the wall deflects less and attracts more load -- also unsafe.
//
// Sizes measured on the verification fixtures are quoted in the messages themselves (KV-STR-006,
// KV-STR-007), because a warning that cannot say how much is a warning the reader must guess at.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <Eigen/Core>

#include <katai/analysis/nonlinear_solver.hpp>      // Structures
#include <katai/analysis/results.hpp>
#include <katai/analysis/structural_diagrams.hpp>   // DiagSpec / IfaceDiag / force_diagram
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace katai::core {

// `disp_total` is the phase increment PLUS the parent's displacement datum: structural elements
// are total-displacement formulated, so a wall installed in an earlier phase carries its force at
// that datum and reporting the increment alone would understate it by exactly the parent's share.
inline void report_coupled_structures(const Structures& structures,
                                      const std::vector<DiagSpec>* diagrams,
                                      const std::vector<IfaceDiag>* iface_diagrams,
                                      const katai::mesh::Mesh& mesh, const DofMap& dofs,
                                      const Eigen::VectorXd& disp_total, SolveResult& R) {
    if (!diagrams && !iface_diagrams) return;
    int over = 0;
    std::string over_note;
    int slack = 0;
    std::string slack_note;
    for (const auto& sp : (diagrams ? *diagrams : std::vector<DiagSpec>{})) {
        StructForce d = force_diagram(sp, structures, mesh, dofs, disp_total, {}, {},
                                      /*elastic=*/true);
        const auto env = force_envelope(d.stations);
        d.max_N = env.max_abs_N; d.max_Q = env.max_abs_Q; d.max_M = env.max_abs_M;
        // THE COST OF AN ELASTIC BRANCH, MADE VISIBLE. Nothing here capped the force at the
        // capacity the engineer entered, so a line that has passed it is reporting a force
        // the element could not carry -- and the redistribution that yielding would have
        // caused did not happen anywhere else either. That is exactly the kind of limit that
        // must not be reachable in silence, so the run says which line and by how much.
        double util = 0.0;
        if ((sp.kind == 0 || sp.kind == 5) && sp.begin < structures.plates.size() &&
            sp.kind == 0) {
            const auto& pp = structures.plates[sp.begin].props;
            const double iN = pp.Np > 0.0 ? 1.0 / pp.Np : 0.0;
            const double iM = pp.Mp > 0.0 ? 1.0 / pp.Mp : 0.0;
            for (const auto& st : d.stations)
                util = std::fmax(util, std::fabs(st.N) * iN + std::fabs(st.M) * iM);
        } else if (sp.kind == 1 && sp.begin < structures.anchors.size()) {
            const auto& an = structures.anchors[sp.begin];
            // The element's cap is per metre of wall and the diagram is per anchor, so the
            // comparison is made in the diagram's units (the audit finding behind DiagSpec).
            const double ft = an.Fmax_tens > 0.0 ? an.Fmax_tens * sp.anchor_spacing : 0.0;
            const double fc = an.Fmax_comp > 0.0 ? an.Fmax_comp * sp.anchor_spacing : 0.0;
            for (const auto& st : d.stations) {
                if (st.N > 0.0 && ft > 0.0) util = std::fmax(util, st.N / ft);
                if (st.N < 0.0 && fc > 0.0) util = std::fmax(util, -st.N / fc);
            }
        } else if (sp.kind == 2 && sp.begin < structures.geogrids.size()) {
            const double np = structures.geogrids[sp.begin].props.Np;
            if (np > 0.0)
                for (const auto& st : d.stations) util = std::fmax(util, st.N / np);
            // A GEOGRID IS NOT AN ELASTIC BAR, AND THE DIFFERENCE IS NOT A CAPACITY. Its
            // defining behaviour is tension-only: in compression it goes slack and carries
            // nothing (reversible). The elastic branch this phase solves with has no such
            // cut, so wherever the sheet is compressed it PUSHES BACK -- it stiffens the
            // ground where the real one would have stopped acting, which errs on the unsafe
            // side. It is exactly measurable, so it is reported rather than declared: on the
            // verified fixture a sheet with compressed ends came out 1.57% off the drained
            // answer, and the same sheet kept wholly in tension came out 0.00000% off
            // (KV-STR-006). Every station in tension = no finding at all.
            double nmax = 0.0;
            for (const auto& st : d.stations) nmax = std::fmax(nmax, std::fabs(st.N));
            int comp = 0;
            for (const auto& st : d.stations) if (st.N < -1e-6 * nmax) ++comp;
            if (comp > 0) {
                ++slack;
                char b[160];
                std::snprintf(b, sizeof(b), "%s%s (%d of %d stations)",
                              slack > 1 ? ", " : "", d.name.c_str(), comp,
                              (int)d.stations.size());
                slack_note += b;
            }
        }
        if (util > 1.0) {
            ++over;
            char b[160];
            std::snprintf(b, sizeof(b), "%s%s at %.0f%% of its capacity",
                          over > 1 ? ", " : "", d.name.c_str(), 100.0 * util);
            over_note += b;
            d.yielded = true;   // it would have, had this phase been able to let it
        }
        R.struct_forces.push_back(std::move(d));
    }
    // The joints, reported the way the dynamic phase reports its own elastic ones: tau,
    // sigma_n and the Coulomb demand/capacity ratio. `slip_checked` is FALSE on purpose --
    // this phase applied no Coulomb return, so a station marked "bonded" means "not checked",
    // not "checked and found bonded", and a reader who cannot tell those apart will read a
    // wall that slipped as a wall that held.
    double iface_util = 0.0;
    std::string iface_note;
    if (iface_diagrams)
        for (const auto& is : *iface_diagrams) {
            InterfaceResult ir = force_diagram(is, structures, mesh, dofs, disp_total,
                                               {}, {}, /*elastic=*/true);
            ir.slip_checked = false;
            const auto& props = (is.order == 15) ? structures.interfaces5[is.begin].props
                                                 : structures.interfaces[is.begin].props;
            int over_st = 0;
            for (auto& st : ir.stations) {
                const double tmax = std::fmax(0.0, props.c_i - st.sigma_n * std::tan(props.phi_i));
                st.utilisation = tmax > 1e-12 ? std::fabs(st.tau) / tmax
                                              : (std::fabs(st.tau) > 1e-12 ? 100.0 : 0.0);
                ir.max_utilisation = std::fmax(ir.max_utilisation, st.utilisation);
                if (st.utilisation > 1.0) ++over_st;
                ir.max_abs_tau = std::fmax(ir.max_abs_tau, std::fabs(st.tau));
                ir.max_abs_sigma_n = std::fmax(ir.max_abs_sigma_n, std::fabs(st.sigma_n));
                ir.max_abs_slip = std::fmax(ir.max_abs_slip, std::fabs(st.slip));
            }
            if (!ir.stations.empty())
                ir.over_fraction = (double)over_st / (double)ir.stations.size();
            if (ir.max_utilisation > 1.0) {
                iface_util = std::fmax(iface_util, ir.max_utilisation);
                char b[200];
                std::snprintf(b, sizeof(b), "%s\"%s\" at %.2fx over %.0f%% of its length",
                              iface_note.empty() ? "" : ", ", ir.name.c_str(),
                              ir.max_utilisation, 100.0 * ir.over_fraction);
                iface_note += b;
            }
            R.interface_forces.push_back(std::move(ir));
        }
    if (iface_util > 1.0)
        add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A016", "interfaces",
                       "an interface in this phase passed its Coulomb capacity: " + iface_note +
                       ". The structural branch of a consolidation phase is ELASTIC, so the "
                       "joint could not slip -- it carried the shear instead. A joint that "
                       "cannot slip is STIFFER than the real one, so the wall it holds "
                       "deflects less and attracts more load than it would: the error is on "
                       "the unsafe side. Check the stage in a Plastic phase, where the "
                       "Coulomb return actually runs.");
    if (over > 0)
        add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A014", "structures",
                       "the structural elements in a consolidation phase are solved ELASTIC: "
                       "no anchor yield, no geogrid tension cut-off, no plate hinge. " +
                       over_note +
                       " -- so those forces are larger than the elements can carry, and the "
                       "redistribution that yielding would have caused is missing from the "
                       "whole result, soil included. Check the capacity in a Plastic phase "
                       "at the same stage before using these numbers.");
    if (slack > 0)
        add_diagnostic(R, DiagnosticSeverity::Warning, "K2D-A015", "structures",
                       "a geogrid in this phase is in COMPRESSION over part of its length: " +
                       slack_note +
                       ". A geogrid carries tension only -- in compression it goes slack and "
                       "carries nothing -- but the structural branch of a consolidation phase "
                       "is elastic and has no such cut, so those stations pushed BACK on the "
                       "soil instead. That stiffens the ground where the real sheet would "
                       "have stopped acting, so the settlement here is on the unsafe side. "
                       "Measured on the verification fixture: 1.57% on the sheet's force "
                       "where it had compressed ends, and 0.00000% where it did not.");
}

} // namespace katai::core
