#include <katai/analysis/strength_reduction.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace katai::core {

// Factor a material's SHEAR STRENGTH by the strength reduction factor srf (phi-c reduction):
//   c_f = c / srf,   phi_f = atan(tan phi / srf)
// applied to BOTH the Mohr-Coulomb fields (used by the MC model) AND the Hardening Soil sub-struct
// (the HS failure surface uses hs.cohesion / hs.friction -- a SEPARATE set of fields, so reducing
// only the MC ones leaves an HS slope at full strength and it never collapses). Dilatancy is clamped
// to the reduced friction afterwards (psi <= phi is required) so the reduced state
// stays physically admissible at high srf.
static void factor_strength(MaterialModel& m, double srf) {
    m.cohesion /= srf;
    m.friction_angle = std::atan(std::tan(m.friction_angle) / srf);
    if (m.dilatancy_angle > m.friction_angle) m.dilatancy_angle = m.friction_angle;
    // Tensile strength is a material strength too: the Safety reduction of the
    // tension cut-off value is defined explicitly for Hoek-Brown; the same rule
    // is applied to the MC cap here (safe direction; the c*cot(phi) clamp inside
    // the return mapping shrinks with c/srf as well).
    m.tensile_strength /= srf;
    m.hs.cohesion /= srf;
    m.hs.friction = std::atan(std::tan(m.hs.friction) / srf);
    if (m.hs.dilatancy > m.hs.friction) m.hs.dilatancy = m.hs.friction;
    // Soft Soil carries its own c/phi fields (ssoil) -- without reducing them Safety would
    // SILENTLY report an inflated FoS. (SS Safety is currently gated off in build_problem;
    // these lines are defensive so it stays correct if the gate opens. The K0nc/M
    // calibration is stiffness structure, not strength: constant.)
    m.ssoil.c /= srf;
    m.ssoil.phi = std::atan(std::tan(m.ssoil.phi) / srf);
    if (m.ssoil.psi > m.ssoil.phi) m.ssoil.psi = m.ssoil.phi;
    // Soft Soil Creep has a block of its own too, and the same reason applies to it.
    m.ssc.c /= srf;
    m.ssc.phi = std::atan(std::tan(m.ssc.phi) / srf);
    if (m.ssc.psi > m.ssc.phi) m.ssc.psi = m.ssc.phi;
}

// An interface's strength is a soil strength -- c_i = R_inter c and tan(phi_i) = R_inter tan(phi) --
// so the reduction is the soil's, term for term, tension cut-off included. The joint has no
// dilatancy to clamp (its Coulomb return is non-dilatant), and its stiffnesses are not strengths.
static void factor_interface_strength(iface::InterfaceProps& p, double srf) {
    p.c_i /= srf;
    p.phi_i = std::atan(std::tan(p.phi_i) / srf);
    p.sigma_t /= srf;
}

double factor_of_safety(const mesh::Mesh& mesh, const DofMap& dofs,
                        const Eigen::VectorXd& gravity_load,
                        const MaterialModel& base, const LinearSolve& linear_solve,
                        const StrengthReductionOptions& options) {
    // Whether the slope reaches equilibrium under gravity with strength factored
    // by srf. Re-solved from the unstressed state each trial.
    auto is_stable = [&](double srf) {
        MaterialModel m = base;
        factor_strength(m, srf);
        const std::vector<MaterialModel> materials = {m};
        try {
            const NewtonResult r = solve_nonlinear(
                mesh, dofs, materials, gravity_load, linear_solve, options.newton);
            return r.converged;
        } catch (...) {
            // A singular/!factorizable tangent at incipient collapse is itself
            // the failure signal -- treat it as an unstable (collapsed) slope.
            return false;
        }
    };

    double lo = options.srf_min, hi = options.srf_max;
    for (int i = 0; i < options.bisection_iterations; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (is_stable(mid))
            lo = mid;  // still stable -> factor of safety is higher
        else
            hi = mid;  // collapsed -> factor of safety is lower
    }
    return 0.5 * (lo + hi);
}

SafetyResult safety_analysis(const mesh::Mesh& mesh, const DofMap& dofs,
                             const Eigen::VectorXd& gravity_load,
                             const std::vector<MaterialModel>& materials,
                             const LinearSolve& linear_solve,
                             const StrengthReductionOptions& options,
                             const std::vector<GaussState>& initial_state,
                             const std::vector<char>& active_element,
                             const std::vector<MaterialProfile>& profile,
                             const Structures& structures, const Eigen::VectorXd& baseline,
                             const StructuralInit& init_struct) {
    // Factor every material's strength by srf and re-solve under gravity from the unstressed state;
    // failure to converge is the standard collapse signal. (This is robust for Mohr-Coulomb. Hardening
    // Soil from a stress-free state over-shoots the cap and is unreliable here -- HS Safety is gated
    // upstream until a path-stable scheme is in place; see build_problem.)
    //
    // The structures enter every trial as they are, except for what a trial changes about them: the
    // interface strengths (reduced below, per trial) and, when the search starts unstressed, the
    // interfaces' seeded normal stress (cleared here, once). See the header for why each.
    Structures base = structures;
    if (initial_state.empty()) {
        for (auto& ie : base.interfaces) ie.sigma_n0.fill(0.0);
        for (auto& ie : base.interfaces5) ie.sigma_n0.fill(0.0);
    }
    auto try_srf = [&](double srf, NewtonResult& out) -> bool {
        std::vector<MaterialModel> m = materials;
        for (auto& mm : m) factor_strength(mm, srf);
        // The COHESION GRADIENT must be reduced by srf too. c(y) = c_ref + c_inc (y_ref - y) is one
        // strength; factoring only c_ref would leave the deep soil at its full gradient strength and
        // report an FoS that is too HIGH -- unconservative, and invisible. (E_inc is a stiffness, not a
        // strength: phi-c reduction never touches it, so it rides through unchanged.)
        std::vector<MaterialProfile> p = profile;
        for (auto& pp : p) pp.c_inc /= srf;
        Structures s = base;
        for (auto& ie : s.interfaces) factor_interface_strength(ie.props, srf);
        for (auto& ie : s.interfaces5) factor_interface_strength(ie.props, srf);
        try {
            out = solve_nonlinear(mesh, dofs, m, gravity_load, linear_solve, options.newton,
                                  initial_state, active_element, s, {}, {}, p);
            return out.converged;
        } catch (...) {
            return false;  // singular tangent at incipient collapse = the failure signal
        }
    };

    double lo = options.srf_min, hi = options.srf_max;
    SafetyResult res;
    NewtonResult trial;

    // --- INCREMENTAL: one gravity solve, then the strength reduced from equilibrium ---------
    // The bisection below re-solves the whole self-weight for every trial, and half of its
    // trials are collapses, each of which only ends after the load stepping has been cut down to
    // its minimum with every retry exhausted: measured on the Griffiths-Lane slope, 7 collapsing
    // trials of 12 took 50 of 58 s (tri6) and most of 673 s (tri15). Reducing the strength from
    // each converged state meets the collapse once, at the end.
    //
    // It is taken only where every piece of state it carries forward is carried: from the
    // unstressed state, no structural element -- their independent trials are what the
    // structures' closed forms were verified on (KV-STR-010) -- and, from a parent phase's state,
    // everything, because that start carries the structures the way a chained phase does.
    //
    // FROM THE PARENT PHASE (initial_state + baseline given). The ground is already in
    // equilibrium at full strength -- the phases before built it -- so the search does not
    // re-solve the self-weight at all: it holds the parent's internal force as the constant
    // load (residual zero by the parent's own equilibrium, the chained-phase rule) and reduces
    // the strength from there. Measured on a geosynthetic-reinforced wall: re-solving the weight
    // from zero failed at full strength on 8 load steps, where the phase before had stood, and
    // the search then reported a factor of 0.80 -- below the 1.0 the model had just carried --
    // and 1.07 on 40 steps. A factor of safety is a property of the state the engineer built,
    // not of how many steps a second construction of it was given.
    const bool has_structures = !structures.plates.empty() || !structures.anchors.empty() ||
                                !structures.geogrids.empty() || !structures.interfaces.empty() ||
                                !structures.plates5.empty() || !structures.interfaces5.empty() ||
                                !structures.embedded_beams.empty();
    const bool from_parent =
        !initial_state.empty() && baseline.size() == (Eigen::Index)dofs.equation_count();
    const bool incremental =
        from_parent || (options.incremental && initial_state.empty() && !has_structures &&
                        std::getenv("KATAI_SRM_BISECTION") == nullptr);
    if (incremental) {
        NewtonResult base_run;
        bool stands = false;
        if (from_parent) {
            // SRF = 1 from the parent is a chained phase of its own: the parent's force held, and
            // whatever this phase's configuration changes (its load vector against the parent's
            // internal force) ramped in at full strength. Nothing changed = a nil step, which
            // converges at once; a change that the ground cannot carry is reported as that, not
            // as a factor of safety below one.
            try {
                const Eigen::VectorXd change = gravity_load - baseline;
                base_run = solve_nonlinear(mesh, dofs, materials, change, linear_solve,
                                           options.newton, initial_state, active_element, base, {},
                                           baseline, profile, init_struct);
                stands = base_run.converged;
            } catch (...) {
                stands = false;
            }
            if (!stands) {
                res.start_failed = true;
                return res;
            }
        } else {
            stands = try_srf(1.0, base_run);
        }
        if (stands) {
            // Equilibrium at full strength: FoS >= 1. Walk the strength down from it.
            res.ok = true;
            res.mechanism = base_run;
            double srf = 1.0;
            std::vector<GaussState> states = base_run.gauss_states;
            Eigen::VectorXd u_total = base_run.displacement;
            const Eigen::VectorXd no_load = Eigen::VectorXd::Zero(gravity_load.size());
            const Eigen::VectorXd& held = gravity_load;   // equilibrated at SRF = 1 either way
            const bool parent_datum =
                from_parent && init_struct.u_datum.size() == (Eigen::Index)dofs.equation_count();
            NewtonOptions step_opt = options.newton;
            step_opt.load_steps = 1;
            step_opt.min_step_fraction = 1.0;   // the load is constant: a failed step is a collapse
            // Every reduction step is held to at least the search's own 1e-3, whatever the phase
            // asks for. A step accepted with a looser residual carries its out-of-balance force
            // into the next one, and along the path that compounds: measured on the
            // Griffiths-Lane slope, 1e-1 on the steps inflated the factor by +26% (1e-2 by +1.7%)
            // where independent trials had given +0.6%. A looser setting still applies to the
            // self-weight solve, and the driver says so (K2D-A006).
            step_opt.tolerance = std::min(step_opt.tolerance, 1.0e-3);
            // A step that fails is not a collapse: a strength drop too large for the iteration
            // budget fails just the same, and measured on a 0.5 m mesh of the Griffiths-Lane slope
            // the first step (1.0 -> 1.1) did exactly that while every smaller step towards 1.1
            // converged in three or four iterations -- a search that kept 1.1 as its upper bound
            // reported 1.100 where the slope stands to 1.361. So a failure only halves the step,
            // the search moves on from every new equilibrium, and the factor is bracketed only by a
            // step no larger than the resolution failing from the last converged state.
            // The structures continue from where the path has taken them: the parent's datum plus
            // what the search has moved so far, with their committed plastic state.
            const auto carry_now = [&]() {
                StructuralInit carry;
                if (has_structures) {
                    const NewtonResult& last = res.mechanism;
                    carry.u_datum = parent_datum ? Eigen::VectorXd(init_struct.u_datum)
                                                 : Eigen::VectorXd::Zero(dofs.equation_count());
                    for (int g = 0; g < dofs.total_dofs(); ++g) {
                        const int eq = dofs.equation(g);
                        if (eq >= 0) carry.u_datum(eq) += u_total(g);
                    }
                    carry.anchor_plastic = last.anchor_plastic;
                    carry.geogrid_plastic = last.geogrid_plastic;
                    carry.interface_slip = last.interface_slip;
                    carry.interface5_slip = last.interface5_slip;
                    carry.embedded_skin_slip = last.embedded_skin_slip;
                    carry.embedded_foot_slip = last.embedded_foot_slip;
                    carry.plate_plastic = last.plate_plastic;
                    carry.plate5_plastic = last.plate5_plastic;
                }
                return carry;
            };
            // THE PATH IS FOLLOWED BY ARC LENGTH (NewtonOptions::srf_reduce; the method is in
            // solve_nonlinear). A search that prescribes the factor at every step asks, near the
            // limit, for equilibria that exist but lie far along a developing mechanism: measured
            // on a c-phi slope (6-noded, 1 m), the last 0.5% of the factor carried the displacement
            // from 9.5 to 78 mm, those steps crawled for hundreds of iterations or were given up,
            // and what the search reported depended on how tightly its steps were solved -- 1.3903
            // at 1e-3 and 1.3838 at 1e-4 with steps of 0.01, 1.3874 with the growing steps -- every
            // one above the converged answer. With the factor an unknown the same slope gives
            // 1.38680, 1.38665 and 1.38680 at 1e-3, 1e-4 and 1e-5. The whole path is followed so:
            // prescribing the first part and continuing from the first failed step was measured
            // too, and with psi < phi the factor then depends on those first coarse steps (1.3820
            // at 1e-5). The prescribed search below remains for a path that cannot be started.
            //
            // WITH STRUCTURES TOO. A plate, an anchor, a geogrid or a pile keeps its own strength
            // while the soil's is reduced, so once the ground has formed a mechanism the structure
            // takes load on and the path goes on rising with the displacement for a long way --
            // which is why the continuation was at first kept from such models: on an anchored
            // excavation (6-noded, 2 m) it walked to 3.8 m for 1.852 where the prescribed search
            // stopped at 1.847 in a quarter of the time. But the prescribed search stops where the
            // iteration gives a step up, and that is not a property of the ground: measured on a
            // dewatered excavation between anchored walls (6-noded, 2 m), 1.555 at 120 mm of
            // movement with one build of the solver and 1.615 at 700 mm with a better one, each
            // step properly converged, on a path that the continuation follows to its plateau at
            // 1.630 -- 1.53 at 95 mm, 1.585 at 235 mm, 1.62 at 800 mm. Both prescribed answers were
            // points on the way. The plateau is the answer that does not move with the solver, so
            // it is the one reported; the prescribed search remains for a path that cannot be
            // started. (Reading the knee of the path instead -- where the ground's mechanism has
            // formed -- was measured and is wrong for a system: 1.726 for a block held by an anchor
            // prestressed to half its capacity, whose factor is 2.0025 once the anchor yields.)
            const auto continue_by_arc = [&](double first_step, NewtonResult& r) {
                NewtonOptions arc_opt = step_opt;
                arc_opt.srf_reduce = [](std::vector<MaterialModel>& mm, Structures& ss,
                                        std::vector<MaterialProfile>& pp, double sf) {
                    for (auto& x : mm) factor_strength(x, sf);
                    for (auto& x : pp) x.c_inc /= sf;
                    for (auto& ie : ss.interfaces) factor_interface_strength(ie.props, sf);
                    for (auto& ie : ss.interfaces5) factor_interface_strength(ie.props, sf);
                };
                arc_opt.srf_start = srf;
                arc_opt.srf_first_step = first_step;
                arc_opt.srf_max = std::max(options.srf_max, options.srf_ceiling);
                arc_opt.srf_resolution = options.fos_resolution;
                try {
                    r = solve_nonlinear(mesh, dofs, materials, no_load, linear_solve, arc_opt,
                                        states, active_element, base, {}, held, profile,
                                        carry_now());
                    // A continuation that stopped without finding the limit (the path could not
                    // be started from here) leaves the search to the prescribed steps.
                    return r.srf_limit || r.srf >= arc_opt.srf_max;
                } catch (...) {
                    return false;
                }
            };
            {
                NewtonResult ra;
                if (continue_by_arc(options.msf_step, ra)) {
                    res.fos = std::min(std::max(ra.srf, srf), std::max(options.srf_max, options.srf_ceiling));
                    res.bracketed = ra.srf_limit;
                    ra.displacement += u_total;
                    res.mechanism = std::move(ra);
                    return res;
                }
            }
            double step = options.msf_step;
            double failed_at = -1.0;            // set only by a step at the resolution failing
            bool ever_failed = false;
            for (int guard = 0; guard < 400; ++guard) {
                double next = srf + step;
                const double top = std::max(options.srf_max, options.srf_ceiling);
                if (next > top) {
                    if (srf >= top) break;
                    next = top;
                }
                std::vector<MaterialModel> m = materials;
                for (auto& mm : m) factor_strength(mm, next);
                std::vector<MaterialProfile> p = profile;
                for (auto& pp : p) pp.c_inc /= next;
                Structures s = base;
                for (auto& ie : s.interfaces) factor_interface_strength(ie.props, next);
                for (auto& ie : s.interfaces5) factor_interface_strength(ie.props, next);
                const StructuralInit carry = carry_now();
                NewtonResult r;
                bool ok = false;
                try {
                    r = solve_nonlinear(mesh, dofs, m, no_load, linear_solve, step_opt, states,
                                        active_element, s, {}, held, p, carry);
                    ok = r.converged;
                } catch (...) {
                    ok = false;
                }
                if (ok) {
                    srf = next;
                    states = r.gauss_states;
                    u_total += r.displacement;
                    r.displacement = u_total;
                    res.mechanism = std::move(r);
                    step *= ever_failed ? 1.25 : 1.5;   // stride out again, more warily once burnt
                } else {
                    ever_failed = true;
                    if (next - srf <= options.fos_resolution * srf) {
                        failed_at = next;               // even the smallest step fails: collapse
                        res.bracketed = true;
                        break;
                    }
                    step = std::max(0.5 * (next - srf), 0.5 * options.fos_resolution * srf);
                }
                if (srf >= std::max(options.srf_max, options.srf_ceiling)) break;
            }
            res.fos = failed_at > 0.0 ? 0.5 * (srf + failed_at) : srf;
            return res;
        }
        // No equilibrium at full strength: the factor is below 1, and the bisection finds it.
        hi = std::min(hi, 1.0);
        res.bracketed = true;
    }

    const auto bisect = [&]() {
        for (int i = 0; i < options.bisection_iterations; ++i) {
            const double mid = 0.5 * (lo + hi);
            if (try_srf(mid, trial)) { lo = mid; res.mechanism = trial; res.ok = true; }
            else { hi = mid; res.bracketed = true; }  // a collapse was seen -> the FoS is finite/bracketed
        }
    };
    bisect();
    // Nothing collapsed anywhere below the upper bracket: [srf_min, srf_max] was a first guess,
    // not a limit on the answer. Widen it, doubling up to the ceiling, and bisect the interval in
    // which a trial first collapses. A search that saw a collapse never comes here, so every factor
    // below srf_max is found exactly as before.
    if (!res.bracketed && res.ok && hi < options.srf_ceiling) {
        double stable = hi;
        for (;;) {
            const double up = std::min(2.0 * stable, options.srf_ceiling);
            if (!try_srf(up, trial)) {
                lo = stable;
                hi = up;
                res.bracketed = true;
                bisect();
                break;
            }
            stable = up;
            res.mechanism = trial;
            lo = hi = up;
            if (up >= options.srf_ceiling) break;
        }
    }
    res.fos = 0.5 * (lo + hi);
    return res;
}

} // namespace katai::core
