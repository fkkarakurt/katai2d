// The analysis-driver bodies, compiled ONCE (section 5.2). As header-inline code every
// driver-including test paid ~220 s of codegen for these functions and the engine
// machinery they pull in; measured in docs/validation/performance-baseline.md.
// Batch 2 moved the ENGINE INCLUDE SET here too: the header now carries only the
// declaration closure, so consumers stop parsing the twenty-odd headers below.
#include <katai/jobs/driver.hpp>

#include <array>
#include <cmath>
#include <exception>
#include <memory>
#include <unordered_map>

#include <katai/analysis/consolidation.hpp>          // time-dependent Biot consolidation phase
#include <katai/analysis/dynamics.hpp>                // DynamicsSolveFactory (injected into the Dynamic strategy)
#include <katai/analysis/transient_flow.hpp>          // transient (un)saturated groundwater flow (W1/W2)
#include <katai/analysis/initial_stress.hpp>
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/seepage.hpp>               // pore/gravity from a flow head field
#include <katai/analysis/general_interface.hpp>    // split_mesh_at_segment, build_soil_interface
#include <katai/analysis/embedded_wall.hpp>      // build_embedded_wall, WallBuild
#include <katai/analysis/staged_construction.hpp> // split_mesh_at_wall, SeamPair
#include <katai/analysis/boundary_conditions.hpp>  // neutral edge BCs (Stage B3: engine-owned)
#include <katai/analysis/design_code.hpp>          // EC7/TBDY partial factors (design approaches)
#include <katai/analysis/hydraulic_boundary.hpp>   // phreatic line + flow-edge services (Stage B4)
#include <katai/analysis/interface_baseline.hpp>   // wished-in-place sigma_n0 baseline (Stage B9)
#include <katai/analysis/phase_solver/consolidation.hpp>  // consolidation phase strategy (Stage B9)
#include <katai/analysis/phase_solver/dynamic.hpp>        // dynamic (seismic) phase strategy (Stage B9)
#include <katai/analysis/phase_solver/fully_coupled.hpp>  // fully-coupled phase strategy (Stage B9)
#include <katai/analysis/phase_solver/safety.hpp>         // safety phase strategy (Stage B9)
#include <katai/analysis/phase_solver/static_phase.hpp>   // K0/gravity/Plastic phase strategy (Stage B9)
#include <katai/analysis/phase_solver/transient_flow.hpp> // transient-flow phase strategy (Stage B9)
#include <katai/analysis/structural_carry.hpp>     // parent structural-state carry (Stage B9)
#include <katai/analysis/structural_diagrams.hpp>  // DiagSpec/IfaceDiag + per-line diagrams (Stage B9)
#include <katai/analysis/structural_dynamics.hpp>  // SSI: structural K + plate mass in a Dynamic phase
#include <katai/analysis/structural_forces.hpp>    // plate/anchor/geogrid force diagrams (Output)
#include <katai/analysis/post/consolidation_recovery.hpp> // committed-stress update (Stage B8)
#include <katai/analysis/post/stress_recovery.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/interface.hpp>           // iface::InterfaceProps, interface_stiffness
#include <katai/fem/elements/plate.hpp>               // PlateProps + set_plate_mass
#include <katai/materials/linear_elastic.hpp>
#include <katai/mesh/boundary_extraction.hpp>      // collect_chain + boundary edges (Stage B5)
#include <katai/linsolve/direct_solver.hpp>

#include <katai/jobs/mesh_builder.hpp>   // point_in_polygon (anchor end-in-soil test)

namespace katai::app {

// Boundary extraction lives in the mesh module (Stage B5); aliases keep the body's spellings.
using katai::mesh::collect_chain;
using katai::mesh::BoundaryEdgeChain;
using katai::mesh::extract_boundary_edges;

// Consolidation stress recovery lives in the engine (Stage B8).
using katai::core::recover_consolidation_stress;

// ---- Body-only seams and schema mappings (moved from the header, section 5.2 batch 2;
// ---- nothing outside the driver uses them -- measured before the move).

// Two solve seams, each stated once instead of at every phase that needs one. Both
// exist for the same reason: the solver must outlive the call. The sparsity pattern
// is fixed by the mesh, so a solver that survives from one call to the next does the
// symbolic analysis once and every later factorization pays only the numerical cost.
// Which backend answers is decided at link time (katai/linsolve); every solve is
// residual-checked, so a factorization that "succeeded" on a singular matrix is
// refused rather than published as a displacement field.

// Newton seam: the tangent changes on every iteration, the pattern does not.
static katai::core::LinearSolve reusing_linear_solve(katai::linsolve::MatrixType type) {
    return [s = std::shared_ptr<katai::linsolve::DirectSolver>(
                katai::linsolve::make_direct_solver(type))](
               const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
        s->factorize(k);
        return s->solve(r);
    };
}

// Time-stepping seam: one factorization, then a back-solve per step.
static std::function<Eigen::VectorXd(const Eigen::VectorXd&)> factorize_once(
    katai::linsolve::MatrixType type, const katai::math::CsrMatrix& a) {
    std::shared_ptr<katai::linsolve::DirectSolver> s = katai::linsolve::make_direct_solver(type);
    s->factorize(a);
    return [s](const Eigen::VectorXd& b) { return s->solve(b); };
}

// Plate inertia lives with the plate element (Stage B6); this wrapper only supplies the
// analysis layer's gravity constant, because fem takes g as an argument rather than
// reaching upward for a constant.
static void set_plate_mass(katai::core::plate::PlateProps& pp, double w) {
    katai::core::plate::set_plate_mass(pp, w, kGravity);
}

// Phreatic-surface elevation at horizontal position x (engine service since Stage B4).
static double water_table_at(const model::Project& pr, double x) {
    if (!pr.has_water) return -1e30;
    return katai::core::phreatic_surface_at(pr.wx, pr.wy, x);
}

// Displacement BC seam (Stage B3): the schema's polygons walked once into the engine's
// neutral edge list. Free edges stay in the list on purpose: they always contributed
// to the nearest-edge tie tolerance.
static katai::core::EdgeFixity to_edge_fixity(model::BCType bc) {
    using M = model::BCType; using E = katai::core::EdgeFixity;
    switch (bc) {
        case M::Free:              return E::Free;
        case M::NormallyFixed:     return E::NormallyFixed;
        case M::HorizontallyFixed: return E::HorizontallyFixed;
        case M::VerticallyFixed:   return E::VerticallyFixed;
        case M::FullyFixed:        return E::FullyFixed;
    }
    return E::Free;
}

static std::vector<katai::core::BcEdge> bc_edges_from(const model::Project& pr) {
    std::vector<katai::core::BcEdge> edges;
    for (const auto& P : pr.polygons) {
        const int n = (int)P.x.size(); if (n < 3 || (int)P.edge_bc.size() != n) continue;
        double cx = 0, cy = 0; for (int k = 0; k < n; ++k) { cx += P.x[k]; cy += P.y[k]; } cx /= n; cy /= n;
        for (int i = 0; i < n; ++i) {
            const double ax = P.x[i], ay = P.y[i], bx = P.x[(i + 1) % n], by = P.y[(i + 1) % n];
            const double ex = bx - ax, ey = by - ay, l2 = ex * ex + ey * ey; if (l2 < 1e-18) continue;
            double nx = ey, ny = -ex;                        // edge normal
            const double mx = 0.5 * (ax + bx), my = 0.5 * (ay + by);
            if (nx * (mx - cx) + ny * (my - cy) < 0) { nx = -nx; ny = -ny; }   // outward
            const double nl = std::hypot(nx, ny);
            edges.push_back({ax, ay, bx, by, nx / nl, ny / nl,
                             to_edge_fixity((model::BCType)P.edge_bc[i])});
        }
    }
    return edges;
}

// Flow-boundary seam (Stage B4), the hydraulic twin of bc_edges_from. Closed edges stay
// in the list; only the apply step filters by kind. Missing head array defaults to 0.
static katai::core::FlowEdgeKind to_flow_kind(model::FlowBCType fb) {
    using M = model::FlowBCType; using F = katai::core::FlowEdgeKind;
    switch (fb) {
        case M::Closed:  return F::Closed;
        case M::Head:    return F::Head;
        case M::Seepage: return F::Seepage;
    }
    return F::Closed;
}

static std::vector<katai::core::FlowEdge> flow_edges_from(const model::Project& pr) {
    std::vector<katai::core::FlowEdge> edges;
    for (const auto& P : pr.polygons) {
        const int n = (int)P.x.size();
        if (n < 3 || (int)P.edge_flow.size() != n) continue;
        for (int i = 0; i < n; ++i) {
            const double h = (int)P.edge_head.size() == n ? P.edge_head[i] : 0.0;
            edges.push_back({P.x[i], P.y[i], P.x[(i + 1) % n], P.y[(i + 1) % n],
                             to_flow_kind((model::FlowBCType)P.edge_flow[i]), h});
        }
    }
    return edges;
}

// The original declaration scan, kept unfiltered on purpose: a flow BC declared on a
// polygon too degenerate to yield a matchable edge must still disable the top-drain
// fallback (the engine header states the same asymmetry from its side).
static bool any_flow_bc_declared(const model::Project& pr) {
    for (const auto& P : pr.polygons)
        for (int fb : P.edge_flow)
            if (fb == (int)model::FlowBCType::Head || fb == (int)model::FlowBCType::Seepage)
                return true;
    return false;
}

// Drainage-boundary mask / prescribed-head boundary for time-dependent flow phases.
static std::vector<char> flow_drained_nodes(const model::Project& pr, const katai::mesh::Mesh& mesh,
                                            const std::vector<char>& act, double yscale) {
    return katai::core::flow_drained_nodes(flow_edges_from(pr), any_flow_bc_declared(pr),
                                           mesh, act, yscale);
}

static void flow_head_nodes(const model::Project& pr, const katai::mesh::Mesh& mesh,
                            std::vector<char>& is_presc, std::vector<double>& head_val) {
    katai::core::flow_head_nodes(flow_edges_from(pr), mesh, is_presc, head_val);
}

// Enum seams: the schema's file-stable enums resolved to the engine's, once, here.
static katai::core::DesignApproach to_core_design_approach(model::DesignApproach da) {
    using M = model::DesignApproach; using C = katai::core::DesignApproach;
    switch (da) {
        case M::EC7_DA1_C1:       return C::EC7_DA1_C1;
        case M::EC7_DA1_C2:       return C::EC7_DA1_C2;
        case M::EC7_DA2:          return C::EC7_DA2;
        case M::EC7_DA3:          return C::EC7_DA3;
        case M::TBDY2018_Static:  return C::TBDY2018_Static;
        case M::TBDY2018_Seismic: return C::TBDY2018_Seismic;
        case M::None:             return C::None;
    }
    return C::None;
}

static katai::core::SeismicWave to_core_seismic_wave(model::SeismicWave w) {
    using M = model::SeismicWave; using C = katai::core::SeismicWave;
    switch (w) {
        case M::Harmonic: return C::Harmonic;
        case M::Ricker:   return C::Ricker;
        case M::Record:   return C::Record;
    }
    return C::Harmonic;
}

// The construction seam between the schema and the constitutive catalogue: schema
// fields ferried into the registry's neutral parameter block (declared in the header).
katai::core::MaterialParams to_material_params(const model::Material& m) {
    katai::core::MaterialParams p;
    const double kPi = std::acos(-1.0);
    p.E = m.E;
    p.nu = m.nu;
    p.c = m.c;
    p.phi_rad = m.phi * kPi / 180.0;
    p.psi_rad = m.psi * kPi / 180.0;
    p.tension_cutoff = m.tension_cutoff;
    p.tensile_strength = m.tensile_strength;
    p.E50_ref = m.E50ref;
    p.Eur_ref = m.Eurref;
    p.Eoed_ref = m.Eoedref;
    p.m = m.m;
    p.p_ref = m.p_ref;
    p.Rf = m.Rf;
    p.nu_ur = m.nu_ur;
    p.G0_ref = m.G0ref;
    p.gamma07 = m.gamma07;
    p.lam_star = m.lam_star;
    p.kap_star = m.kap_star;
    p.mu_star = m.mu_star;
    p.k0nc_auto = m.k0nc_auto;
    p.k0nc = m.k0nc;
    switch (m.drainage) {
        case model::Drainage::Drained:    p.drainage = katai::core::DrainageClass::Drained; break;
        case model::Drainage::Undrained:  p.drainage = katai::core::DrainageClass::UndrainedA; break;
        case model::Drainage::UndrainedB: p.drainage = katai::core::DrainageClass::UndrainedB; break;
        case model::Drainage::NonPorous:  p.drainage = katai::core::DrainageClass::NonPorous; break;
    }
    return p;
}

SolveResult solve_gravity_le(const model::Project& pr, const katai::mesh::Mesh& mesh_in,
                             InitialPhase phase, const Eigen::VectorXd* flow_head,
                             const PhaseIO& io) {
    SolveResult R;
    katai::mesh::Mesh mesh = mesh_in;   // working copy (split for embedded walls below)
    if (mesh.element_count == 0) { R.message = "Empty mesh."; return R; }
    if (pr.materials.empty())    { R.message = "No materials."; return R; }

    // Material tables indexed by material id. The constitutive model follows the GUI Material.model
    // (PLAXIS: the selected model governs), resolved by name in the constitutive catalogue --
    // parameter wiring, K0^NC memory, cap calibration and the Undrained (A/B) machinery all live
    // with the models in katai/materials/registry.hpp now. `mats` (LinearElastic) is kept only as
    // a fallback; the real solve uses `models`. nonlinear_soil drives load-stepping / solver choice.
    std::vector<katai::core::LinearElastic> mats;
    std::vector<katai::core::MaterialModel> models;
    std::vector<double> gamma;      // unsaturated (above water table)
    std::vector<double> gamma_sat;  // saturated (below water table)
    bool nonlinear_soil = false;
    bool has_hardening = false;  // hardening family present (needs PLAXIS-realistic tol)
    bool has_softsoil = false;   // soft-soil family present (FD tangent -> same step/tol class)
    for (const auto& m : pr.materials) {
        katai::core::LinearElastic le; le.youngs_modulus = m.E; le.poisson_ratio = m.nu;
        mats.push_back(le);
        const katai::core::ModelEntry* entry = katai::core::find_model(constitutive_name(m.model));
        if (!entry) {   // unreachable for schema enums; refused by name, never substituted (R2)
            R.message = std::string("Unknown constitutive model '") + constitutive_name(m.model) +
                        "' -- this build cannot solve it.";
            return R;
        }
        models.push_back(entry->build(to_material_params(m)));
        nonlinear_soil |= entry->nonlinear;
        has_hardening |= entry->hardening_family;
        has_softsoil |= entry->softsoil_family;
        gamma.push_back(m.gamma_unsat);
        // NonPorous (PLAXIS): a non-porous material holds NO water -- below the water table its
        // total weight is still gamma_unsat (giving a concrete block gamma_sat would silently
        // saturate it). This one line fixes the gravity_phreatic + gravity_from_head +
        // consolidation dF paths alike.
        gamma_sat.push_back(m.drainage == model::Drainage::NonPorous ? m.gamma_unsat : m.gamma_sat);
    }
    // NonPorous material flags (the shared source for element masks + K0 seeding + flow).
    std::vector<char> mat_nonporous;
    bool any_nonporous = false;
    for (const auto& m : pr.materials) {
        const char np = m.drainage == model::Drainage::NonPorous ? 1 : 0;
        mat_nonporous.push_back(np);
        any_nonporous |= (np != 0);
    }
    for (int e = 0; e < mesh.element_count; ++e)
        if (mesh.element_material[e] < 0 || mesh.element_material[e] >= (int)mats.size()) {
            R.message = "An element has no valid material; assign a material to every soil region.";
            return R;
        }
    // Honest refusals live with the models: each catalogue entry names the parameter/drainage
    // combinations it cannot solve without silently changing meaning (SS/SSC + Undrained (A/B),
    // HS + Undrained (B) -- both audit findings of the silent-wrong class). The first offending
    // material wins and its message reaches the user verbatim.
    for (const auto& m : pr.materials) {
        const std::string why =
            katai::core::find_model(constitutive_name(m.model))->validate(to_material_params(m));
        if (!why.empty()) { R.message = why; return R; }
    }
    // The plate Mp/Np plastic hinge is NOW IN THE CORE (the MMM 18.3 diamond;
    // structural-plate-formulation section 10) -- active in the static family and in nonlinear
    // dynamics. The remaining honest gate: elastoplastic is set but BOTH capacities are 0 (the
    // GUI default) -> the user expects yielding while the plate would solve unbounded-elastic
    // (the silent-wrong audit class). Demand a capacity or have the flag removed.
    for (const auto& s : pr.structs)
        if (s.kind == model::StructKind::Plate && s.material >= 0 &&
            s.material < (int)pr.plates.size() && pr.plates[s.material].elastoplastic &&
            !(pr.plates[s.material].Mp > 0.0) && !(pr.plates[s.material].Np > 0.0)) {
            R.message = "Plate material '" + pr.plates[s.material].name +
                        "' is marked elastoplastic but both Mp and Np are 0 (= unlimited): the "
                        "plate would run elastic while the input promises a plastic hinge. Enter "
                        "Mp and/or Np (> 0), or uncheck 'elastoplastic' for an elastic plate.";
            return R;
        }
    // Depth-varying E'(y) / c'(y), evaluated per stress point by every deformation phase's own assembly.
    const std::vector<katai::core::MaterialProfile> profiles = build_profiles(pr);

    // --- K0 / water infrastructure (used by initial stress, interface seeding and post-processing) -
    const double kPi = std::acos(-1.0);
    const bool water = pr.has_water && pr.wx.size() >= 2;
    const bool axi = pr.axisymmetric;   // axisymmetric (r-z): x = radius, integration is r-weighted
    std::vector<double> k0_by_mat;
    for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
        const auto& m = pr.materials[mi];
        double k0 = m.k0_auto ? (1.0 - std::sin(m.phi * kPi / 180.0)) : m.k0;
        // Overconsolidation raises the automatic K0; the formula and its passive clamp are
        // engine-owned now (Stage B9, initial_stress.hpp: k0_overconsolidated -- audit finding:
        // OCR used to be silently ignored by the MC/LE geostatic field). nu is nu_ur for the
        // advanced families (decided from the engine's models table, not the schema enum), nu
        // for LE/MC; a manually entered K0 (k0_auto = false) always wins.
        if (m.k0_auto && m.oc_mode == 1 && m.OCR > 1.0) {
            const bool adv = models[mi].type != katai::core::MaterialType::LinearElastic &&
                             models[mi].type != katai::core::MaterialType::MohrCoulomb;
            k0 = katai::core::k0_overconsolidated(k0, m.OCR, adv ? m.nu_ur : m.nu,
                                                  std::sin(m.phi * kPi / 180.0));
        }
        k0_by_mat.push_back(k0);
    }
    // Staged construction: active objects of THIS phase (null config = everything active).
    // All the geostatic infrastructure below (material/overburden/surface) sees only the
    // ACTIVE soil -- an excavated region carries no weight and the pit floor is the surface.
    const auto poly_on = [&io](size_t i) { return !io.config || io.config->active_poly(i); };
    const auto material_at = [&pr, &poly_on](double x, double y) {
        int mat = -1;
        for (size_t i = 0; i < pr.polygons.size(); ++i)
            if (poly_on(i) && point_in_polygon(x, y, pr.polygons[i])) mat = pr.polygons[i].material;
        return mat;   // last wins
    };
    const auto eff_unit_weight = [&pr, water, &material_at](double x, double y) {
        const int mat = material_at(x, y);
        if (mat < 0 || mat >= (int)pr.materials.size()) return 0.0;
        const auto& m = pr.materials[mat];
        // SLICE weight, for the sigma'_v integral at a POROUS target: (gamma_total - gamma_w)
        // below water. A NonPorous slice always weighs gamma_unsat (it holds no water); since
        // the pore pressure at a porous point below it still grows hydrostatically through the
        // slice, its effective contribution is gamma_unsat - gamma_w. A NON-POROUS target is
        // itself corrected to TOTAL stress by the pore re-add in the K0 seeding
        // (K0LayeredOptions::nonporous/pore).
        const bool np = m.drainage == model::Drainage::NonPorous;
        const double g_total = (!np && water && y < water_table_at(pr, x)) ? m.gamma_sat
                                                                          : m.gamma_unsat;
        return (water && y < water_table_at(pr, x)) ? (g_total - kGammaWater) : g_total;
    };
    const auto ground_surface = [&pr, &poly_on](double x) {
        double top = -1e30;
        for (size_t pi = 0; pi < pr.polygons.size(); ++pi) {
            if (!poly_on(pi)) continue;
            const auto& P = pr.polygons[pi];
            const int n = (int)P.x.size();
            for (int i = 0; i < n; ++i) {
                const double ax = P.x[i], ay = P.y[i], bx = P.x[(i + 1) % n], by = P.y[(i + 1) % n];
                const double lo = std::min(ax, bx), hi = std::max(ax, bx);
                if (x < lo - 1e-9 || x > hi + 1e-9) continue;
                const double dx = bx - ax;
                const double t = std::fabs(dx) < 1e-12 ? 0.0 : (x - ax) / dx;
                top = std::fmax(top, ay + t * (by - ay));
            }
        }
        return top;
    };
    // Effective vertical stress sigma'_v (<=0) at a point: vertical integral of gamma' above it.
    const auto eff_sigma_v = [&](double x, double y) {
        const double y_surf = ground_surface(x);
        if (y_surf <= y) return 0.0;
        const int ns = 400;
        const double dt = (y_surf - y) / ns;
        double sv = 0.0;
        for (int i = 0; i < ns; ++i) sv += eff_unit_weight(x, y + (i + 0.5) * dt) * dt;
        return -sv;
    };

    // --- Embedded walls: a plate with a positive/negative interface flag becomes a barrier -- the mesh
    // is split along the wall line and the plate sits on independent DOFs, joined to each soil side by a
    // Coulomb interface (PLAXIS plate + interface). K0 seeding keeps the wished-in-place wall in
    // equilibrium (no spurious installation movement). Works on tri6 AND tri15, at ANY orientation: a
    // vertical wall uses the validated x-column split (bit-for-bit), a non-vertical wall the general
    // split_mesh_at_segment (seam keyed by arc length s so the builders' sort stays correct).
    struct WallSpec {
        std::vector<katai::core::SeamPair> seam;
        katai::core::plate::PlateProps pp;
        katai::core::iface::InterfaceProps ip;
        double toe_x, toe_y;     // anchored (deeper) end -> shared toe node
        double nx, ny;           // unit normal (orientation-aware K0 seed: sigma_n0=(K0 nx^2+ny^2) sigma'_v)
        int order;               // 6 or 15 (element order -> build_embedded_wall vs _wall5)
        int soil_mat;
        std::string name;        // drawn element name (for the force-diagram output)
    };
    std::vector<WallSpec> walls;
    std::vector<char> plate_is_wall(pr.structs.size(), 0);
    constexpr bool kEnableEmbeddedWall = true;
    if (kEnableEmbeddedWall) {
        const int order = mesh.nodes_per_element;
        for (size_t si = 0; si < pr.structs.size(); ++si) {
            const auto& s = pr.structs[si];
            if (s.kind != model::StructKind::Plate || !(s.iface_pos || s.iface_neg)) continue;
            // An embedded wall splits the mesh -- the geometry cannot change between phases, so
            // per-phase (de)activation of walls is not supported yet; guard honestly.
            if (io.config && !io.config->active_struct(si)) {
                R.message = "Activating/deactivating an embedded wall (plate with interface) per "
                            "phase is not supported yet -- keep it active in every phase.";
                return R;
            }
            const double dx = s.x2 - s.x1, dy = s.y2 - s.y1;
            const double len = std::hypot(dx, dy);
            if (len < 1e-9) continue;
            WallSpec w;
            w.name = s.name;
            w.order = order;
            // Toe = the DEEPER endpoint (min y; horizontal -> min x); the plate hangs from it (shared node).
            const bool a_is_toe = std::fabs(dy) > 1e-9 ? (s.y1 <= s.y2) : (s.x1 <= s.x2);
            const double tx = a_is_toe ? s.x1 : s.x2, ty = a_is_toe ? s.y1 : s.y2;   // toe end
            const double ux = a_is_toe ? s.x2 : s.x1, uy = a_is_toe ? s.y2 : s.y1;   // top end
            w.toe_x = tx; w.toe_y = ty;
            w.nx = (uy - ty) / len; w.ny = -(ux - tx) / len;   // unit normal of the wall line
            if (std::fabs(dx) <= 1e-6 * len) {                 // VERTICAL: validated x-column split
                w.seam = katai::core::split_mesh_at_wall(mesh, 0.5 * (s.x1 + s.x2),
                                                         std::min(s.y1, s.y2), std::max(s.y1, s.y2));
            } else {                                           // GENERAL: split along the segment, key by s
                const auto seg = katai::core::split_mesh_at_segment(mesh, tx, ty, ux, uy, 1e-3, len + 1.0);
                for (const auto& p : seg) w.seam.push_back({p.orig, p.dup, p.s});
            }
            if ((int)w.seam.size() < (order == 15 ? 4 : 2)) continue;   // line not on mesh edges -> bonded
            // Plate stiffness from the plate material.
            if (s.material >= 0 && s.material < (int)pr.plates.size()) {
                const auto& pm = pr.plates[s.material];
                w.pp.EA = pm.EA; w.pp.EI = pm.EI; w.pp.nu = pm.nu;
                set_plate_mass(w.pp, pm.w);       // inertia (Dynamic phases only; ignored by static)
                if (pm.elastoplastic) {           // M-N hinge caps (<=0 => unlimited, geogrid rule)
                    w.pp.Mp = pm.Mp > 0.0 ? pm.Mp : -1.0;
                    w.pp.Np = pm.Np > 0.0 ? pm.Np : -1.0;
                }
            }
            // Interface strength/stiffness from the adjacent soil + Rinter (PLAXIS strength reduction).
            const double mx = 0.5 * (s.x1 + s.x2), my = 0.5 * (s.y1 + s.y2);
            int smat = s.iface_material >= 0 ? s.iface_material : material_at(mx + 1e-3 * w.nx, my + 1e-3 * w.ny);
            if (smat < 0) smat = material_at(mx - 1e-3 * w.nx, my - 1e-3 * w.ny);
            w.soil_mat = smat;
            const auto& sm = pr.materials[smat >= 0 ? smat : 0];
            const double Rinter = sm.rinter_rigid ? 1.0 : sm.Rinter;
            const double G = sm.E / (2.0 * (1.0 + sm.nu));
            const double avg = len / std::max<size_t>(1, w.seam.size() / (order == 15 ? 4 : 2));
            katai::core::iface::interface_stiffness(Rinter, G, avg, 0.1, w.ip.kn, w.ip.ks);
            w.ip.c_i = Rinter * sm.c;
            w.ip.phi_i = std::atan(Rinter * std::tan(sm.phi * kPi / 180.0));
            // Interface tensile strength = R * sigma_t (the PLAXIS rule) -- only while the
            // material's tension cutoff is on; otherwise 0 (the interface default: carries no
            // tension -- the safe side). Previously the material's sigma_t never reached the
            // interface at all (an audit finding).
            w.ip.sigma_t = sm.tension_cutoff ? Rinter * std::max(0.0, sm.tensile_strength) : 0.0;
            plate_is_wall[si] = 1;
            walls.push_back(std::move(w));
        }
    }

    // --- Standalone interfaces (StructKind::Interface): a Coulomb slip surface in the soil at ANY
    // orientation, on tri6 OR tri15. Mirror the wall: split the mesh along the line (BEFORE the DofMap)
    // and build a soil-soil joint (build_soil_interface, after the DofMap) connecting the two sides.
    // Strength/stiffness come from the adjacent soil + Rinter. Boundary nodes that get duplicated have
    // their fixity propagated (apply_boundary_conditions matches by coordinate). (interface-formulation.md.)
    struct IfaceSpec {
        std::vector<katai::core::SegSeam> seam;
        katai::core::iface::InterfaceProps ip;
        double nx, ny;          // unit normal (orientation-aware K0 seed)
        int soil_mat;
        std::string name;
    };
    std::vector<IfaceSpec> soil_ifaces;
    {
        std::vector<char> is_bnode(mesh.node_count, 0);
        for (int n : mesh.boundary_nodes) is_bnode[n] = 1;
        const int order = mesh.nodes_per_element;
        for (size_t si = 0; si < pr.structs.size(); ++si) {
            const auto& s = pr.structs[si];
            if (s.kind != model::StructKind::Interface) continue;
            if (io.config && !io.config->active_struct(si)) {
                R.message = "Activating/deactivating an interface per phase is not supported yet -- keep "
                            "it active in every phase (it splits the mesh, which is fixed across phases).";
                return R;
            }
            const double dx = s.x2 - s.x1, dy = s.y2 - s.y1;
            const double len = std::hypot(dx, dy);
            if (len < 1e-9) continue;
            IfaceSpec sp;
            sp.name = s.name.empty() ? "Interface" : s.name;
            sp.nx = dy / len; sp.ny = -dx / len;   // unit normal (matches split_mesh_at_segment frame)
            sp.seam = katai::core::split_mesh_at_segment(mesh, s.x1, s.y1, s.x2, s.y2, -1e-6, len + 1e-6);
            const int per_edge = order == 15 ? 4 : 2;
            if ((int)sp.seam.size() < per_edge + 1) continue;   // line not on mesh edges -> leave bonded
            // Propagate boundary fixity to the duplicated boundary nodes (else a side daylighting on a
            // fixed boundary would be left unconstrained -> support hole).
            for (const auto& p : sp.seam) if (is_bnode[p.orig]) mesh.boundary_nodes.push_back(p.dup);
            // Interface strength/stiffness from the adjacent soil + Rinter (PLAXIS strength reduction).
            const double mx = 0.5 * (s.x1 + s.x2), my = 0.5 * (s.y1 + s.y2);
            int smat = s.iface_material >= 0 ? s.iface_material : material_at(mx + 1e-3 * sp.nx, my + 1e-3 * sp.ny);
            if (smat < 0) smat = material_at(mx - 1e-3 * sp.nx, my - 1e-3 * sp.ny);
            sp.soil_mat = smat;
            const auto& sm = pr.materials[smat >= 0 ? smat : 0];
            const double Rinter = sm.rinter_rigid ? 1.0 : sm.Rinter;
            const double G = sm.E / (2.0 * (1.0 + sm.nu));
            const int nedges = std::max<int>(1, ((int)sp.seam.size() - 1) / per_edge);
            katai::core::iface::interface_stiffness(Rinter, G, len / nedges, 0.1, sp.ip.kn, sp.ip.ks);
            sp.ip.c_i = Rinter * sm.c;
            sp.ip.phi_i = std::atan(Rinter * std::tan(sm.phi * kPi / 180.0));
            sp.ip.sigma_t = sm.tension_cutoff ? Rinter * std::max(0.0, sm.tensile_strength) : 0.0;
            soil_ifaces.push_back(std::move(sp));
        }
    }

    // Element activity from the polygon activation of THIS phase (centroid containment, last
    // polygon wins -- matching the material assignment). Empty mask = everything active (the
    // classic single-phase path, bit-identical).
    std::vector<char> act;
    if (io.config) {
        bool any_off = false;
        for (size_t i = 0; i < pr.polygons.size(); ++i)
            if (!io.config->active_poly(i)) any_off = true;
        if (any_off) {
            act.assign(mesh.element_count, 1);
            for (int e = 0; e < mesh.element_count; ++e) {
                const int n0 = mesh.node_of(e, 0), n1 = mesh.node_of(e, 1), n2 = mesh.node_of(e, 2);
                const double cx = (mesh.x[n0] + mesh.x[n1] + mesh.x[n2]) / 3.0;
                const double cy = (mesh.y[n0] + mesh.y[n1] + mesh.y[n2]) / 3.0;
                int owner = -1;
                for (size_t i = 0; i < pr.polygons.size(); ++i)
                    if (point_in_polygon(cx, cy, pr.polygons[i])) owner = (int)i;   // last wins
                act[e] = owner < 0 ? 1 : (io.config->active_poly((size_t)owner) ? 1 : 0);
            }
        }
    }
    const auto struct_on = [&io](size_t i) { return !io.config || io.config->active_struct(i); };
    const auto load_on = [&io](size_t i) { return !io.config || io.config->active_load(i); };

    katai::core::DofMap dofs(mesh.node_count, 2);

    // Structural plates: each plate line is a conforming node chain (corner, mid, ...). Translational
    // DOFs are shared with the soil; one rotation extra-DOF is allocated per chain node.
    // Diagram bookkeeping lives in the engine now (Stage B9, katai/analysis/
    // structural_diagrams.hpp): DiagSpec / IfaceDiag remember which contiguous slice of each
    // structures vector belongs to which drawn line, so force diagrams (PLAXIS Output M/Q/N)
    // and interface results (tau/sigma_n/slip) can be produced per line.
    using katai::core::DiagSpec;
    using katai::core::IfaceDiag;
    std::vector<DiagSpec> diag_specs;
    std::vector<IfaceDiag> iface_diags;
    katai::core::Structures structures;
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Plate || plate_is_wall[si]) continue;   // walls built below
        if (!struct_on(si)) continue;   // not installed in this phase
        const std::vector<int> chain = collect_chain(mesh, s.x1, s.y1, s.x2, s.y2);
        if (chain.size() < 3 || chain.size() % 2 == 0) continue;
        katai::core::plate::PlateProps pp;
        if (s.material >= 0 && s.material < (int)pr.plates.size()) {
            const auto& pm = pr.plates[s.material]; pp.EA = pm.EA; pp.EI = pm.EI; pp.nu = pm.nu;
            set_plate_mass(pp, pm.w);             // inertia (Dynamic phases only; ignored by static)
            if (pm.elastoplastic) {               // M-N hinge caps (<=0 => unlimited, geogrid rule)
                pp.Mp = pm.Mp > 0.0 ? pm.Mp : -1.0;
                pp.Np = pm.Np > 0.0 ? pm.Np : -1.0;
            }
        }
        std::vector<int> rot(chain.size());
        for (size_t i = 0; i < chain.size(); ++i) rot[i] = dofs.add_extra_dof();
        const size_t p0 = structures.plates.size();
        for (size_t e = 0; 2 * e + 2 < chain.size(); ++e) {
            const int A = chain[2 * e], B = chain[2 * e + 2], M = chain[2 * e + 1];
            structures.plates.push_back(katai::core::PlateElement{
                {A, B, M}, {rot[2 * e], rot[2 * e + 2], rot[2 * e + 1]}, pp});
        }
        if (structures.plates.size() > p0)
            diag_specs.push_back({0, s.name, p0, structures.plates.size()});
    }

    // Geogrids: tension-only axial membranes. Like plates, each line is a conforming node chain;
    // translational DOFs are shared with the soil and there is no rotation DOF (no bending).
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Geogrid || !struct_on(si)) continue;
        const std::vector<int> chain = collect_chain(mesh, s.x1, s.y1, s.x2, s.y2);
        if (chain.size() < 3 || chain.size() % 2 == 0) continue;
        katai::core::geogrid::GeogridProps gp;
        if (s.material >= 0 && s.material < (int)pr.geogrids.size()) {
            const auto& gm = pr.geogrids[s.material];
            gp.EA = gm.EA;
            gp.Np = (gm.elastoplastic && gm.Np > 0.0) ? gm.Np : -1.0;   // <=0 => unlimited
        }
        const size_t g0 = structures.geogrids.size();
        for (size_t e = 0; 2 * e + 2 < chain.size(); ++e)
            structures.geogrids.push_back(katai::core::GeogridElement{
                {chain[2 * e], chain[2 * e + 2], chain[2 * e + 1]}, gp});
        if (structures.geogrids.size() > g0)
            diag_specs.push_back({2, s.name, g0, structures.geogrids.size()});
    }

    // Anchors: single axial spring (node-to-node when both ends sit in the soil; fixed-end when one
    // end is drawn outside the soil -> that end is a fixed far anchor point). EA/Fmax are divided by
    // the out-of-plane spacing to get the per-metre (plane-strain) values (PLAXIS MMM §18.1).
    const auto nearest_node = [&](double x, double y) {
        int best = -1; double bd = 1e300;
        for (int n = 0; n < mesh.node_count; ++n) {
            const double d = std::hypot(mesh.x[n] - x, mesh.y[n] - y);
            if (d < bd) { bd = d; best = n; }
        }
        return best;
    };
    const auto in_soil = [&](double x, double y) {
        for (const auto& P : pr.polygons) if (point_in_polygon(x, y, P)) return true;
        return false;
    };
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Anchor || !struct_on(si)) continue;
        katai::core::AnchorElement an;
        double Ls = 1.0;
        if (s.material >= 0 && s.material < (int)pr.anchors.size()) {
            const auto& am = pr.anchors[s.material];
            Ls = am.Lspacing > 1e-9 ? am.Lspacing : 1.0;
            an.EA = am.EA / Ls;
            if (am.elastoplastic) {
                an.Fmax_tens = am.Fmax_tens > 0.0 ? am.Fmax_tens / Ls : -1.0;
                an.Fmax_comp = am.Fmax_comp > 0.0 ? am.Fmax_comp / Ls : -1.0;
            }
        } else {
            an.EA = 1.0e5;
        }
        an.L = 0.0;  // EA/L uses the geometric distance
        const bool in1 = in_soil(s.x1, s.y1), in2 = in_soil(s.x2, s.y2);
        if (in1 && in2) {                          // node-to-node (strut / internal support)
            an.node_a = nearest_node(s.x1, s.y1);
            an.node_b = nearest_node(s.x2, s.y2);
            if (an.node_a < 0 || an.node_b < 0 || an.node_a == an.node_b) continue;
        } else if (in1 || in2) {                   // fixed-end (one end outside the soil)
            const double sx = in1 ? s.x1 : s.x2, sy = in1 ? s.y1 : s.y2;
            const double fx = in1 ? s.x2 : s.x1, fy = in1 ? s.y2 : s.y1;
            an.node_a = nearest_node(sx, sy);
            if (an.node_a < 0) continue;
            an.node_b = -1; an.fixed_point = {fx, fy};
        } else {
            continue;                              // both ends outside the soil — nothing to anchor
        }
        structures.anchors.push_back(an);
        diag_specs.push_back({1, s.name, structures.anchors.size() - 1, structures.anchors.size(), Ls});
    }

    // Embedded walls: build the plate + two interfaces on the split mesh and seed each interface's
    // initial normal stress sigma_n0 = K0*sigma'_v from the layered field (wished-in-place K0).
    for (const auto& w : walls) {
        int toe = -1;
        for (int n = 0; n < mesh.node_count; ++n)
            if (std::fabs(mesh.x[n] - w.toe_x) < 1e-6 && std::fabs(mesh.y[n] - w.toe_y) < 1e-6) { toe = n; break; }
        if (toe < 0) continue;
        const double k0 = (w.soil_mat >= 0 && w.soil_mat < (int)k0_by_mat.size()) ? k0_by_mat[w.soil_mat] : 0.5;
        const double fac = k0 * w.nx * w.nx + w.ny * w.ny;   // sigma_n / sigma'_v (vertical -> k0)
        if (w.order == 15) {
            katai::core::WallBuild5 wb = katai::core::build_embedded_wall5(mesh, w.seam, toe, dofs, w.pp, w.ip);
            const auto ncp = katai::core::iface::nc_points5();
            for (auto& ie : wb.interfaces)
                for (int q = 0; q < 5; ++q) {
                    const int nd = ie.soil_nodes[ncp[q].node];
                    ie.sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
            dofs.fix(wb.dof_phi.front());
            const size_t p0 = structures.plates5.size();
            for (const auto& pe : wb.plates) structures.plates5.push_back(pe);
            const size_t if0 = structures.interfaces5.size();
            for (const auto& ie : wb.interfaces) structures.interfaces5.push_back(ie);
            if (structures.plates5.size() > p0)
                diag_specs.push_back({5, w.name, p0, structures.plates5.size()});   // tri15 wall N/Q/M
            if (structures.interfaces5.size() > if0)
                iface_diags.push_back({w.name, 15, if0, structures.interfaces5.size()});
        } else {
            katai::core::WallBuild wb = katai::core::build_embedded_wall(mesh, w.seam, toe, dofs, w.pp, w.ip);
            const auto ncp = katai::core::iface::nc_points();
            for (auto& ie : wb.interfaces)
                for (int q = 0; q < 3; ++q) {
                    const int nd = ie.soil_nodes[ncp[q].node];
                    ie.sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
            dofs.fix(wb.dof_phi.front());   // remove the toe rotation mode (tied translation, free phi)
            const size_t p0 = structures.plates.size();
            for (const auto& pe : wb.plates) structures.plates.push_back(pe);
            const size_t if0 = structures.interfaces.size();
            for (const auto& ie : wb.interfaces) structures.interfaces.push_back(ie);
            if (structures.plates.size() > p0)
                diag_specs.push_back({0, w.name, p0, structures.plates.size()});
            if (structures.interfaces.size() > if0)
                iface_diags.push_back({w.name, 6, if0, structures.interfaces.size()});
        }
    }

    // Standalone interfaces: build the soil-soil Coulomb joint on the split mesh and seed each NC point's
    // initial normal stress sigma_n0 = (K0*nx^2 + ny^2)*sigma'_v (orientation-aware; vertical -> K0*sigma'_v,
    // horizontal -> sigma'_v) so the wished-in-place interface starts in geostatic equilibrium.
    for (auto& sp : soil_ifaces) {
        const katai::core::InterfaceRange rng =
            katai::core::build_soil_interface(mesh, sp.seam, dofs, sp.ip, mesh.nodes_per_element, structures);
        if (rng.end > rng.begin) iface_diags.push_back({sp.name, rng.order, rng.begin, rng.end});
        const double k0 = (sp.soil_mat >= 0 && sp.soil_mat < (int)k0_by_mat.size()) ? k0_by_mat[sp.soil_mat] : 0.5;
        const double fac = k0 * sp.nx * sp.nx + sp.ny * sp.ny;
        if (rng.order == 15) {
            const auto ncp = katai::core::iface::nc_points5();
            for (std::size_t i = rng.begin; i < rng.end; ++i)
                for (int q = 0; q < 5; ++q) {
                    const int nd = structures.interfaces5[i].soil_nodes[ncp[q].node];
                    structures.interfaces5[i].sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
        } else {
            const auto ncp = katai::core::iface::nc_points();
            for (std::size_t i = rng.begin; i < rng.end; ++i)
                for (int q = 0; q < 3; ++q) {
                    const int nd = structures.interfaces[i].soil_nodes[ncp[q].node];
                    structures.interfaces[i].sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
        }
    }
    const bool has_interfaces = !structures.interfaces.empty() || !structures.interfaces5.empty();

    // Embedded beams (pile rows / nails): a beam line cuts the mesh at any orientation; skin + foot
    // springs couple the beam (on independent DOFs) to the soil AT the beam location (non-conforming,
    // no mesh change). EA/EI and the skin/foot spring stiffnesses come from the embedded-beam material
    // via PLAXIS default interface-stiffness factors (Sluis 2012; PLAXIS 2D Ref §6.6): ISF = 2.5
    // (Lspacing/Deq)^-0.75 for the springs, 25 (...)^-0.75 for the foot, times the soil shear modulus.
    double area_total = 0.0;
    for (int e = 0; e < mesh.element_count; ++e) {
        const int a = mesh.node_of(e, 0), b = mesh.node_of(e, 1), c = mesh.node_of(e, 2);
        area_total += 0.5 * std::fabs((mesh.x[b]-mesh.x[a])*(mesh.y[c]-mesh.y[a]) -
                                      (mesh.x[c]-mesh.x[a])*(mesh.y[b]-mesh.y[a]));
    }
    const double h_char = mesh.element_count > 0 ? std::sqrt(2.0 * area_total / mesh.element_count) : 1.0;
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::EmbeddedBeam || !struct_on(si)) continue;
        if (s.material < 0 || s.material >= (int)pr.embedded.size()) continue;
        const auto& em = pr.embedded[s.material];
        const double L = std::hypot(s.x2 - s.x1, s.y2 - s.y1);
        if (L < 1e-9 || em.diameter <= 0.0) continue;
        const int ne = std::clamp((int)std::lround(L / std::max(1e-6, h_char)), 3, 80);
        const int nbn = 2 * ne + 1;
        std::vector<double> bx(nbn), by(nbn);
        for (int i = 0; i < nbn; ++i) {
            const double t = (double)i / (nbn - 1);
            bx[i] = s.x1 + t * (s.x2 - s.x1); by[i] = s.y1 + t * (s.y2 - s.y1);
        }
        if (by.front() > by.back()) {   // build_embedded_beam: node 0 = toe (foot) -> deeper (lower y) end
            std::reverse(bx.begin(), bx.end()); std::reverse(by.begin(), by.end());
        }
        const double Ls = em.Lspacing > 1e-9 ? em.Lspacing : 1.0;
        const double A = kPi * em.diameter * em.diameter / 4.0;     // circular massive pile
        const double Imom = kPi * std::pow(em.diameter, 4) / 64.0;
        katai::core::plate::PlateProps pp;
        pp.EA = em.E * A / Ls; pp.EI = em.E * Imom / Ls; pp.nu = 0.2;   // per-metre-width (2D row)
        // Pile inertia for a Dynamic phase, on the SAME per-metre-of-wall basis as EA/EI: the row's
        // weight per unit length is gamma * A / L_spacing [kN/m/m], i.e. the plate material's `w`.
        // Static paths never read it; a weightless pile would be a stiffener with no inertia.
        set_plate_mass(pp, em.gamma * A / Ls);
        const int smat = material_at(0.5 * (s.x1 + s.x2), 0.5 * (s.y1 + s.y2));
        const auto& sm = pr.materials[smat >= 0 ? smat : 0];
        const double G = sm.E / (2.0 * (1.0 + sm.nu));
        // Skin/foot spring stiffnesses from the PLAXIS default interface-stiffness factors --
        // engine-owned now (Stage B9, embedded_beam.hpp: default_interface_stiffness).
        double k_axial, k_lateral, D_foot;
        katai::core::ebeam::default_interface_stiffness(Ls, em.diameter, G, k_axial, k_lateral,
                                                        D_foot);
        const double t_max = em.Tskin_max > 0.0 ? em.Tskin_max / Ls : -1.0;
        const double f_max = em.Fmax_base > 0.0 ? em.Fmax_base / Ls : -1.0;
        structures.embedded_beams.push_back(katai::core::ebeam::build_embedded_beam(
            mesh, dofs, bx, by, pp, k_axial, k_lateral, t_max, D_foot, f_max));
        // Force-diagram bookkeeping (kind 3 = embedded beam; produces N/Q/M like a plate).
        diag_specs.push_back({3, s.name, structures.embedded_beams.size() - 1,
                              structures.embedded_beams.size()});
    }
    const bool has_embedded = !structures.embedded_beams.empty();

    // Compliant (absorbing) base: only a Dynamic phase that asked for it frees the base u_x --
    // every static phase and the rigid-base dynamic default keep the user's fixities bit-for-bit.
    const bool compliant_base =
        phase == InitialPhase::Dynamic && io.config && io.config->seismic_compliant_base;
    katai::core::apply_boundary_conditions(bc_edges_from(pr), mesh, dofs, compliant_base);
    // Line prescribed displacements (schema v2): every mesh node ON an active line gets the
    // set components FIXED here and their target values ramped 0 -> u by the static solver
    // (nonzero Dirichlet). Static family only: the validator refuses other phase types at
    // the input contract, and the guard below keeps a direct driver call equally honest.
    const auto disp_on = [&io](size_t i) { return !io.config || io.config->active_disp(i); };
    struct PrescEntry { int node, comp; double value; };
    std::vector<PrescEntry> presc_entries;
    {
        bool any_disp_active = false;
        for (size_t di = 0; di < pr.disps.size(); ++di)
            if (disp_on(di)) any_disp_active = true;
        if (any_disp_active && phase != InitialPhase::GravityLoading) {
            R.message = "Prescribed displacements are only supported in Plastic (staged "
                        "construction) phases and the gravity-loading initial phase in this "
                        "build; deactivate them in this phase.";
            return R;
        }
        for (size_t di = 0; di < pr.disps.size(); ++di) {
            if (!disp_on(di)) continue;
            const auto& D = pr.disps[di];
            const double ex = D.x2 - D.x1, ey = D.y2 - D.y1, l2 = ex * ex + ey * ey;
            if (l2 < 1e-18) continue;   // degenerate: validator reports it; nothing to fix
            for (int n = 0; n < mesh.node_count; ++n) {
                const double t = std::clamp(
                    ((mesh.x[n] - D.x1) * ex + (mesh.y[n] - D.y1) * ey) / l2, 0.0, 1.0);
                const double d = std::hypot(mesh.x[n] - (D.x1 + t * ex),
                                            mesh.y[n] - (D.y1 + t * ey));
                if (d > 1e-6) continue;
                if (D.set_ux) {
                    dofs.fix_node_component(n, 0);
                    presc_entries.push_back({n, 0, D.ux});   // overlapping lines: last wins
                }
                if (D.set_uy) {
                    dofs.fix_node_component(n, 1);
                    presc_entries.push_back({n, 1, D.uy});
                }
            }
        }
    }
    // Staged construction: nodes touched only by passive (excavated / not-yet-filled) elements
    // would be singular -- fix them (their displacement is meaningless this phase).
    if (!act.empty()) katai::core::fix_inactive_nodes(mesh, act, dofs);
    dofs.finalize();
    Eigen::VectorXd presc;   // full-DOF prescribed values; empty = none active
    if (!presc_entries.empty()) {
        presc = Eigen::VectorXd::Zero(dofs.total_dofs());
        for (const PrescEntry& pe : presc_entries)
            presc[dofs.global_dof(pe.node, pe.comp)] = pe.value;
    }
    const int fixed = dofs.total_dofs() - dofs.equation_count();
    if (fixed == 0) {
        R.message = "No boundary conditions. Right-click a soil edge (Soil tab) and set a boundary "
                    "condition before calculating.";
        return R;
    }
    if (dofs.equation_count() == 0) { R.message = "Every DOF is fixed; nothing to solve."; return R; }

    // ---- Parent structural-state carry (Track 1a; consumed by the NONLINEAR dynamic branch and
    // by CHAINED static phases). Rationale and the exact-size contract live with the builder
    // (Stage B9, katai/analysis/structural_carry.hpp); on any mismatch carry_src stays null and
    // the consumers fall back -- honestly, named in their messages -- to increment-from-zero.
    katai::core::StructuralInit carry_init;          // filled only when carry_src != nullptr
    const StructCarryState* carry_src = nullptr;     // parent full_disp (diagram evaluation)
    const bool any_struct_carry =
        !structures.plates.empty() || !structures.plates5.empty() || !structures.anchors.empty() ||
        !structures.geogrids.empty() || !structures.interfaces.empty() ||
        !structures.interfaces5.empty() || !structures.embedded_beams.empty();
    bool static_carry_used = false;      // set by the static tail when it actually consumed carry
    bool static_carry_missing = false;   // chained + structures, but no carriable parent state
    if (any_struct_carry && io.prev && io.prev->ok) {
        const StructCarryState& ps = io.prev->struct_state;
        if (katai::core::build_structural_init(ps, structures, dofs, carry_init)) carry_src = &ps;
    }

    // Axisymmetric (v1) scope: soil-only, dry, no structural elements. Guard the unsupported
    // combinations honestly rather than silently producing a wrong (plane-strain) result for them.
    if (axi) {
        const bool has_structs =
            !structures.plates.empty() || !structures.plates5.empty() || !structures.anchors.empty() ||
            !structures.geogrids.empty() || !structures.interfaces.empty() ||
            !structures.interfaces5.empty() || !structures.embedded_beams.empty();
        if (water) { R.message = "Axisymmetric analysis with a water table is not supported yet "
                                 "(use plane strain, or remove the water table)."; return R; }
        if (has_structs) { R.message = "Structural elements are not supported in axisymmetric mode yet "
                                       "(soil-only). Use plane strain for walls/anchors/plates."; return R; }
    }

    // Staged-phase scope guards (v1): keep the unsupported combinations honest.
    if (io.config && axi) { R.message = "Staged phases are not available in axisymmetric mode yet."; return R; }

    // Groundwater-flow coupling validity: the head field must match the solved mesh (embedded
    // walls SPLIT the mesh -- the flow result was computed on the unsplit one, so reject honestly).
    const bool flow = flow_head != nullptr;
    if (flow) {
        if (axi) { R.message = "Groundwater-flow coupling is not available in axisymmetric mode yet."; return R; }
        if (io.config) { R.message = "Groundwater-flow coupling inside staged phases is not supported yet "
                                     "(use the water level, or single-phase flow coupling)."; return R; }
        if (flow_head->size() != (Eigen::Index)mesh.node_count) {
            R.message = "The groundwater-flow result does not match this mesh (was it regenerated, "
                        "or does the model use an embedded wall?). Recalculate flow on the current "
                        "mesh (Flow conditions > Calculate groundwater flow).";
            return R;
        }
    }

    // NonPorous elements take no pore-pressure LOAD (PLAXIS: a non-porous material carries
    // neither initial nor excess pore pressure -- concrete holds no water of its own;
    // equilibrium in that region is stated in TOTAL stress). The mask overlays activation;
    // without NonPorous it stays empty (the legacy path).
    std::vector<char> act_pore;
    if (any_nonporous) {
        act_pore.assign(mesh.element_count, 1);
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!act.empty() && !act[e]) { act_pore[e] = 0; continue; }
            const int mt = mesh.element_material[e];
            if (mt >= 0 && mt < (int)mat_nonporous.size() && mat_nonporous[mt]) act_pore[e] = 0;
        }
    }
    const std::vector<char>& pore_mask = act_pore.empty() ? act : act_pore;

    // External load vector (free DOFs): gravity body force + point loads (+ pore pressures).
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    if (flow) {
        // Pore pressure + saturation from the steady-state seepage head field (PLAXIS groundwater
        // flow pore generation): total-stress gravity with gamma_sat where psi = h - y >= 0, and
        // the pore load interpolated from the SAME nodal head -> recovered stress is EFFECTIVE.
        katai::core::assemble_gravity_from_head(mesh, dofs, gamma, gamma_sat, *flow_head, f);
        katai::core::assemble_pore_load_from_head(mesh, dofs, *flow_head, kGammaWater, f, pore_mask);
    } else if (water) {
        // Total-stress equilibrium: saturated weight below the phreatic surface, moist above.
        const auto wt = [&pr](double x) { return water_table_at(pr, x); };
        katai::core::assemble_gravity_phreatic(mesh, dofs, gamma, gamma_sat, wt, f, act);
        // Hydrostatic pore-pressure load -> recovered stress is EFFECTIVE (Terzaghi; buoyancy
        // gamma' = gamma_sat - gamma_w emerges naturally). docs/references/effective-stress-formulation.md.
        const auto pore = [&pr](double x, double y) {
            return kGammaWater * std::fmax(0.0, water_table_at(pr, x) - y);
        };
        katai::core::assemble_pore_pressure_load(mesh, dofs, pore, f, pore_mask);
    } else if (axi) {
        katai::core::assemble_axisym_gravity(mesh, dofs, gamma, f);   // r-weighted body force
    } else {
        katai::core::assemble_gravity(mesh, dofs, gamma, f, act);
    }
    // External point loads (kept separate too: for the embedded-wall K0 baseline only loads ramp,
    // while self-weight is carried by the seeded geostatic state).
    Eigen::VectorXd f_loads = Eigen::VectorXd::Zero(dofs.equation_count());
    for (size_t li = 0; li < pr.loads.size(); ++li) {
        const auto& L = pr.loads[li];
        if (L.kind != model::LoadKind::Point || !load_on(li)) continue;
        int best = -1; double bestd = 1e300;
        for (int n = 0; n < mesh.node_count; ++n) {
            const double d = std::hypot(mesh.x[n] - L.x1, mesh.y[n] - L.y1);
            if (d < bestd) { bestd = d; best = n; }
        }
        if (best < 0) continue;
        const int ex = dofs.equation(dofs.global_dof(best, 0));
        const int ey = dofs.equation(dofs.global_dof(best, 1));
        if (ex >= 0) { f[ex] += L.qx1; f_loads[ex] += L.qx1; }
        if (ey >= 0) { f[ey] += L.qy1; f_loads[ey] += L.qy1; }
    }
    // Distributed (line) loads: a linearly-varying traction (qx1,qy1)->(qx2,qy2) along the
    // segment, assembled as CONSISTENT nodal forces over the mesh edge chain it lies on (the
    // common PLAXIS surcharge). The mesh conforms to the load line (build_mesh adds it as a
    // constraint), so collect_chain returns the edge nodes (corner, mid, corner, ...); the
    // varying-traction integral reproduces the q1->q2 ramp exactly. A line that does not resolve
    // to a valid edge chain (non-conforming) is skipped rather than mis-applied.
    {
        const int npe = mesh.nodes_per_element == 15 ? 5 : 3;
        for (size_t li = 0; li < pr.loads.size(); ++li) {
            const auto& L = pr.loads[li];
            if (L.kind != model::LoadKind::Distributed || !load_on(li)) continue;
            const std::vector<int> chain = collect_chain(mesh, L.x1, L.y1, L.x2, L.y2);
            const int cs = (int)chain.size();
            if (cs < npe || (cs - 1) % (npe - 1) != 0) continue;  // needs a conforming edge chain
            const double dx = L.x2 - L.x1, dy = L.y2 - L.y1, L2 = dx * dx + dy * dy;
            std::vector<double> tx(cs), ty(cs);
            for (int i = 0; i < cs; ++i) {
                const double t = L2 > 1e-18
                    ? ((mesh.x[chain[i]] - L.x1) * dx + (mesh.y[chain[i]] - L.y1) * dy) / L2 : 0.0;
                tx[i] = L.qx1 + t * (L.qx2 - L.qx1);
                ty[i] = L.qy1 + t * (L.qy2 - L.qy1);
            }
            if (axi) {
                katai::core::assemble_axisym_traction_varying(mesh, dofs, chain, tx, ty, f);
                katai::core::assemble_axisym_traction_varying(mesh, dofs, chain, tx, ty, f_loads);
            } else {
                katai::core::assemble_surface_traction_varying(mesh, dofs, chain, tx, ty, f);
                katai::core::assemble_surface_traction_varying(mesh, dofs, chain, tx, ty, f_loads);
            }
        }
    }

    R.design_approach = io.config ? to_core_design_approach(io.config->design_approach)
                                  : katai::core::DesignApproach::None;

    // --- EC7 / TBDY 2018 design-code partial factors (v0.3 B1) -------------------------------------
    // Material-factored approaches (EC7 DA1-C2, DA3) reduce c'/tan(phi') on every soil material and
    // scale the VARIABLE (applied) loads by gamma_Q; self-weight is permanent (gamma_G = 1.0 for these
    // approaches) so gravity + pore stay at characteristic values. The geostatic K0 seed below always
    // uses characteristic strength (the in-situ state is real, not factored) -- PLAXIS likewise applies
    // design factors in the calculation phase, not the initial one. None -> no-op (bit-identical). The
    // resistance-factored approaches (DA2, TBDY 2018) do not change the FEM solve; their E_d <= R_d
    // verdict is a report-layer check. Reference: docs/references/design-codes-ec7-tbdy.md.
    if (io.config && io.config->design_approach != model::DesignApproach::None) {
        const katai::core::DesignApproach da = to_core_design_approach(io.config->design_approach);
        if (katai::core::factors_material(da)) {
            const katai::core::PartialFactors pf = katai::core::design_factors(da);
            for (auto& mm : models) katai::core::factor_material_strength(mm, pf);
            if (pf.gamma_Q_unfav != 1.0) f += (pf.gamma_Q_unfav - 1.0) * f_loads;  // variable loads x gamma_Q
        }
    }

    // K0 initial-stress procedure (PLAXIS): the undisturbed ground starts in geostatic equilibrium
    // (sigma'_v = effective overburden, sigma'_h = K0 sigma'_v), so self-weight produces ~zero
    // displacement and the post-processed stress is the real geostatic field. Layered + water-aware
    // via vertical integration; reuses the eff_unit_weight / ground_surface lambdas built above.
    katai::core::K0LayeredOptions k0opt;
    k0opt.k0 = k0_by_mat;
    k0opt.eff_unit_weight = eff_unit_weight;
    k0opt.ground_surface = ground_surface;
    // NonPorous targets are seeded in total stress (header note): u(x,y) hydrostatic.
    if (any_nonporous && water) {
        k0opt.nonporous = mat_nonporous;
        k0opt.pore = [&pr](double x, double y) {
            return kGammaWater * std::fmax(0.0, water_table_at(pr, x) - y);
        };
    }
    double ymin = 1e30, ymax = -1e30;
    for (int n = 0; n < mesh.node_count; ++n) { ymin = std::fmin(ymin, mesh.y[n]); ymax = std::fmax(ymax, mesh.y[n]); }
    k0opt.integration_steps = std::clamp((int)((ymax - ymin) / 0.05), 50, 600);
    // gamma' discontinuity elevations along the column at x: polygon (layer) edge crossings + the
    // water table. With these the overburden integral is exact for piecewise-constant gamma', so the
    // level-ground K0 identity f_int(sigma_K0) = f_gravity holds to round-off even with layers/water.
    // PLAXIS K0 validity condition (Reference Manual, literally): the K0 procedure is correct ONLY
    // when the ground surface, the layer boundaries and the water table are ALL horizontal. Detect
    // violation GEOMETRICALLY (not from the assembled force imbalance: quadrature residuals -- the
    // r-weighted axisymmetric cubic integrand, a water-table kink crossing element interiors -- are
    // of the same order as a gentle slope's genuine imbalance and must NOT trigger a nil-step).
    bool k0_nonlevel = false;
    {
        double smin = 1e300, smax = -1e300;
        for (int n = 0; n < mesh.node_count; ++n) {
            const double s = ground_surface(mesh.x[n]);
            smin = std::fmin(smin, s); smax = std::fmax(smax, s);
        }
        if (smax - smin > 1e-6 * std::fmax(1.0, ymax - ymin)) k0_nonlevel = true;     // sloped surface
        if (water)
            for (size_t i = 1; i < pr.wy.size() && !k0_nonlevel; ++i)
                if (std::fabs(pr.wy[i] - pr.wy[0]) > 1e-9) k0_nonlevel = true;        // sloped water table
        // Non-horizontal boundary between regions whose geostatic properties (gamma, K0) differ:
        // sigma'_v / sigma'_h then jump across a non-horizontal plane -> traction discontinuity.
        for (size_t a = 0; a < pr.polygons.size() && !k0_nonlevel; ++a)
            for (size_t b = 0; b < pr.polygons.size() && !k0_nonlevel; ++b) {
                if (a == b) continue;
                const auto& A = pr.polygons[a]; const auto& Bp = pr.polygons[b];
                const int ma = A.material, mb = Bp.material;
                if (ma == mb) continue;
                if (ma < 0 || mb < 0 || ma >= (int)pr.materials.size() || mb >= (int)pr.materials.size()) continue;
                const bool same_geo =
                    std::fabs(pr.materials[ma].gamma_unsat - pr.materials[mb].gamma_unsat) < 1e-9 &&
                    (!water || std::fabs(pr.materials[ma].gamma_sat - pr.materials[mb].gamma_sat) < 1e-9) &&
                    std::fabs(k0_by_mat[ma] - k0_by_mat[mb]) < 1e-9;
                if (same_geo) continue;
                const int na = (int)A.x.size(), nb = (int)Bp.x.size();
                for (int i = 0; i < na && !k0_nonlevel; ++i) {
                    const double ax = A.x[i], ay = A.y[i], ax2 = A.x[(i + 1) % na], ay2 = A.y[(i + 1) % na];
                    if (std::fabs(ay2 - ay) < 1e-9) continue;   // horizontal edge: always fine
                    for (int j = 0; j < nb; ++j) {
                        const double bx = Bp.x[j], by = Bp.y[j], bx2 = Bp.x[(j + 1) % nb], by2 = Bp.y[(j + 1) % nb];
                        // collinear-overlap test (shared boundary segment between the two regions)
                        const double ex = ax2 - ax, ey = ay2 - ay;
                        const double c1 = ex * (by - ay) - ey * (bx - ax);
                        const double c2 = ex * (by2 - ay) - ey * (bx2 - ax);
                        const double L = std::hypot(ex, ey);
                        if (std::fabs(c1) > 1e-7 * L || std::fabs(c2) > 1e-7 * L) continue;
                        const double t1 = ((bx - ax) * ex + (by - ay) * ey) / (L * L);
                        const double t2 = ((bx2 - ax) * ex + (by2 - ay) * ey) / (L * L);
                        if (std::fmax(std::fmin(t1, t2), 0.0) < std::fmin(std::fmax(t1, t2), 1.0) - 1e-9)
                            k0_nonlevel = true;   // overlapping non-horizontal interface
                        if (k0_nonlevel) break;
                    }
                }
            }
    }
    k0opt.strata_breaks = [&pr, water](double x) {
        std::vector<double> br;
        if (water) br.push_back(water_table_at(pr, x));
        for (const auto& P : pr.polygons) {
            const int n = (int)P.x.size();
            for (int i = 0; i < n; ++i) {
                const double ax = P.x[i], ay = P.y[i], bx = P.x[(i + 1) % n], by = P.y[(i + 1) % n];
                const double lo = std::min(ax, bx), hi = std::max(ax, bx);
                if (x < lo - 1e-9 || x > hi + 1e-9) continue;
                const double dx = bx - ax;
                const double t = std::fabs(dx) < 1e-12 ? 0.0 : std::clamp((x - ax) / dx, 0.0, 1.0);
                br.push_back(ay + t * (by - ay));
            }
        }
        return br;
    };
    // K0 procedure: seed the geostatic stress so undisturbed self-weight gives ~zero displacement.
    // Gravity loading: start from zero stress (self-weight produces settlement). Either way the
    // recovered stress comes from the solver's committed Gauss states below.
    // Embedded walls always need the geostatic seed (interface sigma_n0 carries the lateral earth
    // pressure), so seed K0 whenever there are interfaces, regardless of the chosen phase.
    // Use the K0 procedure (geostatic seed) for the K0 phase and ALWAYS when interfaces are present
    // (interface sigma_n0 carries the lateral earth pressure). Gravity loading starts from zero stress.
    const bool use_k0 = !io.chained && ((phase == InitialPhase::K0Procedure) || has_interfaces);
    std::vector<katai::core::GaussState> init;
    if (io.chained) {
        // Staged phase: start from the previous phase's committed stresses (SumMstage chaining).
        if (!io.init_states || io.init_states->empty() ||
            io.init_states->size() % (size_t)std::max(1, mesh.element_count) != 0) {
            R.message = "Internal: staged phase started without the previous phase's stress state.";
            return R;
        }
        init = *io.init_states;
    } else if (use_k0) {
        init = katai::core::compute_k0_initial_stress_layered(mesh, k0opt);
        const int ng = mesh.element_count > 0 ? (int)(init.size() / mesh.element_count) : 0;
        // Staged initial phase: an INACTIVE region is soil that has not been placed yet -- its
        // seed must be ZERO stress (when a later phase activates it, the fill arrives stress-free
        // and its weight loads the ground; a K0 seed there would smuggle in phantom prestress).
        if (!act.empty())
            for (int e = 0; e < mesh.element_count && ng > 0; ++e)
                if (!act[e])
                    for (int g = 0; g < ng; ++g) init[(size_t)e * ng + g] = katai::core::GaussState{};
        // Cap/hardening preconsolidation seeding (HS pp + gamma_p, SS/SSC pp, OCR/POP
        // equivalence) is engine-owned now (Stage B9, initial_stress.hpp:
        // seed_preconsolidation); this seam only resolves the schema's per-material
        // overconsolidation fields.
        std::vector<katai::core::Overconsolidation> oc(pr.materials.size());
        for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
            oc[mi].mode = pr.materials[mi].oc_mode;
            oc[mi].OCR = pr.materials[mi].OCR;
            oc[mi].POP = pr.materials[mi].POP;
        }
        katai::core::seed_preconsolidation(mesh, models, oc, act, init);
    }

    try {
        // Solver: nonsymmetric tangent for nonlinear soil (non-associated MC/HS), interfaces (Coulomb)
        // or embedded beams (skin/foot plasticity); otherwise SPD. Load steps grow with nonlinearity so
        // plasticity develops gradually (the solver also has adaptive cut-back).
        const bool nonsym = nonlinear_soil || has_interfaces || has_embedded;
        // Hardening Soil uses an analytic CONTINUUM tangent (linear, not quadratic, global
        // convergence), so it needs more, smaller load increments and a PLAXIS-realistic
        // tolerated error (~1%); ramping gravity into the K0 seed keeps each increment small.
        // Mohr-Coulomb has a closed-form CONSISTENT tangent (quadratic) -> tighter 1e-6.
        const int steps = (has_hardening || has_softsoil)
                              ? 40
                              : (nonlinear_soil ? 20 : ((has_interfaces || has_embedded) ? 5 : 1));
        const double tol = (has_hardening || has_softsoil)
                               ? 1e-2
                               : (nonlinear_soil ? 1e-6
                                                 : ((has_interfaces || has_embedded) ? 1e-8 : 1e-10));
        const katai::core::LinearSolve solver = reusing_linear_solve(
            nonsym ? katai::linsolve::MatrixType::RealNonsymmetric
                   : katai::linsolve::MatrixType::RealSymmetricPositiveDefinite);

        // SAFETY analysis (phi-c reduction / SRM): the strategy lives in the engine (Stage B9:
        // katai/analysis/phase_solver/safety.hpp). The seam passes the registry-derived model
        // family flags and the composition root's solver callback; the refusals and the honest
        // lower-bound reporting are engine-owned. On success the phase falls through to the
        // common result tail, exactly as before. Note the Safety branch consumes f BEFORE the
        // structural self-weight below and stays structure-free (safety_analysis takes no
        // structures).
        if (phase == InitialPhase::Safety) {
            katai::core::SafetyPhase sfin;
            sfin.nonlinear_soil = nonlinear_soil;
            sfin.has_hardening = has_hardening;
            sfin.has_softsoil = has_softsoil;
            sfin.axisymmetric = axi;
            sfin.active = act;
            if (!katai::core::solve_safety_phase(mesh, dofs, models, profiles, f, solver, sfin, R))
                return R;
        } else {
        // Structural SELF-WEIGHT (plate w [kN/m/m]; embedded pile gamma*A/L_s) -- the 2026-07 audit
        // fix: the weight fed only the dynamic mass before, so statics silently carried weightless
        // walls/piles. Consistent nodal line load over every ACTIVE chain of THIS phase, added to
        // BOTH f and f_loads: a gravity-start phase ramps it with the body force (ramp = f); a
        // K0/baseline initial phase ramps it as an unbalanced new load (ramp = f_loads -- the
        // wished-in-place wall settles under its own weight, PLAXIS's plastic nil-step equivalent);
        // a chained phase sees it inside ramp = f - B, so a plate already equilibrated by the parent
        // adds NOTHING (nil identity preserved) and a newly activated plate arrives incrementally
        // (SumMstage). Permanent load -> assembled AFTER the gamma_Q variable-load scaling; the
        // Safety branch consumed f above and stays structure-free (consistent: safety_analysis takes
        // no structures). Dynamic reads neither f nor f_loads (static equilibrium comes from the
        // parent baseline). w = 0 (default) contributes nothing -> bit-identical.
        {
            Eigen::VectorXd f_sw = Eigen::VectorXd::Zero(dofs.equation_count());
            katai::core::assemble_structural_weight(mesh, dofs, structures, kGravity, f_sw);
            f += f_sw;
            f_loads += f_sw;
        }
        // K0 procedure: the geostatic state (seeded soil stress + interface sigma_n0) IS the
        // equilibrium; its internal force is held as constant_force B and ONLY external loads ramp, so
        // residual(0) = B - f_int(0) = 0 on any mesh and self-weight is NOT double-counted during the
        // ramp (essential once the soil is plastic -- ramping gravity would yield spuriously). Gravity
        // loading: B = 0, ramp the full body force + loads from a stress-free start.
        // Baseline internal force: the K0 seed (initial phase) or the previous phase's committed
        // stresses (chained phase). Either way residual(0) = 0 and the configuration imbalance
        // (excavation unloading, fill weight, new loads) is what gets ramped.
        const bool baseline = use_k0 || io.chained;
        // Chained static phase with a carriable parent: the structures continue from the parent's
        // displacement datum + committed plastic state (Track 1a generalized to the static chain).
        // The baseline B must then hold the FULL structural internal force at that state -- so
        // residual(0) = 0 by the PARENT'S OWN equilibrium and the ramp f - B carries only the real
        // configuration change (a nil phase becomes a true no-op; before this, the ramp re-applied
        // the parent's structural tractions and the wall moment drifted 32% in a phase where
        // nothing changed).
        const bool static_carry = io.chained && carry_src != nullptr;
        Eigen::VectorXd B = Eigen::VectorXd::Zero(dofs.equation_count());
        if (baseline) {
            if (axi) katai::core::assemble_axisym_internal_force(mesh, dofs, init, B);  // r-weighted baseline
            else katai::core::assemble_internal_force(mesh, dofs, init, B, act);
            if (static_carry) {
                // Full structural baseline at (parent datum + committed plastic state). This
                // SUBSUMES the manual sigma_n0 terms below (coulomb_return adds sigma_n0
                // internally), so those loops must be SKIPPED on this path -- adding both would
                // double-count the wished-in-place lateral earth pressure.
                B += katai::core::structural_internal_force(
                    mesh, dofs, models, structures, carry_init,
                    axi ? katai::core::Kinematics::Axisymmetric
                        : katai::core::Kinematics::PlaneStrain);
            }
            // Wished-in-place interface sigma_n0 baseline is engine-owned now (Stage B9,
            // interface_baseline.hpp). Carry path skips it: sigma_n0 is already inside the
            // structural baseline above (coulomb_return adds it internally) -- adding both
            // would double-count the lateral earth pressure.
            if (!static_carry)
                katai::core::add_interface_sigma_n0_baseline(structures, mesh, dofs, B);
        }

        if (phase == InitialPhase::Consolidation) {
            // --- Time-dependent (Biot) consolidation phase (PLAXIS "Consolidation") --------------
            // The strategy lives in the engine (Stage B9: katai/analysis/phase_solver/
            // consolidation.hpp). This seam resolves the schema materials, reuses the B4 flow-edge
            // vocabulary, and builds the solver factories -- a composition-root decision. On
            // success the phase falls through to the common result tail, exactly as before.
            katai::core::ConsolidationPhase cin;
            cin.materials.resize(pr.materials.size());
            for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
                const auto& Mt = pr.materials[mi];
                katai::core::ConsolidationPhaseMaterial& cm = cin.materials[mi];
                cm.name = Mt.name;
                cm.kx = Mt.kx; cm.ky = Mt.ky;
                cm.porosity = Mt.e_init / (1.0 + Mt.e_init);
                cm.nonporous = mat_nonporous[mi] != 0;
            }
            cin.flow_edges = flow_edges_from(pr);
            cin.have_flow_bcs = any_flow_bc_declared(pr);
            cin.active = act;
            cin.has_structural_elements = has_interfaces || has_embedded ||
                !structures.plates.empty() || !structures.anchors.empty() ||
                !structures.geogrids.empty();
            if (io.config) { cin.duration_day = io.config->duration; cin.time_steps = io.config->time_steps; }
            cin.yscale = ymax - ymin;
            const Eigen::VectorXd dF = f - B;   // configuration imbalance (SumMstage)
            // Sparse Biot solve factories: the coupled saddle-point system A = [K L; Lᵀ -(ΔtH+S)] is
            // symmetric indefinite (nonsymmetric with a non-associated elastoplastic tangent); the
            // solver factors ONCE (constant dt) and back-substitutes every time step. Built here
            // because a test or app executable is the composition root that owns backend choice.
            const katai::core::ConsolidationSolveFactory lef = [](const katai::math::CsrMatrix& Asys) {
                return factorize_once(katai::linsolve::MatrixType::RealSymmetricIndefinite, Asys);
            };
            const katai::core::ConsolidationSolveFactory plf = [](const katai::math::CsrMatrix& Asys) {
                return factorize_once(katai::linsolve::MatrixType::RealNonsymmetric, Asys);
            };
            std::vector<katai::core::GaussState> committed;
            if (!katai::core::solve_consolidation_phase(mesh, dofs, models, profiles, init, dF,
                                                        nonlinear_soil, cin, lef, plf, R, committed))
                return R;
            if (io.out_states) *io.out_states = committed;
        } else if (phase == InitialPhase::FullyCoupled) {
            // --- Fully-coupled flow-deformation (PLAXIS "Fully coupled flow-deformation") ----------
            // The strategy lives in the engine (Stage B9: katai/analysis/phase_solver/
            // fully_coupled.hpp). This seam resolves the schema materials (flow description +
            // van Genuchten retention), reuses the B4 flow-edge vocabulary, and builds the
            // solver factories -- a composition-root decision. On success the phase falls
            // through to the common result tail, exactly as before.
            katai::core::FullyCoupledPhase fin;
            fin.materials.resize(pr.materials.size());
            for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
                const auto& Mt = pr.materials[mi];
                katai::core::FullyCoupledPhaseMaterial& fm = fin.materials[mi];
                fm.name = Mt.name;
                fm.kx = Mt.kx; fm.ky = Mt.ky;
                fm.retention = {Mt.gw_ga, Mt.gw_gn, Mt.gw_gl, Mt.gw_Sres, 1.0};
                fm.porosity = Mt.e_init / (1.0 + Mt.e_init);
                fm.nonporous = mat_nonporous[mi] != 0;
            }
            fin.flow_edges = flow_edges_from(pr);
            fin.have_flow_bcs = any_flow_bc_declared(pr);
            fin.active = act;
            fin.has_structural_elements = has_interfaces || has_embedded ||
                !structures.plates.empty() || !structures.anchors.empty() ||
                !structures.geogrids.empty();
            if (io.config) { fin.duration_day = io.config->duration; fin.time_steps = io.config->time_steps; }
            fin.yscale = ymax - ymin;
            const Eigen::VectorXd dF = f - B;   // configuration imbalance (SumMstage)
            // Sparse coupled solve factories: the saddle-point system is symmetric indefinite
            // (nonsymmetric with a non-associated elastoplastic tangent); the solver factors ONCE
            // (constant dt) and back-substitutes every time step. Built here because a test or app
            // executable is the composition root that owns backend choice.
            const katai::core::ConsolidationSolveFactory lef = [](const katai::math::CsrMatrix& Asys) {
                return factorize_once(katai::linsolve::MatrixType::RealSymmetricIndefinite, Asys);
            };
            const katai::core::ConsolidationSolveFactory plf = [](const katai::math::CsrMatrix& Asys) {
                return factorize_once(katai::linsolve::MatrixType::RealNonsymmetric, Asys);
            };
            std::vector<katai::core::GaussState> committed;
            if (!katai::core::solve_fully_coupled_phase(mesh, dofs, models, profiles, init, dF,
                                                        nonlinear_soil, fin, lef, plf, R, committed))
                return R;
            if (io.out_states) *io.out_states = committed;
        } else if (phase == InitialPhase::TransientFlow) {
            // --- Transient groundwater flow only (PLAXIS "Groundwater flow, transient") ------------
            // The strategy lives in the engine (Stage B9 pilot: katai/analysis/phase_solver/
            // transient_flow.hpp). This seam resolves the schema materials to their flow
            // description and precomputes the initial-head fallback (water table / ground
            // surface); the engine owns the physics, the refusals and the result filling.
            katai::core::TransientFlowPhase tin;
            tin.materials.resize(pr.materials.size());
            for (size_t mi = 0; mi < pr.materials.size(); ++mi) {
                const auto& Mt = pr.materials[mi];
                katai::core::TransientFlowMaterial& fm = tin.materials[mi];
                fm.name = Mt.name;
                fm.kx = Mt.kx; fm.ky = Mt.ky;
                fm.retention = {Mt.gw_ga, Mt.gw_gn, Mt.gw_gl, Mt.gw_Sres, 1.0};
                fm.porosity = Mt.e_init / (1.0 + Mt.e_init);
                fm.nonporous = mat_nonporous[mi] != 0;
            }
            tin.flow_edges = flow_edges_from(pr);
            tin.active = act;
            tin.fallback_head.resize(mesh.node_count);
            for (int n = 0; n < mesh.node_count; ++n)
                tin.fallback_head[n] = pr.has_water ? water_table_at(pr, mesh.x[n]) : mesh.y[n];
            if (io.config) { tin.duration_day = io.config->duration; tin.time_steps = io.config->time_steps; }
            if (!katai::core::solve_transient_flow_phase(mesh, tin, R)) return R;
            R.mesh = std::move(mesh);   // flow-only: skip the deformation post-processing tail
            return R;
        } else if (phase == InitialPhase::Dynamic) {
            // --- Dynamic (seismic) time-history analysis (PLAXIS "Dynamic") -----------------------
            // The strategy lives in the engine (Stage B9: katai/analysis/phase_solver/dynamic.hpp).
            // This seam resolves the schema once -- the seismic configuration and code-spectrum
            // overlays from io.config, the per-material unit weights (drainage-resolved by the
            // common setup) and the water table as a neutral callback -- and injects the solver
            // factories, a composition-root decision. Model-family stiffness selection happens in
            // the ENGINE from its models table (MaterialType), never from the schema enum. The
            // strategy owns the physics, the refusals and the result assembly, and fills
            // R.ok/message itself (TransientFlow-class return discipline); the mesh is moved in
            // here, on success only, exactly as before.
            katai::core::DynamicPhase din;
            din.axisymmetric = axi;
            din.active = act;
            din.compliant_base = compliant_base;
            din.has_water = pr.has_water;
            din.phreatic_mass = water;
            din.water_table_y = [&pr](double x) { return water_table_at(pr, x); };
            din.gamma = gamma;          // gamma_sat is drainage-resolved by the common setup
            din.gamma_sat = gamma_sat;  // (NonPorous holds no water: its saturated weight IS gamma_unsat)
            if (io.config) {
                din.nonlinear = io.config->dynamic_nonlinear;
                din.free_field = io.config->seismic_free_field;
                din.duration = io.config->duration;
                din.time_steps = io.config->time_steps;
                din.rayleigh_f1 = io.config->rayleigh_f1;
                din.rayleigh_f2 = io.config->rayleigh_f2;
                din.damping_ratio = io.config->damping_ratio;
                din.amplitude = io.config->seismic_amp;
                din.frequency = io.config->seismic_freq;
                din.wave = to_core_seismic_wave(io.config->seismic_wave);
                din.accel_record = io.config->accel_record;
                din.record_dt = io.config->record_dt;
                din.tbdy_ss = io.config->tbdy_ss;
                din.tbdy_s1 = io.config->tbdy_s1;
                din.site_class = io.config->site_class;
                din.ec8_enabled = io.config->ec8_enabled;
                din.ec8_gamma = io.config->ec8_gamma;
                din.ec8_agr = io.config->ec8_agr;
                din.ec8_ground = io.config->ec8_ground;
                din.ec8_type = io.config->ec8_type;
            }
            // Solver injection (composition root): K_eff and the f1 estimator's K are SPD and
            // factor once / back-solve many; the nonlinear tangent is nonsymmetric and its solver
            // is built only if the strategy takes the per-step-Newton path.
            const katai::core::DynamicsSolveFactory spdf = [](const katai::math::CsrMatrix& A) {
                return factorize_once(katai::linsolve::MatrixType::RealSymmetricPositiveDefinite, A);
            };
            const auto nonsymf = []() {
                return reusing_linear_solve(katai::linsolve::MatrixType::RealNonsymmetric);
            };
            if (katai::core::solve_dynamic_phase(mesh, dofs, models, mats, profiles, structures,
                                                 diag_specs, iface_diags, carry_init,
                                                 carry_src ? &carry_src->full_disp : nullptr,
                                                 io.prev, io.init_states, din, spdf, nonsymf, R))
                R.mesh = std::move(mesh);
            return R;
        } else {
        // The strategy lives in the engine (Stage B9: katai/analysis/phase_solver/
        // static_phase.hpp) together with the ramp / nil-step semantics (K0 baseline, SumMstage
        // chaining, SSC time apportioning). This seam passes the common setup's neutral products
        // (loads, baseline, activity, Newton class, carry) and the composition root's solver
        // callback; on success the phase falls through to the common result tail, as before.
        katai::core::StaticPhase stin;
        stin.baseline = baseline;
        stin.nil_step = (use_k0 && (k0_nonlevel || flow)) || io.chained;
        stin.axisymmetric = axi;
        stin.load_steps = steps;
        stin.tolerance = tol;
        // Only a CHAINED (staged) phase carries time: in PLAXIS the initial phase is TIMELESS --
        // an unconditional duration here once leaked 1 day of creep into the K0 phase and broke
        // the geostatic identity (measured: K0 max |u| = 26 mm = exactly mu* ln2 H on an SSC
        // column). The initial phase also has io.config (pr.initial, for activation flags), which
        // is why the gate is io.chained and not the config pointer.
        stin.time_interval_day = (io.chained && io.config) ? io.config->duration : 0.0;
        stin.active = act;
        stin.presc = presc;   // active prescribed displacements (empty = none)
        if (!katai::core::solve_static_phase(mesh, dofs, models, profiles, init, f, f_loads, B,
                                             solver, structures, diag_specs, iface_diags,
                                             carry_init,
                                             static_carry ? &carry_src->full_disp : nullptr,
                                             stin, R, io.out_states))
            return R;
        static_carry_used = static_carry;
        static_carry_missing = io.chained && any_struct_carry && !static_carry;
        }  // end ramp/solve (non-consolidation)
        }  // end normal (non-Safety) solve
    } catch (const std::exception& e) {
        R.message = std::string("Calculation failed -- the model is probably insufficiently restrained "
                                "(add boundary conditions so it cannot move or rotate freely). [") + e.what() + "]";
        return R;
    }
    for (int n = 0; n < mesh.node_count; ++n)
        R.max_disp = std::fmax(R.max_disp, std::hypot(R.disp[n * 2], R.disp[n * 2 + 1]));
    // Nodal pore pressure (>= 0), for post-processing display and total-stress recovery: from the
    // seepage head field when flow-coupled, else hydrostatic below the water polyline.
    R.pore.assign(mesh.node_count, 0.0);
    if (flow)
        for (int n = 0; n < mesh.node_count; ++n)
            R.pore[n] = kGammaWater * std::fmax(0.0, (*flow_head)[n] - mesh.y[n]);
    else if (water)
        for (int n = 0; n < mesh.node_count; ++n)
            R.pore[n] = kGammaWater * std::fmax(0.0, water_table_at(pr, mesh.x[n]) - mesh.y[n]);
    R.active = act;             // phase element activity (empty = all active)
    R.mesh = std::move(mesh);   // the (possibly split) mesh the GUI must render
    R.ok = true;
    R.message = "Solved: max |u| = " + std::to_string(R.max_disp);
    if (static_carry_used)
        R.message += " Structural state continued from the parent phase (forces are totals).";
    else if (static_carry_missing)
        R.message += " NOTE: the parent phase supplies no structural state (different structure set, "
                     "or a Safety/consolidation/restored parent), so structural forces re-develop "
                     "from zero this phase.";
    if (!R.consol_time.empty()) {
        const double s_inf = R.consol_settlement.back();
        const double s0 = R.consol_settlement.front();
        const char* label = (phase == InitialPhase::FullyCoupled)
                                ? "Fully-coupled flow-deformation" : "Consolidation";
        R.message = std::string(label) + " solved: settlement " + std::to_string(s0) + " -> " +
                    std::to_string(s_inf) + " m over " + std::to_string(R.consol_time.back()) +
                    " days (excess pore -> " + std::to_string(R.consol_excess_pore.back()) +
                    " kPa). See the settlement-time curve below.";
    }
    if (R.nil_step && !io.chained)
        R.message += "  [K0 on non-level ground (sloped surface / layers / water table): the K0 field "
                     "alone is not in equilibrium there, so an equilibrium (nil) step was applied ---";
    return R;
}

std::vector<SolveResult> solve_phases(const model::Project& pr,
                                      const katai::mesh::Mesh& mesh_in, InitialPhase init_phase,
                                      const PhaseProgress& on_phase,
                                      const std::function<bool()>& cancelled) {
    std::vector<SolveResult> out;
    std::vector<katai::core::GaussState> committed;
    const int total = 1 + (int)pr.phases.size();
    const auto stop = [&cancelled] { return cancelled && cancelled(); };

    PhaseIO io0;
    io0.config = &pr.initial;
    io0.out_states = &committed;
    if (stop()) return out;
    if (on_phase) on_phase(0, total, "Initial phase");
    if (stop()) return out;
    out.push_back(solve_gravity_le(pr, mesh_in, init_phase, nullptr, io0));
    if (!out.back().ok) return out;

    for (size_t pi = 0; pi < pr.phases.size(); ++pi) {
        const auto& ph = pr.phases[pi];
        if (stop()) return out;
        if (on_phase) on_phase((int)pi + 1, total, ph.name);
        if (stop()) return out;
        PhaseIO io;
        io.config = &ph;
        io.chained = true;
        io.init_states = &committed;
        // The parent phase's result: a Dynamic phase superposes its increment onto this static state
        // to report the TOTAL design action. A parent that is itself Dynamic carries an envelope, not
        // a state, so it is not a valid base -- the Dynamic branch checks that and falls back to
        // reporting its own action alone.
        io.prev = out.empty() ? nullptr : &out.back();
        std::vector<katai::core::GaussState> next;
        io.out_states = &next;
        const auto kind = ph.type == model::PhaseType::Safety        ? InitialPhase::Safety
                        : ph.type == model::PhaseType::Consolidation ? InitialPhase::Consolidation
                        : ph.type == model::PhaseType::TransientFlow ? InitialPhase::TransientFlow
                        : ph.type == model::PhaseType::FullyCoupled  ? InitialPhase::FullyCoupled
                        : ph.type == model::PhaseType::Dynamic       ? InitialPhase::Dynamic
                                                                     : InitialPhase::GravityLoading;
        out.push_back(solve_gravity_le(pr, mesh_in, kind, nullptr, io));
        if (!out.back().ok) return out;
        // Safety, TransientFlow and Dynamic leave the committed effective-stress state unchanged (a
        // time-history returns RELATIVE displacements, not a new static state); Plastic, Consolidation
        // and FullyCoupled commit it forward (SumMstage chaining).
        if (ph.type != model::PhaseType::Safety && ph.type != model::PhaseType::TransientFlow &&
            ph.type != model::PhaseType::Dynamic)
            committed = std::move(next);
    }
    return out;
}

}  // namespace katai::app
