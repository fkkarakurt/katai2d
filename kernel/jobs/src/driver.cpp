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
#include <katai/fem/elements/point_location.hpp>      // locate_point (is a point load on the mesh?)
#include <katai/materials/linear_elastic.hpp>
#include <katai/mesh/boundary_extraction.hpp>      // collect_chain + boundary edges (Stage B5)
#include <katai/linsolve/direct_solver.hpp>

#include <katai/jobs/mesh_builder.hpp>   // point_in_polygon (material_at region lookup)

namespace katai::app {

// Boundary extraction lives in the mesh module (Stage B5); aliases keep the body's spellings.
using katai::mesh::collect_chain;
using katai::mesh::BoundaryEdgeChain;
using katai::mesh::extract_boundary_edges;

// Consolidation stress recovery lives in the engine (Stage B8).
using katai::core::recover_consolidation_stress;

// ---- Diagnostics (katai/analysis/results.hpp, Diagnostic). One rule, stated once: an input
// ---- may be used differently from the way it was written, but never in silence.
//
// The distinction the three helpers encode is a judgement about consequence, not about how
// unusual the input is. `warn` is for a run that still answers the user's question -- a line
// clipped to the soil it could reach, a load attached one node away. `refuse` is for a drawn
// object that would contribute NOTHING: a surcharge that lands off the mesh, a wall the mesh
// never sees. Those cannot be warnings, because the analysis would then succeed, report a
// plausible field, and answer a model the engineer did not draw. Every code below is stable:
// tests match on it, users grep logs for it, and it is never reworded or reused.
static void note(SolveResult& R, const char* code, std::string subject, std::string message) {
    katai::core::add_diagnostic(R, katai::core::DiagnosticSeverity::Note, code, std::move(subject),
                                std::move(message));
}
static void warn(SolveResult& R, const char* code, std::string subject, std::string message) {
    katai::core::add_diagnostic(R, katai::core::DiagnosticSeverity::Warning, code,
                                std::move(subject), std::move(message));
}
// Sets the refusal AND records it as a tagged diagnostic, so a front end can act on the code
// instead of matching prose. The caller returns R immediately (R.ok is already false).
static void refuse(SolveResult& R, const char* code, std::string subject, std::string message) {
    katai::core::add_diagnostic(R, katai::core::DiagnosticSeverity::Refusal, code,
                                std::move(subject), std::move(message));
}
// Number for a message, in the file's own spelling.
static std::string dnum(double v) {
    char b[40];
    std::snprintf(b, sizeof(b), "%g", v);
    return b;
}
// A drawn line names itself by its user-given name when it has one, else by its endpoints --
// an unnamed object must still be findable in the file.
static std::string line_subject(const std::string& name, double x1, double y1, double x2, double y2) {
    if (!name.empty()) return name;
    return "(" + dnum(x1) + ", " + dnum(y1) + ") -> (" + dnum(x2) + ", " + dnum(y2) + ")";
}
// How much of a drawn line the mesh actually handed back, as the parameter span [t0, t1] of the
// collected chain along that line (0 = first endpoint, 1 = second). A line drawn past the edge
// of the soil comes back short, and the caller says by how much instead of applying it quietly.
static void chain_span(const katai::mesh::Mesh& mesh, const std::vector<int>& chain,
                       double x1, double y1, double x2, double y2, double& t0, double& t1) {
    t0 = 1.0; t1 = 0.0;
    const double dx = x2 - x1, dy = y2 - y1, L2 = dx * dx + dy * dy;
    if (L2 < 1e-18 || chain.empty()) { t0 = 0.0; t1 = 0.0; return; }
    for (int n : chain) {
        const double t = ((mesh.x[n] - x1) * dx + (mesh.y[n] - y1) * dy) / L2;
        t0 = std::fmin(t0, t);
        t1 = std::fmax(t1, t);
    }
}
// True when the chain covers materially less of the drawn line than the whole of it. The
// threshold is one part in a thousand of the line: below that it is mesh round-off, not clipping.
static bool chain_is_clipped(double t0, double t1) { return t0 > 1e-3 || t1 < 1.0 - 1e-3; }

// The mesh scale where an object sits: the longest corner edge of the element containing (x, y),
// or 0 when nothing contains it. An object that attaches to "the nearest node" has to be judged
// against THIS length rather than the project's target element size, which a region coarseness
// factor may multiply by up to four -- a snap of a fraction of the local element is ordinary
// discretisation, a snap of several is a different model.
static double element_size_at(const katai::mesh::Mesh& mesh, double x, double y) {
    const auto loc = katai::core::ploc::locate_point(mesh, x, y);
    if (!loc.found) return 0.0;
    double h = 0.0;
    for (int k = 0; k < 3; ++k) {
        const int a = mesh.node_of(loc.element, k), b = mesh.node_of(loc.element, (k + 1) % 3);
        h = std::fmax(h, std::hypot(mesh.x[a] - mesh.x[b], mesh.y[a] - mesh.y[b]));
    }
    return h;
}

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
// `config` is the phase being solved: a phase may carry its own phreatic polyline (PLAXIS
// "water conditions per phase"), which is how staged dewatering is expressed. Null config or
// no override falls back to the project's own line, so every existing run is unchanged.
static double water_table_at(const model::Project& pr, double x, const model::Phase* config) {
    if (!pr.has_water) return -1e30;
    if (config && config->water_override && !config->wx.empty())
        return katai::core::phreatic_surface_at(config->wx, config->wy, x);
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
        // silent-drop-ok: an EMPTY edge_bc means the region states no boundary conditions, which
        // is legal input (validate.hpp check_edges returns early on an empty array); any other
        // mismatch, and a region with fewer than three vertices, is an ERROR at its field path,
        // so a project that reaches the driver cannot be losing conditions here.
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
        // silent-drop-ok: as in bc_edges_from -- an empty edge_flow is "no flow conditions
        // stated", and any other mismatch is a validator error at polygons[i].edge_flow.
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
// The mesh nodes of the DRAINS active in this phase, and whether any WELL is active. A drain is
// read by the time-dependent flow solvers (its nodes drain: excess pore pressure zero, PLAXIS Ref
// sec. 5.9.2); a well prescribes a discharge, which those solvers do not take, so an active well
// in such a phase is reported rather than silently ignored.
// Said wherever a pore-pressure field is COMPUTED with walls in the model: an impermeable screen
// is not something this build can express (PLAXIS Sci. Man. sec. 3.4 gives interfaces their own
// flow setting), so the water crosses the wall. It is a warning and not a refusal because the
// deformation answer is still the model's; it is the flow half that is more permeable than drawn.
static void warn_no_flow_barrier(const model::Project& pr, SolveResult& R) {
    for (const auto& st : pr.structs)
        if (st.kind == model::StructKind::Plate || st.kind == model::StructKind::Interface) {
            warn(R, "K2D-A010", st.name,
                 "This phase computes a pore-pressure field with walls or interfaces in the "
                 "model, and a consolidation / fully-coupled phase does not read their cross "
                 "permeability in this build: water crosses the line as if the soil were "
                 "continuous, so a cut-off wall holds back less head here than it would in the "
                 "ground. The steady groundwater-flow calculation does read it.");
            return;
        }
}

static std::vector<int> phase_drain_nodes(const model::Project& pr, const katai::mesh::Mesh& mesh,
                                          const PhaseIO& io, bool& any_well_active) {
    std::vector<int> nodes;
    any_well_active = false;
    for (std::size_t hi = 0; hi < pr.hydros.size(); ++hi) {
        if (io.config && !io.config->active_hydro(hi)) continue;
        const model::HydroLine& H = pr.hydros[hi];
        if (H.kind == model::HydroKind::Well) { any_well_active = true; continue; }
        const auto chain = katai::mesh::collect_chain(mesh, H.x1, H.y1, H.x2, H.y2);
        nodes.insert(nodes.end(), chain.begin(), chain.end());
    }
    return nodes;
}

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
    p.dilatancy_cutoff = m.dilatancy_cutoff;
    p.e_init = m.e_init;
    p.e_max = m.e_max;
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
    p.skempton_mode = (m.und_mode == 1);
    p.nu_u = m.nu_u;
    p.skempton_B = m.skempton_B;
    switch (m.drainage) {
        case model::Drainage::Drained:    p.drainage = katai::core::DrainageClass::Drained; break;
        case model::Drainage::Undrained:  p.drainage = katai::core::DrainageClass::UndrainedA; break;
        case model::Drainage::UndrainedB: p.drainage = katai::core::DrainageClass::UndrainedB; break;
        case model::Drainage::NonPorous:  p.drainage = katai::core::DrainageClass::NonPorous; break;
        case model::Drainage::UndrainedC: p.drainage = katai::core::DrainageClass::UndrainedC; break;
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
        // "Ignore undrained behaviour" (PLAXIS): this phase treats Undrained (A)/(B) soil as
        // drained -- no excess pore pressure is generated. Only the volumetric coupling is
        // switched off; the strength parameters stay exactly as the material declares them, so
        // an Undrained (B) soil keeps its c_u with phi = 0, which is what the option means in
        // PLAXIS too. Reported per material, because a phase that quietly stops being undrained
        // is the difference between a short-term and a long-term answer.
        if (io.config && io.config->ignore_undrained && models.back().undrained) {
            models.back().undrained = false;
            note(R, "K2D-A008", m.name,
                 "This phase ignores undrained behaviour, so material \"" + m.name +
                     "\" is solved DRAINED: no excess pore pressure is generated in it. Its "
                     "strength parameters are unchanged. This is a long-term (or state-setting) "
                     "answer, not a short-term one.");
        }
        // The Rankine tension cut-off is consumed by the Mohr-Coulomb return only: the hardening
        // and soft-soil integrators do not carry the extra planes (material_model.hpp). PLAXIS
        // applies a tension cut-off to these models by DEFAULT, and so does this schema, so a
        // silent omission here is a systematic difference from the reference code in the
        // unsafe direction -- the soil takes tension it should not.
        // CLOSED 2026-08-13: these models now read the cut-off (MMM Eq. 3-11, applied to the
        // principal stresses their own return produced -- materials/mohr_coulomb.hpp,
        // apply_rankine_cap). What remains is a formulation boundary worth stating, and only
        // where it can actually bite: the cap is applied SEQUENTIALLY after the model's own
        // surfaces rather than as one coupled multi-surface solve. The correction is compressive,
        // so it moves away from every tensile surface; the one surface it can in principle
        // disturb is the Hardening Soil family's CAP, and that is not iterated back.
        if (m.tension_cutoff && entry->hardening_family)
            note(R, "K2D-M001", m.name,
                 "Material \"" + m.name + "\" (" + constitutive_name(m.model) +
                     ") applies the tension cut-off at sigma_t = " + dnum(m.tensile_strength) +
                     " kPa after its own return, not as one coupled multi-surface solve. Where "
                     "the cut-off and the volumetric cap are active at the same point, the cap "
                     "is not re-checked against the capped stress.");
        // The model does not take a small-strain stiffness more than 20x its own unload/reload
        // stiffness (MMM sec. 7.5: "Although Alpan suggests that the ratio E0/Eur can exceed 10
        // for very soft clays, the maximum ratio E0/Eur or G0/Gur permitted in the HSsmall model
        // is limited to 20"). The cap is applied in the engine, so the run is the reference
        // code's run rather than a stiffer one -- and it is said out loud, because a G0 that is
        // quietly reduced describes a different soil from the one the file asks for.
        if (m.model == model::SoilModel::HSsmall && m.G0ref > 0.0) {
            const double Gur_ref = m.Eurref / (2.0 * (1.0 + m.nu_ur));
            const double cap = katai::core::HardeningSoilParams::kMaxG0Ratio * Gur_ref;
            if (m.G0ref > cap)
                warn(R, "K2D-M004", m.name,
                     "Material \"" + m.name + "\" asks for G0 = " + dnum(m.G0ref) +
                         " kPa, which is " + dnum(m.G0ref / Gur_ref) +
                         " times its unload/reload shear modulus G_ur = " + dnum(Gur_ref) +
                         " kPa. The model permits at most 20, so it is solved with G0 = " +
                         dnum(cap) + " kPa. Lower G0, or raise Eur, to model the soil as entered.");
        }
        // The pore fluid of a Hardening Soil material is sized at the REFERENCE unload/reload
        // stiffness. That is the model's own elastic pair -- not the unread E box it used to be
        // read from -- but HS stiffness is stress-dependent, so a soil far from p_ref carries a
        // Kw/n that a stress-dependent derivation would not give it. Said out loud rather than
        // left to be discovered from the pore pressures.
        if (models.back().undrained && entry->hardening_family)
            note(R, "K2D-M002", m.name,
                 "The pore-fluid stiffness of \"" + m.name +
                     "\" is derived from the unload/reload pair (Eur_ref = " + dnum(m.Eurref) +
                     " kPa, nu_ur = " + dnum(m.nu_ur) +
                     "), the Hardening Soil model's own elastic constants, evaluated at the "
                     "reference pressure. It does not follow the stress-dependent Eur(sigma3) "
                     "during the run, so far from p_ref the excess pore pressure is that of a "
                     "reference-stiffness pore fluid.");
        // Undrained (C): the results of this material are in a different currency from the rest
        // of the model, and nothing in the output says which is which -- the stress field has one
        // name. Saying it here is the whole reason the diagnostics list exists.
        if (m.drainage == model::Drainage::UndrainedC)
            note(R, "K2D-M003", m.name,
                 "Material \"" + m.name +
                     "\" is Undrained (C), a TOTAL stress analysis: E and nu are read as the "
                     "undrained pair, c as the undrained shear strength with phi = 0, and no pore "
                     "pressure is generated or carried in it. The stresses reported for this "
                     "material are TOTAL stresses, and its K0 refers to total stress -- they "
                     "cannot be compared with the effective stresses of a Drained or Undrained "
                     "(A)/(B) region in the same model.");
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
    // `mat_total_stress` is the wider set: materials whose equilibrium is stated in TOTAL
    // stress and which therefore carry no pore pressure -- Non-porous because it holds no
    // water, Undrained (C) because it declines to separate the water from the skeleton
    // (MMM section 2.7: "all pore pressures are equal to zero"). Everything downstream that
    // asks "does this element get a pore pressure?" asks this mask; the places that ask "does
    // this material hold water at all?" (its saturated weight) keep asking about Non-porous.
    std::vector<char> mat_nonporous, mat_total_stress;
    bool any_nonporous = false, any_total_stress = false;
    for (const auto& m : pr.materials) {
        const char np = m.drainage == model::Drainage::NonPorous ? 1 : 0;
        const char ts = (np || m.drainage == model::Drainage::UndrainedC) ? 1 : 0;
        mat_nonporous.push_back(np);
        mat_total_stress.push_back(ts);
        any_nonporous |= (np != 0);
        any_total_stress |= (ts != 0);
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
    const auto eff_unit_weight = [&pr, &io, water, &material_at](double x, double y) {
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
        const double g_total = (!np && water && y < water_table_at(pr, x, io.config)) ? m.gamma_sat
                                                                          : m.gamma_unsat;
        return (water && y < water_table_at(pr, x, io.config)) ? (g_total - kGammaWater) : g_total;
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
            // silent-drop-ok: a zero-length structural line is an ERROR at structs[i].x2 in the
            // input contract, so it cannot reach a run; the guard is arithmetic self-defence.
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
            // The line is not on mesh edges, so the mesh cannot be split along it and the wall
            // falls back to a plain plate bonded to the soil: no interface, no slip, no gap. That
            // is a stiffer and generally UNCONSERVATIVE model than the one drawn, so it is stated
            // rather than assumed (the plate itself is still built by the plate loop below).
            if ((int)w.seam.size() < (order == 15 ? 4 : 2)) {
                warn(R, "K2D-G009", line_subject(s.name, s.x1, s.y1, s.x2, s.y2),
                     "Wall \"" + s.name +
                         "\" could not be split from the soil along its line, so its interfaces "
                         "were not built: it acts as a plate BONDED to the soil (no slip, no gap). "
                         "Align the wall with mesh edges or refine around it to get the interface.");
                continue;
            }
            // How much of the drawn wall the mesh actually gave: a wall drawn past the top of the
            // soil is split only where soil exists, so the run models a SHORTER wall than the one
            // in the file -- and a retaining wall two metres shorter is a different structure.
            // The toe node is shared by design (split_mesh_at_wall keeps y > y_toe), so the seam
            // legitimately starts one node spacing above the toe; the threshold is twice that
            // spacing, which is the mesh's own scale rather than the project's target size.
            {
                double t0 = 1.0, t1 = 0.0;
                for (const auto& p : w.seam) {
                    const double t = ((mesh.x[p.right] - tx) * (ux - tx) +
                                      (mesh.y[p.right] - ty) * (uy - ty)) / (len * len);
                    t0 = std::fmin(t0, t);
                    t1 = std::fmax(t1, t);
                }
                const double built = (t1 - t0) * len, spacing = len / (double)w.seam.size();
                if (built < len - 2.0 * spacing)
                    warn(R, "K2D-G006", line_subject(s.name, s.x1, s.y1, s.x2, s.y2),
                         "Wall \"" + s.name + "\" is built over " + dnum(built) + " m of the " +
                             dnum(len) +
                             " m drawn: the rest of the line falls outside the soil, and both its "
                             "stiffness and its interfaces are those of the shorter wall.");
            }
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
    // Strength/stiffness come from the adjacent soil + Rinter. Boundary nodes that get duplicated need
    // their support decided rather than copied -- see the seam pass after the loop.
    // (interface-formulation.md.)
    struct IfaceSpec {
        std::vector<katai::core::SegSeam> seam;
        katai::core::iface::InterfaceProps ip;
        double nx, ny;          // unit normal (orientation-aware K0 seed)
        int soil_mat;
        std::string name;
    };
    std::vector<IfaceSpec> soil_ifaces;
    std::vector<char> bc_released;   // nodes the domain boundary must NOT fix (filled by the seam pass)
    {
        std::vector<char> is_bnode(mesh.node_count, 0);
        for (int n : mesh.boundary_nodes) is_bnode[n] = 1;
        const int order = mesh.nodes_per_element;
        std::vector<std::pair<int, int>> boundary_seam;   // (orig, dup) pairs sitting ON the boundary
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
            // silent-drop-ok: a zero-length structural line is an ERROR at structs[i].x2 in the
            // input contract, so it cannot reach a run; the guard is arithmetic self-defence.
            if (len < 1e-9) continue;
            IfaceSpec sp;
            sp.name = s.name.empty() ? "Interface" : s.name;
            sp.nx = dy / len; sp.ny = -dx / len;   // unit normal (matches split_mesh_at_segment frame)
            sp.seam = katai::core::split_mesh_at_segment(mesh, s.x1, s.y1, s.x2, s.y2, -1e-6, len + 1e-6);
            const int per_edge = order == 15 ? 4 : 2;
            // Same fallback as the wall, and the same reason for saying so: with no split there is
            // no joint, and the soil across the drawn slip surface stays fully bonded.
            if ((int)sp.seam.size() < per_edge + 1) {
                warn(R, "K2D-G009", line_subject(s.name, s.x1, s.y1, s.x2, s.y2),
                     "Interface \"" + sp.name +
                         "\" does not lie on mesh edges, so no joint was created: the soil across "
                         "it stays BONDED (no slip, no gap). Align it with mesh edges or refine "
                         "around it.");
                continue;
            }
            // A slip surface drawn past the edge of the soil is split only where soil exists, so
            // the joint is shorter than the drawn line. Here the seam carries arc length directly
            // and both ends are inside the split window, so one node spacing is threshold enough.
            {
                double s0 = len, s1 = 0.0;
                for (const auto& p : sp.seam) { s0 = std::fmin(s0, p.s); s1 = std::fmax(s1, p.s); }
                const double spacing = len / (double)sp.seam.size();
                if (s1 - s0 < len - 2.0 * spacing)
                    warn(R, "K2D-G006", line_subject(s.name, s.x1, s.y1, s.x2, s.y2),
                         "Interface \"" + sp.name + "\" is built over " + dnum(s1 - s0) +
                             " m of the " + dnum(len) +
                             " m drawn: the rest of the line falls outside the soil, and the soil "
                             "there stays bonded.");
            }
            // A seam node that sits on the domain boundary needs its support DECIDED, not copied;
            // the decision needs the finished mesh, so collect the pairs and settle them below.
            // `is_bnode` was sized before the first split, so a seam that catches a node an EARLIER
            // interface created indexes past it -- guard the read rather than the growth: a twin is
            // never a domain-boundary node in its own right, so "not in the snapshot" is the answer.
            for (const auto& p : sp.seam)
                if (p.orig < (int)is_bnode.size() && is_bnode[p.orig])
                    boundary_seam.push_back({p.orig, p.dup});
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
        // Who holds the support on a seam that touches the domain boundary. The two sides of a seam
        // sit at the SAME coordinates and apply_boundary_conditions matches by coordinate, so this
        // cannot be left to it: measured 2026-08-10, an interface drawn ALONG a fixed boundary got
        // both of its sides fixed and became inert -- the mesh split, the joint assembled, and the
        // block welded itself to its own base (PLAXIS Validation Manual 3.3 read 5.4e6 kN/m instead
        // of 60). The discriminator is which side carries soil once the split has re-wired the
        // elements, which is why this runs after every interface has been built.
        if (!boundary_seam.empty()) {
            std::vector<int> refs(mesh.node_count, 0);
            for (int e = 0; e < mesh.element_count; ++e)
                for (int k = 0; k < mesh.nodes_per_element; ++k) ++refs[mesh.node_of(e, k)];
            bc_released.assign(mesh.node_count, 0);
            for (const auto& [orig, dup] : boundary_seam) {
                const bool orig_soil = refs[orig] > 0, dup_soil = refs[dup] > 0;
                if (orig_soil && dup_soil) {
                    // The seam CROSSES the boundary here (a slip surface reaching the fixed base):
                    // both sides carry soil, so both keep the support or the split opens a hole.
                    mesh.boundary_nodes.push_back(dup);
                } else if (dup_soil) {
                    // The seam runs ALONG the boundary and the twin took the soil. The original
                    // carries nothing: it IS the rigid outside world, it is already in the boundary
                    // list, and the twin must stay free to slide against it. Nothing to do.
                } else {
                    // Mirror image -- the drawn line's direction decides which side is duplicated,
                    // and a user cannot be expected to know that. Hand the support to the empty
                    // twin and take it away from the soil-carrying original. (`!orig_soil &&
                    // !dup_soil` lands here too: nothing to hold either way, and the twin taking
                    // the fixity keeps the node count of fixed DOFs unchanged.)
                    mesh.boundary_nodes.push_back(dup);
                    if (orig_soil) bc_released[orig] = 1;
                }
            }
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
        const std::string psub = line_subject(s.name, s.x1, s.y1, s.x2, s.y2);
        // Structural lines are mesh constraints too, so an unusable chain means the plate is not
        // on the soil. Refused rather than skipped: a retaining wall that silently does not exist
        // returns the unretained field and calls it a success.
        if (chain.size() < 3 || chain.size() % 2 == 0) {
            refuse(R, "K2D-G005", psub,
                   "Plate \"" + s.name + "\" from (" + dnum(s.x1) + ", " + dnum(s.y1) + ") to (" +
                       dnum(s.x2) + ", " + dnum(s.y2) +
                       ") does not lie on the mesh, so it would carry nothing. Draw it inside or "
                       "along a soil region.");
            return R;
        }
        {
            double t0 = 0.0, t1 = 1.0;
            chain_span(mesh, chain, s.x1, s.y1, s.x2, s.y2, t0, t1);
            if (chain_is_clipped(t0, t1)) {
                const double drawn = std::hypot(s.x2 - s.x1, s.y2 - s.y1);
                warn(R, "K2D-G006", psub,
                     "Plate \"" + s.name + "\" is built over " + dnum((t1 - t0) * drawn) +
                         " m of the " + dnum(drawn) +
                         " m drawn: the rest of the line falls outside the soil, and its forces "
                         "are those of the shorter element.");
            }
        }
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
        const std::string gsub = line_subject(s.name, s.x1, s.y1, s.x2, s.y2);
        if (chain.size() < 3 || chain.size() % 2 == 0) {
            refuse(R, "K2D-G007", gsub,
                   "Geogrid \"" + s.name + "\" from (" + dnum(s.x1) + ", " + dnum(s.y1) +
                       ") to (" + dnum(s.x2) + ", " + dnum(s.y2) +
                       ") does not lie on the mesh, so it would carry nothing. Draw it inside a "
                       "soil region.");
            return R;
        }
        {
            double t0 = 0.0, t1 = 1.0;
            chain_span(mesh, chain, s.x1, s.y1, s.x2, s.y2, t0, t1);
            if (chain_is_clipped(t0, t1)) {
                const double drawn = std::hypot(s.x2 - s.x1, s.y2 - s.y1);
                warn(R, "K2D-G006", gsub,
                     "Geogrid \"" + s.name + "\" is built over " + dnum((t1 - t0) * drawn) +
                         " m of the " + dnum(drawn) +
                         " m drawn: the rest of the line falls outside the soil.");
            }
        }
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
    // Is this anchor end in the soil -- i.e. is there something there to pull on? The question is
    // answered against the MESH, not against the polygon outline: an even-odd point-in-polygon
    // test reports a point lying exactly ON an edge as outside, and an anchor head placed exactly
    // on the ground surface or on a slope face is ordinary input, not an error. Reading it as
    // "outside" put such an anchor on the both-ends-outside path, where it was dropped without a
    // word (tests/test_anchor_mesh_repro.cpp had a fixture doing exactly that). Locating the
    // containing element also removes a discontinuity: a head one millimetre inside the slope and
    // a head exactly on it now build the same structural model.
    const auto in_soil = [&](double x, double y) {
        return katai::core::ploc::locate_point(mesh, x, y).found;
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
            // Per metre of wall, like the stiffness and the capacities: the schema states the
            // lock-off force of ONE anchor and the spacing it repeats at.
            an.prestress = am.prestress / Ls;
            if (am.elastoplastic) {
                an.Fmax_tens = am.Fmax_tens > 0.0 ? am.Fmax_tens / Ls : -1.0;
                an.Fmax_comp = am.Fmax_comp > 0.0 ? am.Fmax_comp / Ls : -1.0;
            }
        } else {
            an.EA = 1.0e5;
        }
        an.L = 0.0;  // EA/L uses the geometric distance
        const bool in1 = in_soil(s.x1, s.y1), in2 = in_soil(s.x2, s.y2);
        const std::string asub = line_subject(s.name, s.x1, s.y1, s.x2, s.y2);
        // An anchor needs at least one end in the soil to pull on. Both ends outside is not an
        // unusual anchor, it is no anchor at all, and a strutted excavation that quietly loses
        // its strut is exactly the run that must never report success.
        if (!in1 && !in2) {
            refuse(R, "K2D-G008", asub,
                   "Anchor \"" + s.name + "\" has neither end in the soil: (" + dnum(s.x1) + ", " +
                       dnum(s.y1) + ") and (" + dnum(s.x2) + ", " + dnum(s.y2) +
                       ") both fall outside every soil region, so it would carry nothing.");
            return R;
        }
        if (in1 && in2) {                          // node-to-node (strut / internal support)
            an.node_a = nearest_node(s.x1, s.y1);
            an.node_b = nearest_node(s.x2, s.y2);
            if (an.node_a < 0 || an.node_b < 0 || an.node_a == an.node_b) {
                refuse(R, "K2D-G008", asub,
                       "Anchor \"" + s.name +
                           "\" collapses onto a single mesh node, so it would carry nothing. Make "
                           "it longer than one element, or refine the mesh around it.");
                return R;
            }
        } else {                                   // fixed-end (one end outside the soil)
            const double sx = in1 ? s.x1 : s.x2, sy = in1 ? s.y1 : s.y2;
            const double fx = in1 ? s.x2 : s.x1, fy = in1 ? s.y2 : s.y1;
            an.node_a = nearest_node(sx, sy);
            if (an.node_a < 0) {
                refuse(R, "K2D-G008", asub,
                       "Anchor \"" + s.name + "\" has no mesh node at its soil end (" + dnum(sx) +
                           ", " + dnum(sy) + "), so it would carry nothing.");
                return R;
            }
            an.node_b = -1; an.fixed_point = {fx, fy};
        }
        // An anchor end attaches to the NEAREST node, exactly like a point load, so the same
        // discretisation question applies: a strut whose head moved a third of an element is
        // still the drawn strut, one that moved several elements is not. Both ends are measured
        // and the larger move is reported.
        {
            // Each attached node, beside the endpoint it was attached FOR: node-to-node takes
            // both endpoints in order, a fixed-end anchor only its soil end (which may be
            // either endpoint, so it is named rather than assumed).
            const double sx = in1 ? s.x1 : s.x2, sy = in1 ? s.y1 : s.y2;
            const bool fixed_end = an.node_b < 0;
            struct Attached { int node; double x, y; };
            const Attached att[2] = {{an.node_a, fixed_end ? sx : s.x1, fixed_end ? sy : s.y1},
                                     {an.node_b, s.x2, s.y2}};
            double moved = 0.0, at_x = 0.0, at_y = 0.0, h = 0.0;
            for (const Attached& a : att) {
                if (a.node < 0) continue;
                const double d = std::hypot(mesh.x[a.node] - a.x, mesh.y[a.node] - a.y);
                if (d > moved) {
                    moved = d;
                    at_x = mesh.x[a.node]; at_y = mesh.y[a.node];
                    h = element_size_at(mesh, at_x, at_y);
                }
            }
            if (h > 0.0 && moved > 0.25 * h)
                warn(R, "K2D-G002", asub,
                     "Anchor \"" + s.name + "\" attaches at node (" + dnum(at_x) + ", " +
                         dnum(at_y) + "), " + dnum(moved) +
                         " m from the end drawn for it -- the element there is " + dnum(h) +
                         " m. Refine the mesh around the anchor if that matters.");
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
        // The toe is the node the plate hangs from, so without it there is no wall at all --
        // not a shorter one, none: no plate, no interfaces, and a run that returns the field of
        // an unsupported excavation. It goes missing when the deeper end is drawn outside the
        // soil, because the mesher clips a structural line to the soil before it becomes a
        // constraint, and the endpoint is then never inserted as a vertex.
        if (toe < 0) {
            refuse(R, "K2D-G010", w.name.empty() ? "wall" : w.name,
                   "Wall \"" + w.name + "\" has no mesh node at its toe (" + dnum(w.toe_x) +
                       ", " + dnum(w.toe_y) +
                       "): that end is outside the soil, and without the toe neither the plate "
                       "nor its interfaces can be built. Draw the wall inside a soil region.");
            return R;
        }
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
        const std::string esub = line_subject(s.name, s.x1, s.y1, s.x2, s.y2);
        // A pile without a material or with no section has no stiffness to contribute, so it
        // would be drawn in the model and absent from the analysis.
        if (s.material < 0 || s.material >= (int)pr.embedded.size()) {
            refuse(R, "K2D-G012", esub,
                   "Embedded beam \"" + s.name +
                       "\" has no embedded-beam material assigned, so it would carry nothing.");
            return R;
        }
        const auto& em = pr.embedded[s.material];
        const double L = std::hypot(s.x2 - s.x1, s.y2 - s.y1);
        if (L < 1e-9 || em.diameter <= 0.0) {
            refuse(R, "K2D-G012", esub,
                   "Embedded beam \"" + s.name + "\" has length " + dnum(L) +
                       " m and diameter " + dnum(em.diameter) +
                       " m: both must be positive, or it would carry nothing.");
            return R;
        }
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
        // CONNECTION POINT (PLAXIS Ref. sec 5.6.3). Hinged -- PLAXIS's default when no structure
        // shares the point -- ties the beam's top translations to the mesh node there, which the
        // mesher carries as a vertex precisely so this is an exact DOF identity and not an
        // interpolation or a penalty. Free leaves the top coupled through the skin springs alone.
        int conn_beam_node = -1, conn_mesh_node = -1;
        if (s.conn == 0) {
            double ccx = 0.0, ccy = 0.0;
            embedded_connection_point(s, ccx, ccy);
            // Which END of the (possibly reversed) polyline is the connection point?
            conn_beam_node = std::hypot(bx.front() - ccx, by.front() - ccy) <=
                             std::hypot(bx.back() - ccx, by.back() - ccy) ? 0 : nbn - 1;
            double bestd = 1e300;
            for (int n = 0; n < mesh.node_count; ++n) {
                const double d = std::hypot(mesh.x[n] - ccx, mesh.y[n] - ccy);
                if (d < bestd) { bestd = d; conn_mesh_node = n; }
            }
            // The tie lands on the NEAREST node, which is what this tree does with every
            // structural point attachment (a point load, an anchor end) and for the same reason:
            // the mesher is not asked to insert a vertex, because measuring that showed a single
            // isolated interior point costing more refinement than conforming the whole shaft
            // would (mesh_builder.cpp). The two questions are answered separately, as with a
            // point load: "is it on the soil?" is geometry, and "is the snap large?" is
            // discretisation, measured against the element size THERE.
            const double h_here = element_size_at(mesh, ccx, ccy);
            if (conn_mesh_node < 0 || h_here <= 0.0) {
                refuse(R, "K2D-G012", esub,
                       "Embedded beam \"" + s.name + "\" has a hinged connection point at (" +
                           dnum(ccx) + ", " + dnum(ccy) +
                           ") that is not on the mesh, so its head would be tied to nothing. "
                           "Move the pile head onto a soil region, or set the connection to free.");
                return R;
            }
            if (bestd > 0.25 * h_here)
                warn(R, "K2D-G013", esub,
                     "Embedded beam \"" + s.name + "\" is connected to the soil at node (" +
                         dnum(mesh.x[conn_mesh_node]) + ", " + dnum(mesh.y[conn_mesh_node]) +
                         "), " + dnum(bestd) + " m from the head as drawn (" + dnum(ccx) + ", " +
                         dnum(ccy) + ") -- the element there is " + dnum(h_here) +
                         " m. Refine the mesh at the pile head if that distance matters.");
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
        katai::core::ebeam::default_interface_stiffness(Ls, em.diameter, G, pp.EA, pp.EI,
                                                        k_axial, k_lateral, D_foot);
        const double t_max = em.Tskin_max > 0.0 ? em.Tskin_max / Ls : -1.0;
        const double f_max = em.Fmax_base > 0.0 ? em.Fmax_base / Ls : -1.0;
        structures.embedded_beams.push_back(katai::core::ebeam::build_embedded_beam(
            mesh, dofs, bx, by, pp, k_axial, k_lateral, t_max, D_foot, f_max,
            conn_beam_node, conn_mesh_node));
        // Force-diagram bookkeeping (kind 3 = embedded beam; produces N/Q/M like a plate).
        diag_specs.push_back({3, s.name, structures.embedded_beams.size() - 1,
                              structures.embedded_beams.size()});
    }
    const bool has_embedded = !structures.embedded_beams.empty();

    // Compliant (absorbing) base: only a Dynamic phase that asked for it frees the base u_x --
    // every static phase and the rigid-base dynamic default keep the user's fixities bit-for-bit.
    const bool compliant_base =
        phase == InitialPhase::Dynamic && io.config && io.config->seismic_compliant_base;
    katai::core::apply_boundary_conditions(bc_edges_from(pr), mesh, dofs, compliant_base,
                                           bc_released.empty() ? nullptr : &bc_released);
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
            // silent-drop-ok: a zero-length prescribed-displacement line is an ERROR at
            // disps[i].x2 in the input contract, so it cannot reach a run.
            if (l2 < 1e-18) continue;
            int fixed_nodes = 0;
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
                ++fixed_nodes;
            }
            // A prescribed-displacement line that catches no node imposes nothing. The run then
            // reports a converged, unsettled model -- and every displacement-controlled benchmark
            // in the corpus (a rigid footing pushed into the soil, and the reaction read back off
            // it) depends on this line actually reaching the mesh.
            if (fixed_nodes == 0) {
                refuse(R, "K2D-G011", line_subject(D.name, D.x1, D.y1, D.x2, D.y2),
                       "Prescribed displacement \"" + D.name + "\" from (" + dnum(D.x1) + ", " +
                           dnum(D.y1) + ") to (" + dnum(D.x2) + ", " + dnum(D.y2) +
                           ") does not lie on the mesh: no node was constrained, so nothing would "
                           "be imposed. Draw it along a soil boundary or inside a soil region.");
                return R;
            }
        }
        // A structural element standing on a driven node does NOT see that motion: the
        // prescribed-displacement ramp enters the soil element loop alone, and the structural
        // elements read fixed DOFs as zero (internal_forces.hpp states the limit at the seam
        // that has it). Its M / Q / N are then those of an undriven element -- near zero for a
        // plate whose only load is the imposed settlement. The condition is checked rather than
        // assumed, so a model whose structures stand clear of the driven line hears nothing.
        //
        // A component prescribed to ZERO is not a motion, it is a SUPPORT -- the format
        // document's own words: "a set component with value 0 is a rigid support line". There
        // is nothing for the structure to miss, its diagram is right, and warning here would
        // send the user away from the very diagram this schema's only point-support idiom
        // makes correct. Only a nonzero imposed value can be understated.
        if (!presc_entries.empty()) {
            std::vector<char> driven(mesh.node_count, 0);
            for (const auto& e : presc_entries)
                if (e.value != 0.0 && e.node >= 0 && e.node < mesh.node_count) driven[e.node] = 1;
            std::string hit;
            const auto touches = [&](int n) { return n >= 0 && n < mesh.node_count && driven[n]; };
            for (const auto& p : structures.plates)
                for (int n : p.nodes) if (touches(n) && hit.empty()) hit = "a plate";
            for (const auto& p : structures.plates5)
                for (int n : p.nodes) if (touches(n) && hit.empty()) hit = "a plate";
            for (const auto& g : structures.geogrids)
                for (int n : g.nodes) if (touches(n) && hit.empty()) hit = "a geogrid";
            for (const auto& a : structures.anchors)
                if ((touches(a.node_a) || touches(a.node_b)) && hit.empty()) hit = "an anchor";
            if (!hit.empty())
                warn(R, "K2D-A003", "",
                     "A prescribed displacement drives a node that " + hit +
                         " stands on. Structural elements do not receive the imposed motion in "
                         "this build, so that element's internal forces (M, Q, N) understate the "
                         "action -- read the soil results, not the structural diagram, here.");
        }
    }
    // Staged construction: nodes touched only by passive (excavated / not-yet-filled) elements
    // would be singular -- fix them (their displacement is meaningless this phase).
    // A node a PLATE runs through is NOT orphaned -- the beam's own stiffness holds it. Pinning
    // it would weld the structure to the outside world and the phase would converge on a model
    // in which that plate carries nothing, without saying so. The manual's own beam verification
    // (PLAXIS 2D Validation Manual sec. 2.3) is built exactly this way: the soil cluster is
    // deactivated so that only the beams remain, supported at their end points.
    if (!act.empty()) {
        std::vector<char> carried(mesh.node_count, 0);
        const auto carry = [&](int n) { if (n >= 0 && n < mesh.node_count) carried[(size_t)n] = 1; };
        for (const auto& p : structures.plates) for (int n : p.nodes) carry(n);
        for (const auto& p : structures.plates5) for (int n : p.nodes) carry(n);
        katai::core::fix_inactive_nodes(mesh, act, dofs, carried);
    }
    dofs.finalize();
    // Reported support reactions are the SOIL's discrete B^T sigma: a structural element that
    // ends on a support adds its own end force to that support in reality, and this build does
    // not include it (results.hpp states the same limit at the field). It is asked AFTER finalize(),
    // because that is when a DOF's fixity becomes readable. It matters exactly when a
    // structure stands on a fixed node, so that is the condition -- a note rather than a warning,
    // because the number is right for what it says it is, and only incomplete for what a reader
    // may take it to mean.
    {
        std::vector<int> struct_nodes;
        for (const auto& p : structures.plates) for (int n : p.nodes) struct_nodes.push_back(n);
        for (const auto& p : structures.plates5) for (int n : p.nodes) struct_nodes.push_back(n);
        for (const auto& g : structures.geogrids) for (int n : g.nodes) struct_nodes.push_back(n);
        for (const auto& a : structures.anchors) {
            struct_nodes.push_back(a.node_a);
            if (a.node_b >= 0) struct_nodes.push_back(a.node_b);
        }
        bool on_support = false;
        for (int n : struct_nodes) {
            if (n < 0 || n >= mesh.node_count) continue;
            for (int c = 0; c < 2 && !on_support; ++c)
                if (dofs.is_fixed(dofs.global_dof(n, c))) on_support = true;
            if (on_support) break;
        }
        if (on_support)
            note(R, "K2D-A004", "",
                 "A structural element ends on a supported node. The reported reactions are the "
                 "soil's contribution only; the element's own end force at that support is not "
                 "included in this build. Read it from the element's own force diagram.");
    }
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

    // Staged-phase scope guards: the static path below is axisymmetric-aware end to end
    // (r-weighted gravity and internal-force baseline, axisym Newton kinematics, axisym
    // reactions), so a PLASTIC phase that keeps every element active -- load / prescribed-
    // displacement staging, the bearing-capacity workflow -- runs the same verified machinery
    // as the initial phase. What is NOT plumbed for axisymmetry is refused honestly rather
    // than silently integrated as plane strain: the time-dependent phase inputs carry no
    // kinematics flag, and the axisym gravity / internal-force assemblies take no
    // element-activity mask.
    if (io.config && axi) {
        if (io.config->type != model::PhaseType::Plastic) {
            R.message = "Only Plastic phases are available in axisymmetric mode yet "
                        "(consolidation, flow, dynamic and safety phases are plane-strain only).";
            return R;
        }
        if (!act.empty()) {
            R.message = "Excavation / fill (activation changes) is not available in "
                        "axisymmetric mode yet -- keep every soil polygon active.";
            return R;
        }
    }

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
    if (any_total_stress) {
        act_pore.assign(mesh.element_count, 1);
        for (int e = 0; e < mesh.element_count; ++e) {
            if (!act.empty() && !act[e]) { act_pore[e] = 0; continue; }
            const int mt = mesh.element_material[e];
            if (mt >= 0 && mt < (int)mat_total_stress.size() && mat_total_stress[mt]) act_pore[e] = 0;
        }
    }
    const std::vector<char>& pore_mask = act_pore.empty() ? act : act_pore;

    // External load vector (free DOFs): gravity body force + point loads (+ pore pressures).
    Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
    if (flow) {
        // A barrier makes the head DISCONTINUOUS across its line, and this hand-over carries one
        // value per node: the pore load is therefore applied from one side of the wall, so the
        // differential water pressure on it is missing from the structural forces. Said out loud
        // -- it is the load a cut-off wall is designed for.
        for (const auto& st : pr.structs)
            if (st.flow_barrier != 0 && (st.kind == model::StructKind::Plate ||
                                         st.kind == model::StructKind::Interface)) {
                warn(R, "K2D-A011", st.name,
                     "The flow field driving this phase was computed with \"" + st.name +
                         "\" as a barrier, so the head jumps across it -- but the pore load here "
                         "is applied from ONE side. The differential water pressure on that wall "
                         "is not in its structural forces. Read both sides from the flow result "
                         "(head and head_far) until the deformation mesh is split with it.");
                break;
            }
        // Pore pressure + saturation from the steady-state seepage head field (PLAXIS groundwater
        // flow pore generation): total-stress gravity with gamma_sat where psi = h - y >= 0, and
        // the pore load interpolated from the SAME nodal head -> recovered stress is EFFECTIVE.
        katai::core::assemble_gravity_from_head(mesh, dofs, gamma, gamma_sat, *flow_head, f);
        katai::core::assemble_pore_load_from_head(mesh, dofs, *flow_head, kGammaWater, f, pore_mask);
    } else if (water) {
        // Total-stress equilibrium: saturated weight below the phreatic surface, moist above.
        const auto wt = [&pr, &io](double x) { return water_table_at(pr, x, io.config); };
        katai::core::assemble_gravity_phreatic(mesh, dofs, gamma, gamma_sat, wt, f, act);
        // Hydrostatic pore-pressure load -> recovered stress is EFFECTIVE (Terzaghi; buoyancy
        // gamma' = gamma_sat - gamma_w emerges naturally). docs/references/effective-stress-formulation.md.
        const auto pore = [&pr, &io](double x, double y) {
            return kGammaWater * std::fmax(0.0, water_table_at(pr, x, io.config) - y);
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
        // A point load is carried by the nearest node. Unlike a line load, the mesher does not
        // insert the point as a vertex -- it only refines around it (mesh_builder.cpp, SizeSrc)
        // -- so landing between nodes is normal and a snap of a fraction of an element is the
        // expected discretisation. What is NOT expected is a point off the soil altogether: the
        // nearest-node search always succeeds, so such a load used to be relocated to whatever
        // node happened to be closest, however far, and the run reported success for a model the
        // engineer never drew. An input may be used differently from the way it was written, but
        // never in silence -- see docs/diagnostics.md for the codes that rule.
        //
        // The two questions are answered by two different instruments, on purpose. "Is the load
        // on the soil?" is geometry, answered exactly by locating the containing element -- no
        // length scale, so it holds under any mesh density or coarseness factor. "Is the snap
        // large?" is discretisation, so it is measured against that element's OWN size rather
        // than the project's target size, which a region coarseness factor may multiply by up
        // to four.
        const std::string lsub = line_subject(L.name, L.x1, L.y1, L.x1, L.y1);
        const double h_elem = element_size_at(mesh, L.x1, L.y1);
        if (best < 0 || h_elem <= 0.0) {
            refuse(R, "K2D-G001", lsub,
                   "Point load \"" + L.name + "\" at (" + dnum(L.x1) + ", " + dnum(L.y1) +
                       ") is not on the mesh: no element contains it" +
                       (best < 0 ? std::string()
                                 : ", and the nearest node is " + dnum(bestd) + " m away") +
                       ". Place the load on a soil region, or it carries nothing.");
            return R;
        }
        if (bestd > 0.25 * h_elem)
            warn(R, "K2D-G002", lsub,
                 "Point load \"" + L.name + "\" acts at node (" + dnum(mesh.x[best]) + ", " +
                     dnum(mesh.y[best]) + "), " + dnum(bestd) + " m from where it is drawn (" +
                     dnum(L.x1) + ", " + dnum(L.y1) + ") -- the element there is " +
                     dnum(h_elem) + " m. Refine the mesh there if that matters.");
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
            const std::string lsub = line_subject(L.name, L.x1, L.y1, L.x2, L.y2);
            // The mesher adds every distributed load line as a mesh constraint, so a line that
            // does not come back as an edge chain is not on the soil -- it was drawn above the
            // surface, outside the model, or across a hole. Refused: the alternative is a phase
            // that runs unloaded and reports "ok" (docs/diagnostics.md, the refusal rule).
            if (cs < npe || (cs - 1) % (npe - 1) != 0) {
                refuse(R, "K2D-G003", lsub,
                       "Distributed load \"" + L.name + "\" from (" + dnum(L.x1) + ", " +
                           dnum(L.y1) + ") to (" + dnum(L.x2) + ", " + dnum(L.y2) +
                           ") does not lie on the mesh, so it would carry nothing. Draw it along a "
                           "soil boundary or inside a soil region.");
                return R;
            }
            // It resolved, but possibly only in part: a line drawn past the edge of the soil is
            // applied over the stretch the mesh could give it. That is a defensible run and a
            // different load from the one drawn, so it is said out loud.
            double t0 = 0.0, t1 = 1.0;
            chain_span(mesh, chain, L.x1, L.y1, L.x2, L.y2, t0, t1);
            if (chain_is_clipped(t0, t1)) {
                const double drawn = std::hypot(L.x2 - L.x1, L.y2 - L.y1);
                warn(R, "K2D-G004", lsub,
                     "Distributed load \"" + L.name + "\" is applied over " +
                         dnum((t1 - t0) * drawn) + " m of the " + dnum(drawn) +
                         " m drawn: the rest of the line falls outside the soil.");
            }
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
    // Total-stress targets are seeded in total stress (header note): u(x,y) hydrostatic. For
    // Undrained (C) this is the manual's "the K0-value refers to total stresses rather than
    // effective stresses in this case" -- the seed has to be the stress the material is
    // analysed in, or the initial state contradicts the constitutive law from the first step.
    if (any_total_stress && water) {
        k0opt.nonporous = mat_total_stress;
        k0opt.pore = [&pr, &io](double x, double y) {
            return kGammaWater * std::fmax(0.0, water_table_at(pr, x, io.config) - y);
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
        // silent-drop-scope: none -- this pair scan DECIDES something (is the geostatic field
        // level?); it builds nothing from the regions, so skipping a pair discards no input.
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
    k0opt.strata_breaks = [&pr, &io, water](double x) {
        std::vector<double> br;
        if (water) br.push_back(water_table_at(pr, x, io.config));
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
        // io.numeric overrides either of them (0 = keep the derived value), so the same problem
        // can be re-run at other numerics to show whether its answer depends on them.
        const int steps_default =
            (has_hardening || has_softsoil)
                ? 40
                : (nonlinear_soil ? 20 : ((has_interfaces || has_embedded) ? 5 : 1));
        const double tol_default =
            (has_hardening || has_softsoil)
                ? 1e-2
                : (nonlinear_soil ? 1e-6 : ((has_interfaces || has_embedded) ? 1e-8 : 1e-10));
        const int steps = io.numeric.steps > 0 ? io.numeric.steps : steps_default;
        const double tol = io.numeric.tolerance > 0.0 ? io.numeric.tolerance : tol_default;
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
            // A factor of safety from strength reduction with a NON-ASSOCIATED flow rule
            // (psi < phi, and psi = 0 is the usual engineering choice) is mesh-dependent, and
            // the dependence is one-sided: failure localises into a shear band whose width is
            // set by the elements, so refining the mesh narrows the band and LOWERS the
            // computed factor. It is not a defect of this code -- "the result obtained from a
            // phi/c reduction is influenced by the mesh size, element type and convergence
            // tolerances" (Tschuchnigg, Schweiger and Sloan 2015, Computers and Geotechnics,
            // Part I) -- but a number that moves with the mesh must never be handed over as if
            // it did not. Measured on this program for the Griffiths and Lane benchmark:
            // -7.9% across a fourfold refinement (KV-SLP-003).
            // silent-drop-scope: none -- this loop only looks for a material that makes the
            // factor of safety mesh-dependent, so that the run can say so; it builds nothing.
            for (const auto& mm : pr.materials) {
                if (!(mm.phi > 1e-9) || mm.psi >= mm.phi - 1e-9) continue;
                warn(R, "K2D-A005", mm.name,
                     "Strength reduction with a non-associated flow rule (material \"" + mm.name +
                         "\": phi = " + dnum(mm.phi) + " deg, psi = " + dnum(mm.psi) +
                         " deg): the factor of safety depends on the mesh and falls as the mesh "
                         "is refined, because the shear band narrows with the elements. Quote it "
                         "with the mesh it was computed on, and confirm it with a refinement "
                         "study.");
                break;   // one statement per run: the property is the method's, not the material's
            }
            katai::core::SafetyPhase sfin;
            sfin.nonlinear_soil = nonlinear_soil;
            sfin.has_hardening = has_hardening;
            sfin.has_softsoil = has_softsoil;
            sfin.axisymmetric = axi;
            sfin.active = act;
            // The trial solves inside the strength-reduction search read the SAME controls as any
            // other phase. Only an explicit request reaches them (0 = the search's own defaults),
            // so every factor of safety this program has ever reported is unchanged -- but a
            // tolerance study can now actually change the tolerance it claims to be studying.
            sfin.tolerance = io.numeric.tolerance;
            sfin.load_steps = io.numeric.steps;
            sfin.max_iterations = io.numeric.max_iterations;
            // A loose stopping rule does not add scatter to a factor of safety, it adds BIAS,
            // and always the unsafe way: the search reads "this trial converged" as "the slope
            // stands", so a solver allowed to stop early makes it stand at strengths it cannot
            // carry. Measured on the Griffiths and Lane benchmark (KV-NUM-007): +2.0% at 1e-2,
            // +45.6% at 1e-1. Nobody should have to discover that from a manual.
            if (io.numeric.tolerance > 1e-3)
                warn(R, "K2D-A006", "Safety",
                     "This Safety phase is asked for a tolerated error of " +
                         dnum(io.numeric.tolerance) +
                         ", looser than the strength-reduction search's own 1e-3. The factor of "
                         "safety it reports will be too HIGH, not merely less precise: measured "
                         "on the Griffiths and Lane benchmark, +2.0% at 1e-2 and +45.6% at 1e-1.");
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
                cm.total_stress = Mt.drainage == model::Drainage::UndrainedC;
            }
            cin.flow_edges = flow_edges_from(pr);
            cin.have_flow_bcs = any_flow_bc_declared(pr);
            bool well_in_phase = false;
            cin.drain_nodes = phase_drain_nodes(pr, mesh, io, well_in_phase);
            warn_no_flow_barrier(pr, R);
            if (well_in_phase)
                warn(R, "K2D-A009", "wells",
                     "A well is active in this consolidation phase, but a consolidation analysis "
                     "solves the EXCESS pore pressure and takes no prescribed discharge: the "
                     "well's pumping is not applied here. PLAXIS offers wells for groundwater "
                     "flow and fully coupled analyses. Drains, which set the excess pore pressure "
                     "to zero, ARE applied.");
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
                fm.total_stress = Mt.drainage == model::Drainage::UndrainedC;
            }
            fin.flow_edges = flow_edges_from(pr);
            fin.have_flow_bcs = any_flow_bc_declared(pr);
            bool well_in_fc = false;
            fin.drain_nodes = phase_drain_nodes(pr, mesh, io, well_in_fc);
            warn_no_flow_barrier(pr, R);
            if (well_in_fc)
                warn(R, "K2D-A009", "wells",
                     "A well is active in this fully-coupled phase. This build's coupled solver "
                     "takes drainage boundaries and drains, but not a prescribed well discharge, "
                     "so the well's pumping is not applied here. Run the dewatering as a "
                     "groundwater-flow calculation, or model it with a drain at the target head.");
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
                // In flow, "non-porous" means one thing: no water moves through here and none is
                // stored. An Undrained (C) cluster is in the same position -- PLAXIS does not even
                // let its permeability be entered ("the input field for permeabilities are greyed
                // out when the drainage type is either Non-porous or Undrained C") -- so it takes
                // the same impermeable-barrier treatment rather than a permeability nobody set.
                fm.nonporous = mat_total_stress[mi] != 0;
            }
            tin.flow_edges = flow_edges_from(pr);
            tin.active = act;
            // Wells and drains in a TRANSIENT flow phase: refused by name rather than dropped.
            // This solver takes a prescribed-head boundary set and no source term, so a well's
            // discharge has nowhere to enter and a NORMAL drain's one-sided rule ("pore pressures
            // lower than the equivalent head are not affected") cannot be re-evaluated as the
            // head moves through the time steps -- applying it as a plain fixed head would let a
            // drain FEED water into ground that is drier than it, which is the opposite of what
            // a drain does. PLAXIS applies both here; this build does not yet, and says so.
            for (std::size_t hi = 0; hi < pr.hydros.size(); ++hi) {
                if (io.config && !io.config->active_hydro(hi)) continue;
                const model::HydroLine& H = pr.hydros[hi];
                R.message = std::string(model::hydro_kind_name(H.kind)) + " \"" + H.name +
                            "\" is active in a transient flow phase, which this build cannot "
                            "solve with it: the transient solver takes prescribed heads and no "
                            "discharge, and a drain's one-sided rule is not re-evaluated as the "
                            "head moves. Run the dewatering as a steady groundwater-flow "
                            "calculation (both are exact there), or deactivate it in this phase.";
                return R;
            }
            tin.fallback_head.resize(mesh.node_count);
            for (int n = 0; n < mesh.node_count; ++n)
                tin.fallback_head[n] = pr.has_water ? water_table_at(pr, mesh.x[n], io.config) : mesh.y[n];
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
            din.water_table_y = [&pr, &io](double x) { return water_table_at(pr, x, io.config); };
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
        stin.max_iterations = io.numeric.max_iterations;
        // SumMstage: a partial stage is a construction step only where there IS a stage, so the
        // fraction is read on chained phases and left at 1 on the initial one (the validator
        // refuses it there rather than letting a scaled gravity look like a partial excavation).
        stin.stage_fraction = (io.chained && io.config) ? io.config->sum_mstage : 1.0;
        if (stin.stage_fraction != 1.0)
            note(R, "K2D-A007", io.config->name,
                 "This phase applies only " + dnum(100.0 * stin.stage_fraction) +
                     "% of its staged change (SumMstage = " + dnum(stin.stage_fraction) +
                     "), so the configuration it describes is NOT reached: the remainder is "
                     "still carried by the soil. Results belong to the partial stage.");
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
            R.pore[n] = kGammaWater * std::fmax(0.0, water_table_at(pr, mesh.x[n], io.config) - mesh.y[n]);
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
                                      const std::function<bool()>& cancelled,
                                      const NumericalControls& numeric) {
    std::vector<SolveResult> out;
    std::vector<katai::core::GaussState> committed;
    const int total = 1 + (int)pr.phases.size();
    const auto stop = [&cancelled] { return cancelled && cancelled(); };

    // The caller's controls win where they are set, then the phase's own from the file, then
    // (still zero) the material-class default inside the phase solve. The order is not arbitrary:
    // the argument exists so that a FILE can be re-run at other numerics, which it could not do
    // if the file always won.
    const auto controls_for = [&numeric](const model::Phase& ph) {
        NumericalControls n = numeric;
        if (!(n.tolerance > 0.0)) n.tolerance = ph.tolerance;
        if (n.steps <= 0) n.steps = ph.load_steps;
        if (n.max_iterations <= 0) n.max_iterations = ph.max_iterations;
        return n;
    };

    PhaseIO io0;
    io0.config = &pr.initial;
    io0.numeric = controls_for(pr.initial);
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
        io.numeric = controls_for(ph);
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
