// The analysis-driver bodies, compiled ONCE (section 5.2). As header-inline code every
// driver-including test paid ~220 s of codegen for these functions and the engine
// machinery they pull in; measured in docs/validation/performance-baseline.md.
// Batch 2 moved the ENGINE INCLUDE SET here too: the header now carries only the
// declaration closure, so consumers stop parsing the twenty-odd headers below.
#include <katai/jobs/driver.hpp>

#include <array>
#include <cmath>
#include <cstdio>
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
#include <katai/analysis/excess_pore_force.hpp>    // the inherited excess pore pressure in a baseline
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
// `config` is the phase being solved: a phase may carry its own phreatic polyline (per-phase
// water conditions), which is how staged dewatering is expressed. Null config or
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
// read by the time-dependent flow solvers (its nodes drain: excess pore pressure zero); a well
// prescribes a discharge, which those solvers do not take, so an active well
// in such a phase is reported rather than silently ignored.
// Said where a consolidation / fully-coupled phase meets a flow barrier it cannot read. A split
// seam -- a wall with interfaces, or an interface on mesh edges -- IS read: fully permeable ties
// its two sides to one pore equation, impermeable keeps two, semi-permeable is refused by the
// phase. Measured with a surcharge on one side of the line (KV-STR-007's block and test):
// impermeable interface 12.62 kPa across the seam in consolidation and 1.91 kPa fully coupled,
// a wall with interfaces 12.57 / 1.71 kPa, their permeable twins 0. What is NOT read is a
// barrier on a line the mesh was not split along -- a plate without interfaces, or a wall or
// interface that fell back to bonded (K2D-G009): an impermeable plate there gives the permeable
// plate's pore field bit for bit, in both phases. Until 2026-09 this warning fired for every plate or interface
// in such a phase, whatever its barrier setting and whether or not its seam was read. A warning
// and not a refusal, because the deformation answer is still the model's; it is the flow half
// that is more permeable than drawn.
static void warn_no_flow_barrier(const model::Project& pr, const std::vector<char>& unread,
                                 SolveResult& R) {
    for (size_t si = 0; si < pr.structs.size() && si < unread.size(); ++si)
        if (unread[si]) {
            const auto& st = pr.structs[si];
            warn(R, "K2D-A010", st.name,
                 "\"" + st.name + "\" is declared a flow barrier, but the mesh is not split along "
                 "its line -- a plate without interfaces, or a wall or interface that could not be "
                 "split (K2D-G009) -- so this consolidation / fully-coupled phase has no seam to "
                 "hold the water at: it crosses the line as if the soil were continuous, and a "
                 "cut-off holds back less head than it would in the ground. A wall with "
                 "interfaces, or an interface along mesh edges, is split and its cross "
                 "permeability is read. The steady groundwater-flow calculation splits the mesh "
                 "along any barrier line on mesh edges, a plate without interfaces included, and "
                 "reads it there.");
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

static katai::core::ConsolidationStop to_core_consol_stop(model::ConsolStop c) {
    using M = model::ConsolStop; using C = katai::core::ConsolidationStop;
    switch (c) {
        case M::MinExcessPore:         return C::MinExcessPore;
        case M::DegreeOfConsolidation: return C::DegreeOfConsolidation;
        case M::TimeInterval:          return C::TimeInterval;
    }
    return C::TimeInterval;
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
    p.sig_ci = m.sig_ci;
    p.mi = m.mi;
    p.gsi = m.gsi;
    p.hb_D = m.hb_D;
    p.sig_psi = m.sig_psi;
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
    // (the selected model governs), resolved by name in the constitutive catalogue --
    // parameter wiring, K0^NC memory, cap calibration and the Undrained (A/B) machinery all live
    // with the models in katai/materials/registry.hpp now. `mats` (LinearElastic) is kept only as
    // a fallback; the real solve uses `models`. nonlinear_soil drives load-stepping / solver choice.
    std::vector<katai::core::LinearElastic> mats;
    std::vector<katai::core::MaterialModel> models;
    std::vector<double> gamma;      // unsaturated (above water table)
    std::vector<double> gamma_sat;  // saturated (below water table)
    bool nonlinear_soil = false;
    bool has_hardening = false;  // hardening family present (needs a 1% tolerated error)
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
        // "Ignore undrained behaviour" (`ignoreund`): this phase treats Undrained (A)/(B) soil as
        // drained -- no excess pore pressure is generated and the water adds no stiffness. The
        // strength parameters stay exactly as the material declares them, so an Undrained (B) soil
        // keeps its c_u with phi = 0, and the excess pore pressure generated in earlier phases
        // STAYS: the phase has no mechanism that could dissipate it. The flag is the material's
        // own (MaterialModel::ignore_undrained) rather than `undrained = false`, which deleted that
        // pressure and let the phase's drained volume change return in the next undrained phase as
        // a pore pressure nobody had generated.
        // Reported per material, because a phase that quietly stops being undrained
        // is the difference between a short-term and a long-term answer.
        if (io.config && io.config->ignore_undrained && models.back().undrained) {
            models.back().ignore_undrained = true;
            note(R, "K2D-A008", m.name,
                 "This phase ignores undrained behaviour, so material \"" + m.name +
                     "\" is solved DRAINED: no new excess pore pressure is generated in it, and "
                     "the excess pore pressure generated in earlier phases is kept as it is. Its "
                     "strength parameters are unchanged. This is a long-term (or state-setting) "
                     "answer, not a short-term one.");
        }
        // The Rankine tension cut-off is consumed by the Mohr-Coulomb return only: the hardening
        // and soft-soil integrators do not carry the extra planes (material_model.hpp). This
        // schema applies a tension cut-off to these models by DEFAULT, so a silent omission here
        // is a systematic error in the unsafe direction -- the soil takes tension it should not.
        // CLOSED 2026-08-13: these models now read the cut-off (the Rankine condition
        // sigma_i - sigma_t <= 0 on each principal stress, tension positive, applied to the
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
        // stiffness: the ratio E0/Eur or G0/Gur is limited to 20 in the HS small model, although
        // Alpan suggests the ratio can exceed 10 for very soft clays. The cap is applied in the
        // engine, so the run is the model's as defined rather than a stiffer one -- and it is
        // said out loud, because a G0 that is quietly reduced describes a different soil from
        // the one the file asks for.
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
        // NonPorous: a non-porous material holds NO water -- below the water table its
        // total weight is still gamma_unsat (giving a concrete block gamma_sat would silently
        // saturate it). This one line fixes the gravity_phreatic + gravity_from_head +
        // consolidation dF paths alike.
        gamma_sat.push_back(m.drainage == model::Drainage::NonPorous ? m.gamma_unsat : m.gamma_sat);
    }
    // NonPorous material flags (the shared source for element masks + K0 seeding + flow).
    // `mat_total_stress` is the wider set: materials whose equilibrium is stated in TOTAL
    // stress and which therefore carry no pore pressure -- Non-porous because it holds no
    // water, Undrained (C) because it declines to separate the water from the skeleton
    // (all its pore pressures are zero). Everything downstream that
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
    // The plate Mp/Np plastic hinge is NOW IN THE CORE (the M-N interaction diamond;
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
    // Not const: a material-factored design approach divides the cohesion gradient as well (below).
    std::vector<katai::core::MaterialProfile> profiles = build_profiles(pr);

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
        // JAKY'S FORMULA NEEDS A FRICTION ANGLE. Hoek-Brown has none, so 1 - sin(phi') reads the
        // schema's unused default and lands on a number by accident. It is not refused: the
        // value it lands on, K0 = 1.00 at phi' = 0, is the lithostatic state a competent rock
        // mass is usually assumed to be in, so the default is defensible -- but it must not be
        // arrived at silently, because the user could equally have meant the 0.3 a jointed,
        // relaxed mass shows. This is the one place in the run that can say so.
        if (m.model == model::SoilModel::HoekBrown && m.k0_auto)
            warn(R, "K2D-A017", m.name,
                 "\"" + m.name + "\" uses the Hoek-Brown model, which has no phi', so the "
                 "automatic K0 = 1 - sin(phi') falls back on this material's unused friction box "
                 "and gives K0 = " + dnum(k0) +
                     ". That is not a measurement of anything. A competent rock mass is often "
                     "taken as lithostatic (K0 = 1), but a jointed or relaxed one can be far "
                     "lower -- switch the automatic K0 off and enter the value your site "
                     "investigation supports.");
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
    // Coulomb interface (plate + interface). K0 seeding keeps the wished-in-place wall in
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
        int strength_mat = 0;    // the material whose strength the interfaces take (resolved index)
        size_t si = 0;           // index of the drawn structure in the project
        bool active = true;      // installed in this phase; inactive = the seam is tied (continuous soil)
        std::string name;        // drawn element name (for the force-diagram output)
        int flow_barrier = 0;    // 0 permeable / 1 impermeable / 2 semi-permeable (model vocabulary)
        // Which sides carry an interface, in the builder's frame: "right" is the side the unit
        // normal (nx, ny) points into, "left" the twin side. A side without one is bonded.
        bool iface_right = true, iface_left = true;
    };
    std::vector<WallSpec> walls;
    std::vector<char> plate_is_wall(pr.structs.size(), 0);
    constexpr bool kEnableEmbeddedWall = true;
    if (kEnableEmbeddedWall) {
        const int order = mesh.nodes_per_element;
        for (size_t si = 0; si < pr.structs.size(); ++si) {
            const auto& s = pr.structs[si];
            if (s.kind != model::StructKind::Plate || !(s.iface_pos || s.iface_neg)) continue;
            // An embedded wall splits the mesh in EVERY phase, active or not, so that every phase
            // has the same nodes (the carried datum and the displacement fields are node by node).
            // In a phase where it is inactive nothing is built on the seam and its two sides are
            // tied into one unknown each: the ground is continuous there, as if no wall were drawn.
            const double dx = s.x2 - s.x1, dy = s.y2 - s.y1;
            const double len = std::hypot(dx, dy);
            // silent-drop-ok: a zero-length structural line is an ERROR at structs[i].x2 in the
            // input contract, so it cannot reach a run; the guard is arithmetic self-defence.
            if (len < 1e-9) continue;
            WallSpec w;
            w.si = si;
            w.active = !io.config || io.config->active_struct(si);
            w.name = s.name;
            w.order = order;
            w.flow_barrier = s.flow_barrier;
            // Toe = the DEEPER endpoint (min y; horizontal -> min x); the plate hangs from it (shared node).
            const bool a_is_toe = std::fabs(dy) > 1e-9 ? (s.y1 <= s.y2) : (s.x1 <= s.x2);
            const double tx = a_is_toe ? s.x1 : s.x2, ty = a_is_toe ? s.y1 : s.y2;   // toe end
            const double ux = a_is_toe ? s.x2 : s.x1, uy = a_is_toe ? s.y2 : s.y1;   // top end
            w.toe_x = tx; w.toe_y = ty;
            w.nx = (uy - ty) / len; w.ny = -(ux - tx) / len;   // unit normal of the wall line
            // The drawn line's POSITIVE side is the one on the right walking from its first point
            // to its second, normal (dy, -dx)/len -- the side the Studio draws the "+" marks on.
            // The builder's right side is +(nx, ny), the normal of the toe->top direction, so the
            // two frames agree when the first point is the toe and are mirrored otherwise.
            w.iface_right = a_is_toe ? s.iface_pos : s.iface_neg;
            w.iface_left = a_is_toe ? s.iface_neg : s.iface_pos;
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
            // Interface strength/stiffness from the adjacent soil + Rinter (strength reduced as
            // c_i = Rinter c, tan phi_i = Rinter tan phi).
            const double mx = 0.5 * (s.x1 + s.x2), my = 0.5 * (s.y1 + s.y2);
            int smat = s.iface_material >= 0 ? s.iface_material : material_at(mx + 1e-3 * w.nx, my + 1e-3 * w.ny);
            if (smat < 0) smat = material_at(mx - 1e-3 * w.nx, my - 1e-3 * w.ny);
            w.soil_mat = smat;
            const auto& sm = pr.materials[smat >= 0 ? smat : 0];
            const double Rinter = sm.rinter_rigid ? 1.0 : sm.Rinter;
            const double G = sm.E / (2.0 * (1.0 + sm.nu));
            const double avg = len / std::max<size_t>(1, w.seam.size() / (order == 15 ? 4 : 2));
            katai::core::iface::interface_stiffness(Rinter, G, avg, 0.1, w.ip.kn, w.ip.ks);
            if (s.iface_kn > 0.0) w.ip.kn = s.iface_kn;   // entered: the joint's own stiffness
            if (s.iface_ks > 0.0) w.ip.ks = s.iface_ks;
            // A JOINT AGAINST ROCK CANNOT BORROW A STRENGTH THE ROCK DOES NOT HAVE. The two
            // lines below read the material's c' and phi', and a Hoek-Brown material's are the
            // schema DEFAULTS -- 1 kPa and 30 degrees -- because that model never reads them.
            // The joint would then be given a leftover Mohr-Coulomb strength that has nothing
            // to do with the rock it is cut into, and nothing would say so. The schema already
            // carries the remedy: `iface_material` points the interface at a DIFFERENT soil
            // material for its strength, so a rock joint is modelled by naming a Mohr-Coulomb
            // material whose c'/phi' IS the joint strength intended, the interface being
            // Mohr-Coulomb whatever the surrounding model.
            if (sm.model == model::SoilModel::HoekBrown) {
                refuse(R, "K2D-G014", s.name,
                    "The interface on \"" + s.name +
                    "\" takes its strength from the adjacent material \"" + sm.name +
                    "\", which uses the Hoek-Brown model -- and that model has no c' or phi' for "
                    "a joint to be reduced from: the interface would silently be given this "
                    "material's unused c'/phi' boxes instead. An interface is Mohr-Coulomb "
                    "whatever the surrounding model, so set the interface's soil material "
                    "(iface_material) to a Mohr-Coulomb material whose c' and phi' are the joint "
                    "strength you intend -- for a rock joint, typically the discontinuity's own "
                    "friction and cohesion, not the rock mass's.");
                return R;
            }
            // THE JOINT TAKES THE STRENGTH ITS SOIL IS SOLVED WITH, read from the constitutive model
            // rather than from the schema boxes. The two differ exactly where the model overrides a
            // box: Undrained (B) and (C) solve the soil as a Tresca material, c = su with phi = 0,
            // and the input contract says so ("phi is forced to 0 and the entered value is ignored
            // for strength"). Read from the box, the joint beside such a clay kept that ignored
            // angle -- measured on the sliding block: 59.72 kN/m of capacity with the box left at
            // 26.6 deg, against 9.90 with it at 0 and B su = 10 for the Tresca joint.
            w.strength_mat = smat >= 0 ? smat : 0;
            const katai::core::MaterialModel& smodel = models[(size_t)w.strength_mat];
            w.ip.c_i = Rinter * smodel.cohesion;
            w.ip.phi_i = std::atan(Rinter * std::tan(smodel.friction_angle));
            // Interface tensile strength = R * sigma_t (reduced like c) -- only while the
            // material's tension cutoff is on; otherwise 0 (the interface default: carries no
            // tension -- the safe side). Previously the material's sigma_t never reached the
            // interface at all (an audit finding).
            w.ip.sigma_t = smodel.tension_cutoff ? Rinter * smodel.tensile_strength : 0.0;
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
        int strength_mat = 0;   // the material whose strength the joint takes (resolved index)
        size_t si = 0;          // index of the drawn structure in the project
        bool active = true;     // installed in this phase; inactive = the seam is tied (continuous soil)
        std::string name;
        int flow_barrier = 0;   // 0 permeable / 1 impermeable / 2 semi-permeable (model vocabulary)
    };
    std::vector<IfaceSpec> soil_ifaces;
    std::vector<char> iface_split(pr.structs.size(), 0);   // interfaces whose line the mesh was split along
    std::vector<char> bc_released;   // nodes the domain boundary must NOT fix (filled by the seam pass)
    {
        std::vector<char> is_bnode(mesh.node_count, 0);
        for (int n : mesh.boundary_nodes) is_bnode[n] = 1;
        const int order = mesh.nodes_per_element;
        std::vector<std::pair<int, int>> boundary_seam;   // (orig, dup) pairs sitting ON the boundary
        for (size_t si = 0; si < pr.structs.size(); ++si) {
            const auto& s = pr.structs[si];
            if (s.kind != model::StructKind::Interface) continue;
            // Split in every phase, like a wall; inactive = its two sides tied (see the wall loop).
            const double dx = s.x2 - s.x1, dy = s.y2 - s.y1;
            const double len = std::hypot(dx, dy);
            // silent-drop-ok: a zero-length structural line is an ERROR at structs[i].x2 in the
            // input contract, so it cannot reach a run; the guard is arithmetic self-defence.
            if (len < 1e-9) continue;
            IfaceSpec sp;
            sp.si = si;
            sp.active = !io.config || io.config->active_struct(si);
            sp.name = s.name.empty() ? "Interface" : s.name;
            sp.flow_barrier = s.flow_barrier;
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
            // Interface strength/stiffness from the adjacent soil + Rinter (strength reduced as
            // c_i = Rinter c, tan phi_i = Rinter tan phi).
            const double mx = 0.5 * (s.x1 + s.x2), my = 0.5 * (s.y1 + s.y2);
            int smat = s.iface_material >= 0 ? s.iface_material : material_at(mx + 1e-3 * sp.nx, my + 1e-3 * sp.ny);
            if (smat < 0) smat = material_at(mx - 1e-3 * sp.nx, my - 1e-3 * sp.ny);
            sp.soil_mat = smat;
            const auto& sm = pr.materials[smat >= 0 ? smat : 0];
            const double Rinter = sm.rinter_rigid ? 1.0 : sm.Rinter;
            const double G = sm.E / (2.0 * (1.0 + sm.nu));
            const int nedges = std::max<int>(1, ((int)sp.seam.size() - 1) / per_edge);
            katai::core::iface::interface_stiffness(Rinter, G, len / nedges, 0.1, sp.ip.kn, sp.ip.ks);
            if (s.iface_kn > 0.0) sp.ip.kn = s.iface_kn;
            if (s.iface_ks > 0.0) sp.ip.ks = s.iface_ks;
            // A JOINT AGAINST ROCK CANNOT BORROW A STRENGTH THE ROCK DOES NOT HAVE. The two
            // lines below read the material's c' and phi', and a Hoek-Brown material's are the
            // schema DEFAULTS -- 1 kPa and 30 degrees -- because that model never reads them.
            // The joint would then be given a leftover Mohr-Coulomb strength that has nothing
            // to do with the rock it is cut into, and nothing would say so. The schema already
            // carries the remedy: `iface_material` points the interface at a DIFFERENT soil
            // material for its strength, so a rock joint is modelled by naming a Mohr-Coulomb
            // material whose c'/phi' IS the joint strength intended, the interface being
            // Mohr-Coulomb whatever the surrounding model.
            if (sm.model == model::SoilModel::HoekBrown) {
                refuse(R, "K2D-G014", s.name,
                    "The interface on \"" + s.name +
                    "\" takes its strength from the adjacent material \"" + sm.name +
                    "\", which uses the Hoek-Brown model -- and that model has no c' or phi' for "
                    "a joint to be reduced from: the interface would silently be given this "
                    "material's unused c'/phi' boxes instead. An interface is Mohr-Coulomb "
                    "whatever the surrounding model, so set the interface's soil material "
                    "(iface_material) to a Mohr-Coulomb material whose c' and phi' are the joint "
                    "strength you intend -- for a rock joint, typically the discontinuity's own "
                    "friction and cohesion, not the rock mass's.");
                return R;
            }
            // The strength the soil is SOLVED with, as for the wall above (Undrained (B)/(C) is a
            // Tresca material, and its joint is a Tresca joint).
            sp.strength_mat = smat >= 0 ? smat : 0;
            const katai::core::MaterialModel& smodel = models[(size_t)sp.strength_mat];
            sp.ip.c_i = Rinter * smodel.cohesion;
            sp.ip.phi_i = std::atan(Rinter * std::tan(smodel.friction_angle));
            sp.ip.sigma_t = smodel.tension_cutoff ? Rinter * smodel.tensile_strength : 0.0;
            iface_split[si] = 1;
            soil_ifaces.push_back(std::move(sp));
        }
        // Who holds the support on a seam that touches the domain boundary. The two sides of a seam
        // sit at the SAME coordinates and apply_boundary_conditions matches by coordinate, so this
        // cannot be left to it: measured 2026-08-10, an interface drawn ALONG a fixed boundary got
        // both of its sides fixed and became inert -- the mesh split, the joint assembled, and the
        // block welded itself to its own base (the sliding-block benchmark read 5.4e6 kN/m instead
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
    // structures vector belongs to which drawn line, so force diagrams (M/Q/N)
    // and interface results (tau/sigma_n/slip) can be produced per line.
    using katai::core::DiagSpec;
    using katai::core::IfaceDiag;
    std::vector<DiagSpec> diag_specs;
    std::vector<IfaceDiag> iface_diags;
    katai::core::Structures structures;
    // Which material each interface element takes its strength from, parallel to
    // structures.interfaces / interfaces5: a phase that factors ground strength has to factor the
    // joints by the same rule as the soil they came from, after they have been built.
    std::vector<int> iface_strength_mat, iface5_strength_mat;
    // What each drawn structure built, so the phase after can match its state by structure rather
    // than by position (StructCarryRecord, structural_carry.hpp). A snapshot of every Structures
    // vector and of the DOF count before a structure is built, and the record of what it added.
    struct StructSnap { size_t pl, pl5, an, gg, if3, if5, eb, skin; int extra; };
    const auto snap = [&]() {
        size_t skin = 0;
        for (const auto& eb : structures.embedded_beams) skin += eb.skin.size();
        return StructSnap{structures.plates.size(), structures.plates5.size(),
                          structures.anchors.size(), structures.geogrids.size(),
                          structures.interfaces.size(), structures.interfaces5.size(),
                          structures.embedded_beams.size(), skin, dofs.total_dofs()};
    };
    std::vector<katai::core::StructCarryRecord> struct_records;
    const auto record = [&](size_t si, const StructSnap& a) {
        const StructSnap b = snap();
        katai::core::StructCarryRecord r;
        r.si = (int)si;
        r.plates = {a.pl, b.pl}; r.plates5 = {a.pl5, b.pl5}; r.anchors = {a.an, b.an};
        r.geogrids = {a.gg, b.gg}; r.interfaces = {a.if3, b.if3}; r.interfaces5 = {a.if5, b.if5};
        r.embedded = {a.eb, b.eb}; r.skin = {a.skin, b.skin}; r.extra_dof = {a.extra, b.extra};
        const bool built = b.pl > a.pl || b.pl5 > a.pl5 || b.an > a.an || b.gg > a.gg ||
                           b.if3 > a.if3 || b.if5 > a.if5 || b.eb > a.eb || b.extra > a.extra;
        if (built) struct_records.push_back(r);
    };
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Plate || plate_is_wall[si]) continue;   // walls built below
        if (!struct_on(si)) continue;   // not installed in this phase
        const StructSnap snap0 = snap();
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
        record(si, snap0);
    }

    // Geogrids: tension-only axial membranes. Like plates, each line is a conforming node chain;
    // translational DOFs are shared with the soil and there is no rotation DOF (no bending).
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Geogrid || !struct_on(si)) continue;
        const StructSnap snap0 = snap();
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
        record(si, snap0);
    }

    // Anchors: single axial spring (node-to-node when both ends sit in the soil; fixed-end when one
    // end is drawn outside the soil -> that end is a fixed far anchor point). EA/Fmax are divided by
    // the out-of-plane spacing to get the per-metre (plane-strain) values.
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
    std::vector<std::array<double, 4>> anchor_ends;   // drawn end points, per anchor: a (x,y), b (x,y)
    for (size_t si = 0; si < pr.structs.size(); ++si) {
        const auto& s = pr.structs[si];
        if (s.kind != model::StructKind::Anchor || !struct_on(si)) continue;
        const StructSnap snap0 = snap();
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
        {
            const double sx = in1 ? s.x1 : s.x2, sy = in1 ? s.y1 : s.y2;
            anchor_ends.push_back(an.node_b >= 0 ? std::array<double, 4>{s.x1, s.y1, s.x2, s.y2}
                                                 : std::array<double, 4>{sx, sy, 1e300, 1e300});
        }
        diag_specs.push_back({1, s.name, structures.anchors.size() - 1, structures.anchors.size(), Ls});
        record(si, snap0);
    }

    // Embedded walls: build the plate + two interfaces on the split mesh and seed each interface's
    // initial normal stress sigma_n0 = K0*sigma'_v from the layered field (wished-in-place K0).
    // The twins an inactive wall or interface ties back together, for the inactive-node pass below.
    std::vector<std::pair<int, int>> tied_nodes;
    const auto tie_twins = [&](int twin, int original) {
        for (int c = 0; c < 2; ++c) dofs.tie(dofs.global_dof(twin, c), dofs.global_dof(original, c));
        tied_nodes.push_back({twin, original});
    };
    // WHAT IS DRAWN ON A WALL HANGS ON THE WALL. An embedded wall's translations are its own
    // DOFs, joined to the soil only through the interfaces, and the mesh node at a point of the
    // wall line is the SOIL beside it. Anchors, geogrids, plates and point loads find their node
    // by position, so without this an anchor drawn from the wall tensioned the retained ground
    // against itself and the wall never felt it -- measured on an anchored sheet pile: 1.8 kN
    // left of a 20 kN prestress, the free-length wall held by the interfaces alone. Every node of
    // an active wall's line (either twin) is therefore mapped to the wall's own translation.
    std::vector<std::array<int, 2>> wall_dof_at((size_t)mesh.node_count, std::array<int, 2>{-1, -1});
    const size_t plain_plate_count = structures.plates.size();   // the walls' plates follow these
    const auto attach_to_wall = [&wall_dof_at](const std::vector<int>& nr, const std::vector<int>& nl,
                                               const std::vector<int>& dx, const std::vector<int>& dy) {
        for (size_t i = 0; i < nr.size(); ++i) {
            wall_dof_at[(size_t)nr[i]] = {dx[i], dy[i]};
            wall_dof_at[(size_t)nl[i]] = {dx[i], dy[i]};
        }
    };
    for (const auto& w : walls) {
        // silent-drop-ok: an inactive wall is not dropped -- the phase asked for it to be absent, and
        // its seam is tied so that the ground across the line is continuous, exactly as without it.
        if (!w.active) {
            for (const auto& sp : w.seam) tie_twins(sp.left, sp.right);
            continue;
        }
        const StructSnap snap0 = snap();
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
            katai::core::WallBuild5 wb = katai::core::build_embedded_wall5(mesh, w.seam, toe, dofs, w.pp, w.ip,
                                                                           w.iface_right, w.iface_left);
            attach_to_wall(wb.node_r, wb.node_l, wb.dof_x, wb.dof_y);
            const auto ncp = katai::core::iface::nc_points5();
            for (auto& ie : wb.interfaces)
                for (int q = 0; q < 5; ++q) {
                    const int nd = ie.soil_nodes[ncp[q].node];
                    ie.sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
            const size_t p0 = structures.plates5.size();
            for (const auto& pe : wb.plates) structures.plates5.push_back(pe);
            const size_t if0 = structures.interfaces5.size();
            for (const auto& ie : wb.interfaces) structures.interfaces5.push_back(ie);
            iface5_strength_mat.resize(structures.interfaces5.size(), w.strength_mat);
            if (structures.plates5.size() > p0)
                diag_specs.push_back({5, w.name, p0, structures.plates5.size()});   // tri15 wall N/Q/M
            if (structures.interfaces5.size() > if0)
                iface_diags.push_back({w.name, 15, if0, structures.interfaces5.size()});
        } else {
            katai::core::WallBuild wb = katai::core::build_embedded_wall(mesh, w.seam, toe, dofs, w.pp, w.ip,
                                                                         w.iface_right, w.iface_left);
            attach_to_wall(wb.node_r, wb.node_l, wb.dof_x, wb.dof_y);
            const auto ncp = katai::core::iface::nc_points();
            for (auto& ie : wb.interfaces)
                for (int q = 0; q < 3; ++q) {
                    const int nd = ie.soil_nodes[ncp[q].node];
                    ie.sigma_n0[q] = fac * eff_sigma_v(mesh.x[nd], mesh.y[nd]);
                }
            const size_t p0 = structures.plates.size();
            for (const auto& pe : wb.plates) structures.plates.push_back(pe);
            const size_t if0 = structures.interfaces.size();
            for (const auto& ie : wb.interfaces) structures.interfaces.push_back(ie);
            iface_strength_mat.resize(structures.interfaces.size(), w.strength_mat);
            if (structures.plates.size() > p0)
                diag_specs.push_back({0, w.name, p0, structures.plates.size()});
            if (structures.interfaces.size() > if0)
                iface_diags.push_back({w.name, 6, if0, structures.interfaces.size()});
        }
        record(w.si, snap0);
    }

    // Standalone interfaces: build the soil-soil Coulomb joint on the split mesh and seed each NC point's
    // initial normal stress sigma_n0 = (K0*nx^2 + ny^2)*sigma'_v (orientation-aware; vertical -> K0*sigma'_v,
    // horizontal -> sigma'_v) so the wished-in-place interface starts in geostatic equilibrium.
    for (auto& sp : soil_ifaces) {
        // silent-drop-ok: an inactive interface is not dropped -- its seam is tied, so the ground
        // across it is continuous, which is what an interface the phase does not install means.
        if (!sp.active) {
            for (const auto& sg : sp.seam) tie_twins(sg.dup, sg.orig);
            continue;
        }
        const StructSnap snap0 = snap();
        const katai::core::InterfaceRange rng =
            katai::core::build_soil_interface(mesh, sp.seam, dofs, sp.ip, mesh.nodes_per_element, structures);
        iface_strength_mat.resize(structures.interfaces.size(), sp.strength_mat);
        iface5_strength_mat.resize(structures.interfaces5.size(), sp.strength_mat);
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
        record(sp.si, snap0);
    }
    // Rebind the structures built above (they came first, by position) onto the wall DOFs.
    {
        const auto on_wall = [&wall_dof_at](int n) {
            return n >= 0 && n < (int)wall_dof_at.size() && wall_dof_at[(size_t)n][0] >= 0;
        };
        for (auto& an : structures.anchors) {
            if (on_wall(an.node_a)) { an.end_dof[0] = wall_dof_at[(size_t)an.node_a][0]; an.end_dof[1] = wall_dof_at[(size_t)an.node_a][1]; }
            if (on_wall(an.node_b)) { an.end_dof[2] = wall_dof_at[(size_t)an.node_b][0]; an.end_dof[3] = wall_dof_at[(size_t)an.node_b][1]; }
        }
        for (auto& ge : structures.geogrids)
            for (int k = 0; k < 3; ++k)
                if (on_wall(ge.nodes[k])) {
                    ge.trans_dof[2 * k + 0] = wall_dof_at[(size_t)ge.nodes[k]][0];
                    ge.trans_dof[2 * k + 1] = wall_dof_at[(size_t)ge.nodes[k]][1];
                }
        // A plain plate meeting the wall (a strut, a slab) shares its translation there -- a
        // hinged connection, which is what two plates meeting at a node are everywhere else.
        for (size_t pi = 0; pi < plain_plate_count; ++pi) {
            auto& pe = structures.plates[pi];
            for (int k = 0; k < 3; ++k)
                if (pe.trans_dof[2 * k] < 0 && on_wall(pe.nodes[k])) {
                    pe.trans_dof[2 * k + 0] = wall_dof_at[(size_t)pe.nodes[k]][0];
                    pe.trans_dof[2 * k + 1] = wall_dof_at[(size_t)pe.nodes[k]][1];
                }
        }
    }
    // A JOINT NEEDS GROUND ON BOTH SIDES. Where a phase has excavated the soil on one side of a
    // wall (or of a soil-soil interface), the nodes there are orphaned and fixed further down,
    // and the joint would hold the wall against points fixed in space while pressing it with the
    // earth pressure of the removed soil (its sigma_n0) -- the excavation never unloaded the wall.
    // Measured on an anchored sheet pile with interfaces: the lower anchor ended in COMPRESSION
    // (-2.4 kN) because the pit side still pushed. Such a joint is switched off for the phase; it
    // keeps its place in the arrays, so the phase chain's carried state stays aligned.
    if (!act.empty()) {
        const std::vector<char> na = katai::core::active_nodes(mesh, act);
        const auto live = [&](int n) { return n >= 0 && n < mesh.node_count && na[(size_t)n]; };
        // The structure side is a wall's own DOF (always there), or a mesh node's: the other soil
        // of a soil-soil joint, or the bonded side of a one-sided wall -- which the wall holds.
        std::vector<char> held((size_t)mesh.node_count, 0);
        for (const auto& p : structures.plates)
            for (int d : p.trans_dof) if (d >= 0 && d < 2 * mesh.node_count) held[(size_t)(d / 2)] = 1;
        for (const auto& p : structures.plates5)
            for (int d : p.trans_dof) if (d >= 0 && d < 2 * mesh.node_count) held[(size_t)(d / 2)] = 1;
        const auto side_live = [&](int gdof) {
            return gdof >= 2 * mesh.node_count || live(gdof / 2) || held[(size_t)(gdof / 2)];
        };
        for (auto& ie : structures.interfaces) {
            for (int k = 0; k < 3; ++k)
                if (!live(ie.soil_nodes[k]) || !side_live(ie.struct_dof[2 * k])) ie.active = false;
        }
        for (auto& ie : structures.interfaces5) {
            for (int k = 0; k < 5; ++k)
                if (!live(ie.soil_nodes[k]) || !side_live(ie.struct_dof[2 * k])) ie.active = false;
        }
    }
    const bool has_interfaces = !structures.interfaces.empty() || !structures.interfaces5.empty();

    // Embedded beams (pile rows / nails): a beam line cuts the mesh at any orientation; skin + foot
    // springs couple the beam (on independent DOFs) to the soil AT the beam location (non-conforming,
    // no mesh change). EA/EI and the skin/foot spring stiffnesses come from the embedded-beam material
    // via default interface-stiffness factors (Sluis 2012): ISF = 2.5
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
        const StructSnap snap0 = snap();
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
        // build_embedded_beam: node 0 = toe (foot), the end that is NOT the connection point -- the
        // lower end, and for an exactly horizontal pile the higher-x one, since the connection
        // point takes the lower-x end there. Ordering by y alone left a horizontal pile drawn
        // left to right with its foot AT its connection point, and the far end with no base.
        {
            double hx = 0.0, hy = 0.0;
            embedded_connection_point(s, hx, hy);
            if (bx.front() == hx && by.front() == hy) {   // front is s.(x1,y1) exactly (t = 0)
                std::reverse(bx.begin(), bx.end()); std::reverse(by.begin(), by.end());
            }
        }
        // CONNECTION POINT (`conn`). Hinged -- the default when no structure
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
        // Skin/foot spring stiffnesses from the default interface-stiffness factors (Sluis 2012) --
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
        record(si, snap0);
    }
    const bool has_embedded = !structures.embedded_beams.empty();
    // AN ANCHOR THAT ENDS ON A FREE PILE PULLS THE PILE. A free connection is the grout body of
    // a ground anchor: the beam carries the force into the ground through its skin. The anchor
    // finds its end node by position, like every point attachment, and the mesh node there is
    // the SOIL -- so the free length used to pull the ground beside the grout body, which then
    // carried nothing (measured: 0.6-3.8 kN along a grout body behind a 20 kN tieback, whose
    // free length the model had tied to the soil instead). An end drawn exactly on a free
    // beam's node takes that node's translations. A hinged beam needs nothing: its connection
    // node's translations already are the mesh node's.
    for (size_t ai = 0; ai < structures.anchors.size() && ai < anchor_ends.size(); ++ai) {
        auto& an = structures.anchors[ai];
        for (int e = 0; e < 2; ++e) {
            if (an.end_dof[2 * e] >= 0) continue;   // already on a wall
            const double px = anchor_ends[ai][2 * e], py = anchor_ends[ai][2 * e + 1];
            if (px > 1e299) continue;               // a fixed far end
            for (const auto& eb : structures.embedded_beams) {
                if (eb.connection != katai::core::ebeam::Connection::Free) continue;
                for (size_t k = 0; k < eb.node_x.size(); ++k)
                    if (std::hypot(eb.node_x[k] - px, eb.node_y[k] - py) < 1e-6 &&
                        eb.dof_x[k] >= 0 && eb.dof_y[k] >= 0) {
                        an.end_dof[2 * e] = eb.dof_x[k];
                        an.end_dof[2 * e + 1] = eb.dof_y[k];
                    }
            }
        }
    }

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
        // K2D-A003 IS RETIRED, and the reason is two measurements, not one. It warned that a
        // structural element standing on a driven node does not receive the imposed motion. The
        // SOLVER stopped being that way on 2026-08-13 (internal_forces.hpp, u_at), yet the
        // warning kept firing -- including on KV-STR-005's own run, whose driven geogrid force
        // that case asserts at +0.000%. The REPORT was half true. The plate, geogrid, embedded-
        // beam and interface post-processors read the full displacement vector, whose fixed
        // entries hold the prescribed values -- measured for the plate (a two-span plate whose
        // middle support settles returns M_B = 3EI delta / L^2 with its shear term, to 1.4e-12)
        // and for the geogrid (KV-STR-005). The anchor force skipped fixed DOFs and reported
        // N = 0 where the solve carried EA d / L.
        // With that report aligned to the solver (structural_forces.hpp, anchor_force), nothing
        // the warning described remains, and a code is never reused: it is simply not raised.
    }
    // Staged construction: nodes touched only by passive (excavated / not-yet-filled) elements
    // would be singular -- fix them (their displacement is meaningless this phase).
    // A node a PLATE runs through is NOT orphaned -- the beam's own stiffness holds it. Pinning
    // it would weld the structure to the outside world and the phase would converge on a model
    // in which that plate carries nothing, without saying so. A plain beam verification is built
    // exactly this way: the soil cluster is deactivated so that only the beams remain, supported
    // at their end points.
    if (!act.empty()) {
        std::vector<char> carried(mesh.node_count, 0);
        const auto carry = [&](int n) { if (n >= 0 && n < mesh.node_count) carried[(size_t)n] = 1; };
        for (const auto& p : structures.plates) for (int n : p.nodes) carry(n);
        for (const auto& p : structures.plates5) for (int n : p.nodes) carry(n);
        // ... nor a node whose DOFs a plate moves on: the bonded side of a one-sided wall.
        for (const auto& p : structures.plates)
            for (int d : p.trans_dof) if (d >= 0 && d < 2 * mesh.node_count) carry(d / 2);
        for (const auto& p : structures.plates5)
            for (int d : p.trans_dof) if (d >= 0 && d < 2 * mesh.node_count) carry(d / 2);
        // A twin tied to a node that an active element still holds is not orphaned either: the two
        // are one unknown. Only when neither side has an active element is the pair fixed.
        if (!tied_nodes.empty()) {
            const std::vector<char> na = katai::core::active_nodes(mesh, act);
            for (const auto& [a, b] : tied_nodes)
                if (na[(size_t)a] || na[(size_t)b]) { carry(a); carry(b); }
        }
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
    // Matched by drawn structure, so that a phase which activates or removes a structure still
    // carries every other one (build_structural_carry states the three cases). carry_plan.full_datum
    // is the parent's total displacement in THIS phase's DOF numbering -- the datum every consumer
    // below evaluates its structures at.
    katai::core::StructuralCarry carry_plan;
    if (any_struct_carry && io.prev && io.prev->ok) {
        const StructCarryState& ps = io.prev->struct_state;
        carry_plan = katai::core::build_structural_carry(ps, struct_records, structures, dofs,
                                                         2 * mesh.node_count, carry_init);
        if (carry_plan.carried) carry_src = &ps;
    }
    // A wall or interface that was active in the parent and is inactive here leaves its two sides
    // where the joint let them go -- apart, if it slipped or opened -- and ties them for this phase's
    // increments. The ground is fine with that: it carries its stresses, not a displacement datum. A
    // structure CARRIED on such a node is not: its datum is the node's total displacement, and a tied
    // pair is one unknown with one datum while the two twins have two totals. Refused rather than
    // continued from whichever twin happens to be read.
    if (carry_plan.carried && !tied_nodes.empty()) {
        std::vector<char> offset_twin(mesh.node_count, 0);
        for (const auto& [a, b] : tied_nodes)
            for (int c = 0; c < 2; ++c) {
                const double ua = carry_plan.full_datum[2 * a + c], ub = carry_plan.full_datum[2 * b + c];
                if (std::fabs(ua - ub) > 1e-12 * std::max(1.0, std::max(std::fabs(ua), std::fabs(ub))))
                    offset_twin[(size_t)a] = offset_twin[(size_t)b] = 1;
            }
        const auto on_offset = [&](int n) { return n >= 0 && n < mesh.node_count && offset_twin[(size_t)n]; };
        for (size_t c = 0; c < struct_records.size(); ++c) {
            if (c < carry_plan.record_installed.size() && carry_plan.record_installed[c]) continue;
            const auto& r = struct_records[c];
            bool hit = false;
            for (size_t i = r.plates[0]; i < r.plates[1] && !hit; ++i)
                for (int n : structures.plates[i].nodes) hit = hit || on_offset(n);
            for (size_t i = r.plates5[0]; i < r.plates5[1] && !hit; ++i)
                for (int n : structures.plates5[i].nodes) hit = hit || on_offset(n);
            for (size_t i = r.anchors[0]; i < r.anchors[1] && !hit; ++i)
                hit = on_offset(structures.anchors[i].node_a) || on_offset(structures.anchors[i].node_b);
            for (size_t i = r.geogrids[0]; i < r.geogrids[1] && !hit; ++i)
                for (int n : structures.geogrids[i].nodes) hit = hit || on_offset(n);
            for (size_t i = r.interfaces[0]; i < r.interfaces[1] && !hit; ++i)
                for (int n : structures.interfaces[i].soil_nodes) hit = hit || on_offset(n);
            for (size_t i = r.interfaces5[0]; i < r.interfaces5[1] && !hit; ++i)
                for (int n : structures.interfaces5[i].soil_nodes) hit = hit || on_offset(n);
            for (size_t i = r.embedded[0]; i < r.embedded[1] && !hit; ++i) {
                const auto& eb = structures.embedded_beams[i];
                hit = on_offset(eb.conn_mesh_node);
                for (const auto& sp : eb.skin)
                    for (int k = 0; k < mesh.nodes_per_element && !hit && sp.soil_elem >= 0; ++k)
                        hit = on_offset(mesh.node_of(sp.soil_elem, k));
            }
            if (hit) {
                const std::string nm = r.si >= 0 && (size_t)r.si < pr.structs.size() ? pr.structs[(size_t)r.si].name : "";
                R.message = "Structure \"" + nm + "\" is carried into this phase on the line of a wall or "
                            "interface that is deactivated here, after the joint had let its two sides move "
                            "apart. The two sides are tied for this phase, and a structure standing on them "
                            "cannot be continued from both. Deactivate \"" + nm + "\" in this phase as well, "
                            "or keep the wall or interface active.";
                return R;
            }
        }
    }

    // Axisymmetric scope: soil-only. A WATER TABLE is now carried -- the r-weighted phreatic
    // body force and the pore-pressure load with its hoop term (assembler.hpp), on top of a K0
    // seed that was already effective-stress and already set the hoop -- so the tank, the silo
    // and the circular footing on saturated ground are reachable. What is still missing is
    // refused rather than silently integrated as plane strain.
    if (axi) {
        const bool has_structs =
            !structures.plates.empty() || !structures.plates5.empty() || !structures.anchors.empty() ||
            !structures.geogrids.empty() || !structures.interfaces.empty() ||
            !structures.interfaces5.empty() || !structures.embedded_beams.empty();
        if (has_structs) { R.message = "Structural elements are not supported in axisymmetric mode yet "
                                       "(soil-only). A plate in axisymmetry is a shell with a HOOP "
                                       "membrane force and an anchor is a ring, so they are different "
                                       "elements rather than the same ones integrated differently -- "
                                       "which is why this is refused instead of approximated. Use "
                                       "plane strain for walls/anchors/plates."; return R; }
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

    // NonPorous elements take no pore-pressure LOAD (a non-porous material carries
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
        // Pore pressure + saturation from the steady-state seepage head field (pore pressures
        // generated by a groundwater flow calculation): total-stress gravity with gamma_sat where
        // psi = h - y >= 0, and the pore load interpolated from the SAME nodal head -> recovered
        // stress is EFFECTIVE.
        katai::core::assemble_gravity_from_head(mesh, dofs, gamma, gamma_sat, *flow_head, f);
        katai::core::assemble_pore_load_from_head(mesh, dofs, *flow_head, kGammaWater, f, pore_mask);
    } else if (water) {
        // Total-stress equilibrium: saturated weight below the phreatic surface, moist above.
        const auto wt = [&pr, &io](double x) { return water_table_at(pr, x, io.config); };
        // Hydrostatic pore-pressure load -> recovered stress is EFFECTIVE (Terzaghi; buoyancy
        // gamma' = gamma_sat - gamma_w emerges naturally). docs/references/effective-stress-formulation.md.
        const auto pore = [&pr, &io](double x, double y) {
            return kGammaWater * std::fmax(0.0, water_table_at(pr, x, io.config) - y);
        };
        if (axi) {
            // The same two terms, r-weighted, and the pore load carrying the HOOP component the
            // plane-strain one has no equation for (assembler.hpp). Everything the water table
            // decides is decided by the same lambdas, so the two modes cannot drift apart in
            // where the water is -- only in how the integral is weighted.
            katai::core::assemble_axisym_gravity_phreatic(mesh, dofs, gamma, gamma_sat, wt, f, act);
            katai::core::assemble_axisym_pore_pressure_load(mesh, dofs, pore, f, pore_mask);
        } else {
            katai::core::assemble_gravity_phreatic(mesh, dofs, gamma, gamma_sat, wt, f, act);
            katai::core::assemble_pore_pressure_load(mesh, dofs, pore, f, pore_mask);
        }
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
        // A load on the line of an embedded wall acts on the wall, not on the soil beside it.
        const bool on_wall = wall_dof_at[(size_t)best][0] >= 0;
        const int ex = dofs.equation(on_wall ? wall_dof_at[(size_t)best][0] : dofs.global_dof(best, 0));
        const int ey = dofs.equation(on_wall ? wall_dof_at[(size_t)best][1] : dofs.global_dof(best, 1));
        if (ex >= 0) { f[ex] += L.qx1; f_loads[ex] += L.qx1; }
        if (ey >= 0) { f[ey] += L.qy1; f_loads[ey] += L.qy1; }
    }
    // Distributed (line) loads: a linearly-varying traction (qx1,qy1)->(qx2,qy2) along the
    // segment, assembled as CONSISTENT nodal forces over the mesh edge chain it lies on (the
    // common surcharge). The mesh conforms to the load line (build_mesh adds it as a
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
    // uses characteristic strength (the in-situ state is real, not factored) -- design factors
    // belong to the calculation phase, not the initial one. None -> no-op (bit-identical). The
    // resistance-factored approaches (DA2, TBDY 2018) do not change the FEM solve; their E_d <= R_d
    // verdict is a report-layer check. Reference: docs/references/design-codes-ec7-tbdy.md.
    if (io.config && io.config->design_approach != model::DesignApproach::None) {
        const katai::core::DesignApproach da = to_core_design_approach(io.config->design_approach);
        if (katai::core::factors_material(da)) {
            // A MATERIAL FACTOR NEEDS A MATERIAL PARAMETER TO DIVIDE. factor_material_strength
            // divides c' and raises tan(phi'), and Hoek-Brown has neither: the factoring would
            // touch nothing, the rock would run at its CHARACTERISTIC strength, and the report
            // would say the design approach had been applied. A design verification that
            // silently used unfactored strength is the worst outcome this seam can produce, so
            // it is refused. EN 1997-1 factors c' and tan(phi') by name and gives no M-set for a
            // Hoek-Brown envelope; deriving one -- factoring sigma_ci, or m_b, or the equivalent
            // c'/phi' of Hoek, Carranza-Torres & Corkum 2002 -- would be this program inventing a
            // design rule, and the three choices do not agree with each other.
            for (const auto& mm : models)
                if (mm.type == katai::core::MaterialType::HoekBrown) {
                    refuse(R, "K2D-G015", "design approach",
                        "A material-factored design approach (EC7 DA1-C2 / DA3) cannot be applied "
                        "to a Hoek-Brown material: the partial factors divide c' and tan(phi'), "
                        "and this model has neither -- the rock would be solved at its full "
                        "CHARACTERISTIC strength while the report said the design approach had "
                        "been applied. EN 1997-1 gives no partial factor for a Hoek-Brown "
                        "envelope. Convert the envelope to an equivalent c' and phi' over the "
                        "confining range the problem spans (Hoek, Carranza-Torres & Corkum 2002), "
                        "apply the factors to those, and run the design phase "
                        "on a Mohr-Coulomb material -- or use a resistance-factored approach "
                        "(EC7 DA2, TBDY 2018), which does not touch the material at all.");
                    return R;
                }
            const katai::core::PartialFactors pf = katai::core::design_factors(da);
            for (auto& mm : models) katai::core::factor_material_strength(mm, pf);
            // The strengths DERIVED from those materials take the same factors: the depth gradient
            // of cohesion (part of the strength of the material it belongs to) and every interface
            // (a ratio of the strength of the material it was built from). Each is decided by its
            // own material -- an undrained one takes gamma_cu -- and none of them carried a design
            // factor before 2026-09 (design_code.hpp states what was measured).
            if (iface_strength_mat.size() != structures.interfaces.size() ||
                iface5_strength_mat.size() != structures.interfaces5.size()) {
                R.message = "Internal: an interface was built without the material its strength "
                            "comes from, so the design approach could not be applied to it.";
                return R;
            }
            for (size_t mi = 0; mi < profiles.size() && mi < models.size(); ++mi)
                katai::core::factor_profile_strength(profiles[mi], models[mi], pf);
            for (size_t k = 0; k < structures.interfaces.size(); ++k)
                katai::core::factor_interface_strength(structures.interfaces[k].props,
                                                       models[(size_t)iface_strength_mat[k]], pf);
            for (size_t k = 0; k < structures.interfaces5.size(); ++k)
                katai::core::factor_interface_strength(structures.interfaces5[k].props,
                                                       models[(size_t)iface5_strength_mat[k]], pf);
            if (pf.gamma_Q_unfav != 1.0) f += (pf.gamma_Q_unfav - 1.0) * f_loads;  // variable loads x gamma_Q
        }
    }

    // K0 initial-stress procedure: the undisturbed ground starts in geostatic equilibrium
    // (sigma'_v = effective overburden, sigma'_h = K0 sigma'_v), so self-weight produces ~zero
    // displacement and the post-processed stress is the real geostatic field. Layered + water-aware
    // via vertical integration; reuses the eff_unit_weight / ground_surface lambdas built above.
    katai::core::K0LayeredOptions k0opt;
    k0opt.k0 = k0_by_mat;
    k0opt.eff_unit_weight = eff_unit_weight;
    k0opt.ground_surface = ground_surface;
    // Total-stress targets are seeded in total stress (header note): u(x,y) hydrostatic. For
    // Undrained (C) the K0 value refers to total stresses rather than effective stresses -- the
    // seed has to be the stress the material is analysed in, or the initial state contradicts
    // the constitutive law from the first step.
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
    // K0 validity condition: the K0 procedure is correct ONLY
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
        // Staged phase: start from the previous phase's committed stresses (staged chaining).
        if (!io.init_states || io.init_states->empty() ||
            io.init_states->size() % (size_t)std::max(1, mesh.element_count) != 0) {
            R.message = "Internal: staged phase started without the previous phase's stress state.";
            return R;
        }
        init = *io.init_states;
        // RESET SMALL STRAIN (the `resetsmall` phase option). The
        // inherited state carries the small-strain history with it, which is exactly right when
        // the phases continue one loading path and exactly wrong when they do not -- the
        // typical case being a surcharge placed and removed to leave a preconsolidation
        // pressure behind, where the strain history it also leaves is an artefact of how the
        // state was built rather than something the soil would still remember.
        //
        // Only gamma_hist is cleared. Stress, gamma_p and the cap pressure pp are the state of
        // the ground and are NOT the strain history: resetting them would delete the
        // overconsolidation the surcharge was applied to create, which is the opposite of what
        // this option is for.
        if (io.config && io.config->reset_small_strain) {
            int cleared = 0;
            for (auto& gs : init)
                if (gs.gamma_hist != 0.0) { gs.gamma_hist = 0.0; ++cleared; }
            // Whether any material can even feel it is a different question from whether the
            // flag was honoured, and the user asked a question that deserves an answer either
            // way: a reset that met no small-strain material is a no-op, and saying so is
            // cheaper than letting someone conclude their stiffness was restored.
            const bool any_hssmall = std::any_of(
                pr.materials.begin(), pr.materials.end(), [](const model::Material& m) {
                    return m.model == model::SoilModel::HSsmall && m.G0ref > 0.0;
                });
            const std::string& pname = io.config->name;
            if (!any_hssmall)
                note(R, "K2D-M005", pname,
                     "Phase \"" + pname +
                         "\" asks for the small-strain history to be reset, but no material in "
                         "this model is Hardening Soil with small-strain stiffness, so there is "
                         "no history to reset and the phase is unchanged.");
            else
                note(R, "K2D-M005", pname,
                     "Phase \"" + pname + "\" starts with the small-strain history RESET: " +
                         std::to_string(cleared) +
                         " stress points meet it at G0 instead of the stiffness the earlier "
                         "phases had degraded them to. Stress, shear hardening and the "
                         "preconsolidation pressure are carried over untouched.");
        }
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
        // convergence), so it needs more, smaller load increments and a looser
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
        // Iterations per increment. This is a PATIENCE setting, not an accuracy one: an increment
        // that needs more is not solved differently, it is cut back and retried -- so the number
        // decides, on its own, whether a slowly-converging increment is reported as a capacity.
        // It was left at the phase strategy's 80 until KV-STR-004 measured what that costs. On
        // that fixture's refined mesh the increments need up to 96 iterations while their
        // out-of-balance force falls monotonically the whole way; at 80 the run reported "a
        // collapse mechanism after 10% of the load" on the Eigen backend and full convergence on
        // MKL, because which increments survive the budget turns on the linear solver's rounding.
        // At 200 BOTH backends converge to max|u| = 0.1822167 m -- seven figures apart -- and the
        // MKL run is 1.8x FASTER (37 s against 66 s), because the cut-backs a truncated increment
        // triggers cost more than the iterations it was denied.
        // Raising it is free where a real mechanism forms: a genuinely collapsing run ends when
        // the tangent goes singular and the solver refuses (measured: same limit load, same 122
        // iterations, at 80 and at 200), not by exhausting the budget.
        const int iters_default = nonsym ? 200 : 0;   // 0 = the phase strategy's own limit
        const int steps = io.numeric.steps > 0 ? io.numeric.steps : steps_default;
        const double tol = io.numeric.tolerance > 0.0 ? io.numeric.tolerance : tol_default;
        const int iters = io.numeric.max_iterations > 0 ? io.numeric.max_iterations : iters_default;
        const katai::core::LinearSolve solver = reusing_linear_solve(
            nonsym ? katai::linsolve::MatrixType::RealNonsymmetric
                   : katai::linsolve::MatrixType::RealSymmetricPositiveDefinite);

        // Structural SELF-WEIGHT (plate w [kN/m/m]; embedded pile gamma*A/L_s) -- the 2026-07 audit
        // fix: the weight fed only the dynamic mass before, so statics silently carried weightless
        // walls/piles. Consistent nodal line load over every ACTIVE chain of THIS phase, added to
        // BOTH f and f_loads: a gravity-start phase ramps it with the body force (ramp = f); a
        // K0/baseline initial phase ramps it as an unbalanced new load (ramp = f_loads -- the
        // wished-in-place wall settles under its own weight, the equivalent of a plastic nil-step);
        // a chained phase sees it inside ramp = f - B, so a plate already equilibrated by the parent
        // adds NOTHING (nil identity preserved) and a newly activated plate arrives incrementally
        // (the staged change). Permanent load -> assembled AFTER the gamma_Q variable-load scaling.
        // The Safety search reads f too: its structures are solved with the soil, so they bear on
        // the factor of safety with their weight as well as their stiffness. Dynamic reads neither
        // f nor f_loads (static equilibrium comes from the parent baseline). w = 0 (default)
        // contributes nothing -> bit-identical.
        {
            Eigen::VectorXd f_sw = Eigen::VectorXd::Zero(dofs.equation_count());
            katai::core::assemble_structural_weight(mesh, dofs, structures, kGravity, f_sw);
            f += f_sw;
            f_loads += f_sw;
        }
        // SAFETY analysis (phi-c reduction / SRM): the strategy lives in the engine (Stage B9:
        // katai/analysis/phase_solver/safety.hpp). The seam passes the registry-derived model
        // family flags, the phase's active structures and the composition root's solver callback;
        // the refusals and the honest lower-bound reporting are engine-owned. On success the phase
        // falls through to the common result tail, exactly as before.
        if (phase == InitialPhase::Safety) {
            // THE STRUCTURES TAKE PART IN THE SEARCH (since 2026-09). Until then the search
            // was handed none of them, and an active geogrid, anchor, plate carrying 150 kN/m/m or
            // embedded beam returned the factor of safety of the same mesh with the element
            // DEACTIVATED, bit for bit, while an interface left the soil on its two sides
            // unconnected; every one of them was refused. What the reduction does to each is
            // stated with safety_analysis: an interface's strength is reduced with the soil's, a
            // structure's own capacity is not.
            //
            // One element still cannot enter, and it is refused rather than approximated: a
            // PRESTRESSED anchor. The search re-solves the ground from the unstressed state, and a
            // lock-off force is a force applied to a ground that has already moved -- on the
            // unstressed mesh it would pull the wall into soil that carries no stress yet, with
            // nothing in the search to say when the anchor was locked. Starting the search from the
            // parent phase is what gives that force a meaning, and this build does not do that yet.
            // silent-drop-scope: none -- this loop only looks for a prestressed anchor to refuse; it
            // builds nothing, and every element it passes over was built by its own loop above.
            for (size_t si = 0; si < pr.structs.size(); ++si) {
                const auto& s = pr.structs[si];
                if (s.kind != model::StructKind::Anchor || !struct_on(si)) continue;
                if (s.material < 0 || s.material >= (int)pr.anchors.size()) continue;
                const auto& am = pr.anchors[s.material];
                if (!(am.prestress > 0.0)) continue;
                refuse(R, "K2D-G016", line_subject(s.name, s.x1, s.y1, s.x2, s.y2),
                       "Anchor \"" + s.name + "\" is prestressed (lock-off force " +
                           dnum(am.prestress) +
                           " kN) and active in a Safety analysis. The strength-reduction search in "
                           "this build re-solves the ground from the unstressed state, and a "
                           "lock-off force belongs to a ground that has already moved: applied to "
                           "the unstressed mesh it would pull on soil that carries no stress yet, "
                           "and the factor of safety would depend on that. Plates, geogrids, "
                           "embedded beams, interfaces and anchors without prestress do take part "
                           "in the search. Deactivate this anchor in the Safety phase to obtain the "
                           "factor of safety without it.");
                return R;
            }
            // What starting from the unstressed state means once structures are in it, said
            // whenever one is: a structure is present from the first increment of every trial, so
            // it also carries whatever the ground's settlement under its own weight does to it --
            // which, in the phases the user built, may have happened before it was installed.
            {
                size_t n_active = 0;
                for (size_t si = 0; si < pr.structs.size(); ++si) {
                    if (struct_on(si)) ++n_active;
                }
                if (n_active > 0)
                    note(R, "K2D-A018", "Safety",
                         "The factor of safety includes the " + std::to_string(n_active) +
                             " active structural element(s). The strength-reduction search "
                             "re-solves this phase from the unstressed state with them present "
                             "from the start, so each also carries what the ground's settlement "
                             "under its own weight does to it -- including settlement that, in "
                             "the phases built before it, may have happened before the element was "
                             "installed. No structural force is reported for this phase, because "
                             "the state the search stops at is not the state those phases built.");
            }
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
            if (!katai::core::solve_safety_phase(mesh, dofs, models, profiles, f, solver,
                                                 structures, sfin, R))
                return R;
            // A Safety phase leaves the ground as it found it -- solve_phases does not commit its
            // stresses forward -- and it leaves the structures as it found them too. Without this
            // the phase after it met a parent with no structural state and re-developed every
            // structural force from zero, although nothing had changed.
            if (io.prev && io.prev->ok) R.struct_state = io.prev->struct_state;
        } else {
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
        // AN INTERFACE INSTALLED IN THIS PHASE takes the stress the ground has at its line. Until the
        // phase before, its seam was tied and the ground carried that stress across it; untied now,
        // each side stops holding the other, and a joint that did not take the stress over at once
        // would see it applied as a load of the phase. The seed is the committed stress recovered at
        // each Newton-Cotes node: normal n.sigma.n as the joint's sigma_n0, shear t.sigma.n entered as
        // an initial slip offset (capped by the joint's Coulomb strength). It generalises the K0
        // formula the initial phase uses -- after a level K0 phase the recovered stress IS K0 sigma'_v
        // -- so the usual sequence, a wall activated right after the K0 initial phase, installs it
        // stress-free. On a seam whose soil sides differ (a split soil-soil joint) the two sides'
        // stresses are averaged.
        int seeded_interfaces = 0;
        if (io.chained) {
            std::vector<char> new3(structures.interfaces.size(), 0);
            std::vector<char> new5(structures.interfaces5.size(), 0);
            for (size_t c = 0; c < struct_records.size(); ++c) {
                if (c < carry_plan.record_installed.size() && !carry_plan.record_installed[c]) continue;
                for (size_t i = struct_records[c].interfaces[0]; i < struct_records[c].interfaces[1]; ++i)
                    new3[i] = 1;
                for (size_t i = struct_records[c].interfaces5[0]; i < struct_records[c].interfaces5[1]; ++i)
                    new5[i] = 1;
            }
            const bool any_new = std::find(new3.begin(), new3.end(), 1) != new3.end() ||
                                 std::find(new5.begin(), new5.end(), 1) != new5.end();
            if (any_new && phase != InitialPhase::GravityLoading) {
                R.message = std::string("An interface or a wall with interfaces is activated in a ") +
                            (phase == InitialPhase::Consolidation ? "consolidation"
                             : phase == InitialPhase::FullyCoupled ? "fully-coupled"
                                                                   : "dynamic") +
                            " phase. Activating one takes over the ground's stress at its line, and "
                            "this build does that in a Plastic (staged construction) phase only: "
                            "activate it in a Plastic phase first, then run this phase.";
                return R;
            }
            if (any_new) {
                const auto rec = katai::core::recover_nodal_stresses_from_gauss(mesh, init, act);
                const int nn = mesh.node_count;
                if (!static_carry) {
                    carry_init = katai::core::StructuralInit{};
                    carry_plan.full_datum = Eigen::VectorXd::Zero(dofs.total_dofs());
                }
                if (carry_init.interface_slip.size() !=
                    structures.interfaces.size() * (size_t)katai::core::iface::kPointCount)
                    carry_init.interface_slip.assign(
                        structures.interfaces.size() * (size_t)katai::core::iface::kPointCount, 0.0);
                if (carry_init.interface5_slip.size() !=
                    structures.interfaces5.size() * (size_t)katai::core::iface::kPointCount5)
                    carry_init.interface5_slip.assign(
                        structures.interfaces5.size() * (size_t)katai::core::iface::kPointCount5, 0.0);
                // The stress at one joint node: the soil side's, averaged with the other side's when
                // that side is ground too (its structure-side DOF is a mesh node's).
                const auto node_stress = [&](int soil_node, int struct_gdof) {
                    Eigen::Vector3d sig = rec.stress[(size_t)soil_node];
                    if (struct_gdof >= 0 && struct_gdof < 2 * nn)
                        sig = 0.5 * (sig + rec.stress[(size_t)(struct_gdof / 2)]);
                    return sig;
                };
                const auto seed = [&](const Eigen::Vector3d& sig, double c, double s,
                                      const katai::core::iface::InterfaceProps& props,
                                      double& sigma_n0, double& slip0) {
                    const double nx = -s, ny = c, tx = c, ty = s;
                    const double sn = nx * nx * sig(0) + ny * ny * sig(1) + 2.0 * nx * ny * sig(2);
                    double tau = tx * nx * sig(0) + ty * ny * sig(1) + (tx * ny + ty * nx) * sig(2);
                    const double tau_max =
                        std::max(0.0, props.c_i - std::min(sn, props.sigma_t) * std::tan(props.phi_i));
                    if (std::fabs(tau) > tau_max) tau = std::copysign(tau_max, tau);
                    sigma_n0 = sn;
                    slip0 = props.ks > 0.0 ? -tau / props.ks : 0.0;
                };
                const auto ncp = katai::core::iface::nc_points();
                for (size_t i = 0; i < structures.interfaces.size(); ++i) {
                    if (!new3[i]) continue;
                    auto& ie = structures.interfaces[i];
                    katai::core::iface::NodeCoords Xe;
                    for (int k = 0; k < 3; ++k) {
                        Xe(k, 0) = mesh.x[ie.soil_nodes[k]];
                        Xe(k, 1) = mesh.y[ie.soil_nodes[k]];
                    }
                    for (int q = 0; q < katai::core::iface::kPointCount; ++q) {
                        const int nd = ncp[q].node;
                        const auto fr = katai::core::iface::edge_frame(Xe, ncp[q].xi);
                        seed(node_stress(ie.soil_nodes[nd], ie.struct_dof[2 * nd]), fr.c, fr.s, ie.props,
                             ie.sigma_n0[q],
                             carry_init.interface_slip[i * katai::core::iface::kPointCount + q]);
                    }
                    ++seeded_interfaces;
                }
                const auto ncp5 = katai::core::iface::nc_points5();
                for (size_t i = 0; i < structures.interfaces5.size(); ++i) {
                    if (!new5[i]) continue;
                    auto& ie = structures.interfaces5[i];
                    katai::core::iface::NodeCoords5 Xe;
                    for (int k = 0; k < 5; ++k) {
                        Xe(k, 0) = mesh.x[ie.soil_nodes[k]];
                        Xe(k, 1) = mesh.y[ie.soil_nodes[k]];
                    }
                    for (int q = 0; q < katai::core::iface::kPointCount5; ++q) {
                        const int nd = ncp5[q].node;
                        const auto fr = katai::core::iface::edge_frame5(Xe, ncp5[q].xi);
                        seed(node_stress(ie.soil_nodes[nd], ie.struct_dof[2 * nd]), fr.c, fr.s, ie.props,
                             ie.sigma_n0[q],
                             carry_init.interface5_slip[i * katai::core::iface::kPointCount5 + q]);
                    }
                    ++seeded_interfaces;
                }
            }
        }
        // The structural baseline is assembled whenever the structures start from a state -- the
        // carried one, or the stress a newly installed joint took over.
        const bool struct_baseline = static_carry || seeded_interfaces > 0;
        Eigen::VectorXd B = Eigen::VectorXd::Zero(dofs.equation_count());
        if (baseline) {
            if (axi) katai::core::assemble_axisym_internal_force(mesh, dofs, init, B);  // r-weighted baseline
            else katai::core::assemble_internal_force(mesh, dofs, init, B, act);
            if (struct_baseline) {
                // Full structural baseline at (parent datum + committed plastic state). This
                // SUBSUMES the manual sigma_n0 terms below (coulomb_return adds sigma_n0
                // internally), so those loops must be SKIPPED on this path -- adding both would
                // double-count the wished-in-place lateral earth pressure.
                //
                // An anchor INSTALLED in this phase contributes nothing to that baseline -- its
                // installation datum is this phase's datum -- except its lock-off force, which is a
                // constant in its internal force. That force is what this phase applies, so it must
                // stay out of B, exactly as it did when no structure was carried at all.
                katai::core::Structures without_new_lockoff;
                const katai::core::Structures* base_structures = &structures;
                for (size_t ai = 0; ai < structures.anchors.size(); ++ai)
                    if ((!carry_plan.carried ||
                         (ai < carry_plan.new_anchor.size() && carry_plan.new_anchor[ai])) &&
                        structures.anchors[ai].prestress != 0.0) {
                        if (base_structures == &structures) {
                            without_new_lockoff = structures;
                            base_structures = &without_new_lockoff;
                        }
                        without_new_lockoff.anchors[ai].prestress = 0.0;
                    }
                B += katai::core::structural_internal_force(
                    mesh, dofs, models, *base_structures, carry_init,
                    axi ? katai::core::Kinematics::Axisymmetric
                        : katai::core::Kinematics::PlaneStrain);
            }
            // Wished-in-place interface sigma_n0 baseline is engine-owned now (Stage B9,
            // interface_baseline.hpp). Carry path skips it: sigma_n0 is already inside the
            // structural baseline above (coulomb_return adds it internally) -- adding both
            // would double-count the lateral earth pressure.
            if (!struct_baseline)
                katai::core::add_interface_sigma_n0_baseline(structures, mesh, dofs, B);
        }

        // WHAT THE WATER DOES ACROSS EVERY SPLIT SEAM, read from the interface's own cross
        // permeability rather than inherited from the fact that the mesh was split. Fully
        // permeable (the default, and what every model written before that field said) ties
        // the two sides to ONE pore equation; impermeable leaves the two the split produced.
        // The seams are the mesh splitter's own pairs, so nothing here has to re-find them.
        std::vector<int> pore_tie(mesh.node_count, -1);
        bool semi_permeable = false;
        const auto tie = [&](int dup, int orig) {
            if (dup >= 0 && dup < mesh.node_count && orig >= 0 && orig < mesh.node_count)
                pore_tie[dup] = orig;
        };
        for (const auto& w : walls) {
            // silent-drop-ok: nothing is skipped here -- a semi-permeable joint is RECORDED
            // and the phase refuses on it a few lines below (the engine owns that message,
            // which is why it is not raised at this seam).
            if (w.active && w.flow_barrier == 2) { semi_permeable = true; continue; }
            if (!w.active || w.flow_barrier == 0)   // an inactive seam is fully permeable
                for (const auto& sp : w.seam) tie(sp.left, sp.right);
        }
        for (const auto& si : soil_ifaces) {
            // silent-drop-ok: as above -- recorded, then refused by the phase strategy.
            if (si.active && si.flow_barrier == 2) { semi_permeable = true; continue; }
            if (!si.active || si.flow_barrier == 0)   // an inactive seam is fully permeable
                for (const auto& sg : si.seam) tie(sg.dup, sg.orig);
        }
        // The barriers the coupled phases CANNOT read: declared impermeable or semi-permeable, but
        // standing on a line the mesh was not split along (a plate without interfaces, or a wall
        // or interface that fell back to bonded, K2D-G009). With no seam there is no second pore
        // equation for the declaration to act on.
        std::vector<char> barrier_unread(pr.structs.size(), 0);
        for (size_t si = 0; si < pr.structs.size(); ++si) {
            const auto& st = pr.structs[si];
            const bool line_kind =
                st.kind == model::StructKind::Plate || st.kind == model::StructKind::Interface;
            barrier_unread[si] = line_kind && struct_on(si) && st.flow_barrier != 0 &&
                                 !plate_is_wall[si] && !iface_split[si];
        }
        if (phase == InitialPhase::Consolidation) {
            // --- Time-dependent (Biot) consolidation phase ---------------------------------------
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
                cm.eoed = schema_oedometer_modulus(Mt);   // only the automatic first time step reads it
            }
            cin.flow_edges = flow_edges_from(pr);
            cin.have_flow_bcs = any_flow_bc_declared(pr);
            bool well_in_phase = false;
            cin.drain_nodes = phase_drain_nodes(pr, mesh, io, well_in_phase);
            warn_no_flow_barrier(pr, barrier_unread, R);
            if (well_in_phase)
                warn(R, "K2D-A009", "wells",
                     "A well is active in this consolidation phase, but a consolidation analysis "
                     "solves the EXCESS pore pressure and takes no prescribed discharge: the "
                     "well's pumping is not applied here. Wells are applied in a steady "
                     "groundwater-flow calculation. Drains, which set the excess pore pressure "
                     "to zero, ARE applied.");
            cin.active = act;
            // Plates, anchors and geogrids now take part in the coupled solve; interfaces (and
            // the embedded walls built from them, which split the mesh) and embedded beams do
            // not, and the phase says why. The structural system, the lines to report and the
            // parent's displacement datum go in together: a wall installed in an earlier phase
            // carries its force at that datum, and this phase solves the increment from there.
            cin.has_embedded_beams = has_embedded;
            cin.structures = &structures;
            cin.diagrams = &diag_specs;
            cin.iface_diagrams = &iface_diags;
            cin.carry_full = static_carry ? &carry_plan.full_datum : nullptr;
            cin.has_semi_permeable_interface = semi_permeable;
            cin.pore_tie = &pore_tie;
            if (io.config) {
                cin.duration_day = io.config->duration; cin.time_steps = io.config->time_steps;
                // How the phase ends (`cstop`). With a state criterion the two above
                // are not read at all -- the phase marches until the state is reached and its
                // ANSWER is the time that took.
                cin.stop = to_core_consol_stop(io.config->consol_stop);
                cin.stop_min_pore = io.config->consol_min_pore;
                cin.stop_degree = io.config->consol_degree / 100.0;   // file is per cent, engine a ratio
                cin.first_dt = io.config->consol_first_step;
                cin.max_steps = io.config->consol_max_steps;
            }
            cin.yscale = ymax - ymin;
            const Eigen::VectorXd dF = f - B;   // configuration imbalance (the staged change)
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
            // --- Fully-coupled flow-deformation --------------------------------------------------
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
            warn_no_flow_barrier(pr, barrier_unread, R);
            if (well_in_fc)
                warn(R, "K2D-A009", "wells",
                     "A well is active in this fully-coupled phase. This build's coupled solver "
                     "takes drainage boundaries and drains, but not a prescribed well discharge, "
                     "so the well's pumping is not applied here. Run the dewatering as a "
                     "groundwater-flow calculation, or model it with a drain at the target head.");
            fin.active = act;
            fin.has_embedded_beams = has_embedded;
            fin.has_semi_permeable_interface = semi_permeable;
            fin.structures = &structures;
            fin.diagrams = &diag_specs;
            fin.iface_diagrams = &iface_diags;
            fin.carry_full = static_carry ? &carry_plan.full_datum : nullptr;
            fin.pore_tie = &pore_tie;
            if (io.config) { fin.duration_day = io.config->duration; fin.time_steps = io.config->time_steps; }
            fin.yscale = ymax - ymin;
            const Eigen::VectorXd dF = f - B;   // configuration imbalance (the staged change)
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
            // --- Transient groundwater flow only -------------------------------------------------
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
                // stored. An Undrained (C) cluster is in the same position -- its permeability is
                // not an input at all, the permeability fields not applying to the Non-porous and
                // Undrained (C) drainage types -- so it takes
                // the same impermeable-barrier treatment rather than a permeability nobody set.
                fm.nonporous = mat_total_stress[mi] != 0;
            }
            tin.flow_edges = flow_edges_from(pr);
            tin.active = act;
            // Wells and drains in a TRANSIENT flow phase: refused by name rather than dropped.
            // This solver takes a prescribed-head boundary set and no source term, so a well's
            // discharge has nowhere to enter and a NORMAL drain's one-sided rule (pore pressures
            // below the drain's head are left untouched) cannot be re-evaluated as the
            // head moves through the time steps -- applying it as a plain fixed head would let a
            // drain FEED water into ground that is drier than it, which is the opposite of what
            // a drain does. This build does not apply them here yet, and says so.
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
            // --- Dynamic (seismic) time-history analysis -----------------------------------------
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
                                                 carry_src ? &carry_plan.full_datum : nullptr,
                                                 io.prev, io.init_states, din, spdf, nonsymf, R))
                R.mesh = std::move(mesh);
            return R;
        } else {
        // The strategy lives in the engine (Stage B9: katai/analysis/phase_solver/
        // static_phase.hpp) together with the ramp / nil-step semantics (K0 baseline, staged
        // chaining, SSC time apportioning). This seam passes the common setup's neutral products
        // (loads, baseline, activity, Newton class, carry) and the composition root's solver
        // callback; on success the phase falls through to the common result tail, as before.
        katai::core::StaticPhase stin;
        stin.baseline = baseline;
        stin.nil_step = (use_k0 && (k0_nonlevel || flow)) || io.chained;
        stin.axisymmetric = axi;
        stin.load_steps = steps;
        stin.tolerance = tol;
        stin.max_iterations = iters;
        stin.line_search_window = io.numeric.line_search_window;
        stin.enforce_local_criteria = io.numeric.enforce_local_criteria;
        stin.substep_tolerance = io.numeric.substep_tolerance;
        // Stage fraction: a partial stage is a construction step only where there IS a stage, so the
        // fraction is read on chained phases and left at 1 on the initial one (the validator
        // refuses it there rather than letting a scaled gravity look like a partial excavation).
        stin.stage_fraction = (io.chained && io.config) ? io.config->sum_mstage : 1.0;
        if (stin.stage_fraction != 1.0)
            note(R, "K2D-A007", io.config->name,
                 "This phase applies only " + dnum(100.0 * stin.stage_fraction) +
                     "% of its staged change (mstage = " + dnum(stin.stage_fraction) +
                     "), so the configuration it describes is NOT reached: the remainder is "
                     "still carried by the soil. Results belong to the partial stage.");
        // Only a CHAINED (staged) phase carries time: the initial phase is TIMELESS --
        // an unconditional duration here once leaked 1 day of creep into the K0 phase and broke
        // the geostatic identity (measured: K0 max |u| = 26 mm = exactly mu* ln2 H on an SSC
        // column). The initial phase also has io.config (pr.initial, for activation flags), which
        // is why the gate is io.chained and not the config pointer.
        stin.time_interval_day = (io.chained && io.config) ? io.config->duration : 0.0;
        stin.active = act;
        stin.presc = presc;   // active prescribed displacements (empty = none)
        // The excess pore pressure the phase inherits belongs in its baseline too
        // (excess_pore_force.hpp): a chained phase starts in the parent's equilibrium, water
        // included, and ramps only its own change. Only here -- the coupled phases above keep the
        // effective-stress baseline their imbalance re-generates the pressure from.
        if (baseline && io.chained)
            katai::core::add_excess_pore_force(mesh, dofs, models, profiles, init, act, axi, B);
        if (!katai::core::solve_static_phase(mesh, dofs, models, profiles, init, f, f_loads, B,
                                             solver, structures, diag_specs, iface_diags,
                                             carry_init,
                                             struct_baseline ? &carry_plan.full_datum : nullptr,
                                             stin, R, io.out_states))
            return R;
        static_carry_used = static_carry;
        static_carry_missing = io.chained && any_struct_carry && !static_carry;
        // What the next phase needs to find this phase's structures again by structure: the
        // records (with the installation cohort each ended up in), the cohorts' datums in this
        // phase's numbering, and the count of mesh DOFs that keep their numbers.
        R.struct_state.records = struct_records;
        R.struct_state.install_datum = katai::core::install_datum_full(structures, dofs);
        R.struct_state.node_dofs = 2 * mesh.node_count;
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
    if (static_carry_used) {
        R.message += " Structural state continued from the parent phase (forces are totals).";
        if (carry_plan.installed > 0)
            R.message += " " + std::to_string(carry_plan.installed) +
                         " structure(s) installed in this phase on the ground as the parent phase "
                         "left it (their forces start from zero here).";
    }
    else if (static_carry_missing)
        R.message += " NOTE: the parent phase supplies no structural state (no structure in common "
                     "with this phase, or a consolidation/restored parent), so structural forces "
                     "re-develop from zero this phase.";
    if (!R.consol_time.empty()) {
        const double s_inf = R.consol_settlement.back();
        const double s0 = R.consol_settlement.front();
        const char* label = (phase == InitialPhase::FullyCoupled)
                                ? "Fully-coupled flow-deformation" : "Consolidation";
        R.message = std::string(label) + " solved: settlement " + std::to_string(s0) + " -> " +
                    std::to_string(s_inf) + " m over " + std::to_string(R.consol_time.back()) +
                    " days (excess pore -> " + std::to_string(R.consol_excess_pore.back()) +
                    " kPa). See the settlement-time curve below.";
        // A phase that ended on a STATE answers a question the timed one does not ask -- how long
        // -- so the time leads, and the definition of the number travels WITH the number. The
        // degree of consolidation here is a pressure ratio, and saying so is not pedantry:
        // the settlement ratio of the same name reaches 90% about 21% of the time earlier.
        if (R.consol_stop != katai::core::ConsolidationStop::TimeInterval) {
            char cbuf[520];
            if (R.consol_stop == katai::core::ConsolidationStop::MinExcessPore)
                std::snprintf(cbuf, sizeof(cbuf),
                    "%s: the excess pore pressure fell to %.4g kPa after %.6g day(s) and %d time "
                    "step(s) -- %.1f%% of the %.4g kPa this stage generated has dissipated. "
                    "Settlement %.6g -> %.6g m.",
                    label, R.consol_excess_pore.back(), R.consol_time.back(), R.iterations,
                    100.0 * R.consol_degree_reached, R.consol_pore_reference, s0, s_inf);
            else
                std::snprintf(cbuf, sizeof(cbuf),
                    "%s: %.1f%% consolidation reached after %.6g day(s) and %d time step(s). The "
                    "degree of consolidation here is a PRESSURE ratio -- maximum excess pore "
                    "pressure now (%.4g kPa) against the maximum this stage generated (%.4g kPa) -- "
                    "and NOT the settlement ratio of the same name, which reaches the same figure "
                    "earlier. Settlement %.6g -> %.6g m.",
                    label, 100.0 * R.consol_degree_reached, R.consol_time.back(), R.iterations,
                    R.consol_excess_pore.back(), R.consol_pore_reference, s0, s_inf);
            R.message = cbuf;
            R.message += " See the settlement-time curve below.";
        }
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
        if (!(n.substep_tolerance > 0.0)) n.substep_tolerance = ph.substep_tolerance;
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
        // and FullyCoupled commit it forward (staged chaining).
        if (ph.type != model::PhaseType::Safety && ph.type != model::PhaseType::TransientFlow &&
            ph.type != model::PhaseType::Dynamic)
            committed = std::move(next);
    }
    return out;
}

}  // namespace katai::app
