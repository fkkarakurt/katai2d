#pragma once
// KATAI 2D project data model (GUI). The usual geotechnical input workflow: a project holds a
// material database and a soil column defined by a borehole (layers with top/bottom levels + a water
// head). Materials are defined independently and assigned to layers. UI-agnostic plain data; the app
// edits it and the solver consumes it (mesh/solve wiring comes in later GUI steps).
// Reference: docs/references/gui-design.md; the file contract is docs/k2d-format.md.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace katai::model {

// Enum VALUES are file-stable (written as int in the project file) — a new model is appended at the END.
enum class SoilModel { LinearElastic, MohrCoulomb, HardeningSoil, HSsmall, SoftSoil, SoftSoilCreep,
                      HoekBrown };   // values are FILE-STABLE: append only
// Drainage type (material General data). Undrained (A): EFFECTIVE strength c', phi' + the pore-
// fluid bulk stiffness Kw/n (excess pore pressure generated, effective stress path; su is PREDICTED).
// Undrained (B): the same Kw/n machinery but the UNDRAINED strength is entered directly -- c = su,
// phi = 0 (a Tresca su envelope), which the solver enforces. Both share the effective-stress + Kw
// computation. Enum VALUES are kept file-stable (UndrainedB appended as 3); the GUI presents them in
// its own display order via an explicit index<->enum map (drainage_names is indexed by the enum
// value).
// Undrained (C): a TOTAL stress analysis. Stiffness is the undrained pair (E_u, nu_u close to
// 0.5) and strength the undrained shear strength (c = su, phi = 0), entered in the same boxes;
// no pore pressure is generated OR carried, so what the output calls effective stress is total
// stress, and K0 refers to total stress too. It is not a variant of (A)/(B):
// those separate the water from the skeleton, this one declines to.
enum class Drainage { Drained = 0, Undrained = 1, NonPorous = 2, UndrainedB = 3, UndrainedC = 4 };

inline const char* const* soil_model_names() {
    static const char* n[] = {"Linear elastic", "Mohr-Coulomb", "Hardening Soil", "HS small",
                              "Soft Soil", "Soft Soil Creep", "Hoek-Brown (rock)"};
    return n;
}
inline constexpr int kSoilModelCount = 7;   // combo count — updated together with the enum
// BOUNDS-CHECKED display name: the enum int comes from a file (a forward-version project
// can carry an out-of-range value) — raw `names()[i]` or a `& 3` mask silently prints the
// WRONG model name (measured: with the 5th model added, `& 3` showed "Linear elastic").
// All display surfaces use these two.
inline const char* soil_model_name(SoilModel m) {
    const int i = (int)m;
    return (i >= 0 && i < kSoilModelCount) ? soil_model_names()[i] : "Unknown model";
}
inline const char* const* drainage_names() {   // indexed by the enum VALUE (status text)
    static const char* n[] = {"Drained", "Undrained (A)", "Non-porous", "Undrained (B)",
                              "Undrained (C)"};
    return n;
}
inline constexpr int kDrainageCount = 5;
inline const char* drainage_name(Drainage d) {
    const int i = (int)d;
    return (i >= 0 && i < kDrainageCount) ? drainage_names()[i] : "Unknown drainage";
}

// One soil material data set. Parameters cover all supported models; only
// the fields relevant to `model` are used. Stiffness in kN/m^2, unit weight kN/m^3, angles in degrees.
struct Material {
    std::string name = "New material";
    SoilModel model = SoilModel::MohrCoulomb;
    Drainage drainage = Drainage::Drained;
    float color[3] = {0.85f, 0.78f, 0.55f};   // display colour (RGB 0..1)

    // General — unit weights + void ratio.
    double gamma_unsat = 17.0, gamma_sat = 20.0;
    double e_init = 0.5;   // initial void ratio

    // Mechanical — Linear elastic / Mohr-Coulomb.
    double E = 1.3e4;     // E'ref (Young's modulus)
    double nu = 0.3;      // Poisson's ratio
    double c = 1.0;       // c'ref cohesion
    double phi = 30.0;    // phi' friction angle [deg]
    double psi = 0.0;     // psi dilatancy angle [deg]
    double E_inc = 0.0;   // E' increment per metre depth
    double c_inc = 0.0;   // c' increment per metre depth
    double y_ref = 0.0;   // reference level for the increments
    bool tension_cutoff = true;     // tension cut-off (on by default)
    double tensile_strength = 0.0;  // sigma_t
    // Dilatancy cut-off (entered on the material's General tab): a dilating soil arrives at a
    // critical void ratio where dilatancy ends -- once the volume change has taken the soil to
    // e_max, the mobilised dilatancy angle is set to zero. With it OFF -- the default -- a dense
    // sand dilates without limit and its bearing capacity is over-predicted. e_init above is the
    // in-situ void ratio; e_max is the critical one. A minimum void ratio e_min plays no part in
    // the Hardening Soil dilatancy rule, so it is not carried here.
    bool dilatancy_cutoff = false;
    double e_max = 1.0;

    // Groundwater — stiffness of the pore fluid for Undrained (A)/(B), the parameters of the excess
    // pore pressure calculation. The bulk modulus of water here is not a property of the water:
    // it is a numerical value tied to the soil stiffness,
    //   K_w/n = 3(ν_u − ν')/((1 − 2ν_u)(1 + ν')) K',
    // and what the user chooses is how ν_u is arrived at. Two ways of defining the undrained
    // Poisson's ratio:
    //   und_mode = 0  Direct        — ν_u is entered (0.495, nearly incompressible; exactly 0.5 is
    //                                 singular)
    //   und_mode = 1  Skempton-B    — B is entered and, with α_Biot = 1,
    //                                 ν_u = (3ν' + B(1 − 2ν'))/(3 − B(1 − 2ν')) (Skempton 1954)
    // ν_u is NOT ν_ur (the unloading/reloading ratio above) — the two are easily confused.
    // These are per material because the pore fluid's stiffness is: two clays with different ν'
    // and different B do not share a K_w/n, and until now every undrained material in a model was
    // given ν_u = 0.495 whatever its data said.
    // A third definition (the Biot effective stress concept, α_Biot < 1 with K_w entered
    // directly) is deliberately absent: α_Biot enters the effective-stress split itself
    // (σ' = σ − α_Biot p), so honouring it in the undrained corner alone would be a half-truth.
    int und_mode = 0;
    double nu_u = 0.495;       // equivalent undrained Poisson ratio (und_mode = 0)
    double skempton_B = 0.0;   // Skempton's B (und_mode = 1); 0 is refused there, never assumed

    // Mechanical — Hardening Soil (+ HS small) (Schanz, Vermeer & Bonnier 1999; Benz 2007).
    double E50ref = 3.0e4, Eoedref = 3.0e4, Eurref = 9.0e4;
    double m = 0.5, nu_ur = 0.2, p_ref = 100.0;
    double Rf = 0.9;             // failure ratio qf/qa
    bool k0nc_auto = true;       // K0^nc auto = 1 - sin(phi')
    double k0nc = 0.5;
    // HS small.
    double G0ref = 1.2e5, gamma07 = 1.5e-4;

    // Mechanical — Soft Soil. λ*/κ* = the modified compression/swelling
    // indices (strain-based; conversion from Cc/Cs/e₀: λ* = Cc/(2.3(1+e₀)),
    // κ* ≈ 2Cs/(2.3(1+e₀)) — the GUI editor offers a converter, the stored
    // value is ALWAYS λ*/κ*). ν_ur, c, φ, ψ and K0NC are read from the shared fields above
    // (nu_ur, c, phi, psi, k0nc_auto/k0nc); M is derived internally from K0NC (Brinkgreve
    // 1994). OCR/POP scales the initial p_p seed (Initial tab).
    double lam_star = 0.10, kap_star = 0.02;
    // Soft Soil Creep (Vermeer & Neher 1999): modified creep index μ* = Cα/(2.3(1+e₀)); λ*/μ*
    // typically 15-25. Reference time τ = 1 day (a fixed constant; the 24-hour definition
    // of the NC line) — not an input.
    double mu_star = 0.005;

    // Hoek-Brown (the 2002 criterion; Hoek, Carranza-Torres & Corkum 2002). Rock keeps Hooke's
    // law, so the elastic pair is the shared E / nu above -- E is the ROCK MASS modulus E_rm
    // there, not the intact rock's.
    // The strength is these five, and they are the geologist's own vocabulary rather than a fitted
    // c' and phi': the intact rock's uni-axial compressive strength, its material constant, how
    // jointed the mass is, how much the excavation disturbed it, and how it dilates.
    //   sig_ci  |sigma_ci|, uni-axial compressive strength of the INTACT rock [kPa, > 0]
    //   mi      intact rock parameter [-] (about 4 for claystone, ~33 for granite)
    //   gsi     Geological Strength Index [-] (100 = intact, ~10 = crushed)
    //   hb_D    disturbance factor [-] (0 undisturbed, 1 heavily blasted)
    //   sig_psi the confining stress at which dilatancy has died out [kPa]; the mobilised
    //           dilatancy falls linearly from the shared psi (its value at sigma'_3 = 0) to
    //           zero at that confinement
    double sig_ci = 50000.0;
    double mi = 10.0;
    double gsi = 50.0;
    double hb_D = 0.0;
    double sig_psi = 0.0;

    // Groundwater — permeability (Groundwater tab; used by seepage/consolidation/flow).
    // Default ≈ a fine-medium sand, ~1 m/day ≈ 1.2e-5 m/s (Das, Principles of Geotech. Eng., Table 7.1).
    double kx = 1.0, ky = 1.0;   // [m/day]
    // Change of permeability with void ratio in a consolidation phase: log10(k / k0) = de / ck,
    // de the change in void ratio since the initial state (e_init). 0 = the permeability does
    // not change. Generally of the order of the compression index.
    double ck = 0.0;
    // Groundwater — unsaturated water-retention (van Genuchten 1980 + Mualem 1976; Groundwater tab
    // "van Genuchten"). Used by transient + fully-coupled flow in the unsaturated zone (suction ψ>0).
    // In the saturated regime (ψ≤0) these do NOT activate (S_e=1, k_rel=1) → saturated behaviour
    // reduces to classic seepage/consolidation. Defaults = USDA "sand" (Carsel & Parrish 1988, WRR
    // 24:755): α=0.145/cm=14.5/m, n=2.68, θr/θs=0.045/0.43≈0.10.
    double gw_ga = 14.5;    // g_a [1/m]  (van Genuchten α — inverse air-entry)
    double gw_gn = 2.68;    // g_n        (van Genuchten n, > 1)
    double gw_gl = 0.5;     // g_l        (Mualem pore-connectivity)
    double gw_Sres = 0.10;  // residual saturation S_res (S_sat = 1)

    // Interfaces — strength reduction factor Rinter (rigid = 1.0).
    bool rinter_rigid = true;
    double Rinter = 1.0;

    // Initial — lateral earth pressure coefficient (K0). Auto = 1 - sin(phi').
    bool k0_auto = true;
    double k0 = 0.5;
    // Initial — stress history / pre-overburden: 0 = none, 1 = OCR, 2 = POP.
    int oc_mode = 0;
    double OCR = 1.0;   // overconsolidation ratio
    double POP = 0.0;   // pre-overburden pressure [kN/m2]
};

// --- Structural material data sets (set types: Plates / Anchors / Geogrids / Embedded beams).
//     Each structural element type has its own parameter set, separate from soil materials.

// Plate: EA, EI, weight, Poisson, optional plastic Mp/Np.
struct PlateMaterial {
    std::string name = "Plate";
    float color[3] = {0.20f, 0.42f, 0.80f};
    bool elastoplastic = false;
    double EA = 5.0e6;    // axial stiffness [kN/m]
    double EI = 8.5e3;    // bending stiffness [kN m2/m]
    double w = 0.0;       // weight [kN/m/m]
    double nu = 0.0;      // Poisson's ratio
    double Mp = 0.0;      // plastic moment [kN m/m]  (0 = elastic)
    double Np = 0.0;      // plastic axial force [kN/m] (0 = elastic)
    // End bearing: a plate has no thickness and so no tip. With this set, soil plasticity is
    // excluded in the ground within the plate's equivalent thickness d = sqrt(12 EI / EA) of
    // its lower end, which keeps a loaded wall from punching through the ground under it.
    // Not the wall's true tip resistance, and not for sheet piles, whose tip is far thinner
    // than d.
    bool prevent_punching = false;
    double d() const { return EA > 0.0 ? std::sqrt(12.0 * EI / EA) : 0.0; }  // equivalent thickness
};

// Anchor (node-to-node / fixed-end): axial stiffness + capacities + out-of-plane spacing.
struct AnchorMaterial {
    std::string name = "Anchor";
    float color[3] = {0.70f, 0.35f, 0.12f};
    bool elastoplastic = false;
    double EA = 1.0e5;        // axial stiffness [kN]
    double Fmax_tens = 0.0;   // |Fmax,tension|   (0 = unlimited)
    double Fmax_comp = 0.0;   // |Fmax,compression| (0 = unlimited)
    double Lspacing = 1.0;    // out-of-plane spacing [m]
    // Lock-off force per anchor [kN], tension-positive. An anchor or a strut is tensioned
    // against the wall when it is installed; that force holds the excavation before any
    // further movement, and it is what every anchored-excavation benchmark specifies (the
    // DGGT / Schweiger triple-anchored wall locks off 768 / 945 / 980 kN). The anchor is an
    // elastic spring from that state afterwards, so the force follows the wall rather than
    // staying constant. 0 = installed slack.
    double prestress = 0.0;
};

// Geogrid (tension-only membrane): axial stiffness + optional tension cap.
struct GeogridMaterial {
    std::string name = "Geogrid";
    float color[3] = {0.10f, 0.60f, 0.32f};
    bool elastoplastic = false;
    double EA = 1.0e3;   // axial stiffness [kN/m]
    double Np = 0.0;     // tension cap [kN/m] (0 = unlimited)
};

// Embedded beam / pile row: material + geometry + skin & base resistance + spacing.
struct EmbeddedBeamMaterial {
    std::string name = "Embedded beam";
    float color[3] = {0.40f, 0.40f, 0.42f};
    double E = 3.0e7;        // Young's modulus [kN/m2]
    double gamma = 24.0;     // unit weight [kN/m3]
    double diameter = 0.40;  // pile diameter (circular) [m]
    double Lspacing = 2.5;   // out-of-plane spacing [m]
    double Tskin_max = 0.0;  // max skin resistance [kN/m] (0 = unlimited)
    // Skin resistance at the END of the beam (the second point it is drawn to) [kN/m]: the
    // resistance then varies linearly along the beam from Tskin_max at its first point to
    // this. 0 = the same along the whole beam.
    double Tskin_end = 0.0;
    double Fmax_base = 0.0;  // base resistance [kN]
};

// --- Structural elements + loads (Structures mode) ----------------------------------------------
enum class StructKind { Plate, Anchor, Geogrid, EmbeddedBeam, Interface };
inline const char* const* struct_kind_names() {
    static const char* n[] = {"Plate", "Anchor", "Geogrid", "Embedded beam", "Interface"};
    return n;
}
// A structural element drawn as a line (two endpoints). `material` indexes the material list for its
// kind (plates/anchors/geogrids/embedded); interfaces use the adjacent soil (material = -1).
struct StructElement {
    StructKind kind = StructKind::Plate;
    std::string name = "Element";
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    int material = -1;
    // Local mesh density (coarseness factor): target element size near this line is
    // multiplied by it (1 = global, 0.5 = twice as fine, ...). Refine/coarsen halve/double it.
    double coarseness = 1.0;
    // Interfaces attached to the element side (right-click -> positive / negative interface).
    // They use the adjacent soil's Rinter unless iface_material >= 0 (a soil material override).
    bool iface_pos = false, iface_neg = false;
    int iface_material = -1;
    // Interface elastic stiffnesses [kN/m3 = kPa per metre of relative movement]: normal kn and
    // shear ks. 0 = derived from the adjacent soil (R_inter^2 G over a virtual thickness of 0.1
    // element lengths), which is what every file before v19 means. A joint whose stiffness is
    // known -- measured, or given by a design guide or a reference solution -- is entered here,
    // since the derived value follows the mesh and the soil, not the joint.
    double iface_kn = 0.0, iface_ks = 0.0;
    // Virtual thickness factor of the interfaces: the derived stiffnesses are G_i / t and
    // E_oed,i / t with t = this factor times the model's element size (mesh.elem_size). The
    // thicker, the more elastic movement the contact has before it slips; 0.1 by default.
    double iface_vt = 0.1;
    // CROSS PERMEABILITY in a groundwater calculation (docs/k2d-format.md, `flow_barrier`). A wall
    // or interface is a line the water either crosses freely, cannot cross at all, or crosses
    // against a resistance:
    //   0 fully permeable   — no effect on flow. The DEFAULT, and what every project written
    //                         before this field said: the flow net runs through the line.
    //   1 impermeable       — a screen. The pore-pressure degrees of freedom on the two sides are
    //                         fully separated.
    //   2 semi-permeable    — the line passes q_n = dh / R per unit area, R = d/k the HYDRAULIC
    //                         RESISTANCE in units of TIME: for a wall of thickness d and
    //                         permeability k, R = d/k, and it is measured as the head difference
    //                         dh over the average discharge q through the wall per unit area,
    //                         d/k = dh/q. Neither d nor k separately matters.
    // The end points of an interface stay permeable, and that comes for free: the mesh splitter
    // keeps the end nodes shared, so water goes around the ends.
    // Longitudinal drainage conductivity (dk, the interface acting as a drain along its length) is
    // NOT carried yet — a declared gap rather than a silent zero.
    int flow_barrier = 0;
    double hydraulic_resistance = 0.0;   // d/k [day], read when flow_barrier = 2
    // EMBEDDED BEAM ONLY -- how the beam's connection point is attached (docs/k2d-format.md,
    // `conn`). The connection point is the pile's top: the end with the highest y, or for an
    // exactly horizontal pile the end with the lowest x.
    //   0 hinged  — the beam's top translations ARE the soil's there: the two undergo the same
    //               displacement, but not necessarily the same rotation. This is the DEFAULT
    //               when no structure shares the point, and it is what makes a pile loadable at
    //               its head at all.
    //   1 free    — the beam top is not coupled directly to the soil element it lies in; it
    //               still interacts with the soil through the skin and foot springs. This is
    //               the natural choice for a grout body, so that a ground anchor does not lose
    //               axial force into the soil at the connection.
    // A third type, Rigid (rotation coupled too, to a plate), is out of scope and
    // declared in docs/references/embedded-beam-formulation.md sec 8.
    int conn = 0;
};

enum class LoadKind { Point, Distributed };
inline const char* const* load_kind_names() {
    static const char* n[] = {"Point load", "Distributed load"};
    return n;
}
// Point load: applied at (x1,y1) with components (qx1,qy1). Distributed load: along (x1,y1)-(x2,y2),
// components vary linearly from (qx1,qy1) at the start to (qx2,qy2) at the end.
struct Load {
    LoadKind kind = LoadKind::Point;
    std::string name = "Load";
    double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
    double qx1 = 0, qy1 = -10, qx2 = 0, qy2 = -10;
    double coarseness = 1.0;   // local mesh density near the load (see StructElement)
};

// Deformation boundary conditions (model conditions, deformations). Automatic by default
// (sides on rollers, base fully fixed, top free); the user may override each edge.
enum class BCType { Free, NormallyFixed, HorizontallyFixed, VerticallyFixed, FullyFixed };
inline const char* const* bc_type_names() {
    static const char* n[] = {"Free", "Normally fixed", "Horizontally fixed", "Vertically fixed", "Fully fixed"};
    return n;
}
struct BoundaryConditions {   // default: all Free (the user applies BC explicitly)
    int xmin = (int)BCType::Free;   // left
    int xmax = (int)BCType::Free;   // right
    int ymin = (int)BCType::Free;   // base
    int ymax = (int)BCType::Free;   // top
};

// Groundwater flow boundary conditions (per polygon edge). Closed
// (impermeable, q_n = 0) is the natural FE condition and the default; Head prescribes the hydraulic
// head h [m elevation] (reservoir / far-field water level); Seepage marks a free-drainage face
// (downstream dam face, excavation wall) where water may exit at atmospheric pressure -- the exit
// point is found by the solver's active-set iteration.
// Appended, file-stable: an older build meets the value 3 outside its range and REFUSES the
// file, which is what the enum bound is for -- no version bump needed to stay honest.
enum class FlowBCType { Closed, Head, Seepage, Flux };
inline const char* const* flow_bc_names() {
    static const char* n[] = {"Closed (impermeable)", "Prescribed head", "Seepage face",
                              "Prescribed flux (inflow +)"};
    return n;
}

// A construction phase (staged construction). The INITIAL phase
// establishes the initial stress state (K0 procedure / gravity loading, selected separately);
// each USER phase changes the active configuration -- excavate (deactivate a soil region),
// fill (activate one), install or remove structures and loads -- and re-equilibrates starting
// from the previous phase's committed stresses (the staged-construction multiplier).
// A Safety phase computes the factor of safety of the configuration instead.
// Consolidation is a time-dependent (Biot) calculation phase: the configuration's load increment is
// applied at t=0 (undrained response generates excess pore pressure), which then dissipates over the
// phase's time interval, developing settlement (classic Terzaghi U-t).
// TransientFlow: time-dependent groundwater flow only (no deformation) -- the pore/head field evolves
// under time-varying hydraulic BCs (storage + backward-Euler).
// FullyCoupled: simultaneous flow + deformation (the most
// general analysis) -- Biot consolidation generalized with unsaturated van Genuchten/Mualem retention
// and Bishop effective stress (chi = S_eff). Both are time-dependent (reuse duration / time_steps).
// Dynamic: time-history (seismic) analysis -- M u'' + C u' + K u = -M r a_g(t) integrated by Newmark-β
// with Rayleigh damping. Linear-elastic skeleton (v1); horizontal base acceleration
// (site response / SSI). Boundaries: absorbing base + free-field sides (docs/references/dynamic-seismic-
// formulation.md). Enum VALUES are file-stable (appended); names array is indexed by the enum value.
enum class PhaseType { Plastic, Safety, Consolidation, TransientFlow, FullyCoupled, Dynamic };
inline constexpr int kPhaseTypeCount = 6;
inline const char* const* phase_type_names() {
    static const char* n[] = {"Plastic (staged construction)", "Safety (phi-c reduction)",
                              "Consolidation (time-dependent)", "Transient groundwater flow",
                              "Fully coupled flow-deformation", "Dynamic (seismic time-history)"};
    static_assert(sizeof(n) / sizeof(n[0]) == kPhaseTypeCount, "one name per PhaseType");
    return n;
}

// How a Consolidation phase ENDS (docs/k2d-format.md, `cstop`). Model-level enum,
// mapped to katai::core::ConsolidationStop in build_problem (like SoilModel/Drainage/PhaseType).
// TimeInterval is the classic one: the phase's own duration, split into equal steps. The other two
// end it on a STATE instead, and for them a time interval is not applicable at all -- the phase
// runs until the excess pore pressure has fallen far enough and REPORTS how long that took, which
// is the form the design question has ("how long until it has settled out?").
// DegreeOfConsolidation is a PRESSURE ratio, not Terzaghi's settlement ratio -- it is defined from
// the current maximum excess pore pressure over the maximum initial excess pore pressure, and that
// is stated explicitly because the name suggests the other definition. The two differ by ~20% in
// time on the 1-D column. katai::core::ConsolidationStop (analysis/results.hpp) carries the full
// statement.
// Values file-stable (append only).
enum class ConsolStop { TimeInterval, MinExcessPore, DegreeOfConsolidation };
inline const char* const* consol_stop_names() {
    static const char* n[] = {"Time interval", "Minimum excess pore pressure",
                              "Degree of consolidation"};
    return n;
}

// Base-motion waveform for a Dynamic phase (v1 synthetic input; a loaded record is a later phase).
// Record = a user accelerogram (real earthquake record): Phase.accel_record [m/s^2] +
// Phase.record_dt sampling interval; seismic_amp becomes a SCALE factor (1 = record as-is).
enum class SeismicWave { Harmonic, Ricker, Record };
inline const char* const* seismic_wave_names() {
    static const char* n[] = {"Harmonic (sine)", "Ricker pulse"};
    return n;
}

// EC7 (EN 1997-1) / TBDY 2018 design approach applied to a phase (v0.3 workstream B). Model-level
// enum, mapped to katai::core::DesignApproach in build_problem (like SoilModel/Drainage/PhaseType).
// None = characteristic values (no partial factors). Values file-stable (append only).
enum class DesignApproach { None, EC7_DA1_C1, EC7_DA1_C2, EC7_DA2, EC7_DA3, TBDY2018_Static, TBDY2018_Seismic };
inline const char* const* design_approach_names() {
    static const char* n[] = {"None (characteristic)", "Eurocode 7 - DA1 Comb.1",
                              "Eurocode 7 - DA1 Comb.2", "Eurocode 7 - DA2", "Eurocode 7 - DA3",
                              "TBDY 2018 (static)", "TBDY 2018 (seismic)"};
    return n;
}
struct Phase {
    std::string name = "Phase";
    PhaseType type = PhaseType::Plastic;
    // Groundwater conditions for THIS phase (per-phase water conditions). When
    // water_override is set, the phreatic surface is this polyline instead of the project's,
    // and the change from the previous phase is carried by the staged-construction imbalance
    // like any other change: the pore pressure and the effective weight both follow it.
    //
    // Staged dewatering is not a refinement. Lowering the water inside an excavation before
    // digging is a step in nearly every real excavation, and it is what loads the wall first:
    // the benchmark this was written for (the DGGT / Schweiger triple-anchored wall in Berlin
    // sand) lowers the level from -3 m to -17.9 m inside the pit before a single cubic metre
    // is dug. Without it that problem cannot be built at all.
    bool water_override = false;
    std::vector<double> wx, wy;   // the phase's own phreatic polyline (used when overriding)
    // Consolidation (PhaseType::Consolidation) time control: total time interval [day] and the
    // number of equal time steps used to integrate the dissipation. Ignored by other phase types.
    // Dynamic (PhaseType::Dynamic) reuses these: duration = total shaking time [s], time_steps = the
    // number of Newmark steps over it (so the step dt = duration/time_steps [s]).
    double duration = 1.0;
    int time_steps = 25;
    // Consolidation stop criterion. With consol_stop != TimeInterval the two
    // fields above are NOT used -- there is no interval to divide -- and the phase instead marches
    // until |p|max <= consol_min_pore [kPa], or until the degree of consolidation (as a PRESSURE
    // ratio, see ConsolStop) reaches consol_degree [%]. consol_first_step [day] is the first time
    // step; 0 = automatic (the Vermeer-Verruijt critical step).
    // consol_max_steps caps the march: reaching it REFUSES the phase, it does not report a time as
    // if the target had been met. Ignored by every other phase type.
    ConsolStop consol_stop = ConsolStop::TimeInterval;
    double consol_min_pore = 1.0;    // [kPa] (1 stress unit)
    double consol_degree = 90.0;     // [%]   (target degree of consolidation, a pressure ratio)
    double consol_first_step = 0.0;  // [day] 0 = automatic
    int consol_max_steps = 1000;
    // Dynamic (seismic) time-history input (PhaseType::Dynamic; ignored otherwise). Horizontal base
    // acceleration a_g(t): waveform seismic_wave, peak amplitude seismic_amp [m/s^2], dominant
    // frequency seismic_freq [Hz]. Rayleigh damping: ratio damping_ratio at the two target frequencies
    // rayleigh_f1, rayleigh_f2 [Hz] (the band where xi ~ target; outside it damping rises).
    SeismicWave seismic_wave = SeismicWave::Harmonic;
    double seismic_amp = 1.0;
    double seismic_freq = 2.0;
    // Accelerogram input (SeismicWave::Record): uniformly sampled base acceleration a_g [m/s^2]
    // at record_dt [s]. seismic_amp scales it (1 = as recorded). Beyond the record the input is
    // ZERO (free vibration). Stored IN the project file -- a run is reproducible without loose
    // side files (provenance belongs in the record's own documentation).
    std::vector<double> accel_record;
    double record_dt = 0.02;
    double damping_ratio = 0.05;
    double rayleigh_f1 = 1.0;
    double rayleigh_f2 = 8.0;
    // Lateral boundary treatment (Dynamic phase). false = free/roller sides (site-response columns
    // and idealized 1D checks); true = Lysmer FREE-FIELD sides (a boundary dashpot + a 1D free-field
    // driving force) so the sides follow the free field and do not reflect the interior's scattered
    // waves back in -- required for realistic 2D soil-structure models. The base stays rigid.
    bool seismic_free_field = false;
    // COMPLIANT (absorbing) BASE (Dynamic phase; Joyner & Chen 1975, input force with factor 2).
    // false = rigid base (the user's fixed bottom; the historical bit-identical path; the base
    // reflects all downgoing energy). true = the bottom boundary absorbs outgoing waves through
    // Lysmer dashpots (rho Vs of the deepest layer = the halfspace continues it) and the input is
    // applied THERE as the upward-propagating wave: KATAI takes the phase's a_g(t) as the bedrock
    // (within) motion and applies HALF of it internally (the upgoing incident wave). The solve is
    // then in TOTAL motion (the base moves). Combinable with free-field sides (the 1D side columns
    // solve on a compliant base too) and with dynamic_nonlinear (linear-limit identity verified).
    // Scope: horizontal SH only (base u_y stays fixed; disclosed in the phase message).
    bool seismic_compliant_base = false;
    // NONLINEAR dynamic solve (Dynamic phase): false = the LINEAR Newmark solve (default; the internal
    // force is K u, so the soil stays elastic and interfaces/geogrids/anchors do not slip/yield during
    // shaking -- fast, and the historical bit-identical path). true = the fully NONLINEAR Newmark solve
    // (per-step Newton on f_int(u) from the constitutive return mapping, so the soil can plastify and the
    // structures can slip/yield DURING the earthquake). Much slower (the tangent is
    // refactored every iteration), so it is opt-in. Needs a parent phase for the initial stress state.
    bool dynamic_nonlinear = false;
    // TBDY 2018 design spectrum overlay (Dynamic phase): map spectral acceleration coefficients S_S
    // (short period) and S_1 (1.0 s) from the AFAD hazard map, and the local site class (0=ZA..4=ZE).
    // The design spectrum S_ae(T) is compared with the computed surface response spectrum.
    double tbdy_ss = 1.0;
    double tbdy_s1 = 0.4;
    int site_class = 2;   // 0=ZA, 1=ZB, 2=ZC, 3=ZD, 4=ZE (TBDY Table 16.1)
    // EC8 (EN 1998-1) elastic design-spectrum overlay (Dynamic phase; optional, next to TBDY --
    // the Eurocode user compares the computed surface spectrum against THEIR code's target).
    // a_gR = reference peak ground acceleration on type-A ground [g]; gamma_I = importance factor
    // (a_g = gamma_I * a_gR); ground 0..4 = A..E (EN Table 3.1); type 0 = Type 1 (M_s > 5.5),
    // 1 = Type 2. Recommended EN parameter values (national annexes may differ; stated in output).
    bool ec8_enabled = false;
    double ec8_agr = 0.2;
    double ec8_gamma = 1.0;
    int ec8_ground = 2;   // 0=A, 1=B, 2=C, 3=D, 4=E
    int ec8_type = 0;     // 0=Type 1, 1=Type 2
    // EC7 / TBDY 2018 design approach for this phase (None = characteristic values). The material-
    // factored approaches (DA1-C2, DA3) reduce c'/tan(phi') and scale variable loads before the solve.
    DesignApproach design_approach = DesignApproach::None;
    // Stage-fraction TARGET (docs/k2d-format.md `mstage`): the fraction of this phase's
    // staged change that is actually applied. 1 = the whole stage, which is what staged
    // construction normally means and what every project written before this field said. A
    // smaller value applies part of it and leaves the rest -- half an excavation lift, a partly
    // built embankment -- and the next phase continues from that partly-changed state.
    //
    // It applies to a STAGED (chained) phase only. The initial phase establishes the in-situ
    // state; scaling gravity there would not be a partial construction step, it would be a
    // different planet, so the value is refused rather than quietly obeyed.
    double sum_mstage = 1.0;
    // IGNORE UNDRAINED BEHAVIOUR (docs/k2d-format.md `ignoreund`): for this phase, materials
    // whose drainage is Undrained (A) or (B) are treated as drained -- no new excess pore pressure is
    // generated, and the excess pore pressure generated in earlier phases is kept as it is.
    // Strength parameters are untouched, so an Undrained (B) material still carries its
    // c_u with phi = 0. The standard use is a phase where the undrained response is not the
    // question being asked (establishing an initial state, or a long-term stage).
    bool ignore_undrained = false;
    // RESET SMALL STRAIN (docs/k2d-format.md `resetsmall`): start
    // this phase with the small-strain history cleared, so a Hardening Soil small material meets
    // it at its maximum stiffness G0 instead of wherever the previous phases left it.
    //
    // The reason for the option is that strain history outlives its cause. The typical case is
    // overconsolidation modelled by placing and removing a surcharge: the loading is there to
    // leave a preconsolidation pressure behind, but it also leaves a strain history that -- in
    // the real soil -- ageing erased long before the analysis begins. Placing and removing the
    // surcharge has already triggered that history; this option clears it.
    //
    // The strain history is a tensor that a reversal resets in the direction it happens in
    // (Benz 2007), so an unloading that follows a loading recovers its stiffness on its own.
    // What this flag adds is the engineer's statement that the history up to here -- reversal
    // or not -- is not the history that matters now.
    bool reset_small_strain = false;
    // Reset displacements to zero: the total displacement reported from this phase on counts
    // from the START of this phase -- what the phases before it moved is left out. The state is
    // untouched (stresses, structural forces); it is the read-out that starts again, which is
    // what "the settlement caused by the excavation" or "by the tunnel" means when the building
    // stood there before. One part of the state does start again with it: the small-strain
    // history of Hardening Soil with small-strain stiffness (see reset_small_strain), because the
    // strain the earlier phases left is not the history a new event meets.
    bool reset_displacements = false;
    // NUMERICAL CONTROLS for this phase. Zero means
    // "let the program choose", which is what every project that never touches them says, so
    // the defaults stay exactly where they are: derived from the material class (Hardening Soil
    // at a tolerated error of 1%, Mohr-Coulomb at 1e-6, a linear problem at 1e-10) or, in a
    // Safety phase, from the strength-reduction search's own trial settings.
    //
    // They belong in the FILE and not only in the program because a published number has to be
    // reproducible by someone else: "the settlement is 42 mm" is a claim about a model AND about
    // the numerics it was solved with, and if the second half cannot be written down, nobody can
    // re-run it. Reviewers may also legitimately ask whether an answer moved because the physics
    // moved or because the stopping rule did -- a question only a file that carries the stopping
    // rule can answer.
    //
    // `load_steps` is the number of increments the load is split into. This is NOT a cap on an
    // automatically chosen step count ("max steps"): KATAI splits the load into a fixed
    // number of increments (the solver still cuts back adaptively when one will not converge).
    // The field is named for what it does rather than for the box it resembles.
    double tolerance = 0.0;     // tolerated relative force residual ("tolerated error")
    int load_steps = 0;         // load increments for this phase
    int max_iterations = 0;     // Newton iterations per increment ("max iterations")
    // The FOURTH control, and the one that took longest to earn its place here. The three above
    // govern the EQUILIBRIUM iteration; this one governs the CONSTITUTIVE integration -- the
    // error tolerance of the substepping that walks a stress point along its material law inside
    // one increment. It has to be a phase control for the same reason the other three are: an
    // answer is a claim about the numerics as well as the model, and until 2026-08-25 this
    // particular number could only be set through an environment variable, which is exactly the
    // half of the claim that cannot be written down or re-run by someone else.
    //
    // It also is not interchangeable with `tolerance`. Tightening the equilibrium residual below
    // the integration error does not buy accuracy, it buys iterations grinding against noise --
    // measured on the Hardening Soil oedometer, where a LOOSER integration costs MORE Newton
    // iterations, not fewer (0.9.0 N-1). A run whose two tolerances are set against each other is
    // a run nobody can interpret; carrying both in the file is what makes that visible.
    //
    // 0 = the material class's own default (Hardening Soil: 1e-5, set by measurement). Only the
    // models with an error-controlled integrator read it; every other model ignores it, so a file
    // that sets it on a Mohr-Coulomb phase is not wrong, it is simply not consulted.
    double substep_tolerance = 0.0;  // constitutive integration error tolerance (STOL)
    // Active flags per project object, by index (kept in sync by the GUI; an entry missing
    // because the vector is short counts as ACTIVE -- new objects default to active).
    std::vector<char> poly_active, struct_active, load_active, disp_active, hydro_active;
    bool active_poly(size_t i) const { return i >= poly_active.size() || poly_active[i]; }
    bool active_struct(size_t i) const { return i >= struct_active.size() || struct_active[i]; }
    bool active_load(size_t i) const { return i >= load_active.size() || load_active[i]; }
    bool active_disp(size_t i) const { return i >= disp_active.size() || disp_active[i]; }
    bool active_hydro(size_t i) const { return i >= hydro_active.size() || hydro_active[i]; }
    // DRY REGIONS: a soil region switched off in this phase whose space holds no free water -- a
    // tunnel kept dry, a pit pumped dry. By default an excavated region below the water table
    // holds the water that stands in it, and that water presses on the ground around it; a region
    // flagged here does not. A MISSING entry counts as not dry, so every file written before the
    // field existed reads as it did. It means nothing for an active region.
    std::vector<char> poly_dry;
    bool dry_poly(size_t i) const { return i < poly_dry.size() && poly_dry[i]; }
    // REGION WATER: the pore pressure of one soil region in this phase, where the water table does
    // not describe it. kRegionWaterTable (the default, and every missing entry) takes the phase's
    // water table; kRegionWaterLevel puts a hydrostatic level of the region's own at
    // poly_water_level[i] -- the pit lowered to its floor while the ground around it keeps the
    // table; kRegionWaterInterpolate makes the pore pressure vary linearly with depth between the
    // values just above and just below the region along each vertical -- the ground between a pit
    // floor pumped dry and the toe of its wall, where the water arriving from outside meets the
    // pumped level. That last is the classic way a design analysis states the water under an
    // excavation without solving its flow; a benchmark specified that way (ERTC7) cannot be
    // reproduced without it: hydrostatic water from the floor instead left the passive ground
    // with twice the effective stress at the toe.
    static constexpr int kRegionWaterTable = 0, kRegionWaterLevel = 1, kRegionWaterInterpolate = 2;
    std::vector<int> poly_water;
    std::vector<double> poly_water_level;
    int water_mode(size_t i) const { return i < poly_water.size() ? poly_water[i] : 0; }
    double water_level(size_t i) const {
        return i < poly_water_level.size() ? poly_water_level[i] : 0.0;
    }
    bool any_region_water() const {
        for (int m : poly_water)
            if (m != kRegionWaterTable) return true;
        return false;
    }
    // THE MATERIAL OF A REGION IN THIS PHASE, by soil-region index: -1 (and every missing entry)
    // = the material the region was drawn with. A region can change what it is made of between
    // phases -- ground that is improved, a fill that sets, a clay that is undrained while it is
    // loaded and drained in the long term -- and the phase says which material it is in THIS
    // phase. The stresses the region carries are kept; what is re-read from them is the state
    // only the new model has (its cap and its shear hardening).
    std::vector<int> poly_material;
    int material_of(size_t i, int drawn) const {
        return i < poly_material.size() && poly_material[i] >= 0 ? poly_material[i] : drawn;
    }
    bool any_material_change() const {
        for (int m : poly_material)
            if (m >= 0) return true;
        return false;
    }
    // THE SECTION OF A PLATE IN THIS PHASE, by structure index: an index into the project's
    // plate materials, -1 (and every missing entry) = the one the plate was drawn with. A
    // lining that replaces a shield, a wall that is thickened: the plate keeps the forces it
    // carries, and what it does from here on it does with the new stiffness and weight.
    std::vector<int> struct_material;
    int plate_of(size_t si, int drawn) const {
        return si < struct_material.size() && struct_material[si] >= 0 ? struct_material[si] : drawn;
    }
    bool any_plate_change() const {
        for (int m : struct_material)
            if (m >= 0) return true;
        return false;
    }
    // LINE CONTRACTION of a plate, by structure index: the area loss of the closed line the
    // plate is part of, as a fraction (0.005 = 0.5%), TOTAL at the end of this phase. The plate
    // is given the shrinkage strain that goes with it -- half the area loss, since the area of
    // a closed line goes with the square of its length -- which by itself changes no stress; it
    // is the ground around the plate that decides how much of it is realised. Missing entries
    // are 0. A later phase repeats the value to keep it, and 0 takes it back.
    std::vector<double> struct_contraction;
    double contraction_of(size_t si) const {
        return si < struct_contraction.size() ? struct_contraction[si] : 0.0;
    }
    bool any_contraction() const {
        for (double c : struct_contraction)
            if (c != 0.0) return true;
        return false;
    }
    // RE-TENSIONING of an anchor, by structure index [kN per anchor, tension-positive]: the anchor
    // -- already in place, carried from the phase before -- is jacked to this force in this phase
    // and locked again where the wall comes to. It is an event of this phase only: the next phase
    // continues from the lock (a spring from there), and a missing entry or 0 is no jack.
    std::vector<double> struct_retension;
    double retension_of(size_t si) const {
        return si < struct_retension.size() ? struct_retension[si] : 0.0;
    }
    bool any_retension() const {
        for (double p : struct_retension)
            if (p != 0.0) return true;
        return false;
    }
    // PORE PRESSURES FROM A STEADY-STATE GROUNDWATER FLOW CALCULATION instead of the phreatic
    // line. The flow is solved on the ground active in this phase: below the water table (this
    // phase's, else the project's) every boundary of that ground -- the model's sides and the faces
    // an excavation exposes -- is held at the head of the table above it, above the table it is a
    // seepage face, the base is closed, and a wall or interface marked impermeable closes the line
    // it runs along. So water flows under a cut-off wall towards a pit pumped down to its floor,
    // and the pore pressures in front of the wall are those of that flow, not of two separate
    // hydrostatic levels. Off = the phreatic line, as every file written before it.
    bool pore_from_flow = false;
    // A consolidation phase with a duration places its load over that duration, linearly. On =
    // the load goes on at the phase's first instant instead, and the whole duration is left for
    // the pore pressure it generates to dissipate: the column of the one-dimensional theory,
    // which is loaded at once. Read by a consolidation phase that ends on a time only.
    bool consol_load_at_start = false;
    // TRANSIENT FLOW ONLY (`headramp`): the Head edges reach the heads of this phase linearly over
    // its duration, from the heads the phase starts from -- a reservoir lowered or raised at a
    // rate. Off (the default), they change at the phase's first instant: an instantaneous
    // drawdown. Read when the phase starts from the head a groundwater-flow phase before it ended
    // with; from the water table there is no head to move from.
    bool head_ramp = false;
    // PSEUDO-STATIC ACCELERATION (`kh`, `kv`): the inertia of the ground in the equivalent static
    // method of a seismic check, as body forces kh W horizontally (toward +x) and kv W vertically
    // (upward positive: it takes kv W off the weight, as (1 - kv) does in TBDY 2018 Eq. 16.31), W
    // the total weight of the active ground -- the pore water moves with the soil, its pressure
    // stays the static one -- and of the structures that weigh something. Fractions of g. Each
    // Plastic phase states its own (one that states none has none); a Safety phase reduces the
    // strength of the state its parent left under its parent's coefficients.
    double pseudo_kh = 0.0, pseudo_kv = 0.0;
};

// A soil region: an explicit polygon (drawn by the user) with a material. Classic, rock-solid:
// what you draw is exactly the polygon -- no auto-split, no clusters, no surprises.
struct SoilPolygon {
    std::string name = "Soil";
    int material = -1;
    std::vector<double> x, y;       // polygon vertices (implicitly closed)
    std::vector<int> edge_bc;       // per-edge boundary condition (BCType; size == #vertices), 0 = Free
    // Groundwater flow BC per edge (FlowBCType; empty or size == #vertices, 0 = Closed) and the
    // prescribed head [m] for FlowBCType::Head edges (same size; ignored otherwise).
    std::vector<int> edge_flow;
    std::vector<double> edge_head;
    // Prescribed boundary flux [m/day] for FlowBCType::Flux edges (same size; ignored otherwise),
    // INFLOW POSITIVE -- the same sign convention as a source term, which is positive for a
    // recharge well. It is a specific discharge normal to the boundary: the flow solve integrates
    // it over the edge as the boundary flux term q of the weak form. Flow calculations only -- a
    // CONSOLIDATION analysis does not carry a non-zero prescribed outflow.
    std::vector<double> edge_flux;
    // Local mesh density of the region (coarseness factor): the target element size
    // inside this polygon is multiplied by it (1 = global, 0.5 = twice as fine, 2 = coarser).
    double coarseness = 1.0;
};

// ---------------------------------------------------------------- stratigraphy from boreholes --
//
// How ground data actually arrives: a log per borehole, giving the level of each layer boundary
// and the water table at one x. Until this existed the engineer had to turn those logs into
// polygons by hand, which is both the slowest part of setting a model up and the easiest place to
// put a number in the wrong row.
//
// The contract (docs/k2d-format.md, `strata` / `boreholes`) is the borehole workflow every user of
// this kind of program already knows:
//
//   A borehole is a location in the drawing area at which the position of the soil layers and the
//   water table is given. With several boreholes the levels are interpolated between them. Every
//   soil layer is used throughout the whole model contour, so all layers appear in all boreholes;
//   their top and bottom boundaries may vary from borehole to borehole, which gives
//   non-horizontal layers of non-uniform thickness as well as layers that locally have a zero
//   thickness.
//
// Two consequences are load-bearing and are why the types look like this:
//
//   * THE LAYER LIST IS GLOBAL, THE LEVELS ARE PER BOREHOLE. A stratum carries the name and the
//     material; a borehole carries where its boundaries are. A layer missing at one location is
//     not a missing row, it is a ZERO THICKNESS there (bottom level equal to top level in that
//     borehole).
//   * ONE BOREHOLE MEANS HORIZONTAL, AND THE LEVELS REACH THE MODEL EDGE. A single borehole gives
//     a horizontal water surface that extends to the model boundaries, and the layers, used
//     throughout the whole model contour, behave the same way. So outside the outermost borehole
//     the nearest one's levels continue unchanged; nothing is extrapolated on a slope.
//
// WHAT THIS DOES NOT DO, and the boundary is the point: boreholes GENERATE polygons, they do not
// replace them. `polygons` remains the model -- what the mesher meshes and what a phase activates
// -- and the boreholes are the record of where it came from. Generation is an action the user
// takes, not something the solver does on the way past, because a geometry with two sources of
// truth is a geometry that can disagree with itself.
struct Stratum {
    std::string name = "Layer";
    int material = -1;                // index into Project::materials
};

struct Borehole {
    std::string name = "BH";
    double x = 0.0;                   // location along the model [m]
    // Levels of the layer boundaries at this borehole, top down: level[j] is the TOP of stratum j
    // and level[j+1] its bottom, so the size is strata.size() + 1 and the stack is continuous by
    // construction (the top boundary of an underlying layer is the lower boundary of the overlying
    // layer). Equal consecutive values are a layer that pinches out
    // here, which is legal and meant.
    std::vector<double> level;
    // The phreatic level at this borehole (a "Head" water condition). Heads from several
    // boreholes combine into a non-horizontal water surface.
    bool has_head = true;
    double head = 0.0;                // [m]
};

// A line prescribed displacement: every mesh node on
// the segment gets the set components IMPOSED (ramped 0 -> value over the phase that activates
// it, like a load; a following phase in which it stays active holds the nodes where they are, as
// an unchanged load adds nothing), the unset components stay free. A set component with value 0 is a rigid support
// line. Activated per phase exactly like loads (Phase::disp_active); the classic use is a
// displacement-controlled rigid footing, whose force is then read from the reactions.
struct PrescribedDisp {
    std::string name = "Displacement";
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;   // the line [m]
    bool set_ux = false; double ux = 0.0;             // horizontal component [m]
    bool set_uy = true;  double uy = 0.0;             // vertical component [m]
    double coarseness = 1.0;                          // local mesh density factor (like loads)
};

// A hydraulic condition drawn INSIDE the model (docs/k2d-format.md, `hydros`): a line along
// which water is taken out of or put into the ground. Two kinds,
// because the ground is told two different things:
//
//   Well  — a DISCHARGE is prescribed: |Q| per unit time per unit width out of plane. The head
//           is whatever that extraction produces, down to h_min, the lowest head the well can
//           reach: once the groundwater head has fallen to h_min, no further extraction occurs.
//   Drain — a HEAD is prescribed: the pore pressure in every node of the drain is reduced to
//           the value that corresponds to the given head. A Normal drain only takes water away
//           (pore pressures already below the drain's head are untouched); a Vacuum drain holds
//           the head at its value in both directions, which is what vacuum consolidation does.
//
// Both are activated per phase (Phase::hydro_active), exactly like loads: dewatering is a
// construction step, and a well that cannot be switched off after the pit is built is not one.
// Enum values are file-stable.
enum class HydroKind { Well = 0, Drain = 1 };
enum class WellBehaviour { Extraction = 0, Infiltration = 1 };
enum class DrainBehaviour { Normal = 0, Vacuum = 1 };

inline const char* const* hydro_kind_names() {
    static const char* n[] = {"Well", "Drain"};
    return n;
}
inline constexpr int kHydroKindCount = 2;
inline const char* hydro_kind_name(HydroKind k) {
    const int i = (int)k;
    return (i >= 0 && i < kHydroKindCount) ? hydro_kind_names()[i] : "Unknown hydraulic condition";
}

struct HydroLine {
    std::string name = "Well";
    HydroKind kind = HydroKind::Well;
    int behaviour = 0;                                // WellBehaviour or DrainBehaviour by kind
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;    // the line [m]
    // Well: the magnitude of the discharge [m3/day per m out of plane] and the lowest head the
    // well can draw the ground down to [m]. Drain: the head it holds [m].
    double q = 0.0;
    double h_min = 0.0;
    double head = 0.0;
    double coarseness = 1.0;   // local mesh density factor (like loads and displacements)
};

// How the initial phase establishes the in-situ stress state (the "Initial stress"
// selector): the K0 procedure (geostatic equilibrium, the default), gravity loading
// (switch-on self-weight, for ground where the K0 field does not equilibrate), or an
// immediate phi-c reduction of the gravity state (single-phase slope factor of safety).
// File-stable append-only enum (.k2d `initial_procedure`).
enum class InitialProcedure { K0Procedure, GravityLoading, Safety };

// Global mesh-generation settings (Mesh mode). The mesher's target triangle area is
// 0.5 * elem_size^2; the per-object `coarseness` factors scale the target locally.
// Part of the file contract: a .k2d determines its own discretization.
struct MeshSettings {
    double elem_size = 2.0;   // target element edge length [m]
    int order = 6;            // element type: 6 = quadratic tri6, 15 = quartic tri15
    bool auto_refine = true;  // ~2x finer near structural elements and loads
    // RELATIVE ELEMENT SIZE: 0 = elem_size above is the size; above 0 the size is this factor times
    // a twentieth of the model's diagonal (element_size below), so a model drawn wider keeps the
    // same density of elements. 2 / 1.33 / 1 / 0.67 / 0.5 are very coarse ... very fine.
    double rel = 0.0;
};
// A twentieth of the model's diagonal: the element size of a relative factor of 1.
inline constexpr double kRelativeElementBase = 1.0 / 20.0;

struct Project {
    std::string name = "Untitled project";
    std::vector<Material> materials;            // soil & interfaces
    std::vector<PlateMaterial> plates;          // plate material sets
    std::vector<AnchorMaterial> anchors;        // anchor material sets
    std::vector<GeogridMaterial> geogrids;      // geogrid material sets
    std::vector<EmbeddedBeamMaterial> embedded; // embedded beam material sets
    std::vector<SoilPolygon> polygons;          // soil regions (drawn, or generated from boreholes)
    // The ground data the polygons above may have been generated FROM. Empty on a model drawn by
    // hand, which is still the ordinary way to build one. See Stratum/Borehole for the contract.
    std::vector<Stratum> strata;                // the global layer list (name + material)
    std::vector<Borehole> boreholes;            // one log per location
    std::vector<StructElement> structs;         // structural elements (Structures mode)
    std::vector<Load> loads;                    // external loads
    std::vector<PrescribedDisp> disps;          // line prescribed displacements
    std::vector<HydroLine> hydros;              // wells and drains (hydraulic conditions)
    // Staged construction phases AFTER the initial phase (empty = classic single-phase solve).
    // phases[0] runs from the initial state, phases[k] from phases[k-1]'s committed stresses.
    std::vector<Phase> phases;
    // Initial-phase activation (the initial phase is soil only by default; structures and
    // loads are installed in later phases). Same index semantics as Phase.
    Phase initial;
    BoundaryConditions bc;                      // deformation boundary conditions (auto default)
    bool has_water = true;                      // groundwater level (a polyline, may be sloped)
    std::vector<double> wx, wy;                 // water polyline
    double x_min = 0.0, x_max = 40.0;           // model horizontal extent (BC edges)
    double y_min = 0.0, y_max = 20.0;           // model vertical extent (BC edges)
    // Analysis mode: plane strain (default) or axisymmetric. In axisymmetric mode the
    // x-coordinate is the RADIUS r (x=0 is the symmetry axis) and the integration is r-weighted; used
    // for circular footings, piles, shafts and cylinders: soil, water, excavation and fill, and
    // interfaces (a pile is a soil region with interfaces); plates, anchors, geogrids and embedded
    // beams are refused there, and only Plastic and Safety phases follow the initial one.
    bool axisymmetric = false;
    // How the initial phase establishes the in-situ state, and the global mesh settings.
    // Both belong to the file contract: a .k2d must determine the whole run, or no case
    // is reproducible from the file alone (docs/k2d-format.md).
    InitialProcedure initial_procedure = InitialProcedure::K0Procedure;
    MeshSettings mesh;

    double top() const { return y_max; }
    double bottom() const { return y_min; }

    // A new project starts EMPTY (blank drawing area). One "Sand" material is in the database so the
    // user has something to assign; no geometry, no water until the user draws them.
    static Project demo() {
        Project p;
        Material sand;
        sand.name = "Sand"; sand.model = SoilModel::MohrCoulomb;
        sand.color[0] = 0.93f; sand.color[1] = 0.86f; sand.color[2] = 0.55f;
        sand.gamma_unsat = 17.0; sand.gamma_sat = 20.0; sand.e_init = 0.5;
        sand.E = 1.3e4; sand.nu = 0.3; sand.c = 1.0; sand.phi = 30.0; sand.psi = 0.0;
        sand.kx = sand.ky = 1.0;   // fine-medium sand ~ 1 m/day (van Genuchten "sand" defaults apply)
        p.materials.push_back(sand);
        p.has_water = false;
        return p;
    }
};

// EDITING THE OBJECT LISTS. Every phase -- the initial one included -- holds what it does to the
// project's objects BY INDEX: which regions, structures, loads, prescribed displacements, wells
// and drains are active; which regions are dry, have water of their own or change material;
// which plates change section or contract. An edit that shifts an index has to shift every one
// of those with it: deleting structure i and leaving the flags where they were hands structure
// i's activation to the structure after it, silently, in every later phase -- the model solved
// is then another model than the one on the screen. So the lists are edited ONLY through the
// functions below (the Studio's delete commands call them), and test_project_edit checks each
// against the project built without the deleted object.
namespace detail {
template <class V>
inline void erase_index(V& v, size_t i) {
    if (i < v.size()) v.erase(v.begin() + (long)i);
}
// Insert n copies of `value` before position i of a per-object vector whose missing entries
// read as `missing`: a vector that is short of i only needs padding when the inserted value is
// not what a missing entry already says.
template <class V, class T>
inline void insert_index(V& v, size_t i, size_t n, const T& value, const T& missing) {
    if (i < v.size()) {
        v.insert(v.begin() + (long)i, n, value);
    } else if (!(value == missing)) {
        v.resize(i, missing);
        v.insert(v.end(), n, value);
    }
}
template <class F>
inline void for_each_phase(Project& p, F f) {
    f(p.initial);
    for (auto& ph : p.phases) f(ph);
}
// A reference into a list from which entry `gone` was removed: the entries after it move down
// one place; a reference to the removed entry itself becomes `orphan`.
inline void remap_reference(int& ref, size_t gone, int orphan = -1) {
    if (ref < 0) return;
    if (ref == (int)gone) ref = orphan;
    else if (ref > (int)gone) --ref;
}
}  // namespace detail

inline void erase_polygon(Project& p, size_t i) {
    if (i >= p.polygons.size()) return;
    p.polygons.erase(p.polygons.begin() + (long)i);
    detail::for_each_phase(p, [i](Phase& ph) {
        detail::erase_index(ph.poly_active, i);
        detail::erase_index(ph.poly_dry, i);
        detail::erase_index(ph.poly_water, i);
        detail::erase_index(ph.poly_water_level, i);
        detail::erase_index(ph.poly_material, i);
    });
}

// Insert `pieces` before position i. In every phase each piece takes what polygon `like` had
// there (an index into the list BEFORE the insertion; -1 = the defaults) -- its activation, its
// dry flag, its water and its material -- which is what a polygon cut into pieces needs: the
// pieces are excavated, kept dry, flooded and improved wherever the whole was.
inline void insert_polygons(Project& p, size_t i, const std::vector<SoilPolygon>& pieces, int like) {
    if (pieces.empty()) return;
    i = std::min(i, p.polygons.size());
    const size_t n = pieces.size();
    detail::for_each_phase(p, [&](Phase& ph) {
        const bool has = like >= 0;
        const size_t l = has ? (size_t)like : 0;
        const char on = (!has || ph.active_poly(l)) ? 1 : 0;
        const char dry = (has && ph.dry_poly(l)) ? 1 : 0;
        const int wm = has ? ph.water_mode(l) : Phase::kRegionWaterTable;
        const double wl = has ? ph.water_level(l) : 0.0;
        const int mat = (has && l < ph.poly_material.size()) ? ph.poly_material[l] : -1;
        detail::insert_index(ph.poly_active, i, n, on, (char)1);
        detail::insert_index(ph.poly_dry, i, n, dry, (char)0);
        detail::insert_index(ph.poly_water, i, n, wm, Phase::kRegionWaterTable);
        detail::insert_index(ph.poly_water_level, i, n, wl, 0.0);
        detail::insert_index(ph.poly_material, i, n, mat, -1);
    });
    p.polygons.insert(p.polygons.begin() + (long)i, pieces.begin(), pieces.end());
}

// The target element edge length the mesh is generated with [m]: elem_size, or under a relative
// size (mesh.rel > 0) that factor times a twentieth of the model's diagonal. Everything that reads
// the element size -- the mesher and the virtual thickness of an interface -- reads it here.
inline double element_size(const Project& p) {
    if (p.mesh.rel > 0.0) {
        const double d = std::hypot(p.x_max - p.x_min, p.y_max - p.y_min);
        if (d > 0.0) return p.mesh.rel * kRelativeElementBase * d;
    }
    return p.mesh.elem_size;
}

inline void erase_struct(Project& p, size_t i) {
    if (i >= p.structs.size()) return;
    p.structs.erase(p.structs.begin() + (long)i);
    detail::for_each_phase(p, [i](Phase& ph) {
        detail::erase_index(ph.struct_active, i);
        detail::erase_index(ph.struct_material, i);
        detail::erase_index(ph.struct_contraction, i);
        detail::erase_index(ph.struct_retension, i);
    });
}
inline void erase_load(Project& p, size_t i) {
    if (i >= p.loads.size()) return;
    p.loads.erase(p.loads.begin() + (long)i);
    detail::for_each_phase(p, [i](Phase& ph) { detail::erase_index(ph.load_active, i); });
}
inline void erase_disp(Project& p, size_t i) {
    if (i >= p.disps.size()) return;
    p.disps.erase(p.disps.begin() + (long)i);
    detail::for_each_phase(p, [i](Phase& ph) { detail::erase_index(ph.disp_active, i); });
}
inline void erase_hydro(Project& p, size_t i) {
    if (i >= p.hydros.size()) return;
    p.hydros.erase(p.hydros.begin() + (long)i);
    detail::for_each_phase(p, [i](Phase& ph) { detail::erase_index(ph.hydro_active, i); });
}

// THE MATERIAL LISTS. A material is referred to by index from the regions, the borehole layers,
// the interfaces that name a soil of their own, the structures, and the phases that change a
// region's material or a plate's section. Deleting one moves every later index down; what
// pointed at the deleted material itself is left WITHOUT one: a region or a layer or a
// structure without a material is refused when the model is calculated, so the gap shows, and
// an interface or a phase change falls back to what it is when nothing is named (the soil beside
// the joint; the material the region or plate was drawn with). That fallback is a change of
// model, which is why material_uses() exists: the Studio says what will change and asks first.
enum class MaterialList { Soil, Plate, Anchor, Geogrid, EmbeddedBeam };
inline MaterialList material_list_of(StructKind k) {
    switch (k) {
        case StructKind::Plate:   return MaterialList::Plate;
        case StructKind::Anchor:  return MaterialList::Anchor;
        case StructKind::Geogrid: return MaterialList::Geogrid;
        default:                  return MaterialList::EmbeddedBeam;   // Interface has no set
    }
}
struct MaterialUses {
    int regions = 0, layers = 0, interfaces = 0, structures = 0, phase_changes = 0;
    int total() const { return regions + layers + interfaces + structures + phase_changes; }
};
inline MaterialUses material_uses(const Project& p, MaterialList list, size_t m) {
    MaterialUses u;
    const int mi = (int)m;
    auto count_phase = [&](const Phase& ph) {
        if (list == MaterialList::Soil)
            for (int v : ph.poly_material) u.phase_changes += (v == mi);
        if (list == MaterialList::Plate)
            for (size_t s = 0; s < ph.struct_material.size() && s < p.structs.size(); ++s)
                u.phase_changes += (ph.struct_material[s] == mi);
    };
    if (list == MaterialList::Soil) {
        for (const auto& P : p.polygons) u.regions += (P.material == mi);
        for (const auto& L : p.strata) u.layers += (L.material == mi);
        for (const auto& s : p.structs) u.interfaces += (s.iface_material == mi);
    } else {
        for (const auto& s : p.structs)
            if (s.kind != StructKind::Interface && material_list_of(s.kind) == list)
                u.structures += (s.material == mi);
    }
    count_phase(p.initial);
    for (const auto& ph : p.phases) count_phase(ph);
    return u;
}
inline void erase_material(Project& p, MaterialList list, size_t m) {
    switch (list) {
        case MaterialList::Soil:
            if (m >= p.materials.size()) return;
            p.materials.erase(p.materials.begin() + (long)m);
            for (auto& P : p.polygons) detail::remap_reference(P.material, m);
            for (auto& L : p.strata) detail::remap_reference(L.material, m);
            for (auto& s : p.structs) detail::remap_reference(s.iface_material, m);
            detail::for_each_phase(p, [m](Phase& ph) {
                for (int& v : ph.poly_material) detail::remap_reference(v, m);
            });
            return;
        case MaterialList::Plate:
            if (m >= p.plates.size()) return;
            p.plates.erase(p.plates.begin() + (long)m);
            for (auto& s : p.structs)
                if (s.kind == StructKind::Plate) detail::remap_reference(s.material, m);
            detail::for_each_phase(p, [&p, m](Phase& ph) {
                for (size_t s = 0; s < ph.struct_material.size(); ++s)
                    if (s >= p.structs.size() || p.structs[s].kind == StructKind::Plate)
                        detail::remap_reference(ph.struct_material[s], m);
            });
            return;
        case MaterialList::Anchor:
            if (m >= p.anchors.size()) return;
            p.anchors.erase(p.anchors.begin() + (long)m);
            break;
        case MaterialList::Geogrid:
            if (m >= p.geogrids.size()) return;
            p.geogrids.erase(p.geogrids.begin() + (long)m);
            break;
        case MaterialList::EmbeddedBeam:
            if (m >= p.embedded.size()) return;
            p.embedded.erase(p.embedded.begin() + (long)m);
            break;
    }
    for (auto& s : p.structs)
        if (s.kind != StructKind::Interface && material_list_of(s.kind) == list)
            detail::remap_reference(s.material, m);
}

// THE VERTICES OF A REGION. Each edge carries its own boundary conditions (edge_bc, edge_flow,
// edge_head, edge_flux; one entry per edge, edge k running from vertex k to vertex k+1), and a
// list may be empty when nothing is stated for that kind. Splitting edge k at its midpoint gives
// both halves edge k's conditions; deleting vertex k merges edges k-1 and k into the edge that
// arrives at k, which keeps its own.
inline void split_edge(SoilPolygon& P, size_t k) {
    const size_t n = P.x.size();
    if (n < 2 || k >= n) return;
    const size_t j = (k + 1) % n;
    const double xm = 0.5 * (P.x[k] + P.x[j]), ym = 0.5 * (P.y[k] + P.y[j]);
    P.x.insert(P.x.begin() + (long)k + 1, xm);
    P.y.insert(P.y.begin() + (long)k + 1, ym);
    auto dup = [k](auto& v) {
        if (k < v.size()) v.insert(v.begin() + (long)k + 1, v[k]);
    };
    dup(P.edge_bc); dup(P.edge_flow); dup(P.edge_head); dup(P.edge_flux);
}
inline void erase_vertex(SoilPolygon& P, size_t k) {
    if (P.x.size() <= 3 || k >= P.x.size()) return;
    P.x.erase(P.x.begin() + (long)k);
    P.y.erase(P.y.begin() + (long)k);
    detail::erase_index(P.edge_bc, k);
    detail::erase_index(P.edge_flow, k);
    detail::erase_index(P.edge_head, k);
    detail::erase_index(P.edge_flux, k);
}

} // namespace katai::model
