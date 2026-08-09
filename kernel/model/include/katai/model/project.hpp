#pragma once
// KATAI 2D project data model (GUI). Mirrors the PLAXIS 2D Input working logic: a project holds a
// material database and a soil column defined by a borehole (layers with top/bottom levels + a water
// head). Materials are defined independently and assigned to layers. UI-agnostic plain data; the app
// edits it and the solver consumes it (mesh/solve wiring comes in later GUI steps).
// Reference: PLAXIS 2D 2025.1 Reference Manual §3-§4 (Soil mode); docs/references/gui-design.md.

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace katai::model {

// Enum VALUES are file-stable (written as int in the project file) — a new model is appended at the END.
enum class SoilModel { LinearElastic, MohrCoulomb, HardeningSoil, HSsmall, SoftSoil, SoftSoilCreep };
// Drainage type (PLAXIS Material > General). Undrained (A): EFFECTIVE strength c', phi' + the pore-
// fluid bulk stiffness Kw/n (excess pore pressure generated, effective stress path; su is PREDICTED).
// Undrained (B): the same Kw/n machinery but the UNDRAINED strength is entered directly -- c = su,
// phi = 0 (a Tresca su envelope), which the solver enforces. Both share the effective-stress + Kw
// computation. Enum VALUES are kept file-stable (UndrainedB appended as 3); the GUI presents them in
// PLAXIS order via an explicit index<->enum map (drainage_names is indexed by the enum value).
// Undrained (C): a TOTAL stress analysis. Stiffness is the undrained pair (E_u, nu_u close to
// 0.5) and strength the undrained shear strength (c = su, phi = 0), entered in the same boxes;
// no pore pressure is generated OR carried, so what the output calls effective stress is total
// stress, and K0 refers to total stress too (MMM section 2.7). It is not a variant of (A)/(B):
// those separate the water from the skeleton, this one declines to.
enum class Drainage { Drained = 0, Undrained = 1, NonPorous = 2, UndrainedB = 3, UndrainedC = 4 };

inline const char* const* soil_model_names() {
    static const char* n[] = {"Linear elastic", "Mohr-Coulomb", "Hardening Soil", "HS small",
                              "Soft Soil", "Soft Soil Creep"};
    return n;
}
inline constexpr int kSoilModelCount = 6;   // combo count — updated together with the enum
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

// One soil material data set (PLAXIS "Soil" material). Parameters cover all supported models; only
// the fields relevant to `model` are used. Stiffness in kN/m^2, unit weight kN/m^3, angles in degrees.
struct Material {
    std::string name = "New material";
    SoilModel model = SoilModel::MohrCoulomb;
    Drainage drainage = Drainage::Drained;
    float color[3] = {0.85f, 0.78f, 0.55f};   // display colour (RGB 0..1)

    // General — unit weights + void ratio.
    double gamma_unsat = 17.0, gamma_sat = 20.0;
    double e_init = 0.5;   // initial void ratio

    // Mechanical — Linear elastic / Mohr-Coulomb (PLAXIS MMM §3).
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
    // Dilatancy cut-off (PLAXIS MMM Eq. 5.16b / Fig. 5.6, entered on the material's General
    // tab): a dilating soil arrives at a critical void ratio where dilatancy ends. With it OFF
    // -- the default, as in PLAXIS -- a dense sand dilates without limit and its bearing
    // capacity is over-predicted. e_init above is the in-situ void ratio; e_max is the critical
    // one. PLAXIS also stores e_min, and its own manual says e_min "is not used within the
    // context of the Hardening-Soil model", so it is not carried here.
    bool dilatancy_cutoff = false;
    double e_max = 1.0;

    // Groundwater — stiffness of the pore fluid for Undrained (A)/(B) (PLAXIS Reference §6.1.2.17
    // "Parameters for excess pore pressure calculation"; MMM §2.4). The bulk modulus of water is
    // not a property of the water in PLAXIS: it is a numerical value tied to the soil stiffness,
    // K_w/n = 3(ν_u − ν')/((1 − 2ν_u)(1 + ν')) K' (Eq. 2-50), and what the user chooses is how ν_u
    // is arrived at. Two ways, PLAXIS's two suboptions of the ν-undrained definition:
    //   und_mode = 0  Direct        — ν_u is entered (PLAXIS default 0.495; exactly 0.5 is singular)
    //   und_mode = 1  Skempton-B    — B is entered and ν_u follows from Eq. 2-55 with α_Biot = 1
    // ν_u is NOT ν_ur (the unloading/reloading ratio above) — the manual flags that confusion too.
    // These are per material because the pore fluid's stiffness is: two clays with different ν'
    // and different B do not share a K_w/n, and until now every undrained material in a model was
    // given ν_u = 0.495 whatever its data said.
    // PLAXIS's third option (the Biot effective stress concept, α_Biot < 1 with K_w entered
    // directly) is deliberately absent: α_Biot enters the effective-stress split itself
    // (Eq. 2-60), so honouring it in the undrained corner alone would be a half-truth.
    int und_mode = 0;
    double nu_u = 0.495;       // equivalent undrained Poisson ratio (und_mode = 0)
    double skempton_B = 0.0;   // Skempton's B (und_mode = 1); 0 is refused there, never assumed

    // Mechanical — Hardening Soil (+ HS small) (PLAXIS MMM §6.4 / §7).
    double E50ref = 3.0e4, Eoedref = 3.0e4, Eurref = 9.0e4;
    double m = 0.5, nu_ur = 0.2, p_ref = 100.0;
    double Rf = 0.9;             // failure ratio qf/qa
    bool k0nc_auto = true;       // K0^nc auto = 1 - sin(phi')
    double k0nc = 0.5;
    // HS small.
    double G0ref = 1.2e5, gamma07 = 1.5e-4;

    // Mechanical — Soft Soil (PLAXIS MMM §10). λ*/κ* = the modified compression/swelling
    // indices (strain-based; conversion from Cc/Cs/e₀: λ* = Cc/(2.3(1+e₀)),
    // κ* ≈ 2Cs/(2.3(1+e₀)), Table 10-2 — the GUI editor offers a converter, the stored
    // value is ALWAYS λ*/κ*). ν_ur, c, φ, ψ and K0NC are read from the shared fields above
    // (nu_ur, c, phi, psi, k0nc_auto/k0nc); M is derived internally from K0NC (Brinkgreve
    // 1994). OCR/POP scales the initial p_p seed (Initial tab).
    double lam_star = 0.10, kap_star = 0.02;
    // Soft Soil Creep (PLAXIS MMM §11): modified creep index μ* = Cα/(2.3(1+e₀)); λ*/μ*
    // typically 15-25. Reference time τ = 1 day (a PLAXIS constant; the 24-hour definition
    // of the NC line) — not an input.
    double mu_star = 0.005;

    // Groundwater — permeability (PLAXIS Groundwater tab; used by seepage/consolidation/flow).
    // Default ≈ a fine-medium sand, ~1 m/day ≈ 1.2e-5 m/s (Das, Principles of Geotech. Eng., Table 7.1).
    double kx = 1.0, ky = 1.0;   // [m/day]
    // Groundwater — unsaturated water-retention (van Genuchten 1980 + Mualem 1976; PLAXIS Groundwater
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
    // Initial — stress history / pre-overburden (PLAXIS): 0 = none, 1 = OCR, 2 = POP.
    int oc_mode = 0;
    double OCR = 1.0;   // overconsolidation ratio
    double POP = 0.0;   // pre-overburden pressure [kN/m2]
};

// --- Structural material data sets (PLAXIS set types: Plates / Anchors / Geogrids / Embedded beams).
//     Each structural element type has its own parameter set, separate from soil materials.

// Plate (PLAXIS Ref §5.6 / MMM §18.3): EA, EI, weight, Poisson, optional plastic Mp/Np.
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
    double d() const { return EA > 0.0 ? std::sqrt(12.0 * EI / EA) : 0.0; }  // equivalent thickness
};

// Anchor (node-to-node / fixed-end): axial stiffness + capacities + out-of-plane spacing (MMM §18.1).
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

// Geogrid (tension-only membrane): axial stiffness + optional tension cap (MMM §18.2).
struct GeogridMaterial {
    std::string name = "Geogrid";
    float color[3] = {0.10f, 0.60f, 0.32f};
    bool elastoplastic = false;
    double EA = 1.0e3;   // axial stiffness [kN/m]
    double Np = 0.0;     // tension cap [kN/m] (0 = unlimited)
};

// Embedded beam / pile row (MMM §18.4): material + geometry + skin & base resistance + spacing.
struct EmbeddedBeamMaterial {
    std::string name = "Embedded beam";
    float color[3] = {0.40f, 0.40f, 0.42f};
    double E = 3.0e7;        // Young's modulus [kN/m2]
    double gamma = 24.0;     // unit weight [kN/m3]
    double diameter = 0.40;  // pile diameter (circular) [m]
    double Lspacing = 2.5;   // out-of-plane spacing [m]
    double Tskin_max = 0.0;  // max skin resistance [kN/m] (0 = unlimited)
    double Fmax_base = 0.0;  // base resistance [kN]
};

// --- Structural elements + loads (Structures mode; PLAXIS Reference Manual §5) ------------------
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
    // Local mesh density (PLAXIS Coarseness factor): target element size near this line is
    // multiplied by it (1 = global, 0.5 = twice as fine, ...). Refine/coarsen halve/double it.
    double coarseness = 1.0;
    // Interfaces attached to the element side (PLAXIS: right-click -> positive / negative interface).
    // They use the adjacent soil's Rinter unless iface_material >= 0 (a soil material override).
    bool iface_pos = false, iface_neg = false;
    int iface_material = -1;
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

// Deformation boundary conditions (PLAXIS "Model conditions > Deformations"). Automatic by default
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

// Groundwater flow boundary conditions (PLAXIS GroundwaterFlow BCs, Reference Manual). Closed
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

// A construction phase (PLAXIS Phases explorer / staged construction). The INITIAL phase
// establishes the initial stress state (K0 procedure / gravity loading, selected separately);
// each USER phase changes the active configuration -- excavate (deactivate a soil region),
// fill (activate one), install or remove structures and loads -- and re-equilibrates starting
// from the previous phase's committed stresses (PLAXIS SumMstage). A Safety phase computes the
// factor of safety of the configuration instead.
// Consolidation is a time-dependent (Biot) calculation phase: the configuration's load increment is
// applied at t=0 (undrained response generates excess pore pressure), which then dissipates over the
// phase's time interval, developing settlement (PLAXIS "Consolidation" phase, classic Terzaghi U-t).
// TransientFlow: time-dependent groundwater flow only (no deformation) -- the pore/head field evolves
// under time-varying hydraulic BCs (storage + backward-Euler; PLAXIS "Groundwater flow, transient").
// FullyCoupled: simultaneous flow + deformation (PLAXIS "Fully coupled flow-deformation", the most
// general analysis) -- Biot consolidation generalized with unsaturated van Genuchten/Mualem retention
// and Bishop effective stress (chi = S_eff). Both are time-dependent (reuse duration / time_steps).
// Dynamic: time-history (seismic) analysis -- M u'' + C u' + K u = -M r a_g(t) integrated by Newmark-β
// with Rayleigh damping (PLAXIS "Dynamic"). Linear-elastic skeleton (v1); horizontal base acceleration
// (site response / SSI). Boundaries: absorbing base + free-field sides (docs/references/dynamic-seismic-
// formulation.md). Enum VALUES are file-stable (appended); names array is indexed by the enum value.
enum class PhaseType { Plastic, Safety, Consolidation, TransientFlow, FullyCoupled, Dynamic };
inline const char* const* phase_type_names() {
    static const char* n[] = {"Plastic (staged construction)", "Safety (phi-c reduction)",
                              "Consolidation (time-dependent)", "Transient groundwater flow",
                              "Fully coupled flow-deformation", "Dynamic (seismic time-history)"};
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
    // Groundwater conditions for THIS phase (PLAXIS "water conditions per phase"). When
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
    // COMPLIANT (absorbing) BASE (Dynamic phase; Joyner & Chen 1975 -- PLAXIS Sci 6.3.2, factor 2).
    // false = rigid base (the user's fixed bottom; the historical bit-identical path; the base
    // reflects all downgoing energy). true = the bottom boundary absorbs outgoing waves through
    // Lysmer dashpots (rho Vs of the deepest layer = the halfspace continues it) and the input is
    // applied THERE as the upward-propagating wave: KATAI takes the phase's a_g(t) as the bedrock
    // (within) motion and applies HALF of it internally (Tutorial 17.8.5 convention). The solve is
    // then in TOTAL motion (the base moves). Combinable with free-field sides (the 1D side columns
    // solve on a compliant base too) and with dynamic_nonlinear (linear-limit identity verified).
    // Scope: horizontal SH only (base u_y stays fixed; disclosed in the phase message).
    bool seismic_compliant_base = false;
    // NONLINEAR dynamic solve (Dynamic phase): false = the LINEAR Newmark solve (default; the internal
    // force is K u, so the soil stays elastic and interfaces/geogrids/anchors do not slip/yield during
    // shaking -- fast, and the historical bit-identical path). true = the fully NONLINEAR Newmark solve
    // (per-step Newton on f_int(u) from the constitutive return mapping, so the soil can plastify and the
    // structures can slip/yield DURING the earthquake, like PLAXIS Dynamics). Much slower (the tangent is
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
    // SumMstage TARGET (PLAXIS "Σ Mstage"): the fraction of this phase's staged change that is
    // actually applied. 1 = the whole stage, which is what staged construction normally means
    // and what every project written before this field said. A smaller value applies part of it
    // and leaves the rest -- half an excavation lift, a partly built embankment -- and the next
    // phase continues from that partly-changed state.
    //
    // It applies to a STAGED (chained) phase only. The initial phase establishes the in-situ
    // state; scaling gravity there would not be a partial construction step, it would be a
    // different planet, so the value is refused rather than quietly obeyed.
    double sum_mstage = 1.0;
    // IGNORE UNDRAINED BEHAVIOUR (PLAXIS "Ignore und. behaviour (A,B)"): for this phase, materials
    // whose drainage is Undrained (A) or (B) are treated as drained -- no excess pore pressure is
    // generated. Strength parameters are untouched, so an Undrained (B) material still carries its
    // c_u with phi = 0. The standard use is a phase where the undrained response is not the
    // question being asked (establishing an initial state, or a long-term stage), which PLAXIS
    // supports for exactly the same reason.
    bool ignore_undrained = false;
    // NUMERICAL CONTROLS for this phase (PLAXIS "Numerical control parameters"). Zero means
    // "let the program choose", which is what every project that never touches them says, so
    // the defaults stay exactly where they are: derived from the material class (Hardening Soil
    // at a PLAXIS-realistic 1%, Mohr-Coulomb at 1e-6, a linear problem at 1e-10) or, in a
    // Safety phase, from the strength-reduction search's own trial settings.
    //
    // They belong in the FILE and not only in the program because a published number has to be
    // reproducible by someone else: "the settlement is 42 mm" is a claim about a model AND about
    // the numerics it was solved with, and if the second half cannot be written down, nobody can
    // re-run it. Reviewers may also legitimately ask whether an answer moved because the physics
    // moved or because the stopping rule did -- a question only a file that carries the stopping
    // rule can answer.
    //
    // `load_steps` is the number of increments the load is split into. This is NOT PLAXIS's "Max
    // steps": PLAXIS steps automatically and caps the count, KATAI splits the load into a fixed
    // number of increments (the solver still cuts back adaptively when one will not converge).
    // The field is named for what it does rather than for the box it resembles.
    double tolerance = 0.0;     // tolerated relative force residual (PLAXIS "Tolerated error")
    int load_steps = 0;         // load increments for this phase
    int max_iterations = 0;     // Newton iterations per increment (PLAXIS "Max iterations")
    // Active flags per project object, by index (kept in sync by the GUI; an entry missing
    // because the vector is short counts as ACTIVE -- new objects default to active).
    std::vector<char> poly_active, struct_active, load_active, disp_active;
    bool active_poly(size_t i) const { return i >= poly_active.size() || poly_active[i]; }
    bool active_struct(size_t i) const { return i >= struct_active.size() || struct_active[i]; }
    bool active_load(size_t i) const { return i >= load_active.size() || load_active[i]; }
    bool active_disp(size_t i) const { return i >= disp_active.size() || disp_active[i]; }
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
    // INFLOW POSITIVE -- the sign the PLAXIS Scientific Manual gives its wells in the same
    // chapter ("the source term is positive for a recharge well", §3.2.7). It is a specific
    // discharge normal to the boundary: the flow solve integrates it over the edge as the
    // manual's boundary term q (Eqs. 3-31, 3-34). Flow calculations only -- the manual states
    // that a CONSOLIDATION analysis cannot carry a non-zero prescribed outflow (Ch. 4).
    std::vector<double> edge_flux;
    // Local mesh density of the region (PLAXIS Coarseness factor): the target element size
    // inside this polygon is multiplied by it (1 = global, 0.5 = twice as fine, 2 = coarser).
    double coarseness = 1.0;
};

// A line prescribed displacement (PLAXIS "Prescribed displacement"): every mesh node on
// the segment gets the set components IMPOSED (ramped 0 -> value over the phase, like a
// load), the unset components stay free. A set component with value 0 is a rigid support
// line. Activated per phase exactly like loads (Phase::disp_active); the classic use is a
// displacement-controlled rigid footing, whose force is then read from the reactions.
struct PrescribedDisp {
    std::string name = "Displacement";
    double x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;   // the line [m]
    bool set_ux = false; double ux = 0.0;             // horizontal component [m]
    bool set_uy = true;  double uy = 0.0;             // vertical component [m]
    double coarseness = 1.0;                          // local mesh density factor (like loads)
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
};

struct Project {
    std::string name = "Untitled project";
    std::vector<Material> materials;            // soil & interfaces
    std::vector<PlateMaterial> plates;          // plate material sets
    std::vector<AnchorMaterial> anchors;        // anchor material sets
    std::vector<GeogridMaterial> geogrids;      // geogrid material sets
    std::vector<EmbeddedBeamMaterial> embedded; // embedded beam material sets
    std::vector<SoilPolygon> polygons;          // explicit soil regions (user-drawn)
    std::vector<StructElement> structs;         // structural elements (Structures mode)
    std::vector<Load> loads;                    // external loads
    std::vector<PrescribedDisp> disps;          // line prescribed displacements
    // Staged construction phases AFTER the initial phase (empty = classic single-phase solve).
    // phases[0] runs from the initial state, phases[k] from phases[k-1]'s committed stresses.
    std::vector<Phase> phases;
    // Initial-phase activation (PLAXIS: initial phase = soil only by default; structures and
    // loads are installed in later phases). Same index semantics as Phase.
    Phase initial;
    BoundaryConditions bc;                      // deformation boundary conditions (auto default)
    bool has_water = true;                      // groundwater level (a polyline, may be sloped)
    std::vector<double> wx, wy;                 // water polyline
    double x_min = 0.0, x_max = 40.0;           // model horizontal extent (BC edges)
    double y_min = 0.0, y_max = 20.0;           // model vertical extent (BC edges)
    // Analysis mode (PLAXIS "Model"): plane strain (default) or axisymmetric. In axisymmetric mode the
    // x-coordinate is the RADIUS r (x=0 is the symmetry axis) and the integration is r-weighted; used
    // for circular footings, piles, shafts and cylinders. (v1: soil-only, dry, no structural elements.)
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

} // namespace katai::model
