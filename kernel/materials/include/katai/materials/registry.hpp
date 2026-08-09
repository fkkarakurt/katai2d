#pragma once
// The constitutive catalogue: every material model this program can construct,
// registered against one interface and resolved by canonical name.
//
// This is the construction/validation seam only. The hot path is untouched on
// purpose: integration stays a free function over the tagged MaterialModel
// (integrate_point and its switch -- see material_model.hpp for why there is
// no vtable there). What the registry owns is everything that happens *before*
// the first Newton iteration: turning neutral input parameters into a
// MaterialModel, and refusing the combinations a model cannot honestly solve.
//
// Two rules shape the interface:
//   - A model is resolved by name, and an unknown name resolves to nothing.
//     The caller must refuse it by name -- never substitute a different model,
//     because composition must not be able to change an answer (R2).
//   - A registered entry is never replaced. An external model (the planned
//     user-material ABI implements exactly this interface from outside the
//     tree) cannot masquerade as a shipped one.
//
// The materials module depends on nothing else in the tree, so the input here
// is a neutral parameter block, not the project schema: the caller that owns
// the schema (the katai/jobs driver) copies its fields in.

#include <string>
#include <string_view>
#include <vector>

#include <katai/materials/material_model.hpp>

namespace katai::core {

// Drainage condition of the constitutive input, independent of the project
// schema's enum (the materials module cannot see katai/model). The caller maps
// its schema value; the mapping is a plain switch at the seam.
// UndrainedC = a TOTAL stress analysis (MMM section 2.7): undrained stiffness and undrained
// strength, no pore pressure at all. It shares nothing with (A)/(B) but its name -- there is
// no Kw/n, because there is no separation of water from skeleton to make one for.
enum class DrainageClass { Drained, UndrainedA, UndrainedB, NonPorous, UndrainedC };

// Neutral constitutive input. Angles are in radians -- the caller converts
// once, at the seam, so every entry sees the same convention. Fields a model
// does not read are simply ignored by its build function.
struct MaterialParams {
    double E = 0.0;              // Young's modulus E' [kN/m2]
    double nu = 0.0;             // Poisson's ratio nu'
    double c = 0.0;              // cohesion c' (Undrained (B): su)
    double phi_rad = 0.0;        // friction angle phi' [rad]
    double psi_rad = 0.0;        // dilatancy angle psi [rad]
    bool tension_cutoff = false; // Rankine tension cut-off (MMM Eq 3-11)
    double tensile_strength = 0.0;  // sigma_t [kN/m2]; negative input clamps to 0
    // Dilatancy cut-off (MMM Eq 5.16b): a dilating soil arrives at a critical void ratio where
    // dilatancy ends. e_init is the in-situ void ratio, e_max the critical one; when the volume
    // change has taken the soil to e_max the mobilised dilatancy angle is set to zero.
    bool dilatancy_cutoff = false;
    double e_init = 0.5, e_max = 1.0;

    // Hardening Soil (+ HSsmall) stiffness law (MMM sections 6-7).
    double E50_ref = 0.0, Eur_ref = 0.0, Eoed_ref = 0.0;
    double m = 0.0, p_ref = 0.0, Rf = 0.0, nu_ur = 0.0;
    double G0_ref = 0.0, gamma07 = 0.0;   // HSsmall only

    // Soft Soil (Creep) modified indices (MMM sections 10-11).
    double lam_star = 0.0, kap_star = 0.0, mu_star = 0.0;

    // K0^NC memory shared by the HS cap and the SS/SSC preconsolidation law:
    // auto = Jaky 1 - sin(phi'), otherwise the given value.
    bool k0nc_auto = true;
    double k0nc = 0.0;

    // Pore-fluid stiffness for Undrained (A)/(B) (MMM section 2.4). Either the equivalent
    // undrained Poisson ratio is given directly (skempton_mode = false, PLAXIS default 0.495)
    // or Skempton's B is, and nu_u follows from Eq. 2-55. Kw/n comes from Eq. 2-50 either way.
    bool skempton_mode = false;
    double nu_u = 0.495;
    double skempton_B = 0.0;

    DrainageClass drainage = DrainageClass::Drained;
};

// One registered constitutive model.
struct ModelEntry {
    std::string name;      // canonical name, e.g. "HardeningSoil"
    MaterialType type;     // the tag integrate_point dispatches on
    bool nonlinear;        // drives load stepping and solver selection
    // Solver step/tolerance family: the hardening family (HS, HSsmall) needs
    // PLAXIS-realistic tolerances, the soft-soil family (SS, SSC) additionally
    // the FD-tangent step class. Read by the driver exactly like its former
    // has_hardening / has_softsoil flags.
    bool hardening_family;
    bool softsoil_family;
    // Which depth-gradient inputs the model actually reads (PLAXIS E'_inc /
    // c'_inc about y_ref): E'_inc only means anything where E' itself is an
    // input (Linear elastic, Mohr-Coulomb -- Hardening Soil derives stiffness
    // from E50/Eoed/Eur plus its own stress dependency), and c'_inc only where
    // the Mohr-Coulomb strength is. profile_builder.hpp gates on these, so a
    // gradient a model cannot honour never silently reaches its profile.
    bool profile_E = false;
    bool profile_c = false;
    // Refuse a parameter/drainage combination this model cannot honestly
    // solve. Empty string = acceptable; anything else is shown to the user
    // verbatim as the reason the whole solve is refused.
    std::string (*validate)(const MaterialParams&);
    // Construct the MaterialModel. Cannot fail: everything refusable is
    // refused by validate(), and construction is plain arithmetic.
    MaterialModel (*build)(const MaterialParams&);
};

// Resolve a model by canonical name; nullptr when nothing of that name is
// registered. The returned pointer stays valid for the life of the process --
// later registrations never move existing entries.
const ModelEntry* find_model(std::string_view name);

// Registration order, built-ins first. A snapshot: names registered after the
// call do not appear in it.
std::vector<std::string_view> model_names();

// The seam an external model enters through. Returns false -- and registers
// nothing -- when the name is empty, already taken, or either function is
// missing. Not thread-safe by design: registration belongs to program setup,
// before any solve starts.
bool register_model(ModelEntry entry);

} // namespace katai::core
