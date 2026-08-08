// Built-in constitutive catalogue. The build functions reproduce, arithmetic
// operation for arithmetic operation, the construction the analysis driver
// performed inline before the registry existed -- the suite's bit-identical
// gate rests on that. Change a formula here and a validated number moves.
#include <katai/materials/registry.hpp>

#include <algorithm>
#include <cmath>
#include <deque>

#include <katai/materials/hardening_soil.hpp>

namespace katai::core {
namespace {

bool is_undrained(const MaterialParams& p) {
    return p.drainage == DrainageClass::UndrainedA ||
           p.drainage == DrainageClass::UndrainedB;
}

// Fields every model shares. The undrained flag wires the pore fluid's bulk
// stiffness Kw/n into the global tangent (effective stress path); nu_u = 0.495
// is the PLAXIS default (exactly 0.5 is singular).
MaterialModel common_fields(const MaterialParams& p) {
    MaterialModel mm;
    mm.youngs_modulus = p.E;
    mm.poisson_ratio = p.nu;
    mm.cohesion = p.c;
    mm.friction_angle = p.phi_rad;
    mm.dilatancy_angle = p.psi_rad;
    // Rankine tension cut-off (MMM Eq 3-11; default ON in the schema,
    // sigma_t = 0). Applied in the Mohr-Coulomb return mapping; the HS/SS
    // integrators do not read it yet, and the editor says so honestly.
    mm.tension_cutoff = p.tension_cutoff;
    mm.tensile_strength = std::max(0.0, p.tensile_strength);
    // Dilatancy cut-off (MMM Eq. 5.16b): read by the Mohr-Coulomb and Hardening Soil return
    // mappings, which are the models that HAVE a dilatancy angle to switch off.
    mm.dilatancy_cutoff = p.dilatancy_cutoff;
    mm.e_init = p.e_init;
    mm.e_max = p.e_max;
    if (is_undrained(p)) {
        mm.undrained = true;
        mm.undrained_poisson = 0.495;
    }
    return mm;
}

// Jaky's K0^NC = 1 - sin(phi') unless the user pinned a value. Shared by the
// HS cap calibration and the SS/SSC preconsolidation law.
double k0nc_of(const MaterialParams& p) {
    return p.k0nc_auto ? (1.0 - std::sin(p.phi_rad)) : p.k0nc;
}

std::string accept(const MaterialParams&) { return {}; }

// --- Linear elastic ---------------------------------------------------------

MaterialModel build_le(const MaterialParams& p) {
    MaterialModel mm = common_fields(p);
    mm.type = MaterialType::LinearElastic;
    return mm;
}

// --- Mohr-Coulomb -----------------------------------------------------------

MaterialModel build_mc(const MaterialParams& p) {
    MaterialModel mm = common_fields(p);
    mm.type = MaterialType::MohrCoulomb;
    // Undrained (B): c already holds su; a Tresca envelope is phi = psi = 0.
    // This override exists only for Mohr-Coulomb -- see validate_hs below for
    // why HS refuses the same request.
    if (p.drainage == DrainageClass::UndrainedB) {
        mm.friction_angle = 0.0;
        mm.dilatancy_angle = 0.0;
    }
    return mm;
}

// --- Hardening Soil (+ HSsmall) ----------------------------------------------

std::string validate_hs(const MaterialParams& p) {
    if (p.drainage == DrainageClass::UndrainedB)
        return "Hardening Soil with Undrained (B) is not supported yet (the su/Tresca "
               "override is only implemented for Mohr-Coulomb; HS would silently follow "
               "its effective envelope like Undrained (A)). Use Undrained (A) with "
               "effective strength, or a Mohr-Coulomb set with c = su for Undrained (B).";
    return {};
}

MaterialModel build_hs_core(const MaterialParams& p) {
    MaterialModel mm = common_fields(p);
    mm.type = MaterialType::HardeningSoil;
    auto& h = mm.hs;
    h.E50_ref = p.E50_ref;
    h.Eur_ref = p.Eur_ref;
    h.Eoed_ref = p.Eoed_ref;
    h.m = p.m;
    h.p_ref = p.p_ref;
    h.cohesion = p.c;
    h.friction = p.phi_rad;
    h.dilatancy = p.psi_rad;
    h.Rf = p.Rf;
    h.nu_ur = p.nu_ur;
    return mm;
}

MaterialModel build_hs(const MaterialParams& p) {
    MaterialModel mm = build_hs_core(p);
    hs_calibrate_cap(mm.hs, k0nc_of(p));   // cap_alpha/beta from K0^NC + Eoed
    return mm;
}

MaterialModel build_hss(const MaterialParams& p) {
    MaterialModel mm = build_hs_core(p);
    mm.hs.G0_ref = p.G0_ref;               // small-strain overlay (MMM ch. 7)
    mm.hs.gamma07 = p.gamma07;
    hs_calibrate_cap(mm.hs, k0nc_of(p));
    return mm;
}

// --- Soft Soil (Creep) --------------------------------------------------------

std::string validate_ss(const MaterialParams& p) {
    // Honest refusal: the pore-fluid stiffness plumbing derives Kw/n from a
    // constant E, and Soft Soil has no E input (stiffness is the ln law
    // K = p'/kappa*) -- wiring it through would set Kw/n silently to zero and
    // an "undrained" analysis would behave drained (the silent-wrong class).
    if (is_undrained(p))
        return "Soft Soil (Creep) with Undrained (A/B) drainage is not supported yet: the "
               "pore-fluid stiffness plumbing derives Kw/n from a constant E, and Soft Soil's "
               "stiffness is stress-dependent (K = p'/kappa*) -- it would silently behave "
               "drained. Use Drained (or model the short term with a Consolidation phase).";
    return {};
}

MaterialModel build_ss(const MaterialParams& p) {
    // Stiffness comes from the ln law (K_ur = p'/kappa*); youngs_modulus is
    // NOT read by this model (it survives only in the driver's linear-elastic
    // fallback table). M is derived internally from K0^NC (Brinkgreve 1994).
    MaterialModel mm = common_fields(p);
    mm.type = MaterialType::SoftSoil;
    auto& s = mm.ssoil;
    s.lam_star = p.lam_star;
    s.kap_star = p.kap_star;
    s.nu_ur = p.nu_ur;
    s.c = p.c;
    s.phi = p.phi_rad;
    s.psi = p.psi_rad;
    s.K0nc = k0nc_of(p);
    return mm;
}

MaterialModel build_ssc(const MaterialParams& p) {
    // Soft Soil Creep (MMM section 11): SS parameters plus mu*; time enters
    // integration through integrate_point's dt_day tail, tau = 1 day (PLAXIS
    // constant). Same step/tolerance family and Safety gate as SS.
    MaterialModel mm = common_fields(p);
    mm.type = MaterialType::SoftSoilCreep;
    auto& s = mm.ssc;
    s.lam_star = p.lam_star;
    s.kap_star = p.kap_star;
    s.mu_star = p.mu_star;
    s.nu_ur = p.nu_ur;
    s.c = p.c;
    s.phi = p.phi_rad;
    s.psi = p.psi_rad;
    s.K0nc = k0nc_of(p);
    return mm;
}

// A deque so registration never moves an existing entry: find_model hands out
// pointers, and they stay valid for the life of the process.
std::deque<ModelEntry>& table() {
    static std::deque<ModelEntry> entries = {
        {"LinearElastic", MaterialType::LinearElastic, false, false, false, true, false, accept, build_le},
        {"MohrCoulomb", MaterialType::MohrCoulomb, true, false, false, true, true, accept, build_mc},
        {"HardeningSoil", MaterialType::HardeningSoil, true, true, false, false, false, validate_hs, build_hs},
        {"HSsmall", MaterialType::HardeningSoil, true, true, false, false, false, validate_hs, build_hss},
        {"SoftSoil", MaterialType::SoftSoil, true, false, true, false, false, validate_ss, build_ss},
        {"SoftSoilCreep", MaterialType::SoftSoilCreep, true, false, true, false, false, validate_ss, build_ssc},
    };
    return entries;
}

} // namespace

const ModelEntry* find_model(std::string_view name) {
    for (const ModelEntry& entry : table())
        if (entry.name == name) return &entry;
    return nullptr;
}

std::vector<std::string_view> model_names() {
    std::vector<std::string_view> names;
    names.reserve(table().size());
    for (const ModelEntry& entry : table()) names.push_back(entry.name);
    return names;
}

bool register_model(ModelEntry entry) {
    if (entry.name.empty() || !entry.validate || !entry.build) return false;
    if (find_model(entry.name)) return false;
    table().push_back(std::move(entry));
    return true;
}

} // namespace katai::core
