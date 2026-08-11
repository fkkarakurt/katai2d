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

// The two effective-stress undrained types: the ones that add a pore fluid to the skeleton.
// Undrained (C) is deliberately NOT one of them -- it has no pore fluid to add.
bool is_undrained(const MaterialParams& p) {
    return p.drainage == DrainageClass::UndrainedA ||
           p.drainage == DrainageClass::UndrainedB;
}

// Undrained strength entered directly: a Tresca envelope, phi = psi = 0 with c = su. (B) reaches
// it through an effective-stress analysis, (C) through a total-stress one; the envelope is the
// same and so is this substitution.
bool is_su_entered(const MaterialParams& p) {
    return p.drainage == DrainageClass::UndrainedB || p.drainage == DrainageClass::UndrainedC;
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
        // The equivalent undrained Poisson ratio: entered directly, or derived from Skempton's
        // B (MMM Eq. 2-55). It used to be the constant 0.495 for every material in every model,
        // which is PLAXIS's DEFAULT rather than its rule -- a soil with a measured B of 0.90 was
        // silently solved at B = 0.978, and the two differ by more than the excess pore pressure
        // a designer would call negligible.
        mm.undrained_poisson = p.skempton_mode
                                   ? undrained_poisson_from_skempton(p.skempton_B, p.nu)
                                   : p.nu_u;
    }
    // Undrained (C): total stress throughout. E and nu are already the undrained pair -- the
    // user entered them in those boxes, as PLAXIS asks -- so nothing about the elasticity
    // changes here; what the flag carries is that the result is total stress and there are no
    // pore pressures to separate out.
    mm.total_stress = p.drainage == DrainageClass::UndrainedC;
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
    // Undrained (B) and (C): c already holds su; a Tresca envelope is phi = psi = 0.
    // This override exists only for Mohr-Coulomb -- see validate_hs below for
    // why HS refuses the same request.
    if (is_su_entered(p)) {
        mm.friction_angle = 0.0;
        mm.dilatancy_angle = 0.0;
    }
    return mm;
}

// --- Hardening Soil (+ HSsmall) ----------------------------------------------

std::string validate_hs(const MaterialParams& p) {
    // Undrained (C) is offered by PLAXIS for the Linear Elastic and Mohr-Coulomb models (and
    // for NGI-ADP / UDCAM-S, which this build does not have) -- MMM section 2.7.1. A total
    // stress analysis with a hardening model would run: the model would happily take undrained
    // parameters as effective ones and harden along the wrong stress path.
    if (p.drainage == DrainageClass::UndrainedC)
        return "Hardening Soil with Undrained (C) is not a supported combination: a total "
               "stress analysis needs undrained stiffness and strength, and the hardening laws "
               "are written for effective stress (PLAXIS offers Undrained (C) for the Linear "
               "Elastic and Mohr-Coulomb models). Use Undrained (A) or (B) with this model, or "
               "Mohr-Coulomb with Eu, nu_u and su for a total stress analysis.";
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
    // The pore fluid follows THIS model's elasticity. Hardening Soil's elastic pair is
    // (Eur_ref, nu_ur); E and nu belong to the Linear-elastic/Mohr-Coulomb boxes, which an HS
    // data set never fills and the HS integrator never reads. Until this line existed, an
    // undrained HS material was given Kw/n = f(E_box, nu_box) -- a pore fluid sized by a
    // default, on the model most likely to be used for soft clay. nu' also enters Skempton's
    // conversion, so a B entered by the user has to be re-resolved against nu_ur.
    if (mm.undrained) {
        mm.undrained_E_ref = p.Eur_ref;
        mm.undrained_nu_ref = p.nu_ur;
        if (p.skempton_mode)
            mm.undrained_poisson = undrained_poisson_from_skempton(p.skempton_B, p.nu_ur);
    }
    return mm;
}

MaterialModel build_hs(const MaterialParams& p) {
    MaterialModel mm = build_hs_core(p);
    hs_calibrate_cap(mm.hs, k0nc_of(p));   // cap_alpha/beta from K0^NC + Eoed
    return mm;
}

MaterialModel build_hss(const MaterialParams& p) {
    MaterialModel mm = build_hs_core(p);
    mm.hs.gamma07 = p.gamma07;
    // Small-strain overlay (MMM ch. 7), capped at the ratio the model permits (sec. 7.5:
    // G0/Gur <= 20). The driver says so when the cap bites (K2D-M004) -- a G0 quietly reduced
    // is a different soil from the one the file asked for.
    mm.hs.G0_ref = std::min(p.G0_ref, mm.hs.G0_ref_cap());
    hs_calibrate_cap(mm.hs, k0nc_of(p));
    return mm;
}

// --- Soft Soil (Creep) --------------------------------------------------------

std::string validate_ss(const MaterialParams& p) {
    if (p.drainage == DrainageClass::UndrainedC)
        return "Soft Soil (Creep) with Undrained (C) is not a supported combination: a total "
               "stress analysis needs an undrained stiffness, and this model's stiffness is the "
               "stress-dependent ln law (K = p'/kappa*) written for EFFECTIVE stress (PLAXIS "
               "offers Undrained (C) for the Linear Elastic and Mohr-Coulomb models). Use "
               "Drained, or Mohr-Coulomb with Eu, nu_u and su.";
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
