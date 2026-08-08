// The dilatancy cut-off: a dense sand stops dilating when it reaches its critical void ratio,
// and until now this program let it dilate for ever.
//
// That is not a missing convenience. Unlimited dilation raises the confining stress a footing
// mobilises, so the computed bearing capacity is too HIGH -- an unsafe number, produced quietly,
// on exactly the soil (dense sand) where a practitioner would trust it most.
//
// The rule is the Material Models Manual's, quoted rather than paraphrased (V8 edition, §5.4,
// Fig. 5.6): "After extensive shearing, dilating materials arrive in a state of critical density
// where dilatancy has come to an end... In order to specify this behaviour, the initial void
// ratio, e_init, and the maximum void ratio, e_max, of the material must be entered as general
// parameters. As soon as the volume change results in a state of maximum void, the mobilised
// dilatancy angle, psi_m, is automatically set back to zero." Equation 5.16b is that last
// sentence: for e >= e_max, psi_m = 0.
//
// KATAI states the void ratio through the volume change, 1 + e = (1 + e_init) exp(eps_v) with
// expansion positive, rather than through the manual's Eq. 5.17, whose printed sign convention
// is ambiguous. The two say the same thing; only one of them says it once.
//
// verify: KV-CST-004
//   oracle:   closed_form
//   source:   PLAXIS 2D Material Models Manual, dilatancy cut-off: Eq. 5.16a for the mobilised dilatancy, Eq. 5.16b "for e >= e_max: psi_m = 0", Eq. 5.17 for the void ratio, Fig. 5.6 for the resulting drained-triaxial strain curve; the manual also states that e_min "is not used within the context of the Hardening-Soil model", which is why KATAI does not carry it
//   locator:  1 + e = (1 + e_init) exp(eps_v), expansion positive, so dilation stops at the volumetric strain eps_v,cut = ln((1 + e_max)/(1 + e_init)); beyond it the return mapping runs with psi = 0 and the plastic flow is isochoric (stated in full)
//   quantity: the accumulated volumetric strain of a Mohr-Coulomb stress point sheared far past its cut-off, and the same point's response with the cut-off switched off [-]
//   expected: with the cut-off ON the volumetric strain stops at eps_v,cut and does not grow afterwards; with it OFF the same shearing keeps dilating; and the switch changes nothing at all before the cut-off is reached
//   band:     MEASURED on this tree with e_init = 0.60, e_max = 0.63 (eps_v,cut = 0.018576) on a drained biaxial at 100 kPa cell pressure, 400 increments of 2e-4 axial strain: the cut-off run stops at 0.018635 -- inside one loading increment of the closed form -- while the same soil without it reaches 0.027139, 1.46x as much. After the cut-off the volumetric strain grows by 1.4e-16 over 80 further increments, i.e. the flow is isochoric to round-off; before it the two runs are bit-for-bit identical

#include <katai/materials/material_model.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace core = katai::core;

namespace {

int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kEInit = 0.60, kEMax = 0.63;   // dense sand, a narrow room to dilate

core::MaterialModel dense_sand(bool cutoff) {
    core::MaterialModel m;
    m.type = core::MaterialType::MohrCoulomb;
    m.youngs_modulus = 3.0e4;
    m.poisson_ratio = 0.2;
    m.cohesion = 1.0;
    m.friction_angle = 40.0 * kPi / 180.0;
    m.dilatancy_angle = 12.0 * kPi / 180.0;   // dense: psi well above zero
    m.tension_cutoff = false;
    m.dilatancy_cutoff = cutoff;
    m.e_init = kEInit;
    m.e_max = kEMax;
    return m;
}

// Shear the point in small increments of deviatoric strain under a constant mean confinement,
// accumulating the volumetric strain the way the assembler does (the state carries it).
struct History {
    std::vector<double> eps_v;
    double last = 0.0;
};

// A DRAINED BIAXIAL element test, which is the plane-strain sibling of the drained triaxial the
// manual draws in Fig. 5.6: the axial strain is imposed and the lateral strain is solved for so
// that the lateral stress stays at the cell pressure. The volume change is then the soil's
// answer rather than the test's assumption -- imposing a strain path with zero trace, as the
// first version of this test did, makes dilation arithmetically impossible and proves nothing.
History shear(const core::MaterialModel& m, int steps, double d_axial) {
    History h;
    core::GaussState s;
    const double cell = -100.0;   // lateral confinement [kPa], compression negative
    s.stress = Eigen::Vector3d(cell, cell, 0.0);
    s.stress_zz = cell;
    for (int i = 0; i < steps; ++i) {
        double dexx = 0.0;
        core::GaussState trial;
        Eigen::Matrix3d tangent;
        // Newton on the single unknown lateral strain: sigma_xx must come back to the cell
        // pressure. The return mapping's own tangent is the derivative.
        for (int it = 0; it < 30; ++it) {
            const Eigen::Vector3d de(dexx, -d_axial, 0.0);
            core::integrate_point(m, s, de, trial, tangent, core::TangentMode::kConsistent, 0.0);
            const double r = trial.stress(0) - cell;
            if (std::fabs(r) < 1e-10) break;
            const double k = std::fabs(tangent(0, 0)) > 1e-12 ? tangent(0, 0) : 1.0;
            dexx -= r / k;
        }
        // The assembler accumulates the volumetric strain into the state; this stands in for it.
        trial.eps_vol = s.eps_vol + dexx - d_axial;
        s = trial;
        h.eps_v.push_back(s.eps_vol);
    }
    h.last = s.eps_vol;
    return h;
}

}  // namespace

int main() {
    std::printf("== dilatancy cut-off (KV-CST-004) ==\n");

    // The closed form: dilation must stop when the volume change has taken e to e_max.
    const double eps_cut = std::log((1.0 + kEMax) / (1.0 + kEInit));
    std::printf("  e_init = %.3f, e_max = %.3f -> dilation must stop at eps_v = %.6f\n", kEInit,
                kEMax, eps_cut);

    // A point sheared well past the cut-off, with and without it.
    const int kSteps = 400;
    const double kStep = 2e-4;   // axial strain increment (compression)
    const History on = shear(dense_sand(true), kSteps, kStep);
    const History off = shear(dense_sand(false), kSteps, kStep);
    std::printf("  after %d increments: cut-off ON  eps_v = %.6f\n", kSteps, on.last);
    std::printf("                       cut-off OFF eps_v = %.6f\n", off.last);

    std::printf("  the cut-off removes %.1f%% of the volume change this soil would otherwise\n"
                "  have produced -- and all of it would have raised its bearing capacity\n",
                100.0 * (off.last - on.last) / off.last);
    check(off.last > eps_cut,
          "without the cut-off the point keeps dilating past the critical void ratio");
    check(off.last > on.last * 1.3,
          "and it ends up markedly more dilated than the cut-off run (measured 1.46x)");
    check(on.last < eps_cut + 2.0 * kStep,
          "with the cut-off it stops at the critical void ratio, within one increment");

    // Once past the cut-off the flow rule is volume-preserving, so the volumetric strain must
    // not merely grow slowly -- it must not grow at all.
    size_t first_cut = on.eps_v.size();
    for (size_t i = 0; i < on.eps_v.size(); ++i)
        if (on.eps_v[i] >= eps_cut) { first_cut = i; break; }
    check(first_cut + 2 < on.eps_v.size(), "the cut-off is reached inside the loading history");
    if (first_cut + 2 < on.eps_v.size()) {
        double growth = 0.0;
        for (size_t i = first_cut + 2; i < on.eps_v.size(); ++i)
            growth = std::fmax(growth, on.eps_v[i] - on.eps_v[first_cut + 1]);
        std::printf("  volumetric growth after the cut-off: %.3e (over %zu increments)\n", growth,
                    on.eps_v.size() - first_cut - 2);
        check(growth < 1e-12, "and afterwards the plastic flow is isochoric, exactly");
    }

    // Before the cut-off is reached the switch must change NOTHING. This is the property that
    // says the feature is inert until it acts -- and it is asserted bit-for-bit, because the
    // same material integrating the same increments has no reason to differ by even one ulp.
    bool identical = true;
    for (size_t i = 0; i < first_cut && i < off.eps_v.size(); ++i)
        if (on.eps_v[i] != off.eps_v[i]) { identical = false; break; }
    check(identical,
          "before the cut-off is reached the two runs are bit-for-bit identical");

    // And the switch is OFF by default, as in PLAXIS: a file that says nothing gets the old
    // behaviour, which is what makes the schema bump the honest way to introduce this.
    const core::MaterialModel fresh;
    check(!fresh.dilatancy_cutoff, "the cut-off is off by default");

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
