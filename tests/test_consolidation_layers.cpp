// A layered column consolidates the same way whatever order its materials are listed in.
//
// The consolidation and fully-coupled phases carry a storage term, n/Kw: the volume of water a
// soil stores per kPa of pore pressure, because water is compressible and the pores hold it. The
// porosity n belongs to each material, so the storage does too. Until 2026-09-13 both phases took
// ONE value, from the first material in the project's list, and gave it to every element. For a
// soft soil that is harmless (m_v dwarfs n/Kw); for a stiff layer it is not, and the tell was that
// the answer depended on the order of a list:
//
//   two 6 m layers, porosity 0.2 over 0.6, E_oed = 500 MPa, k = 0.1 m/day, 10 kPa surcharge --
//   the degree of consolidation moved by up to 6.8 percentage points, and the undrained pressure
//   in the lower layer from 8.69 to 9.52 kPa, when the two materials swapped places in the list.
//   At E_oed = 50 MPa the same swap moved U by 0.8 points; at 1 MPa by 0.02.
//
// Three things are checked, and they are three different kinds of evidence.
//
//  (1) THE CLOSED FORM, PER LAYER. A confined layer loaded undrained takes u0 = q / (1 + n Eoed/Kw).
//      After a first step far shorter than the drainage time, mid-height of each layer is still
//      undrained, so each layer must show its OWN u0 -- two different numbers in one column.
//  (2) THE INVARIANCE. The same model with the materials listed the other way round: the same
//      answer. A property the physics demands, and the one the defect broke.
//  (3) AN INDEPENDENT PATH for the time history. Two layers of different consolidation coefficient
//      have no short closed form, so this test solves the 1-D equation itself -- a finite-volume
//      column written below, sharing no code with the engine -- and the engine's settlement history
//      must follow it.
//
// verify: KV-CON-004
//   oracle:   independent_path
//   source:   one-dimensional consolidation with a compressible pore fluid: (m_v + n/Kw) du/dt = (k/gamma_w) d2u/dz2 under constant total stress (Terzaghi 1943; Biot 1941 for the storage term), and the undrained response of a laterally confined layer to a surcharge q, u0 = q / (1 + n Eoed / Kw), both stated in full; the two-layer time history from a cell-centred finite-volume solution of that equation (2400 cells, backward Euler at 1e-6 day) written in the test and sharing no code with the engine; KATAI 2D input contract (docs/k2d-format.md, e_init for the porosity n = e/(1+e))
//   locator:  a laterally confined weightless column 1 m wide and 12 m tall, drained at the top only, E_oed = 500 MPa (E with nu = 0) and k = 0.1 m/day in both layers; the upper 6 m has e = 0.25 (n = 0.2), the lower 6 m e = 1.5 (n = 0.6); Kw = 2e6 kPa (the engine's water) and gamma_w = 9.81 kN/m3; a 10 kPa surcharge applied at t = 0+ in a Consolidation phase and in a FullyCoupled phase; each run twice, with the two materials listed in either order
//   quantity: the excess pore pressure at mid-height of each layer after one 2e-5 day step [kPa]; the surface settlement at t = 0.002, 0.005, 0.01, 0.02 and 0.04 day of a 0.05 day, 1000-step phase [m]; and the largest difference between the two listing orders [relative]
//   expected: u0 = 10 / (1 + 0.2 x 500000 / 2e6) = 9.52381 kPa in the upper layer and 10 / (1 + 0.6 x 500000 / 2e6) = 8.69565 kPa in the lower, in both phase types; the settlement history of the finite-volume solution; the same answer, to round-off, whichever material is listed first
//   band:     0.2% on each layer's u0, as asserted below -- measured -0.013% (upper) and +0.000% (lower), identical in both phase types and both listing orders; 0.3% of the final settlement q H / Eoed on the settlement history -- measured -0.090% at 0.002 day falling to -0.019% at 0.04 day, where the same finite-volume column with ONE porosity for both layers sits +0.97% (n = 0.6) or -3.76% (n = 0.2) away -- computed in the test, which asserts that either single value would fall outside the band, so the band can say no; 1e-12 relative between the listing orders -- measured 0.0 (bit-identical) in both phase types. The first step is kept short because a single backward-Euler step drains the column to a depth of about sqrt(cv dt): at 1e-4 day it took 1.4% off the upper layer's mid-height, which is the step, not the storage
#include <katai/analysis/constants.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;
using katai::app::InitialPhase;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

constexpr double kW = 1.0, kH = 12.0, kHalf = 6.0;
constexpr double kEoed = 5.0e5;          // [kPa]; nu = 0 makes E the constrained modulus
constexpr double kPerm = 0.1;            // [m/day]
constexpr double kQ = 10.0;              // [kPa]
constexpr double kKw = 2.0e6;            // the engine's water bulk modulus [kPa]
constexpr double kUpperE = 0.25;         // void ratio of the upper layer -> n = 0.2
constexpr double kLowerE = 1.5;          // void ratio of the lower layer -> n = 0.6

double porosity(double e) { return e / (1.0 + e); }

// The column, with the materials listed dense-first or loose-first. Everything else is identical.
m::Project column(bool upper_listed_first, m::PhaseType type, double duration, int steps) {
    m::Project pr;
    const auto layer = [](const char* name, double e) {
        m::Material s;
        s.name = name;
        s.model = m::SoilModel::LinearElastic;
        s.E = kEoed; s.nu = 0.0; s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = e;
        s.kx = kPerm; s.ky = kPerm;
        return s;
    };
    const int upper = upper_listed_first ? 0 : 1;
    const int lower = 1 - upper;
    pr.materials.resize(2);
    pr.materials[upper] = layer("Upper (n = 0.2)", kUpperE);
    pr.materials[lower] = layer("Lower (n = 0.6)", kLowerE);

    const auto poly = [](int material, double y0, double y1, int bottom_bc, int top_bc,
                         int top_flow) {
        m::SoilPolygon P;
        P.material = material;
        P.x = {0, kW, kW, 0};
        P.y = {y0, y0, y1, y1};
        P.edge_bc = {bottom_bc, (int)m::BCType::HorizontallyFixed, top_bc,
                     (int)m::BCType::HorizontallyFixed};
        P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed, top_flow,
                       (int)m::FlowBCType::Closed};
        P.edge_head = {0.0, 0.0, kH, 0.0};
        return P;
    };
    pr.polygons.push_back(poly(lower, 0.0, kHalf, (int)m::BCType::FullyFixed,
                               (int)m::BCType::Free, (int)m::FlowBCType::Closed));
    pr.polygons.push_back(poly(upper, kHalf, kH, (int)m::BCType::Free, (int)m::BCType::Free,
                               (int)m::FlowBCType::Head));   // the top drains (p = 0)
    pr.has_water = false;   // pore = excess only

    m::Load L; L.kind = m::LoadKind::Distributed; L.name = "Surcharge";
    L.x1 = 0; L.y1 = kH; L.x2 = kW; L.y2 = kH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kQ;
    pr.loads.push_back(L);
    pr.initial.load_active = {0};           // off in the initial phase

    m::Phase ph; ph.name = "Consolidate";
    ph.type = type;
    ph.load_active = {1};                   // on at t = 0+
    ph.duration = duration;
    ph.time_steps = steps;
    pr.phases.push_back(ph);
    return pr;
}

struct Run {
    bool ok = false;
    std::string msg;
    std::vector<double> time, settlement, excess;
    katai::mesh::Mesh mesh;
};

Run solve(const m::Project& pr) {
    Run out;
    const auto M = katai::app::mesh_from_project(pr, 0.4, 6);
    if (!M.ok) { out.msg = "mesh"; return out; }
    const auto res = katai::app::solve_phases(pr, M.mesh, InitialPhase::K0Procedure);
    if (res.size() != 2 || !res[1].ok) {
        out.msg = res.size() == 2 ? res[1].message : std::string("phases");
        return out;
    }
    out.ok = true;
    out.time = res[1].consol_time;
    out.settlement = res[1].consol_settlement;
    out.excess = res[1].excess_pore;
    out.mesh = M.mesh;
    return out;
}

double pore_near(const Run& r, double x, double y) {
    int best = -1;
    double d2 = 1e300;
    for (int n = 0; n < r.mesh.node_count; ++n) {
        const double d = (r.mesh.x[n] - x) * (r.mesh.x[n] - x) + (r.mesh.y[n] - y) * (r.mesh.y[n] - y);
        if (d < d2) { d2 = d; best = n; }
    }
    return best >= 0 && best < (int)r.excess.size() ? std::fabs(r.excess[best]) : -1.0;
}

// ---------------------------------------------------------------------------------------------
// The independent path: cell-centred finite volumes on z in [0, H], z = 0 the closed base, z = H
// the drained top (u = 0 on the face). S_i du_i/dt = (k/gamma_w) (flux differences) / dz with
// S_i = m_v + n_i / Kw, backward Euler, the tridiagonal system solved by the Thomas algorithm.
// Settlement is the integral of m_v (q - u) over the column: the effective stress each cell has
// gained, times its compressibility.
// ---------------------------------------------------------------------------------------------
class FiniteVolumeColumn {
public:
    // n_lower / n_upper: each layer's porosity. The per-material column passes the two real ones;
    // passing one value for both builds the single-porosity model the engine used to solve.
    FiniteVolumeColumn(int cells, double dt_max, double n_lower, double n_upper)
        : n_(cells), dt_max_(dt_max) {
        dz_ = kH / n_;
        u_.resize(n_);
        S_.resize(n_);
        const double mv = 1.0 / kEoed;
        for (int i = 0; i < n_; ++i) {
            const double zc = (i + 0.5) * dz_;
            const double nporo = zc < kHalf ? n_lower : n_upper;
            S_[i] = mv + nporo / kKw;
            u_[i] = kQ * mv / S_[i];   // the undrained response, the same closed form as (1)
        }
    }
    void advance_to(double t_target) {
        while (t_ < t_target - 1e-15) {
            const double dt = std::min(dt_max_, t_target - t_);
            step(dt);
            t_ += dt;
        }
    }
    double settlement() const {
        double s = 0.0;
        for (int i = 0; i < n_; ++i) s += (kQ - u_[i]) / kEoed * dz_;
        return s;
    }

private:
    void step(double dt) {
        const double c = kPerm / katai::core::kGammaWater / (dz_ * dz_);
        std::vector<double> a(n_), b(n_), cc(n_), d(n_);
        for (int i = 0; i < n_; ++i) {
            const double r = dt * c / S_[i];
            double diag = 1.0;
            a[i] = 0.0; cc[i] = 0.0;
            if (i > 0) { a[i] = -r; diag += r; }
            if (i < n_ - 1) { cc[i] = -r; diag += r; }
            else diag += 2.0 * r;   // the drained face half a cell above the last centre
            b[i] = diag;
            d[i] = u_[i];
        }
        for (int i = 1; i < n_; ++i) {
            const double w = a[i] / b[i - 1];
            b[i] -= w * cc[i - 1];
            d[i] -= w * d[i - 1];
        }
        u_[n_ - 1] = d[n_ - 1] / b[n_ - 1];
        for (int i = n_ - 2; i >= 0; --i) u_[i] = (d[i] - cc[i] * u_[i + 1]) / b[i];
    }
    int n_;
    double dz_, dt_max_, t_ = 0.0;
    std::vector<double> u_, S_;
};

void test_layer_closed_form(m::PhaseType type, const char* label) {
    std::printf("\n-- (1) %s: each layer takes its OWN undrained pressure --\n", label);
    const double u_upper = kQ / (1.0 + porosity(kUpperE) * kEoed / kKw);
    const double u_lower = kQ / (1.0 + porosity(kLowerE) * kEoed / kKw);
    for (bool upper_first : {true, false}) {
        const Run r = solve(column(upper_first, type, 2.0e-5, 1));
        check(r.ok, "the one-step phase ran");
        if (!r.ok) { std::printf("   (%s)\n", r.msg.c_str()); continue; }
        const double pu = pore_near(r, 0.5 * kW, kHalf + 0.5 * kHalf);
        const double pl = pore_near(r, 0.5 * kW, 0.5 * kHalf);
        std::printf("   %s listed first: upper %.5f kPa (closed form %.5f, %+.3f%%) | lower %.5f "
                    "kPa (closed form %.5f, %+.3f%%)\n",
                    upper_first ? "upper" : "lower", pu, u_upper, 100.0 * (pu / u_upper - 1.0),
                    pl, u_lower, 100.0 * (pl / u_lower - 1.0));
        check(std::fabs(pu / u_upper - 1.0) < 0.002,
              "the upper layer (n = 0.2) holds q / (1 + n Eoed / Kw) within 0.2%");
        check(std::fabs(pl / u_lower - 1.0) < 0.002,
              "the lower layer (n = 0.6) holds its own, different, value within 0.2%");
    }
}

void test_order_and_history(m::PhaseType type, const char* label, bool with_history) {
    std::printf("\n-- (2)%s %s: the listing order does not matter%s --\n",
                with_history ? " + (3)" : "", label,
                with_history ? ", and the history follows the independent solution" : "");
    const Run a = solve(column(true, type, 0.05, 1000));
    const Run b = solve(column(false, type, 0.05, 1000));
    check(a.ok && b.ok, "both listing orders ran");
    if (!a.ok || !b.ok) { std::printf("   (%s | %s)\n", a.msg.c_str(), b.msg.c_str()); return; }
    check(a.settlement.size() == b.settlement.size() && !a.settlement.empty(),
          "both report the same number of steps");
    if (a.settlement.size() != b.settlement.size() || a.settlement.empty()) return;
    const double s_inf = kQ * kH / kEoed;
    double worst = 0.0;
    for (size_t k = 0; k < a.settlement.size(); ++k)
        worst = std::max(worst, std::fabs(a.settlement[k] - b.settlement[k]) / s_inf);
    std::printf("   largest settlement difference between the two orders: %.3e of q H / Eoed\n", worst);
    check(worst < 1e-12, "the same settlement history, to round-off, whichever material is first");

    if (!with_history) return;
    const double n_up = porosity(kUpperE), n_lo = porosity(kLowerE);
    FiniteVolumeColumn fv(2400, 1.0e-6, n_lo, n_up);
    // The two single-porosity columns the engine used to solve, depending on the list order. They
    // are here so the band below is shown to be able to say NO: a band both of them would also
    // pass would verify nothing about the storage term.
    FiniteVolumeColumn fv_up(2400, 1.0e-6, n_up, n_up), fv_lo(2400, 1.0e-6, n_lo, n_lo);
    double worst_fv = 0.0, worst_up = 0.0, worst_lo = 0.0;
    std::printf("   t [day]     FE settlement   finite-volume   diff / s_inf   one n = 0.2   one n = 0.6\n");
    for (double t_sample : {0.002, 0.005, 0.01, 0.02, 0.04}) {
        size_t k = 0;
        for (size_t j = 0; j < a.time.size(); ++j)
            if (std::fabs(a.time[j] - t_sample) < std::fabs(a.time[k] - t_sample)) k = j;
        fv.advance_to(a.time[k]);
        fv_up.advance_to(a.time[k]);
        fv_lo.advance_to(a.time[k]);
        const double s_fv = fv.settlement();
        const double diff = (a.settlement[k] - s_fv) / s_inf;
        const double d_up = (fv_up.settlement() - s_fv) / s_inf;
        const double d_lo = (fv_lo.settlement() - s_fv) / s_inf;
        worst_fv = std::max(worst_fv, std::fabs(diff));
        worst_up = std::max(worst_up, std::fabs(d_up));
        worst_lo = std::max(worst_lo, std::fabs(d_lo));
        std::printf("   %.4f      %.6e    %.6e    %+.4f%%      %+.3f%%       %+.3f%%\n", a.time[k],
                    a.settlement[k], s_fv, 100.0 * diff, 100.0 * d_up, 100.0 * d_lo);
    }
    check(worst_fv < 0.003,
          "the settlement history follows the independent finite-volume solution within 0.3%");
    check(worst_up > 0.003 + worst_fv && worst_lo > 0.003 + worst_fv,
          "and either single porosity would fall outside that band -- the band can say no");
}

}  // namespace

int main() {
    std::printf("A layered column: the storage term n/Kw belongs to each material\n");
    test_layer_closed_form(m::PhaseType::Consolidation, "consolidation");
    test_layer_closed_form(m::PhaseType::FullyCoupled, "fully coupled");
    test_order_and_history(m::PhaseType::Consolidation, "consolidation", true);
    test_order_and_history(m::PhaseType::FullyCoupled, "fully coupled", false);
    if (g_failures == 0) {
        std::printf("\nOK: each layer stores its own water, whatever order the list is in\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
