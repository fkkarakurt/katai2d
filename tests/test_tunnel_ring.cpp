// A lining drawn as the chords of a circle, around a deep tunnel in elastic ground under an
// isotropic in-situ stress (K0 = 1): unlined, lined as it is dug, lined with an interface (stiff and
// soft), and lined after it is dug.
//
// THE JOINT WHERE TWO SPLIT LINES MEET. A plate with an interface is a wall: the mesh is split along
// it, and the two sides meet the plate only through the interfaces. At the crown of a ring two
// chords end at one point, and the second chord split again the two nodes the first had left there:
// four nodes, a wall position that belonged to no element, and a ring cut at its crown. Measured
// before 2026-10-26: the lined tunnel with an interface carried 23 kN/m of the 2690 below, its two
// halves 0.55 m apart at the crown. And where the walls were not yet built, the twins they tie at
// that point formed a cycle, which left the point without an equation -- fixed: the crown of an
// unlined tunnel did not move at all. And where the deeper ends of two chords meet -- the invert --
// both hung from the ground's node there: the ring was welded to the ground at that point.
//
// THE JOINT A DIG LINES. An interface installed in the phase that digs out the side it is bonded
// across has no ground on that side: it starts as a new contact, sigma_n0 = 0. Measured before
// 2026-10-26: it started at half the ground's stress (-500 of -1000 kPa, the empty side's nil
// stress averaged in), and with kn = 1e5 the ring carried 5.9% more than the closed form, its
// joint's normal stress 500 kPa beyond kn times its closure -- a joint that read 2.3 times as stiff.
//
// THE CLOSED FORMS. A circular hole of radius R in an infinite elastic plane under an isotropic
// stress p0 at infinity: the wall moves in by u = p0 R (1 + nu) / E (Kirsch). A thin ring of axial
// stiffness EA built as the ground is dug, held to the ground at its outside, takes a pressure p_i,
// and three movements must agree -- the ground's, (p0 - p_i) R (1 + nu) / E; the interface's, p_i / kn
// (none without one); and the ring's, p_i R^2 / EA. The ring's hoop force is N = p_i R:
//   N = p0 R a / (a + 1/kn + R^2/EA),   a = R (1 + nu) / E.
// The ground here has weight, so p0 is the stress at the tunnel's centre and the comparison is made on
// the mean over the ring, to which the stress gradient across the tunnel contributes nothing.
//
// verify: KV-STR-024
//   oracle:   closed_form
//   source:   Kirsch, G. (1898). Die Theorie der Elastizitaet und die Beduerfnisse der Festigkeitslehre. Zeitschrift des Vereines deutscher Ingenieure 42, 797-807 -- the circular hole in an elastic plane; Einstein, H.H. & Schwartz, C.W. (1979). Simplified analysis for tunnel supports. Journal of the Geotechnical Engineering Division, ASCE 105(GT4), 499-518 -- the support ring bonded to the ground, here its axisymmetric (isotropic-stress) part, where the ring's flexural stiffness does not enter
//   locator:  a deep tunnel, R = 3 m, centre 50 m below the surface of a 120 x 100 m block; linear elastic ground E = 50000 kPa, nu = 0.3, gamma = 20 kN/m3, K0 = 1 (p0 = 1000 kPa at the centre); the outline 48 chords; the ring EA = 1e6 kN/m, EI = 1e3 kNm2/m; the interface kn = 1e6 kPa/m, and 1e5. Phases: K0 procedure, then the tunnel dug (and lined); or dug, then lined
//   quantity: the mean inward movement of the tunnel wall [m]; the mean hoop force of the ring [kN/m]; the displacement of the two chords at the crown; the interface's normal stress and closure; the movement a later-built ring causes
//   expected: u = p0 R (1 + nu) / E = 0.078 m; N = p0 R a / (a + 1/kn + R^2/EA) = 2689.7 kN/m bonded, 2659.1 with the interface, 2412.4 with the soft one; the chords at the crown move together; the interface closes at the invert as around the ring; the soft interface's normal stress is kn times its closure and N/R; a ring built after the dig takes nothing and moves nothing
//   band:     2% on u (measured -0.79%: the chords' polygon in place of the circle, and the 20-radius block in place of the plane); 3% on N (measured -0.49% bonded, -0.50% with the interface, -0.52% with the soft one); the interface closure at the invert at least half the ring's mean (measured 0.903 against 0.884 mm; 0 when the point was welded); 2% between the soft joint's normal stress and kn times its closure (measured 0.00%), 3% between it and N/R (+0.05%); 1e-9 m between the chords at the crown; the later-built ring: |N| below 0.5% of the bonded value and the largest movement below 1e-3 of u
#include <katai/jobs/driver.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/model/project.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace m = katai::model;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

constexpr double kPi = 3.14159265358979323846;
constexpr double kE = 50000.0, kNu = 0.3, kGamma = 20.0, kR = 3.0, kYc = 50.0;
constexpr int kChords = 48;
constexpr double kEA = 1.0e6, kKn = 1.0e6;

struct Case { bool lining; bool interface; bool line_later; double kn = kKn; };

m::Project tunnel(const Case& c) {
    m::Project pr;
    pr.x_min = -60; pr.x_max = 60; pr.y_min = 0; pr.y_max = 100;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 6.0;
    pr.mesh.auto_refine = true;   // the ring's own factor (0.08 below) applies
    m::Material s;
    s.name = "Ground"; s.model = m::SoilModel::LinearElastic; s.drainage = m::Drainage::Drained;
    s.E = kE; s.nu = kNu; s.gamma_unsat = s.gamma_sat = kGamma; s.k0_auto = false; s.k0 = 1.0;
    pr.materials.push_back(s);
    m::SoilPolygon g;
    g.name = "Ground"; g.material = 0; g.x = {-60, 60, 60, -60}; g.y = {0, 0, 100, 100};
    g.edge_bc.assign(4, 0);
    pr.polygons.push_back(g);
    // The tunnel region, counter-clockwise from the invert: the crown is a vertex.
    m::SoilPolygon t;
    t.name = "Tunnel"; t.material = 0; t.coarseness = 0.08;
    for (int k = 0; k < kChords; ++k) {
        const double a = -0.5 * kPi + 2.0 * kPi * k / kChords;
        t.x.push_back(kR * std::cos(a));
        t.y.push_back(kYc + kR * std::sin(a));
    }
    t.edge_bc.assign(kChords, 0);
    pr.polygons.push_back(t);
    if (c.lining) {
        m::PlateMaterial ring;
        ring.name = "Ring"; ring.EA = kEA; ring.EI = 1.0e3; ring.nu = 0.0;
        pr.plates.push_back(ring);
        for (int k = 0; k < kChords; ++k) {
            m::StructElement e;
            e.kind = m::StructKind::Plate;
            e.name = "Chord " + std::to_string(k + 1);
            e.x1 = t.x[(size_t)k]; e.y1 = t.y[(size_t)k];
            e.x2 = t.x[(size_t)(k + 1) % kChords]; e.y2 = t.y[(size_t)(k + 1) % kChords];
            e.material = 0;
            e.coarseness = 0.08;
            e.iface_pos = c.interface;   // the ground side of a counter-clockwise outline
            e.iface_kn = c.interface ? c.kn : 0.0;
            pr.structs.push_back(e);
        }
        pr.initial.struct_active.assign((size_t)kChords, 0);
    }
    m::Phase dig = pr.initial;
    dig.name = "Dig";
    dig.poly_active = {1, 0};
    dig.poly_dry = {0, 1};
    if (c.lining) dig.struct_active.assign((size_t)kChords, c.line_later ? 0 : 1);
    pr.phases.push_back(dig);
    if (c.lining && c.line_later) {
        m::Phase line = dig;
        line.name = "Line";
        line.struct_active.assign((size_t)kChords, 1);
        pr.phases.push_back(line);
    }
    return pr;
}

std::vector<katai::app::SolveResult> run(const Case& c, const char* what) {
    const m::Project pr = tunnel(c);
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) {
        check(false, std::string(what) + ": meshes (" + M.message + ")");
        return {};
    }
    auto res = katai::app::solve_phases(pr, M.mesh, katai::app::initial_phase_from(pr.initial_procedure));
    bool ok = res.size() == pr.phases.size() + 1;
    for (const auto& r : res) ok = ok && r.ok;
    check(ok, std::string(what) + ": every phase solves");
    if (!ok) res.clear();
    return res;
}

// The mean inward movement of the tunnel wall, over the mesh nodes on the circle, and the fewest a
// node of it moves.
void wall_movement(const katai::app::SolveResult& r, double& mean, double& least) {
    double sum = 0.0;
    int n = 0;
    least = 1e300;
    for (int k = 0; k < r.mesh.node_count; ++k) {
        const double dx = r.mesh.x[k], dy = r.mesh.y[k] - kYc, rr = std::hypot(dx, dy);
        if (std::fabs(rr - kR) > 1e-6) continue;
        const double ur = -(r.disp[2 * k] * dx + r.disp[2 * k + 1] * dy) / rr;
        sum += ur;
        least = std::min(least, ur);
        ++n;
    }
    mean = n ? sum / n : 0.0;
}

double hoop_force(const katai::app::SolveResult& r) {
    double sum = 0.0;
    int n = 0;
    for (const auto& f : r.struct_forces)
        for (const auto& st : f.stations) { sum += st.N; ++n; }
    return n ? -sum / n : 0.0;
}

}  // namespace

int main() {
    std::printf("== a lined circular tunnel against its closed forms (KV-STR-024) ==\n");
    const double p0 = kGamma * (100.0 - kYc);
    const double u_kirsch = p0 * kR * (1.0 + kNu) / kE;
    const double a = kR * (1.0 + kNu) / kE;
    const double n_bonded = p0 * kR * a / (a + kR * kR / kEA);
    const double n_iface = p0 * kR * a / (a + 1.0 / kKn + kR * kR / kEA);

    {
        const auto res = run({false, false, false}, "unlined");
        if (!res.empty()) {
            double u = 0.0, least = 0.0;
            wall_movement(res[1], u, least);
            std::printf("      wall moves in %.6f m, Kirsch %.6f (%+.2f%%)\n", u, u_kirsch, 100.0 * (u / u_kirsch - 1.0));
            check(std::fabs(u / u_kirsch - 1.0) < 0.02, "unlined: the wall moves in by p0 R (1 + nu) / E, within 2%");
        }
    }
    {
        const auto res = run({true, false, false}, "bonded ring");
        if (!res.empty()) {
            const double N = hoop_force(res[1]);
            std::printf("      hoop force %.2f kN/m, closed form %.2f (%+.2f%%)\n", N, n_bonded, 100.0 * (N / n_bonded - 1.0));
            check(std::fabs(N / n_bonded - 1.0) < 0.03, "a ring bonded to the ground carries N = p0 R a / (a + R^2/EA), within 3%");
        }
    }
    {
        const auto res = run({true, true, false}, "ring with an interface");
        if (!res.empty()) {
            const auto& r = res[1];
            const double N = hoop_force(r);
            std::printf("      hoop force %.2f kN/m, closed form %.2f (%+.2f%%)\n", N, n_iface, 100.0 * (N / n_iface - 1.0));
            check(std::fabs(N / n_iface - 1.0) < 0.03,
                  "a ring with an interface carries N = p0 R a / (a + 1/kn + R^2/EA), within 3%");
            // The two chords that end at the crown: the end of the one and the start of the other.
            const katai::core::ForceStation* end_of = nullptr;
            const katai::core::ForceStation* start_of = nullptr;
            for (const auto& f : r.struct_forces)
                for (const auto& st : f.stations)
                    if (std::fabs(st.x) < 1e-9 && std::fabs(st.y - (kYc + kR)) < 1e-9)
                        (end_of ? start_of : end_of) = &st;
            const bool both = end_of && start_of;
            check(both && std::hypot(end_of->ux - start_of->ux, end_of->uy - start_of->uy) < 1e-9,
                  "...and its two chords meet at the crown: one point, one displacement");
            // The invert, where the deeper ends of two chords meet: the interface is there as
            // everywhere else, compressed by the ring's pressure -- not a point welded shut.
            double gap_invert = 0.0, gap_mean = 0.0, sn_mean = 0.0;
            int n_invert = 0, n_all = 0;
            for (const auto& f : r.interface_forces)
                for (const auto& st : f.stations) {
                    gap_mean += st.gap;
                    sn_mean += st.sigma_n;
                    ++n_all;
                    if (std::fabs(st.x) < 1e-9 && std::fabs(st.y - (kYc - kR)) < 1e-9) { gap_invert += st.gap; ++n_invert; }
                }
            gap_mean /= std::max(n_all, 1);
            sn_mean /= std::max(n_all, 1);
            std::printf("      interface normal stress around the ring %.2f kPa (ring pressure N/R %.2f)%c", sn_mean, N / kR, 10);
            if (n_invert) gap_invert /= n_invert;
            std::printf("      interface closure at the invert %.4e m, around the ring %.4e m\n", gap_invert, gap_mean);
            check(n_invert > 0 && gap_invert < 0.5 * gap_mean,
                  "...and at the invert the interface closes as it does around the ring");
        }
    }
    {
        // A softer joint, where 1/kn is a fifth of the ground's compliance. The lining is built in
        // the phase that digs it out, so nothing is across the joint to have held the ground: the
        // joint starts as a new contact (sigma_n0 = 0) and its normal stress is its stiffness
        // times its closure, equal to the ring's pressure. Seeded with half the ground's stress
        // (the bonded side's empty node averaged in) it read as a joint 2.3 times as stiff.
        const double kn = 1.0e5;
        const double n_soft = p0 * kR * a / (a + 1.0 / kn + kR * kR / kEA);
        const auto res = run({true, true, false, kn}, "ring with a soft interface");
        if (!res.empty()) {
            const auto& r = res[1];
            const double N = hoop_force(r);
            double gap_mean = 0.0, sn_mean = 0.0;
            int n_all = 0;
            for (const auto& f : r.interface_forces)
                for (const auto& st : f.stations) { gap_mean += st.gap; sn_mean += st.sigma_n; ++n_all; }
            gap_mean /= std::max(n_all, 1);
            sn_mean /= std::max(n_all, 1);
            std::printf("      hoop force %.2f kN/m, closed form %.2f (%+.2f%%)\n", N, n_soft, 100.0 * (N / n_soft - 1.0));
            std::printf("      interface normal stress %.2f kPa, kn x closure %.2f, ring pressure N/R %.2f\n",
                        sn_mean, kn * gap_mean, N / kR);
            check(std::fabs(N / n_soft - 1.0) < 0.03,
                  "a ring with a soft interface carries N = p0 R a / (a + 1/kn + R^2/EA), within 3%");
            check(std::fabs(std::fabs(sn_mean) / (kn * std::fabs(gap_mean)) - 1.0) < 0.02,
                  "...the joint built with the dig starts as a new contact: its normal stress is kn x its closure");
            check(std::fabs(std::fabs(sn_mean) / (N / kR) - 1.0) < 0.03,
                  "...and is the pressure the ring carries, N/R");
        }
    }
    {
        const auto res = run({true, true, true}, "ring built after the dig");
        if (!res.empty()) {
            double u = 0.0, least = 0.0;
            wall_movement(res[1], u, least);
            std::printf("      the dig: wall moves in %.6f m (Kirsch %.6f), the least of it %.6f m\n", u, u_kirsch, least);
            check(std::fabs(u / u_kirsch - 1.0) < 0.02 && least > 0.5 * u_kirsch,
                  "the dig, with the walls drawn but not built: every point of the wall moves, as unlined");
            double most = 0.0;
            for (int k = 0; k < res[2].mesh.node_count; ++k)
                most = std::max(most, std::hypot(res[2].disp[2 * k], res[2].disp[2 * k + 1]));
            const double N = hoop_force(res[2]);
            std::printf("      the lining: largest movement %.3e m, hoop force %.3f kN/m\n", most, N);
            check(most < 1e-3 * u_kirsch && std::fabs(N) < 0.005 * n_bonded,
                  "a ring built after the dig takes nothing and moves nothing");
        }
    }

    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
