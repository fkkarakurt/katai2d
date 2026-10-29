// Jointed rock: the material point against closed forms that share no code with it.
//
// The model (jointed_rock.hpp) is elastic rock -- isotropic or transversely isotropic about the
// normal to its first joint set -- that fails only by slip or opening on up to three sets of
// joints. Its two halves are checked separately, each against an answer that is older than the
// code:
//   - the elasticity against the directional Young's modulus of a transversely isotropic body
//     (Lekhnitskii): a uni-axial stress along a direction at angle theta to the axis of symmetry
//     strains that direction by sigma / E(theta),
//         1/E(theta) = sin^4/E1 + cos^4/E2 + (1/G2 - 2 nu2/E2) sin^2 cos^2;
//   - the strength against Jaeger's single plane of weakness: a specimen under a confining stress
//     sigma_3 whose major principal stress makes an angle beta with the normal to the joints slips
//     on them when
//         sigma_1 - sigma_3 = 2 (c + sigma_3 tan phi) / ((1 - tan phi cot beta) sin 2 beta),
//     and does not slip at all for beta <= phi or beta = 90 degrees (here: stays elastic, since the
//     intact rock has no strength criterion of its own);
//   - the consistent tangent against a finite difference of the return itself, on one surface, on
//     the slip-opening corner and on two crossing joint sets;
//   - the return itself, over random states: it always finds an admissible active set, ends
//     inside every surface (a joint without a cut-off included: its apex must be reachable with
//     psi = 0), and is continuous in the trial.
//   - why non-associated slip is not verified in whole problems (K2D-M008): with psi < phi the
//     acoustic tensor of the slipping point has a negative determinant for some band (Rudnicki &
//     Rice 1975) -- ellipticity is lost at the onset of slip; with psi = phi it stays neutral.
// The program half (a block through the whole program; a rock slope's factor of safety on a
// daylighting joint set) is in test_jointed_rock_fe.
//
// verify: KV-CST-025
//   oracle:   closed_form
//   source:   Lekhnitskii, S.G. (1963). Theory of Elasticity of an Anisotropic Elastic Body. Holden-Day, San Francisco -- the directional modulus of a transversely isotropic body; Jaeger, J.C. (1960). Shear failure of anisotropic rocks. Geological Magazine 97(1), 65-72 -- the single plane of weakness; Rudnicki, J.W. & Rice, J.R. (1975). Conditions for the localization of deformation in pressure-sensitive dilatant materials. J. Mech. Phys. Solids 23(6), 371-394 -- the acoustic tensor
//   locator:  E1 = 2 GPa, E2 = 0.5 GPa, nu1 = 0.25, nu2 = 0.2, G2 = 0.3 GPa, joints dipping 30 deg; uni-axial stress along 13 directions; a plane strain specimen under sigma_3 = 0 and 500 kPa, sigma_xy held at 0, compressed vertically in 2e-6 steps, joints c = 100 kPa, phi = 30 deg, psi = 0 and 10 deg, dipping 0..90 deg in 5 deg steps; tangents at three active sets, central differences of 1e-8
//   quantity: the strain along the stressed direction [-]; the major principal stress at the end of the path [kPa]; the tangent [kPa]
//   expected: the directional strain = sigma / E(theta); sigma_1 = Jaeger's for phi < beta < 90 deg, and no plastic step outside that range; the tangent = the finite difference of the return
//   band:     1e-10 relative on the strain (round-off); 1e-8 relative on sigma_1 (the return is exact, so only the path's own equilibrium iteration remains); 1e-5 relative on the tangent (the difference step)
#include <katai/materials/jointed_rock.hpp>

#include <cmath>
#include <cstdio>
#include <string>

namespace jr = katai::core::jr;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}
double rel(double a, double b) { return std::fabs(a - b) / std::max(std::fabs(b), 1e-300); }
constexpr double PI = 3.14159265358979323846;
double rad(double d) { return d * PI / 180.0; }

jr::Params rock(double alpha_deg, double c, double phi_deg, double psi_deg, bool aniso) {
    jr::Params P;
    P.E1 = 2.0e6; P.nu1 = 0.25;
    P.aniso = aniso;
    P.E2 = 0.5e6; P.nu2 = 0.2; P.G2 = 0.3e6;
    P.n = 1;
    P.p[0].alpha = rad(alpha_deg);
    P.p[0].c = c; P.p[0].phi = rad(phi_deg); P.p[0].psi = rad(psi_deg);
    P.p[0].tcut = true; P.p[0].sigt = 0.0;
    return P;
}

struct Path { double s1 = 0.0; bool plastic = false; bool ok = true; };

// A plane strain specimen: sigma_xx = -s3 and sigma_xy = 0 held, eps_yy driven in compression.
// The two held components are found by Newton on the return's own tangent.
Path compress(const jr::Params& P, double s3, int steps, double h) {
    Path out;
    const jr::Mat4 D = jr::elastic(P);
    jr::Vec4 sig(-s3, -s3, 0.0, 0.0);
    // the initial isotropic in-plane state, sigma_zz from nothing: take it as it is
    for (int i = 0; i < steps; ++i) {
        jr::Vec4 de(0.0, -h, 0.0, 0.0);
        jr::Return r;
        for (int it = 0; it < 50; ++it) {
            r = jr::return_map(P, D, sig + D * de);
            if (!r.ok) { out.ok = false; return out; }
            const Eigen::Vector2d res(r.sig(0) + s3, r.sig(3));
            if (res.norm() < 1e-9 * (1.0 + s3)) break;
            Eigen::Matrix2d J;
            J << r.T(0, 0), r.T(0, 3), r.T(3, 0), r.T(3, 3);
            const Eigen::Vector2d dx = J.fullPivLu().solve(res);
            de(0) -= dx(0);
            de(3) -= dx(1);
        }
        sig = r.sig;
        out.plastic = out.plastic || r.plastic;
    }
    out.s1 = -sig(1);
    return out;
}

double jaeger(double s3, double c, double phi, double beta) {
    const double t = std::tan(phi);
    return s3 + 2.0 * (c + s3 * t) / ((1.0 - t / std::tan(beta)) * std::sin(2.0 * beta));
}

}  // namespace

int main() {
    std::printf("== Jointed rock at the material point (KV-CST-025) ==\n");
    {
        // (a) the directional modulus
        const jr::Params P = rock(30.0, 100.0, 30.0, 0.0, true);
        const jr::Mat4 D = jr::elastic(P);
        check((D - D.transpose()).cwiseAbs().maxCoeff() < 1e-9 * D.cwiseAbs().maxCoeff(), "the elastic operator is symmetric");
        check(Eigen::SelfAdjointEigenSolver<jr::Mat4>(D).eigenvalues().minCoeff() > 0.0, "...and positive definite");
        check(jr::elastic_admissible(P) == nullptr, "the constants are admissible");
        const jr::Mat4 C = D.inverse();
        const double nx = -std::sin(P.p[0].alpha), ny = std::cos(P.p[0].alpha);
        double worst = 0.0;
        for (int k = 0; k < 13; ++k) {
            const double th = rad(15.0 * k);
            const double dx = std::cos(th), dy = std::sin(th);
            const jr::Vec4 s(dx * dx, dy * dy, 0.0, dx * dy);   // unit uni-axial stress along d
            const jr::Vec4 e = C * s;
            const double edd = dx * dx * e(0) + dy * dy * e(1) + dx * dy * e(3);
            const double cs = dx * nx + dy * ny, sn2 = 1.0 - cs * cs;
            const double inv = sn2 * sn2 / P.E1 + cs * cs * cs * cs / P.E2 +
                               (1.0 / P.G2 - 2.0 * P.nu2 / P.E2) * sn2 * cs * cs;
            worst = std::max(worst, rel(edd, inv));
        }
        std::printf("      directional strain against 1/E(theta): worst %.2e over 13 directions\n", worst);
        check(worst < 1e-10, "a uni-axial stress strains its own direction by sigma / E(theta) (Lekhnitskii)");
        // isotropic: the transversely isotropic operator with E2 = E1, nu2 = nu1, G2 = E1/(2(1+nu1))
        jr::Params Q = P;
        Q.E2 = Q.E1; Q.nu2 = Q.nu1; Q.G2 = Q.E1 / (2.0 * (1.0 + Q.nu1));
        jr::Params I = P; I.aniso = false;
        check((jr::elastic(Q) - jr::elastic(I)).cwiseAbs().maxCoeff() < 1e-8 * jr::elastic(I).cwiseAbs().maxCoeff(),
              "with E2 = E1, nu2 = nu1 and G2 = E1/(2(1 + nu1)) it is the isotropic operator, whatever the dip");
        jr::Params B = P; B.nu2 = 0.9;
        check(jr::elastic_admissible(B) != nullptr, "a set of constants whose energy is not positive is refused");
    }
    {
        // (b) Jaeger's plane of weakness
        const double c = 100.0, phi = rad(30.0);
        for (double psi : {0.0, 10.0})
            for (double s3 : {0.0, 500.0}) {
                double worst = 0.0;
                bool elastic_outside = true, all_ok = true, slipped_inside = true;
                int inside = 0;
                for (int a = 0; a <= 90; a += 5) {
                    const jr::Params P = rock((double)a, c, 30.0, psi, true);
                    const double beta = rad((double)a);
                    const bool in = beta > phi + 1e-9 && a < 90;
                    const double want = in ? jaeger(s3, c, phi, beta) : 0.0;
                    // drive well past the limit: to twice the Jaeger stress, or 40 MPa outside the range
                    const double target = in ? 2.0 * want : 4.0e4;
                    const int steps = 400;
                    const double h = target / 0.3e6 / steps;
                    const Path p = compress(P, s3, steps, h);
                    all_ok = all_ok && p.ok;
                    if (in) {
                        ++inside;
                        slipped_inside = slipped_inside && p.plastic;
                        worst = std::max(worst, rel(p.s1, want));
                    } else {
                        elastic_outside = elastic_outside && !p.plastic;
                    }
                }
                std::printf("      psi %.0f deg, sigma_3 %.0f kPa: worst |sigma_1 - Jaeger| %.2e over %d dips; outside phi < beta < 90 deg %s\n",
                            psi, s3, worst, inside, elastic_outside ? "elastic" : "SLIPPED");
                char lbl[200];
                std::snprintf(lbl, sizeof lbl, "psi %.0f, sigma_3 %.0f: every path integrates", psi, s3);
                check(all_ok, lbl);
                std::snprintf(lbl, sizeof lbl, "psi %.0f, sigma_3 %.0f: the specimen slips at Jaeger's sigma_1 for every dip phi < beta < 90 deg", psi, s3);
                check(slipped_inside && worst < 1e-8, lbl);
                std::snprintf(lbl, sizeof lbl, "psi %.0f, sigma_3 %.0f: and stays elastic for beta <= phi and beta = 90 deg", psi, s3);
                check(elastic_outside, lbl);
            }
    }
    {
        // (c) the consistent tangent against a difference of the return
        const auto fd_check = [](const jr::Params& P, const jr::Vec4& sig0, const jr::Vec4& de, int want_active,
                                 const char* what) {
            const jr::Mat4 D = jr::elastic(P);
            const jr::Return r = jr::return_map(P, D, sig0 + D * de);
            jr::Mat4 F;
            for (int j = 0; j < 4; ++j) {
                jr::Vec4 d1 = de, d2 = de;
                const double h = 1e-8;
                d1(j) += h; d2(j) -= h;
                F.col(j) = (jr::return_map(P, D, sig0 + D * d1).sig - jr::return_map(P, D, sig0 + D * d2).sig) / (2.0 * h);
            }
            const double err = (F - r.T).cwiseAbs().maxCoeff() / D.cwiseAbs().maxCoeff();
            std::printf("      %s: %d active, tangent - difference %.2e\n", what, r.n_active, err);
            check(r.ok && r.plastic && r.n_active == want_active && err < 1e-5, std::string(what) + ": the consistent tangent is the derivative of the return");
        };
        // increments written in the axes of the first set ([ss, nn, zz, sn]) and turned to global
        const auto along = [](double alpha_deg, double enn, double gsn) {
            return jr::Vec4(jr::strain_rotation(rad(alpha_deg)).inverse() * jr::Vec4(0.0, enn, 0.0, gsn));
        };
        const jr::Vec4 s0(-300.0, -300.0, -300.0, 0.0);
        fd_check(rock(50.0, 100.0, 30.0, 10.0, true), s0, along(50.0, 0.0, 2e-3), 1, "slip on one set");
        jr::Params open = rock(50.0, 20.0, 30.0, 10.0, true);
        open.p[0].sigt = 10.0;
        fd_check(open, s0, along(50.0, 2e-3, 2e-3), 2, "slip and opening (the corner)");
        // Two crossing sets, both slipping. (Two sets placed symmetrically about the loading are no
        // test of this: there the two-set step has a negative multiplier -- slip on one set loads
        // the other more than it relieves itself -- and the answer is slip on either set alone, a
        // genuine bifurcation at which the return has no derivative.) Sets at 50 and -30 degrees,
        // the second the weaker, compressed vertically and sheared.
        jr::Params two = rock(50.0, 100.0, 30.0, 10.0, true);
        two.n = 2;
        two.p[1] = two.p[0];
        two.p[1].alpha = rad(-30.0);
        two.p[1].c = 40.0;
        // the first increment of a fan, in growing size, at which both sets slip
        jr::Vec4 de2 = jr::Vec4::Zero();
        bool found = false;
        const jr::Mat4 D2 = jr::elastic(two);
        for (int k = 0; k < 72 && !found; ++k)
            for (double size = 1e-3; size < 2e-2 && !found; size *= 1.5) {
                const double th = rad(5.0 * k);
                const jr::Vec4 d(size * std::cos(th), -size * std::cos(th), 0.0, 2.0 * size * std::sin(th));
                const jr::Return r = jr::return_map(two, D2, s0 + D2 * d);
                if (r.n_active == 2) {
                    jr::Surface S[9];
                    jr::surfaces(two, S);
                    if (S[r.active[0]].plane != S[r.active[1]].plane) { de2 = d; found = true; }
                }
            }
        check(found, "an increment at which both sets slip exists");
        if (found) fd_check(two, s0, de2, 2, "slip on two crossing sets");
    }
    {
        // (d) opening: pulled across its joints the rock carries the tensile strength, no more
        jr::Params P = rock(0.0, 50.0, 30.0, 0.0, false);
        P.p[0].sigt = 20.0;
        const jr::Mat4 D = jr::elastic(P);
        const jr::Return r = jr::return_map(P, D, D * jr::Vec4(0.0, 1e-3, 0.0, 0.0));
        std::printf("      pulled across horizontal joints: sigma_yy %.6f (sigma_t 20)\n", r.sig(1));
        check(r.plastic && std::fabs(r.sig(1) - 20.0) < 1e-9, "a horizontal joint pulled open holds sigma_yy = sigma_t");
        P.p[0].sigt = 1e3;   // above the apex c / tan phi = 86.6
        const jr::Return a = jr::return_map(P, D, D * jr::Vec4(0.0, 1e-2, 0.0, 0.0));
        check(a.plastic && std::fabs(a.sig(1) - 50.0 / std::tan(rad(30.0))) < 1e-8, "a tensile strength above the apex is held at c / tan phi");
    }
    {
        // (e) every return finds its active set. Random states and increments -- from the apex
        // of a cohesionless set at zero stress up to deep inside, one to three sets, isotropic and
        // anisotropic -- and the returned stress must satisfy every surface with every multiplier
        // non-negative. A return that finds no admissible set hands back the trial, and an
        // iteration then sees a stress outside the strength.
        unsigned seed = 12345u;
        const auto rnd = [&seed]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0 / 16777216.0); };
        int failed = 0, outside = 0, plastic = 0;
        for (int k = 0; k < 4000; ++k) {
            jr::Params P = rock(180.0 * rnd() - 90.0, k % 4 == 0 ? 0.0 : 100.0 * rnd(), 20.0 + 25.0 * rnd(), 0.0, k % 2 == 0);
            P.p[0].psi = P.p[0].phi * rnd();
            P.p[0].sigt = k % 3 == 0 ? 0.0 : 20.0 * rnd();
            P.p[0].tcut = k % 5 != 1;   // and without a cut-off: the apex must still be reachable
            P.n = 1 + k % 3;
            for (int i = 1; i < P.n; ++i) {
                P.p[i] = P.p[0];
                P.p[i].alpha = rad(180.0 * rnd() - 90.0);
                P.p[i].c = k % 5 == 0 ? 0.0 : 80.0 * rnd();
            }
            const jr::Mat4 D = jr::elastic(P);
            const double s = k % 7 == 0 ? 0.0 : 500.0 * rnd();
            const jr::Vec4 s0(-s * rnd(), -s * rnd(), -s * rnd(), 0.0);
            const jr::Vec4 de(2e-3 * (rnd() - 0.5), 2e-3 * (rnd() - 0.5), 0.0, 2e-3 * (rnd() - 0.5));
            const jr::Return r = jr::return_map(P, D, s0 + D * de);
            if (!r.ok) { ++failed; continue; }
            plastic += r.plastic;
            jr::Surface S[9];
            const int ns = jr::surfaces(P, S);
            double scale = 1.0 + r.sig.cwiseAbs().maxCoeff();
            for (int q = 0; q < ns; ++q)
                if (S[q].a.dot(r.sig) - S[q].r > 1e-8 * scale) { ++outside; break; }
        }
        std::printf("      4000 random returns: %d plastic, %d without an admissible set, %d outside a surface\n", plastic, failed, outside);
        check(failed == 0 && outside == 0 && plastic > 1000, "(e) every return finds an admissible set and ends inside every surface");
    }
    {
        // (f) the return is continuous in the trial. Newton can only converge on a stress-strain
        // map without jumps; with non-associated flow more than one active set can satisfy the
        // Kuhn-Tucker conditions, and a rule that picks between them by which surfaces the trial
        // violates can jump where the trial crosses a surface. Pairs of trials 1e-9 apart, across
        // the corner of slip and opening of a cohesionless set and through its apex.
        unsigned seed = 777u;
        const auto rnd = [&seed]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0 / 16777216.0); };
        double worst = 0.0;
        jr::Vec4 at_w = jr::Vec4::Zero();
        for (int k = 0; k < 20000; ++k) {
            jr::Params P = rock(35.0, k % 2 ? 0.0 : 20.0, 45.0, 0.0, false);
            P.p[0].sigt = 0.0;
            const jr::Mat4 D = jr::elastic(P);
            const jr::Vec4 t(400.0 * (rnd() - 0.5), 400.0 * (rnd() - 0.5), 400.0 * (rnd() - 0.5), 400.0 * (rnd() - 0.5));
            jr::Vec4 d(rnd() - 0.5, rnd() - 0.5, 0.0, rnd() - 0.5);
            d *= 1e-9 / d.norm();
            const jr::Return a = jr::return_map(P, D, t), b = jr::return_map(P, D, t + D * d);
            const double jump = (a.sig - b.sig).norm() / (D * d).norm();   // a jump is >> 1
            if (jump > worst) { worst = jump; at_w = t; }
        }
        std::printf("      continuity: worst |d sigma| / |d trial| %.2f over 20000 trial pairs 1e-9 apart (at [%.1f %.1f %.1f %.1f])\n",
                    worst, at_w(0), at_w(1), at_w(2), at_w(3));
        check(worst < 10.0, "(f) the return is continuous in the trial (no jump between admissible active sets)");
    }
    {
        // (g) simple shear across horizontal joints: sigma_yy held at -sv, eps_xx = 0, gamma_xy
        // driven. The principal axes rotate as the test goes on, which Jaeger's test (b) never
        // does; the joints slip at tau = sv tan(phi), whatever the dilatancy.
        for (double psi : {0.0, 5.0, 30.0}) {
            jr::Params P = rock(0.0, 0.0, 30.0, psi, true);
            const jr::Mat4 D = jr::elastic(P);
            const double sv = 200.0;
            jr::Vec4 sig(-sv * 0.5, -sv, -sv * 0.5, 0.0);
            for (int i = 0; i < 400; ++i) {
                jr::Vec4 de(0.0, 0.0, 0.0, 2e-5);
                jr::Return r;
                for (int it = 0; it < 50; ++it) {
                    r = jr::return_map(P, D, sig + D * de);
                    const double res = r.sig(1) + sv;
                    if (std::fabs(res) < 1e-9 * sv) break;
                    de(1) -= res / r.T(1, 1);
                }
                sig = r.sig;
            }
            std::printf("      simple shear on horizontal joints, psi %.0f: tau %.6f kPa, sv tan(phi) %.6f kPa\n",
                        psi, sig(3), sv * std::tan(rad(30.0)));
            char lbl[160];
            std::snprintf(lbl, sizeof lbl, "(g) simple shear, psi %.0f: the joints slip at sv tan(phi)", psi);
            check(rel(std::fabs(sig(3)), sv * std::tan(rad(30.0))) < 1e-8, lbl);
        }
    }
    {
        // (h) WHY NON-ASSOCIATED SLIP IS NOT VERIFIED IN WHOLE PROBLEMS (K2D-M008): the loss of
        // ellipticity at the onset of slip. A band of normal m can form where the acoustic tensor
        // Q(m) = N(m)^T D_ep N(m) is singular (Rudnicki, J.W. & Rice, J.R. (1975). Conditions for
        // the localization of deformation in pressure-sensitive dilatant materials. J. Mech. Phys.
        // Solids 23(6), 371-394). With associated flow the tangent is symmetric and positive
        // semi-definite, so det Q >= 0 for every m; with psi < phi it is not, and det Q < 0 for
        // some band: the incremental problem has lost ellipticity at the yield point itself, which
        // is the mesh dependence and the stalls the slope showed. Measured at a point slipping on
        // a set dipping 50 degrees under sigma_3 = 200 kPa.
        for (double psi : {30.0, 10.0, 0.0}) {
            jr::Params P = rock(50.0, 100.0, 30.0, psi, false);
            const double s3 = 200.0;
            const double s1 = jaeger(s3, 100.0, rad(30.0), rad(50.0));
            const jr::Mat4 D = jr::elastic(P);
            // a trial a hair beyond the yield state, pushed along the major stress
            const jr::Vec4 at_yield(-s3, -s1, -0.5 * (s1 + s3) * 0.5, 0.0);
            const jr::Return r = jr::return_map(P, D, at_yield + jr::Vec4(0.0, -1e-3 * s1, 0.0, 0.0));
            Eigen::Matrix3d T;
            const int idx[3] = {0, 1, 3};
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) T(i, j) = r.T(idx[i], idx[j]);
            double worst = 1e300;
            for (int k = 0; k < 3600; ++k) {
                const double a = rad(0.05 * k);
                const double mx = std::cos(a), my = std::sin(a);
                Eigen::Matrix<double, 3, 2> N;
                N << mx, 0.0, 0.0, my, my, mx;
                const Eigen::Matrix2d Q = N.transpose() * T * N;
                worst = std::min(worst, Q.determinant());
            }
            const Eigen::Matrix2d Qe = [&] {
                Eigen::Matrix3d De;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) De(i, j) = D(idx[i], idx[j]);
                Eigen::Matrix<double, 3, 2> N;
                N << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0;
                return Eigen::Matrix2d(N.transpose() * De * N);
            }();
            const double scale = Qe.determinant();
            std::printf("      slipping, psi %.0f deg: min over bands of det Q / det Q_elastic = %+.3e\n", psi, worst / scale);
            char lbl[200];
            if (psi == 30.0) {
                std::snprintf(lbl, sizeof lbl, "(h) associated slip keeps the acoustic tensor non-negative (no band below the limit)");
                check(r.plastic && worst / scale > -1e-10, lbl);
            } else {
                std::snprintf(lbl, sizeof lbl, "(h) psi %.0f < phi: a band has det Q < 0 at the onset of slip (ellipticity lost)", psi);
                check(r.plastic && worst / scale < -1e-6, lbl);
            }
        }
    }
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
