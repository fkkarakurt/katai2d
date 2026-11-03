// A shell of revolution with a capacity: the plate's M-N diamond in its meridian and a cap of its
// own on the hoop force.
//
// In axial symmetry a plate is a shell of revolution (KV-STR-029) with four generalized stresses,
// s = (N_s, N_t, M_s, M_t), coupled through the plate's Poisson's ratio. With a capacity the
// meridional force and moment keep the diamond of the plane strain plate, |N_s|/Np + |M_s|/Mp <= 1;
// the hoop force has a cap of its own, |N_t| <= Np, which takes no part in the diamond (an
// isotropic section has the same Np in both directions); the hoop moment is not capped. The return
// is the closest admissible point in the energy norm of the elastic law (plate::shell_return).
//   (a) The return at a point, over 4000 trial states drawn at random round the admissible set:
//       the result is admissible; the plastic strain lies in the normal cone (the principle of
//       maximum plastic dissipation against every vertex of the set, with any hoop moment, and no
//       plastic hoop curvature); the tangent is the derivative of the return (central differences).
//   (b) With nu = 0 the hoop is decoupled from the meridian: (N_s, M_s) is the plane strain
//       diamond's return bit for bit and N_t the trial hoop force clipped at +-Np -- the caps are
//       independent (a state at the hoop cap and on the diamond at once is admissible).
//   (c) A cylinder of radius R, free at its ends, expanded radially by u_r (all its nodes) through
//       and back from its hoop capacity, 3- and 5-noded, by Newton on the shell alone: N_t = Np once
//       it yields and N_s = 0 throughout (statics), the radial reaction per radian Np H; elastic, the
//       reaction is EA (1 - nu^2) (u_r / R) H; unloading by du_r lowers N_t by EA (1 - nu^2) du_r / R.
//
// verify: KV-STR-034
//   oracle:   closed_form
//   source:   Simo, J.C. & Hughes, T.J.R. (1998). Computational Inelasticity. Springer, New York -- multisurface plasticity: the closest point projection, its Kuhn-Tucker conditions and the consistent tangent; Timoshenko, S.P. & Woinowsky-Krieger, S. (1959). Theory of Plates and Shells, 2nd ed. McGraw-Hill, New York, ch. 15 -- the free cylinder: N_t = p R, N_s = 0
//   locator:  EA = 1e6 kN/m, EI = 1e3 kNm/m, nu = 0.3 (0 in (b)), Np = 100 kN/m, Mp = 10 kNm/m; (a) 4000 trial states, seed 20261030; (c) R = 2 m, height 4 m, 16 elements (3-noded) / 8 (5-noded), u_r to 3 times the yield displacement and back by one tenth of it
//   quantity: (a) the largest plane value, the dissipation residual, the tangent's difference; (b) N_s, M_s, N_t; (c) N_t, N_s at the Gauss points [kN/m] and the radial reaction [kN/rad]
//   expected: (a) admissible to 1e-9, dissipation <= 1e-9, tangent to 1e-5; (b) bit for bit, and N_t = clip(EA eps_t, +-Np); (c) N_t = Np, N_s = 0, reaction Np H; elastic and unloading as stated
//   band:     1e-9 relative on (a)-(c) (exact states of the return and of the element), 1e-5 on the tangent (central differences)
#include <katai/fem/elements/plate.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace pl = katai::core::plate;

namespace {
int g_failures = 0;
void check(bool ok, const std::string& what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what.c_str());
    if (!ok) ++g_failures;
}

pl::PlateProps shell(double nu) {
    pl::PlateProps P;
    P.EA = 1.0e6; P.EI = 1.0e3; P.nu = nu; P.Np = 100.0; P.Mp = 10.0; P.axisym = true;
    return P;
}

double largest_plane(const pl::PlateProps& P, const Eigen::Vector4d& s) {
    return std::max(std::fabs(s(0)) / P.Np + std::fabs(s(2)) / P.Mp, std::fabs(s(1)) / P.Np);
}

void test_return() {
    std::printf("-- (a) the return at a point --\n");
    const pl::PlateProps P = shell(0.3);
    const Eigen::Matrix4d D = pl::shell_elasticity(P);
    std::mt19937 rng(20261030);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    double worst_f = 0.0, worst_diss = 0.0, worst_tan = 0.0, worst_law = 0.0, max_kt = 0.0;
    int yielded = 0;
    for (int n = 0; n < 4000; ++n) {
        // a committed plastic state and a total strain whose trial lies within ~2.5 times the set
        const Eigen::Vector4d ep_c(2e-4 * U(rng), 2e-4 * U(rng), 2e-2 * U(rng), 0.0);
        const Eigen::Vector4d s_tr(2.5 * P.Np * U(rng), 2.5 * P.Np * U(rng), 2.5 * P.Mp * U(rng),
                                   2.5 * P.Mp * U(rng));
        const Eigen::Vector4d e = ep_c + D.inverse() * s_tr;
        const auto r = pl::shell_return(P, e, ep_c);
        yielded += r.yielded ? 1 : 0;
        worst_f = std::max(worst_f, largest_plane(P, r.s) - 1.0);
        worst_law = std::max(worst_law, (r.s - D * (e - r.ep)).norm() / std::max(r.s.norm(), 1.0));
        const Eigen::Vector4d dp = r.ep - ep_c;
        max_kt = std::max(max_kt, std::fabs(dp(3)));
        // maximum dissipation: (s* - s) . dp <= 0 for every admissible s*, here every vertex of the
        // set (diamond corners x hoop extremes) with hoop moments of either sign
        for (double sn : {1.0, -1.0})
            for (int corner = 0; corner < 2; ++corner)
                for (double st : {1.0, -1.0})
                    for (double mt : {-1.0e3, 1.0e3}) {
                        const Eigen::Vector4d star = corner == 0
                            ? Eigen::Vector4d(sn * P.Np, st * P.Np, 0.0, mt)
                            : Eigen::Vector4d(0.0, st * P.Np, sn * P.Mp, mt);
                        worst_diss = std::max(worst_diss, (star - r.s).dot(dp) /
                                                              std::max(1e-300, r.s.norm() * dp.norm()));
                    }
        // the tangent against central differences of the return
        Eigen::Matrix4d fd;
        for (int c = 0; c < 4; ++c) {
            const double h = 1e-7 * std::max(std::fabs(e(c)), 1e-4);
            Eigen::Vector4d ep = e, em = e;
            ep(c) += h; em(c) -= h;
            fd.col(c) = (pl::shell_return(P, ep, ep_c).s - pl::shell_return(P, em, ep_c).s) / (2.0 * h);
        }
        worst_tan = std::max(worst_tan, (fd - r.D).norm() / D.norm());
    }
    std::printf("   %d of 4000 yielded; largest plane - 1 %.2e, law %.2e, dissipation %.2e, "
                "plastic hoop curvature %.1e, tangent %.2e\n",
                yielded, worst_f, worst_law, worst_diss, max_kt, worst_tan);
    check(yielded > 1000 && worst_f < 1e-9, "(a) the return is admissible");
    check(worst_law < 1e-12, "(a) the stresses are the elastic law on the new plastic strain");
    check(worst_diss < 1e-9 && max_kt == 0.0,
          "(a) the plastic strain lies in the normal cone (maximum dissipation), none in the hoop curvature");
    check(worst_tan < 1e-5, "(a) the tangent is the derivative of the return");
}

void test_independent_caps() {
    std::printf("-- (b) nu = 0: the hoop cap apart from the diamond --\n");
    const pl::PlateProps P = shell(0.0);
    pl::PlateProps plane = P;
    plane.axisym = false;
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(-3.0, 3.0);
    bool same = true, clipped = true;
    for (int n = 0; n < 2000; ++n) {
        const Eigen::Vector4d ep_c(1e-5 * U(rng), 1e-5 * U(rng), 1e-3 * U(rng), 0.0);
        const Eigen::Vector4d e = ep_c + Eigen::Vector4d(P.Np / P.EA * U(rng), P.Np / P.EA * U(rng),
                                                         P.Mp / P.EI * U(rng), P.Mp / P.EI * U(rng));
        const auto r = pl::shell_return(P, e, ep_c);
        const auto m = pl::mn_return(plane, e(0), e(2), ep_c(0), ep_c(2));
        const double nt = std::clamp(P.EA * (e(1) - ep_c(1)), -P.Np, P.Np);
        same = same && std::fabs(r.s(0) - m.N) <= 1e-9 * P.Np && std::fabs(r.s(2) - m.M) <= 1e-9 * P.Mp;
        clipped = clipped && std::fabs(r.s(1) - nt) <= 1e-9 * P.Np;
    }
    check(same, "(b) the meridional force and moment are the plane strain diamond's return");
    check(clipped, "(b) the hoop force is its trial value clipped at +-Np");
    // a state at the hoop cap and on the diamond at once
    const Eigen::Vector4d e(2.0 * P.Np / P.EA, 2.0 * P.Np / P.EA, 0.5 * P.Mp / P.EI, 0.0);
    const auto r = pl::shell_return(P, e, Eigen::Vector4d::Zero());
    std::printf("   both caps: N_s/Np + M_s/Mp = %.12f, N_t/Np = %.12f\n",
                std::fabs(r.s(0)) / P.Np + std::fabs(r.s(2)) / P.Mp, r.s(1) / P.Np);
    check(std::fabs(std::fabs(r.s(0)) / P.Np + std::fabs(r.s(2)) / P.Mp - 1.0) < 1e-12 &&
              std::fabs(r.s(1) / P.Np - 1.0) < 1e-12,
          "(b) the caps are independent: on the diamond and at the hoop cap at the same time");
}

// (c) the free cylinder expanded radially, Newton on the shell alone
struct Cyl {
    int order, ne, nn;
    std::vector<double> z;
    std::vector<double> state_c, state_t;
};

double expand(int order, const std::vector<double>& path, std::vector<double>& reaction,
              std::vector<double>& nt_min, std::vector<double>& nt_max, std::vector<double>& ns_max) {
    constexpr double R = 2.0, H = 4.0;
    const pl::PlateProps P = shell(0.3);
    const int ne = order == 3 ? 16 : 8, per = order - 1, nn = ne * per + 1, nd = 3 * nn;
    const int S = order == 3 ? pl::kPlasticStateSize : pl::kPlasticStateSize5;
    std::vector<double> state_c((size_t)(ne * S), 0.0), state_t = state_c;
    Eigen::VectorXd u = Eigen::VectorXd::Zero(nd);
    const auto nodes_of = [&](int e) {
        std::vector<int> nodes;
        if (order == 3) nodes = {e * 2, e * 2 + 2, e * 2 + 1};
        else nodes = {e * 4, e * 4 + 1, e * 4 + 2, e * 4 + 3, e * 4 + 4};
        return nodes;
    };
    const auto zof = [&](int n) { return H * n / (nn - 1); };
    // assemble f_int, K at u from state_c into state_t; per element
    const auto assemble = [&](Eigen::VectorXd& f, Eigen::MatrixXd& K, std::vector<double>& nts,
                              std::vector<double>& nss) {
        f = Eigen::VectorXd::Zero(nd);
        K = Eigen::MatrixXd::Zero(nd, nd);
        nts.clear(); nss.clear();
        for (int e = 0; e < ne; ++e) {
            const auto nodes = nodes_of(e);
            const double* sc = state_c.data() + (size_t)(e * S);
            double* st = state_t.data() + (size_t)(e * S);
            if (order == 3) {
                pl::NodeCoords X; pl::Dof ue, fe; pl::ElementMatrix Ke;
                for (int k = 0; k < 3; ++k) {
                    X(k, 0) = R; X(k, 1) = zof(nodes[k]);
                    for (int c = 0; c < 3; ++c) ue(3 * k + c) = u(3 * nodes[k] + c);
                }
                pl::internal_force_plastic(X, P, ue, sc, st, fe, &Ke);
                std::array<Eigen::Vector4d, 3> sg;
                pl::gauss_forces_shell(X, P, ue, sc, true, sg);
                for (const auto& v : sg) { nts.push_back(v(1)); nss.push_back(v(0)); }
                for (int a = 0; a < 9; ++a) {
                    f(3 * nodes[a / 3] + a % 3) += fe(a);
                    for (int b = 0; b < 9; ++b) K(3 * nodes[a / 3] + a % 3, 3 * nodes[b / 3] + b % 3) += Ke(a, b);
                }
            } else {
                pl::NodeCoords5 X; pl::Dof5 ue, fe; pl::ElementMatrix5 Ke;
                for (int k = 0; k < 5; ++k) {
                    X(k, 0) = R; X(k, 1) = zof(nodes[k]);
                    for (int c = 0; c < 3; ++c) ue(3 * k + c) = u(3 * nodes[k] + c);
                }
                pl::internal_force_plastic5(X, P, ue, sc, st, fe, &Ke);
                std::array<Eigen::Vector4d, 5> sg;
                pl::gauss_forces_shell5(X, P, ue, sc, true, sg);
                for (const auto& v : sg) { nts.push_back(v(1)); nss.push_back(v(0)); }
                for (int a = 0; a < 15; ++a) {
                    f(3 * nodes[a / 3] + a % 3) += fe(a);
                    for (int b = 0; b < 15; ++b) K(3 * nodes[a / 3] + a % 3, 3 * nodes[b / 3] + b % 3) += Ke(a, b);
                }
            }
        }
    };
    // free: u_z (but the bottom node's) and every rotation; u_r prescribed at every node
    std::vector<int> fr;
    for (int n = 0; n < nn; ++n) {
        if (n > 0) fr.push_back(3 * n + 1);
        fr.push_back(3 * n + 2);
    }
    double worst_resid = 0.0;
    for (double ur : path) {
        for (int n = 0; n < nn; ++n) u(3 * n) = ur;
        Eigen::VectorXd f;
        Eigen::MatrixXd K;
        std::vector<double> nts, nss;
        for (int it = 0; it < 30; ++it) {
            assemble(f, K, nts, nss);
            Eigen::VectorXd rf((int)fr.size());
            Eigen::MatrixXd Kf((int)fr.size(), (int)fr.size());
            for (size_t i = 0; i < fr.size(); ++i) {
                rf((int)i) = f(fr[i]);
                for (size_t j = 0; j < fr.size(); ++j) Kf((int)i, (int)j) = K(fr[i], fr[j]);
            }
            if (rf.norm() < 1e-11 * P.Np) break;
            const Eigen::VectorXd du = Kf.fullPivLu().solve(-rf);
            for (size_t i = 0; i < fr.size(); ++i) u(fr[i]) += du((int)i);
        }
        assemble(f, K, nts, nss);
        double rfree = 0.0, react = 0.0;
        for (int d : fr) rfree = std::max(rfree, std::fabs(f(d)));
        for (int n = 0; n < nn; ++n) react += f(3 * n);
        worst_resid = std::max(worst_resid, rfree);
        reaction.push_back(react);
        nt_min.push_back(*std::min_element(nts.begin(), nts.end()));
        nt_max.push_back(*std::max_element(nts.begin(), nts.end()));
        double ns = 0.0;
        for (double v : nss) ns = std::max(ns, std::fabs(v));
        ns_max.push_back(ns);
        state_c = state_t;   // commit
    }
    return worst_resid;
}

void test_cylinder() {
    std::printf("-- (c) a free cylinder expanded through its hoop capacity --\n");
    constexpr double R = 2.0, H = 4.0, EA = 1.0e6, nu = 0.3, Np = 100.0;
    const double k_el = EA * (1.0 - nu * nu) / R;   // N_t per unit u_r with N_s = 0
    const double uy = Np / k_el;                     // the yield displacement
    const std::vector<double> path = {0.5 * uy, 1.5 * uy, 2.25 * uy, 3.0 * uy, 2.9 * uy};
    for (int order : {3, 5}) {
        std::vector<double> react, ntmin, ntmax, nsmax;
        const double resid = expand(order, path, react, ntmin, ntmax, nsmax);
        std::printf("   %d-noded: reaction %.9f / %.9f / %.9f / %.9f / %.9f kN/rad (Np H %.1f), residual %.1e\n",
                    order, react[0], react[1], react[2], react[3], react[4], Np * H, resid);
        std::printf("   %d-noded: N_t at 3 uy %.12f..%.12f, max |N_s| %.2e; after unloading %.9f (want %.9f)\n",
                    order, ntmin[3], ntmax[3], nsmax[3], ntmax[4], Np - k_el * 0.1 * uy);
        const auto rel = [](double a, double b) { return std::fabs(a - b) / std::fabs(b); };
        char lbl[200];
        std::snprintf(lbl, sizeof lbl, "(c) %d-noded: elastic, the reaction is EA (1 - nu^2) (u_r / R) H", order);
        check(rel(react[0], k_el * 0.5 * uy * H) < 1e-9 && rel(ntmax[0], 0.5 * Np) < 1e-9, lbl);
        std::snprintf(lbl, sizeof lbl, "(c) %d-noded: yielded, N_t = Np at every point, N_s = 0, reaction Np H", order);
        bool at_cap = true;
        for (int k = 1; k <= 3; ++k)
            at_cap = at_cap && rel(ntmin[(size_t)k], Np) < 1e-9 && rel(ntmax[(size_t)k], Np) < 1e-9 &&
                     nsmax[(size_t)k] < 1e-9 * Np && rel(react[(size_t)k], Np * H) < 1e-9;
        check(at_cap && resid < 1e-9 * Np, lbl);
        std::snprintf(lbl, sizeof lbl, "(c) %d-noded: unloading by du_r lowers N_t by EA (1 - nu^2) du_r / R", order);
        check(rel(ntmax[4], Np - k_el * 0.1 * uy) < 1e-9 && rel(ntmin[4], Np - k_el * 0.1 * uy) < 1e-9, lbl);
    }
}

}  // namespace

int main() {
    std::printf("== A shell of revolution with a capacity (KV-STR-034) ==\n");
    test_return();
    test_independent_caps();
    test_cylinder();
    std::printf(g_failures ? "\n%d CHECK(S) FAILED\n" : "\nall checks passed\n", g_failures);
    return g_failures ? 1 : 0;
}
