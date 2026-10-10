// HSsmall -- the tangent is the derivative of the stress, and the stiffness of a step can be read
// from another increment than the one the stress is integrated over.
//
// In this model the stiffness of a step is set by the direction of the step's strain increment
// against the strain history (Benz 2007). The stress is then a function of the increment whose
// slope contains the slope of the stiffness, of order 1 / |d eps| where the increment is small
// beside the history. Two statements are checked at random stress points with a strain history,
// for increments from 1e-7 to 3e-4, in plane strain and in axisymmetry (below 1e-8 the slope of the
// stiffness is left out of the tangent on purpose: it is of order 1 / |d eps| and multiplies the
// rounding of the stress there):
//
//   (a) the tangent handed back agrees with the central difference of the stress, taken with a
//       step that is a fraction of the increment. Without the slope of the stiffness -- and with
//       a difference taken at a fixed 1.5e-8, as large as the increments it differences -- a
//       quarter of the points were more than 10% off in plane strain (748 of 3000), and 59% were
//       within 1%; it is 99.4% and 98.9% now. What is left are the points where the walk meets
//       a surface inside the increment or the increment's principal directions are not unique.
//   (b) with the stiffness read from a FIXED increment (small_strain_increment) the stress is
//       the same when that increment is the increment itself, the history is advanced by the
//       increment the stress was integrated over whichever one the stiffness was read from, and
//       the stress no longer follows the direction of the increment: the tangent without any
//       slope of the stiffness is then the derivative (99.3% of the points within 1%). That
//       needs the difference of the plastic walk taken with a step below the increment as well:
//       with the fixed step it was 90%, every miss at an increment below 1e-7.
// verify: KV-NUM-014
//   oracle:   independent_path
//   source:   the central difference of the returned stress with respect to the strain increment, taken with a step of 1e-5 of the increment's norm along a random direction, against the tangent the same call hands back (Hardening Soil with small-strain stiffness after Benz, T. (2007), Small-Strain Stiffness of Soils and its Numerical Consequences, Mitteilung 55, Institut fuer Geotechnik, Universitaet Stuttgart, Sec. 6.1 and App. A: the stiffness of a step depends on the direction of the step's strain increment)
//   locator:  integrate_point / integrate_point_axisym at 2000 random stress points each (sigma_v 50-300 kPa, K0 0.4-0.9, a strain history of 1e-6 to 1e-3, increments of 1e-7 to 3e-4; below 1e-8 the slope is left out on purpose), with the stiffness read from the increment itself and from a second, fixed increment
//   quantity: share of the points at which |D d - (sigma(de + h d) - sigma(de - h d)) / 2h| < 1% of the difference [-]
//   expected: all of them, less the points where the walk meets a yield surface inside the increment or the increment's principal directions are not unique
//   band:     more than 97% with the stiffness read from the increment (measured 99.4% in plane strain, 98.9% in axisymmetry; 59% before the tangent carried the slope of the stiffness) and more than 99% with it read from a fixed increment (measured 99.3% and 99.4%; 90% while the plastic walk was differenced with a fixed step of 1.5e-8)
#include <katai/materials/hardening_soil.hpp>
#include <katai/materials/material_model.hpp>

#include <cmath>
#include <cstdio>
#include <random>

#include <Eigen/Dense>

using namespace katai::core;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf("%s %s\n", ok ? "ok:  " : "FAIL:", what);
    if (!ok) ++g_failures;
}

MaterialModel sand() {
    HardeningSoilParams p;
    p.E50_ref = 25e3; p.Eoed_ref = 25e3; p.Eur_ref = 75e3; p.m = 0.5; p.p_ref = 100; p.nu_ur = 0.2;
    p.Rf = 0.9; p.G0_ref = 80e3; p.gamma07 = 0.2e-3;
    p.friction = 31 * 3.14159265358979 / 180; p.dilatancy = 1 * 3.14159265358979 / 180;
    p.cohesion = 0.0; p.K0nc = 1 - std::sin(p.friction);
    MaterialModel m{MaterialType::HardeningSoil, 0, 0};
    m.hs = p;
    return m;
}

GaussState state(std::mt19937& rng, std::uniform_real_distribution<double>& U) {
    GaussState gs;
    const double sv = -(50.0 + 250.0 * std::fabs(U(rng)));
    const double k0 = 0.4 + 0.5 * std::fabs(U(rng));
    gs.stress << k0 * sv, sv, 0.15 * sv * U(rng);
    gs.stress_zz = k0 * sv;
    gs.pp = 1.3 * std::fabs(sv);
    const double hmag = std::pow(10.0, -6.0 + 3.0 * std::fabs(U(rng)));
    gs.hs_hist << hmag * U(rng), hmag * U(rng), hmag * U(rng), hmag * U(rng);
    gs.gamma_hist = hmag;
    return gs;
}

// The share of points at which the tangent is within 1% of the central difference of the stress
// along a random direction; `frozen` reads the stiffness from a second, fixed increment.
template <int N>
double share_within(bool axisym, bool frozen, unsigned seed) {
    using Vec = Eigen::Matrix<double, N, 1>;
    using Mat = Eigen::Matrix<double, N, N>;
    const MaterialModel m = sand();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    int good = 0;
    const int n = 2000;
    const auto run = [&](const GaussState& gs, const Vec& de, const Vec* small, Mat& D,
                         TangentMode mode) {
        GaussState tr;
        Vec s;
        if constexpr (N == 3) {
            integrate_point(m, gs, de, tr, D, mode, 0.0, nullptr, 0.0, small);
            s = tr.stress;
        } else {
            integrate_point_axisym(m, gs, de, tr, D, mode, 0.0, nullptr, 0.0, small);
            s << tr.stress, tr.stress_zz;
        }
        return s;
    };
    (void)axisym;
    for (int c = 0; c < n; ++c) {
        const GaussState gs = state(rng, U);
        const double dmag = std::pow(10.0, -7.0 + 3.5 * std::fabs(U(rng)));
        Vec de, dir, fixed;
        for (int k = 0; k < N; ++k) { de(k) = dmag * U(rng); dir(k) = U(rng); fixed(k) = dmag * U(rng); }
        dir.normalize();
        const Vec* small = frozen ? &fixed : nullptr;
        Mat D, Dd;
        run(gs, de, small, D, TangentMode::kConsistent);
        const double h = 1.0e-5 * de.norm();
        const Vec fd = (run(gs, de + h * dir, small, Dd, TangentMode::kNone) -
                        run(gs, de - h * dir, small, Dd, TangentMode::kNone)) / (2.0 * h);
        if ((D * dir - fd).norm() < 1.0e-2 * std::max(fd.norm(), 1e-12)) ++good;
    }
    return (double)good / n;
}
}  // namespace

int main() {
    char line[256];

    // (a) the tangent against the difference of the stress.
    const double ps = share_within<3>(false, false, 11), ax = share_within<4>(true, false, 12);
    std::snprintf(line, sizeof(line),
                  "(a) the tangent is within 1%% of the stress difference at %.1f%% of the points in "
                  "plane strain and %.1f%% in axisymmetry (> 97%%)", 100 * ps, 100 * ax);
    check(ps > 0.97 && ax > 0.97, line);

    // (b) the stiffness read from a fixed increment.
    const MaterialModel m = sand();
    std::mt19937 rng(21);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    double worst_same = 0.0, worst_hist = 0.0;
    int differs = 0;
    for (int c = 0; c < 400; ++c) {
        const GaussState gs = state(rng, U);
        const double dmag = std::pow(10.0, -7.0 + 3.0 * std::fabs(U(rng)));
        const Eigen::Vector3d de(dmag * U(rng), dmag * U(rng), dmag * U(rng));
        const Eigen::Vector3d other(dmag * U(rng), dmag * U(rng), dmag * U(rng));
        GaussState t0, t1, t2;
        Eigen::Matrix3d D;
        integrate_point(m, gs, de, t0, D, TangentMode::kNone, 0.0);
        integrate_point(m, gs, de, t1, D, TangentMode::kNone, 0.0, nullptr, 0.0, &de);
        integrate_point(m, gs, de, t2, D, TangentMode::kNone, 0.0, nullptr, 0.0, &other);
        worst_same = std::max(worst_same, (t1.stress - t0.stress).norm() /
                                              std::max((t0.stress - gs.stress).norm(), 1e-300));
        worst_hist = std::max(worst_hist, (t2.hs_hist - t0.hs_hist).norm() /
                                              std::max(t0.hs_hist.norm(), 1e-300));
        if ((t2.stress - t0.stress).norm() > 1e-3 * (t0.stress - gs.stress).norm()) ++differs;
    }
    std::snprintf(line, sizeof(line),
                  "(b) read from the increment itself the stress is the same (worst %.1e of the "
                  "stress increment)", worst_same);
    check(worst_same < 1e-12, line);
    std::snprintf(line, sizeof(line),
                  "(b) the history is advanced by the increment the stress was integrated over, "
                  "whichever one the stiffness was read from (worst %.1e)", worst_hist);
    check(worst_hist < 1e-12, line);
    std::snprintf(line, sizeof(line),
                  "(b) read from another increment the stress is another one at %d of 400 points "
                  "(the stiffness does depend on the direction)", differs);
    check(differs > 100, line);
    const double fps = share_within<3>(false, true, 31), fax = share_within<4>(true, true, 32);
    std::snprintf(line, sizeof(line),
                  "(b) and the tangent, with no slope of the stiffness in it, is within 1%% of the "
                  "stress difference at %.1f%% and %.1f%% of the points (> 99%%)", 100 * fps, 100 * fax);
    check(fps > 0.99 && fax > 0.99, line);

    return g_failures == 0 ? 0 : 1;
}
