// Axisymmetric Mohr-Coulomb at the material point (P1.6): the return mapping is
// reused from plane strain (same principal structure: r-z block + hoop sigma_theta),
// and the 4x4 consistent tangent is D_T = Psi * D_e. We verify (i) inadmissible
// trials return to the yield surface (f = 0) and (ii) the analytic 4x4 tangent
// equals a central finite difference over the 4-strain space (associated and
// non-associated), the decisive check of the hoop coupling in the tangent.
#include <katai/materials/material_model.hpp>
#include <katai/materials/mohr_coulomb.hpp>

#include <cmath>
#include <cstdio>
#include <random>

using katai::core::GaussState;
using katai::core::integrate_point_axisym;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::mc_yield;
using katai::core::PlaneStrainStress;
using katai::core::principal_stresses;

namespace {

int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
double deg(double d) { return d * 3.14159265358979323846 / 180.0; }

// sigma (4-vector [r,z,rz,theta]) from the axisymmetric material update.
Eigen::Vector4d stress_of(const MaterialModel& m, const GaussState& committed,
                          const Eigen::Vector4d& de) {
    GaussState t; Eigen::Matrix4d dummy;
    integrate_point_axisym(m, committed, de, t, dummy);
    Eigen::Vector4d s; s << t.stress, t.stress_zz; return s;
}

void test_return_to_surface() {
    const double c = 10.0, phi = deg(30.0);
    const MaterialModel m{MaterialType::MohrCoulomb, 20000.0, 0.3, c, phi, deg(10.0)};
    GaussState committed;  // virgin
    // A strongly deviatoric strain increment (with r-z shear) must yield.
    const Eigen::Vector4d de(3e-3, -3e-3, 1e-3, 0.0);
    GaussState trial; Eigen::Matrix4d T;
    integrate_point_axisym(m, committed, de, trial, T);
    PlaneStrainStress s; s.in_plane = trial.stress; s.zz = trial.stress_zz;
    check(std::fabs(mc_yield(principal_stresses(s), c, phi)) < 1e-7,
          "axisym MC: inadmissible trial returns to yield surface (f=0)");
}

// Analytic 4x4 tangent vs central finite difference, over many random plastic
// states (associated and non-associated), filtering samples whose stencil straddles
// a region boundary (forward and backward slopes disagree).
void test_tangent_vs_fd() {
    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    const double phi = deg(25.0);
    int nplastic = 0, nfail = 0;
    double worst = 0.0;
    for (double psi : {phi, deg(5.0)}) {
        const MaterialModel m{MaterialType::MohrCoulomb, 1.0e5, 0.3, 5.0, phi, psi};
        const MaterialModel le{MaterialType::LinearElastic, 1.0e5, 0.3, 0, 0, 0};
        for (int t = 0; t < 8000; ++t) {
            GaussState committed;
            committed.stress = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 40.0;
            committed.stress_zz = u(rng) * 40.0;
            const double scale = std::pow(10.0, u(rng) * 2.0 - 2.5);
            const Eigen::Vector4d de(u(rng), u(rng), u(rng), u(rng));
            const Eigen::Vector4d dstrain = de * scale;

            GaussState trial; Eigen::Matrix4d T;
            integrate_point_axisym(m, committed, dstrain, trial, T);
            GaussState el; Eigen::Matrix4d dummy;
            integrate_point_axisym(le, committed, dstrain, el, dummy);
            if ((trial.stress - el.stress).cwiseAbs().maxCoeff() +
                    std::fabs(trial.stress_zz - el.stress_zz) < 1e-9)
                continue;  // elastic
            ++nplastic;

            constexpr double h = 5e-7;
            Eigen::Matrix4d fwd, bwd;
            const Eigen::Vector4d s0 = stress_of(m, committed, dstrain);
            for (int j = 0; j < 4; ++j) {
                Eigen::Vector4d ep = dstrain, em = dstrain; ep(j) += h; em(j) -= h;
                fwd.col(j) = (stress_of(m, committed, ep) - s0) / h;
                bwd.col(j) = (s0 - stress_of(m, committed, em)) / h;
            }
            const double sc = T.cwiseAbs().maxCoeff() + 1.0;
            if ((fwd - bwd).cwiseAbs().maxCoeff() > 1e-3 * sc) continue;  // straddle
            const double rel = (T - 0.5 * (fwd + bwd)).cwiseAbs().maxCoeff() / sc;
            worst = std::max(worst, rel);
            if (rel > 1e-3) ++nfail;
        }
    }
    std::printf("  axisym MC tangent: %d plastic samples, %d fail, worst rel=%.2e\n",
                nplastic, nfail, worst);
    check(nfail == 0, "axisym MC: analytic 4x4 tangent matches central FD");
}

} // namespace

int main() {
    test_return_to_surface();
    test_tangent_vs_fd();
    if (g_failures == 0) {
        std::printf("OK: axisymmetric Mohr-Coulomb material point (P1.6) verified\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
