// van Genuchten (1980) SWCC + Mualem (1976) k_rel — the unsaturated-flow water-retention
// core (W2a). Verified: (1) saturated limit ψ≤0 → S_e=1,k_rel=1,dS/dψ=0; (2) monotonicity;
// (3) asymptotic S_e ~ (g_a·ψ)^{−(g_n−1)} (large suction); (4) Mualem-vG special case n=2
// closed-form cross-check; (5) CRITICAL: analytic moisture capacity dS/dψ = central finite
// difference (algorithmic correctness).
// Kaynaklar: van Genuchten (1980) SSSAJ 44:892; Mualem (1976) WRR 12:513; PLAXIS MMM.
#include <katai/materials/water_retention.hpp>

#include <cmath>
#include <cstdio>

using katai::core::effective_saturation;
using katai::core::moisture_capacity;
using katai::core::relative_permeability;
using katai::core::relative_permeability_psi;
using katai::core::saturation;
using katai::core::WaterRetention;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}

void test_saturated_limit() {
    WaterRetention w;  // g_a=2, g_n=2, g_l=0.5
    for (double psi : {-1.0, -0.001, 0.0}) {
        check(effective_saturation(w, psi) == 1.0, "S_e=1 for psi<=0 (saturated)");
        check(saturation(w, psi) == w.S_sat, "S=S_sat for psi<=0");
        check(moisture_capacity(w, psi) == 0.0, "dS/dpsi=0 for psi<=0");
        check(relative_permeability_psi(w, psi) == 1.0, "k_rel=1 for psi<=0");
    }
    std::printf("  (1) saturated limit (psi<=0): S_e=1, k_rel=1, dS/dpsi=0  OK\n");
}

void test_monotonic_limits() {
    WaterRetention w{1.5, 2.5, 0.5, 0.05, 1.0};
    double prevSe = 2.0, prevkr = 2.0;
    for (double psi = 0.01; psi < 100.0; psi *= 1.3) {
        const double Se = effective_saturation(w, psi);
        const double kr = relative_permeability_psi(w, psi);
        check(Se > 0.0 && Se <= 1.0, "S_e in (0,1]");
        check(kr >= w.k_rel_min - 1e-15 && kr <= 1.0 + 1e-15, "k_rel in [k_rel_min,1]");
        check(Se <= prevSe + 1e-12, "S_e monotone decreasing in psi");
        check(kr <= prevkr + 1e-12, "k_rel monotone decreasing in psi");
        prevSe = Se; prevkr = kr;
    }
    // k_rel(saturated)=1 exactly at S_e=1.
    check(std::fabs(relative_permeability(w, 1.0) - 1.0) < 1e-12, "k_rel(S_e=1)=1");
    std::printf("  (2) monotonic + bounded S_e/k_rel  OK\n");
}

void test_asymptotic_slope() {
    // Large suction: S_e ≈ (g_a·ψ)^{−g_n·g_m} = (g_a·ψ)^{−(g_n−1)}; log-log slope −(g_n−1).
    WaterRetention w{2.0, 3.0, 0.5, 0.0, 1.0};  // n=3 → slope −2
    const double psi1 = 50.0, psi2 = 100.0;
    const double s1 = effective_saturation(w, psi1), s2 = effective_saturation(w, psi2);
    const double slope = std::log(s2 / s1) / std::log(psi2 / psi1);
    const double expected = -(w.g_n - 1.0);
    std::printf("  (3) asymptotic log-log slope = %.4f (expected %.4f)\n", slope, expected);
    check(std::fabs(slope - expected) < 0.02, "large-suction S_e slope = -(g_n-1)");
}

void test_mualem_special_case() {
    // n=2 → m=1/2: k_rel = S_e^l·[1 − (1 − S_e²)^{1/2}]².  S_e=0.5, l=0.5 → 0.0126918 (elle).
    WaterRetention w{2.0, 2.0, 0.5, 0.0, 1.0};
    const double Se = 0.5;
    const double s2 = Se * Se;                       // 0.25
    const double hand = std::pow(Se, 0.5) * std::pow(1.0 - std::sqrt(1.0 - s2), 2.0);
    const double kr = relative_permeability(w, Se);
    std::printf("  (4) Mualem n=2 special case: k_rel(S_e=0.5)=%.7f (hand %.7f)\n", kr, hand);
    check(std::fabs(kr - hand) < 1e-12, "Mualem n=2 closed form matches");
    check(std::fabs(kr - 0.0126918) < 1e-6, "Mualem n=2 known value 0.0126918");
}

void test_capacity_vs_fd() {
    // CRITICAL: analytic dS/dψ = central finite difference of saturation(ψ).
    const WaterRetention sets[] = {
        {2.0, 2.0, 0.5, 0.0, 1.0},
        {1.0, 1.8, 0.5, 0.10, 0.95},
        {4.0, 3.5, 0.4, 0.05, 1.0},
    };
    std::printf("  (5) moisture capacity vs central finite difference:\n");
    int checked = 0;
    for (const auto& w : sets)
        for (double psi : {0.05, 0.2, 0.5, 1.0, 3.0, 10.0}) {
            const double h = 1.0e-6 * psi;
            const double fd = (saturation(w, psi + h) - saturation(w, psi - h)) / (2.0 * h);
            const double an = moisture_capacity(w, psi);
            const double rel = std::fabs(an - fd) / (std::fabs(fd) + 1e-14);
            if (psi == 1.0)
                std::printf("    g_a=%.1f g_n=%.1f psi=%.2f  analytic=%.6e  FD=%.6e  rel=%.2e\n",
                            w.g_a, w.g_n, psi, an, fd, rel);
            check(rel < 1e-5, "analytic dS/dpsi matches central FD within 1e-5");
            ++checked;
        }
    check(checked >= 18, "(5) checked several params/suctions");
}

}  // namespace

int main() {
    test_saturated_limit();
    test_monotonic_limits();
    test_asymptotic_slope();
    test_mualem_special_case();
    test_capacity_vs_fd();
    if (g_failures == 0) {
        std::printf("OK: van Genuchten/Mualem water retention verified (limits, monotonicity, "
                    "asymptote, Mualem special case, analytic capacity = FD)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
