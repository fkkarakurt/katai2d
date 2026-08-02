// Dynamics core (seismic track D1): Newmark-beta time integration + Rayleigh damping, verified
// against closed-form structural dynamics and 1D free-field site-response theory
// (docs/references/dynamic-seismic-formulation.md; Chopra, Dynamics of Structures; Kramer,
// Geotechnical Earthquake Engineering ch.7). This is the reusable engine that will later drive the
// 2D FE system; here it is pinned on exact 1D theory:
//   (a) SDOF oscillator: free-vibration period, damped log-decrement, forced steady-state amplitude.
//   (b) 1D shear column: natural frequencies f_n = (2n-1) Vs/(4H) via K phi = w^2 M phi (checks M, K).
//   (c) base-harmonic transfer function |u_t,surface / a_g|: sub-resonant = 1/cos(wH/Vs) and the
//       fundamental resonance peak = 2/(pi xi) -- the whole Newmark + Rayleigh + base-excitation chain.
#include <katai/analysis/dynamics.hpp>
#include <katai/analysis/free_field.hpp>
#include <katai/fem/assembly/assembler.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/materials/linear_elastic.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/math/sparse_matrix.hpp>
#include <katai/mesh/mesh.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::NewmarkResult;
using katai::core::rayleigh_from_modes;
using katai::core::solve_newmark;
using katai::math::CsrMatrix;
using katai::math::SparseMatrixBuilder;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

CsrMatrix scalar_mat(double v) { SparseMatrixBuilder b(1, 1); b.add_entry(0, 0, v); return b.build(); }
CsrMatrix zero_mat(int n) {
    SparseMatrixBuilder b(n, n);
    for (int i = 0; i < n; ++i) b.add_entry(i, i, 0.0);  // structural zeros -> valid n x n, C*v = 0
    return b.build();
}
Eigen::MatrixXd dense(const CsrMatrix& A) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(A.rows, A.cols);
    for (int r = 0; r < A.rows; ++r)
        for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p) D(r, A.col_indices[p]) = A.values[p];
    return D;
}
// alpha*X + beta*Y assembled as a CsrMatrix (both n x n, same free-dof layout).
CsrMatrix combine(const CsrMatrix& X, double a, const CsrMatrix& Y, double b) {
    SparseMatrixBuilder B(X.rows, X.cols);
    auto add = [&](const CsrMatrix& A, double s) {
        for (int r = 0; r < A.rows; ++r)
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p) B.add_entry(r, A.col_indices[p], s * A.values[p]);
    };
    add(X, a); add(Y, b);
    return B.build();
}

// (a) Single-degree-of-freedom oscillator: the simplest exact oracle for Newmark + Rayleigh.
void test_sdof() {
    std::printf("-- (a) SDOF oscillator (Newmark vs exact) --\n");
    const double m = 1.0, fn = 2.0, wn = 2 * kPi * fn, k = wn * wn * m, xi = 0.05;
    const double c = 2 * xi * std::sqrt(k * m), T = 1.0 / fn;
    const CsrMatrix M = scalar_mat(m), K = scalar_mat(k);
    Eigen::VectorXd one(1); one << 1.0;
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(1);
    auto no_force = [&](int) { return z; };

    // Free UNDAMPED vibration: period = T, amplitude conserved (gamma=1/2 -> no algorithmic damping).
    {
        const double dt = T / 200; const int ns = 200 * 6;
        const auto R = solve_newmark(M, zero_mat(1), K, no_force, dt, ns, one, z);
        std::vector<double> peaks_t;
        for (size_t i = 1; i + 1 < R.u.size(); ++i)
            if (R.u[i][0] > R.u[i - 1][0] && R.u[i][0] >= R.u[i + 1][0]) peaks_t.push_back(R.t[i]);
        const double Tmeas = peaks_t.size() >= 2 ? (peaks_t.back() - peaks_t.front()) / (peaks_t.size() - 1) : 0.0;
        double amp_end = 0.0;
        for (size_t i = R.u.size() * 3 / 4; i < R.u.size(); ++i) amp_end = std::max(amp_end, std::fabs(R.u[i][0]));
        std::printf("   free:   T_meas=%.5f  T_exact=%.5f (err %.2f%%)  amp_end=%.4f (u0=1)\n",
                    Tmeas, T, 100 * (Tmeas - T) / T, amp_end);
        check(std::fabs(Tmeas - T) < 0.01 * T, "free-vibration period = 2 pi sqrt(m/k)");
        check(std::fabs(amp_end - 1.0) < 0.02, "no algorithmic damping (amplitude conserved, gamma=1/2)");
    }
    // Damped free vibration: log-decrement delta = ln(peak_i/peak_{i+1}) = 2 pi xi / sqrt(1-xi^2).
    {
        const CsrMatrix C = scalar_mat(c);
        const double dt = T / 400; const int ns = 400 * 8;
        const auto R = solve_newmark(M, C, K, no_force, dt, ns, one, z);
        std::vector<double> pk;
        for (size_t i = 1; i + 1 < R.u.size(); ++i)
            if (R.u[i][0] > R.u[i - 1][0] && R.u[i][0] >= R.u[i + 1][0] && R.u[i][0] > 0) pk.push_back(R.u[i][0]);
        double dmeas = 0.0; int cnt = 0;
        for (size_t i = 0; i + 1 < pk.size(); ++i) { dmeas += std::log(pk[i] / pk[i + 1]); ++cnt; }
        dmeas = cnt ? dmeas / cnt : 0.0;
        const double dexact = 2 * kPi * xi / std::sqrt(1 - xi * xi);
        std::printf("   damped: log-decrement meas=%.5f  exact=%.5f (err %.2f%%)\n",
                    dmeas, dexact, 100 * (dmeas - dexact) / dexact);
        check(std::fabs(dmeas - dexact) < 0.03 * dexact, "damped log-decrement = 2 pi xi / sqrt(1-xi^2)");
    }
    // Forced steady-state: F = F0 sin(w t), r = w/wn = 0.8. amplitude = (F0/k) Rd.
    {
        const CsrMatrix C = scalar_mat(c);
        const double r = 0.8, w = r * wn, F0 = 10.0, dt = (2 * kPi / w) / 200; const int ns = 200 * 40;
        auto force = [&](int step) { Eigen::VectorXd f(1); f << F0 * std::sin(w * step * dt); return f; };
        const auto R = solve_newmark(M, C, K, force, dt, ns, z, z);
        double amp = 0.0;
        for (size_t i = R.u.size() * 3 / 4; i < R.u.size(); ++i) amp = std::max(amp, std::fabs(R.u[i][0]));
        const double Rd = 1.0 / std::sqrt((1 - r * r) * (1 - r * r) + (2 * xi * r) * (2 * xi * r));
        const double amp_exact = (F0 / k) * Rd;
        std::printf("   forced: amp_meas=%.5e  amp_exact=%.5e  Rd=%.3f (err %.2f%%)\n",
                    amp, amp_exact, Rd, 100 * (amp - amp_exact) / amp_exact);
        check(std::fabs(amp - amp_exact) < 0.02 * amp_exact, "steady-state amplitude = (F0/k) Rd(r,xi)");
    }
}

// 1D shear column: nodes 0..ne up the depth, base (node 0) fixed. Free dofs = nodes 1..ne (index
// node-1). 2-node linear shear element: k = (G/L)[[1,-1],[-1,1]], consistent mass m = (rho L/6)[[2,1],[1,2]].
struct Column { CsrMatrix M, K; int nd; double H, Vs; };
Column shear_column(double H, double Vs, double rho, int ne) {
    const double L = H / ne, G = rho * Vs * Vs, ke = G / L, me = rho * L / 6.0;
    SparseMatrixBuilder bM(ne, ne), bK(ne, ne);
    const double Ke[2][2] = {{ke, -ke}, {-ke, ke}};
    const double Me[2][2] = {{2 * me, me}, {me, 2 * me}};
    for (int e = 0; e < ne; ++e) {
        const int node[2] = {e, e + 1};
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) {
                const int di = node[i] - 1, dj = node[j] - 1;
                if (di < 0 || dj < 0) continue;
                bK.add_entry(di, dj, Ke[i][j]);
                bM.add_entry(di, dj, Me[i][j]);
            }
    }
    return {bM.build(), bK.build(), ne, H, Vs};
}

// (b) Natural frequencies from the generalized eigenproblem K phi = w^2 M phi (checks M and K only).
void test_column_frequencies() {
    std::printf("-- (b) 1D shear column natural frequencies (M, K assembly) --\n");
    const double H = 30.0, Vs = 200.0, rho = 2.0; const int ne = 40;
    const Column col = shear_column(H, Vs, rho, ne);
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(col.K), dense(col.M));
    const Eigen::VectorXd lam = es.eigenvalues();  // ascending
    for (int n = 1; n <= 3; ++n) {
        const double f_fe = std::sqrt(lam[n - 1]) / (2 * kPi);
        const double f_th = (2 * n - 1) * Vs / (4 * H);
        std::printf("   f_%d: FE=%.4f Hz  theory (2n-1)Vs/4H=%.4f Hz  (err %.2f%%)\n",
                    n, f_fe, f_th, 100 * (f_fe - f_th) / f_th);
        const double tol = n == 1 ? 0.01 : (n == 2 ? 0.02 : 0.05);  // consistent-mass error grows with n
        check(std::fabs(f_fe - f_th) < tol * f_th, "natural frequency f_n = (2n-1) Vs/(4H)");
    }
}

// Steady-state surface amplification |u_t,surface / a_g| for a rigid-base harmonic acceleration
// a_g(t) = A sin(w t): effective earthquake force F = -M r a_g (r = 1 on all free dofs); Rayleigh
// damping with ratio xi targeted at (fr1, fr2). Total surface accel = relative surface accel + a_g.
double steady_amplification(const Column& col, double f_drive, double xi, double fr1, double fr2,
                            int periods, int steps_per) {
    const double w = 2 * kPi * f_drive, A = 1.0, dt = (1.0 / f_drive) / steps_per;
    const int ns = periods * steps_per;
    const auto ray = rayleigh_from_modes(fr1, xi, fr2, xi);
    const CsrMatrix C = combine(col.M, ray.alpha, col.K, ray.beta);
    const Eigen::VectorXd Mr = col.M * Eigen::VectorXd::Ones(col.nd);  // M r,  r = 1
    auto force = [&](int step) {
        Eigen::VectorXd f = -(A * std::sin(w * step * dt)) * Mr;  // -M r a_g
        return f;
    };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(col.nd);
    const auto R = solve_newmark(col.M, C, col.K, force, dt, ns, z, z);
    double amp = 0.0;
    for (size_t i = R.u.size() * 3 / 4; i < R.u.size(); ++i) {
        const double ag = A * std::sin(w * R.t[i]);
        amp = std::max(amp, std::fabs(R.a[i][col.nd - 1] + ag));  // total surface accel; surface dof = nd-1
    }
    return amp / A;
}

// (c) Base-harmonic transfer function: exercises the whole chain (Newmark + Rayleigh + base motion).
void test_column_transfer() {
    std::printf("-- (c) 1D column base-harmonic transfer function (whole chain) --\n");
    const double H = 30.0, Vs = 200.0, rho = 2.0; const int ne = 40, xi_pct = 5;
    const double xi = xi_pct / 100.0;
    const Column col = shear_column(H, Vs, rho, ne);
    // FE fundamental / 2nd frequency (drive + Rayleigh target at the actual FE resonance).
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(col.K), dense(col.M));
    const double f1 = std::sqrt(es.eigenvalues()[0]) / (2 * kPi);
    const double f2 = std::sqrt(es.eigenvalues()[1]) / (2 * kPi);

    // Sub-resonant f = Vs/(8H) -> wH/Vs = pi/4 exactly; |F| = 1/cos = 1.414 (damping negligible here).
    {
        const double fd = Vs / (8 * H);
        const double amp = steady_amplification(col, fd, xi, f1, f2, 80, 60);
        const double kH = 2 * kPi * fd * H / Vs, th = 1.0 / std::cos(kH);
        std::printf("   sub-resonant f=%.3f Hz: amp_FE=%.3f  1/cos(wH/Vs)=%.3f  (err %.1f%%)\n",
                    fd, amp, th, 100 * (amp - th) / th);
        check(std::fabs(amp - th) < 0.05 * th, "sub-resonant amplification = 1/cos(wH/Vs)");
    }
    // Resonance at the FE fundamental: peak amplification = 2/(pi xi) (Kramer Eq 7.30). Rayleigh xi is
    // exactly the target at f1, so the fundamental-mode damping is xi.
    {
        const double amp = steady_amplification(col, f1, xi, f1, f2, 160, 48);
        const double th = 2.0 / (kPi * xi);
        std::printf("   resonant  f=f1=%.3f Hz: amp_FE=%.2f  2/(pi xi)=%.2f  (err %.1f%%)\n",
                    f1, amp, th, 100 * (amp - th) / th);
        check(std::fabs(amp - th) < 0.10 * th, "fundamental resonance amplification = 2/(pi xi)");
    }
}

// (d) 2D consistent mass (tri6/tri15): partition of unity conserves total mass exactly.
void test_mass_conservation() {
    std::printf("-- (d) 2D consistent mass: 1^T M 1 = 2 rho Area (partition of unity) --\n");
    const double W = 4.0, H = 6.0, rho = 2.0;
    const katai::geometry::RectangularDomain dom{0.0, 0.0, W, H, 0};
    for (int order : {6, 15}) {
        const auto mesh = order == 15 ? katai::mesh::generate_structured_tri15(dom, 3, 5)
                                      : katai::mesh::generate_structured_tri6(dom, 3, 5);
        katai::core::DofMap dofs(mesh.node_count, 2);
        dofs.finalize();  // all free
        SparseMatrixBuilder b(dofs.equation_count());
        katai::core::assemble_mass(mesh, dofs, {rho}, b);
        const CsrMatrix M = b.build();
        double total = 0.0; for (double v : M.values) total += v;  // 1^T M 1
        const double exact = 2.0 * rho * W * H;  // x and y blocks each sum to rho * Area
        std::printf("   %-5s: total=%.6f  2 rho A=%.6f  (err %.2e)\n",
                    order == 15 ? "tri15" : "tri6", total, exact, std::fabs(total - exact) / exact);
        check(std::fabs(total - exact) < 1e-9 * exact, "consistent mass conserves total mass");
    }
}

// (e) 2D shear column == 1D site response: uy = 0 everywhere (pure SH shear), rigid base ux = 0, free
// sides + top, horizontal base excitation. The exact solution is the 1D shear field ux = ux(y), so the
// same Newmark integrator on the assembled 2D M/K/C must reproduce f_n = (2n-1) Vs/(4H) and the
// transfer function -- proving the 2D consistent mass + the dimension-agnostic integrator together.
void test_column_2d() {
    std::printf("-- (e) 2D shear column (SH free-field) == 1D site response --\n");
    const double W = 4.0, H = 30.0, E = 208000.0, nu = 0.3, rho = 2.0;
    const double G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);  // Vs = 200
    const katai::geometry::RectangularDomain dom{0.0, 0.0, W, H, 0};
    const auto mesh = katai::mesh::generate_structured_tri6(dom, 2, 20);
    katai::core::DofMap dofs(mesh.node_count, 2);
    for (int n = 0; n < mesh.node_count; ++n) dofs.fix_node_component(n, 1);  // uy = 0 (SH)
    for (int n : mesh.bottom_nodes) dofs.fix_node_component(n, 0);            // rigid base ux = 0
    dofs.finalize();
    const int neq = dofs.equation_count();
    SparseMatrixBuilder bK(neq), bM(neq);
    katai::core::assemble_stiffness(mesh, dofs, {{E, nu}}, bK);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    const CsrMatrix K = bK.build(), M = bM.build();

    // Natural frequencies K phi = w^2 M phi (lowest modes = the 1D shear modes; x-varying modes are
    // far stiffer). Validates the 2D consistent mass against theory.
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(dense(K), dense(M));
    const double f1 = std::sqrt(es.eigenvalues()[0]) / (2 * kPi);
    const double f2 = std::sqrt(es.eigenvalues()[1]) / (2 * kPi);
    for (int n = 1; n <= 2; ++n) {
        const double f_fe = std::sqrt(es.eigenvalues()[n - 1]) / (2 * kPi);
        const double f_th = (2 * n - 1) * Vs / (4 * H);
        std::printf("   f_%d: FE=%.4f Hz  theory=%.4f Hz  (err %.2f%%)\n",
                    n, f_fe, f_th, 100 * (f_fe - f_th) / f_th);
        check(std::fabs(f_fe - f_th) < (n == 1 ? 0.02 : 0.04) * f_th, "2D SH column frequency f_n = (2n-1) Vs/(4H)");
    }

    // Surface (top-mid) ux equation + base-harmonic transfer function via the SAME solve_newmark.
    int surf = -1; double best = 1e30;
    for (int n : mesh.top_nodes) { const double d = std::fabs(mesh.x[n] - W / 2); if (d < best) { best = d; surf = n; } }
    const int surf_eq = dofs.equation(dofs.global_dof(surf, 0));
    const double xi = 0.05;
    const CsrMatrix C = [&] { const auto r = rayleigh_from_modes(f1, xi, f2, xi); return combine(M, r.alpha, K, r.beta); }();
    const Eigen::VectorXd Mr = M * Eigen::VectorXd::Ones(neq);  // M r, r = 1 (all free dofs horizontal)
    auto amp_at = [&](double fd, int periods, int steps_per) {
        const double w = 2 * kPi * fd, A = 1.0, dt = (1.0 / fd) / steps_per;
        auto force = [&](int step) { Eigen::VectorXd f = -(A * std::sin(w * step * dt)) * Mr; return f; };
        const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
        const auto R = solve_newmark(M, C, K, force, dt, periods * steps_per, z, z);
        double a = 0.0;
        for (size_t i = R.u.size() * 3 / 4; i < R.u.size(); ++i)
            a = std::max(a, std::fabs(R.a[i][surf_eq] + A * std::sin(w * R.t[i])));
        return a / A;
    };
    {
        const double fd = Vs / (8 * H), amp = amp_at(fd, 80, 60), kH = 2 * kPi * fd * H / Vs, th = 1.0 / std::cos(kH);
        std::printf("   sub-resonant f=%.3f Hz: amp_FE=%.3f  1/cos(wH/Vs)=%.3f  (err %.1f%%)\n",
                    fd, amp, th, 100 * (amp - th) / th);
        check(std::fabs(amp - th) < 0.06 * th, "2D sub-resonant amplification = 1/cos(wH/Vs)");
    }
    {
        const double amp = amp_at(f1, 160, 48), th = 2.0 / (kPi * xi);
        std::printf("   resonant  f=f1=%.3f Hz: amp_FE=%.2f  2/(pi xi)=%.2f  (err %.1f%%)\n",
                    f1, amp, th, 100 * (amp - th) / th);
        check(std::fabs(amp - th) < 0.12 * th, "2D fundamental resonance = 2/(pi xi) (whole 2D chain)");
    }
}

// (f) Lysmer-Kuhlemeyer absorbing boundary (1D): a d'Alembert one-way downward pulse hits an end
// dashpot c; the reflected energy fraction is exactly R^2 = ((Z-c)/(Z+c))^2, Z = rho Vs A (impedance).
// c = Z -> R = 0 (perfect absorption); c = 0 (free) -> R = 1; c = Z/2 or 2Z -> R = +/-1/3 -> R^2 = 1/9.
void test_absorbing_1d() {
    std::printf("-- (f) 1D absorbing boundary: reflected energy E/E0 = R^2 = ((Z-c)/(Z+c))^2 --\n");
    const double L = 100.0, Vs = 200.0, rho = 2.0, A = 1.0;
    const int ne = 200;
    const double Lh = L / ne, G = rho * Vs * Vs, ke = G * A / Lh, me = rho * A * Lh / 6.0;
    const double Z = rho * Vs * A, z0 = 0.7 * L, w = 3.0;  // Gaussian pulse, width w
    const int nn = ne + 1;  // free-free rod, node index = dof index
    SparseMatrixBuilder bM(nn), bK(nn);
    const double Ke[2][2] = {{ke, -ke}, {-ke, ke}}, Me[2][2] = {{2 * me, me}, {me, 2 * me}};
    for (int e = 0; e < ne; ++e)
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j) { bK.add_entry(e + i, e + j, Ke[i][j]); bM.add_entry(e + i, e + j, Me[i][j]); }
    const CsrMatrix M = bM.build(), K = bK.build();
    // d'Alembert downward pulse: u0(z) = exp(-((z-z0)/w)^2), v0(z) = Vs u0'(z) -> travels toward z=0.
    Eigen::VectorXd u0(nn), v0(nn);
    for (int i = 0; i < nn; ++i) {
        const double z = i * Lh, g = std::exp(-std::pow((z - z0) / w, 2));
        u0[i] = g; v0[i] = Vs * (-2 * (z - z0) / (w * w)) * g;
    }
    const double dt = (Lh / Vs) * 0.25, tm = 1.15 * L / Vs;
    const int ns = static_cast<int>(std::round(tm / dt));
    auto energy = [&](const Eigen::VectorXd& u, const Eigen::VectorXd& v) {
        return 0.5 * v.dot(M * v) + 0.5 * u.dot(K * u);
    };
    const double E0 = energy(u0, v0);
    auto no_force = [&](int) { return Eigen::VectorXd::Zero(nn); };
    const double cvals[4] = {0.0, Z / 2, Z, 2 * Z};
    const double R2[4] = {1.0, 1.0 / 9, 0.0, 1.0 / 9};
    for (int k = 0; k < 4; ++k) {
        SparseMatrixBuilder bC(nn); bC.add_entry(0, 0, cvals[k]);  // nodal dashpot at z=0
        const CsrMatrix C = bC.build();
        const auto Rn = solve_newmark(M, C, K, no_force, dt, ns, u0, v0);
        const double ratio = energy(Rn.u.back(), Rn.v.back()) / E0;
        std::printf("   c=%.3g Z: E/E0=%.4f  R^2=((Z-c)/(Z+c))^2=%.4f  (diff %.3f)\n",
                    cvals[k] / Z, ratio, R2[k], std::fabs(ratio - R2[k]));
        check(std::fabs(ratio - R2[k]) < 0.04, "reflected energy = R^2, R=(Z-c)/(Z+c)");
    }
}

// (g) 2D absorbing base via assemble_boundary_dashpot: a downward SH pulse must radiate out through
// the dashpot base (energy -> 0), whereas a free base traps it (energy conserved). Validates the 2D
// boundary-integral C_b = integral(N^T D_c N ds) assembly.
void test_absorbing_2d() {
    std::printf("-- (g) 2D SH column absorbing base (assemble_boundary_dashpot) radiates energy --\n");
    const double W = 4.0, H = 40.0, E = 208000.0, nu = 0.3, rho = 2.0;
    const double G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), Vp = std::sqrt((lam + 2 * G) / rho);
    const auto mesh = katai::mesh::generate_structured_tri6({0.0, 0.0, W, H, 0}, 2, 40);
    katai::core::DofMap dofs(mesh.node_count, 2);
    for (int n = 0; n < mesh.node_count; ++n) dofs.fix_node_component(n, 1);  // uy = 0 (SH); base ux FREE
    dofs.finalize();
    const int neq = dofs.equation_count();
    SparseMatrixBuilder bK(neq), bM(neq);
    katai::core::assemble_stiffness(mesh, dofs, {{E, nu}}, bK);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    const CsrMatrix K = bK.build(), M = bM.build();
    // d'Alembert downward SH pulse (uniform in x): ux0(y) = exp(-((y-y0)/w)^2), vx0 = Vs ux0'.
    const double y0 = 0.65 * H, w = 3.0;
    Eigen::VectorXd u0 = Eigen::VectorXd::Zero(neq), v0 = Eigen::VectorXd::Zero(neq);
    for (int n = 0; n < mesh.node_count; ++n) {
        const int eq = dofs.equation(dofs.global_dof(n, 0));
        if (eq < 0) continue;
        const double y = mesh.y[n], g = std::exp(-std::pow((y - y0) / w, 2));
        u0[eq] = g; v0[eq] = Vs * (-2 * (y - y0) / (w * w)) * g;
    }
    auto energy = [&](const Eigen::VectorXd& u, const Eigen::VectorXd& v) {
        return 0.5 * v.dot(M * v) + 0.5 * u.dot(K * u);
    };
    const double E0 = energy(u0, v0);
    const double dt = (H / (2.0 * 40) / Vs) * 0.3, tm = 2.5 * H / Vs;
    const int ns = static_cast<int>(std::round(tm / dt));
    auto no_force = [&](int) { return Eigen::VectorXd::Zero(neq); };
    // Free base (no dashpot): energy conserved.
    { SparseMatrixBuilder bC(neq); for (int i = 0; i < neq; ++i) bC.add_entry(i, i, 0.0);
      const auto Rn = solve_newmark(M, bC.build(), K, no_force, dt, ns, u0, v0);
      const double ratio = energy(Rn.u.back(), Rn.v.back()) / E0;
      std::printf("   free base:      E_final/E0=%.4f (energy conserved, wave trapped)\n", ratio);
      check(ratio > 0.9, "free base traps the wave (energy conserved)"); }
    // Absorbing base (Lysmer dashpot, c_t = rho Vs absorbs the tangential SH motion): energy radiates out.
    { SparseMatrixBuilder bC(neq);
      katai::core::assemble_boundary_dashpot(mesh, dofs, mesh.bottom_nodes, rho * Vp, rho * Vs, bC);
      const auto Rn = solve_newmark(M, bC.build(), K, no_force, dt, ns, u0, v0);
      const double ratio = energy(Rn.u.back(), Rn.v.back()) / E0;
      std::printf("   absorbing base: E_final/E0=%.4f (radiated out through the dashpot)\n", ratio);
      check(ratio < 0.05, "absorbing base radiates the wave out (energy -> 0)"); }
}

// (h) Free-field lateral boundary (D3b): under seismic base input the sides must FOLLOW the 1D free-
// field (site response), not be absorbed to zero, while still absorbing waves scattered from the
// interior. Lysmer free-field = boundary dashpot C_b + a driving force C_b*v_ff (net boundary force
// C_b*(v_ff - v_2D)). Solve A (free sides) is the reference free-field (laterally uniform, ux=ux(y));
// solve B (C=C_b, force += C_b*v_A) must reproduce A; solve C (absorbing sides, NO driving) wrongly
// damps the free-field motion at the sides. B == A, C deviates -> the driving force is essential.
void test_freefield_lateral() {
    std::printf("-- (h) free-field lateral boundary preserves the free-field (vs absorbing sides) --\n");
    const double W = 8.0, H = 30.0, E = 208000.0, nu = 0.3, rho = 2.0;
    const double G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double lam = E * nu / ((1 + nu) * (1 - 2 * nu)), Vp = std::sqrt((lam + 2 * G) / rho);
    const auto mesh = katai::mesh::generate_structured_tri6({0.0, 0.0, W, H, 0}, 4, 20);
    katai::core::DofMap dofs(mesh.node_count, 2);
    for (int n = 0; n < mesh.node_count; ++n) dofs.fix_node_component(n, 1);  // uy = 0 (SH)
    for (int n : mesh.bottom_nodes) dofs.fix_node_component(n, 0);            // rigid base ux = 0
    dofs.finalize();
    const int neq = dofs.equation_count();
    SparseMatrixBuilder bK(neq), bM(neq);
    katai::core::assemble_stiffness(mesh, dofs, {{E, nu}}, bK);
    katai::core::assemble_mass(mesh, dofs, {rho}, bM);
    const CsrMatrix K = bK.build(), M = bM.build();
    // Lateral (left + right) boundary dashpot C_b: a horizontally-propagating ux disturbance leaves at
    // Vp (constrained), so the normal coefficient c_n = rho*Vp absorbs it (uy=0 -> c_t unused here).
    SparseMatrixBuilder bC(neq);
    katai::core::assemble_boundary_dashpot(mesh, dofs, mesh.left_nodes, rho * Vp, rho * Vs, bC);
    katai::core::assemble_boundary_dashpot(mesh, dofs, mesh.right_nodes, rho * Vp, rho * Vs, bC);
    const CsrMatrix Cb = bC.build();
    SparseMatrixBuilder bZ(neq); for (int i = 0; i < neq; ++i) bZ.add_entry(i, i, 0.0);
    const CsrMatrix Zero = bZ.build();
    // Rigid-base seismic input: Ricker pulse a_g(t); relative-formulation force -M r a_g.
    const Eigen::VectorXd Mr = M * Eigen::VectorXd::Ones(neq);
    const double f1 = Vs / (4 * H), fr = 3.0 * f1, t0 = 1.2 / fr;
    auto ag = [&](double t) { const double a = kPi * fr * (t - t0), a2 = a * a; return (1 - 2 * a2) * std::exp(-a2); };
    const double dt = (1.0 / f1) / 120, tend = 3.0 / f1;
    const int ns = static_cast<int>(std::round(tend / dt));
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
    // A: free sides -> the reference (laterally uniform) free-field.
    auto fA = [&](int s) { Eigen::VectorXd f = -ag(s * dt) * Mr; return f; };
    const auto RA = solve_newmark(M, Zero, K, fA, dt, ns, z, z);
    // B: free-field boundary = dashpot C_b + driving force C_b * v_A(step).
    auto fB = [&](int s) { Eigen::VectorXd f = -ag(s * dt) * Mr + Cb * RA.v[s]; return f; };
    const auto RB = solve_newmark(M, Cb, K, fB, dt, ns, z, z);
    // C: absorbing sides only (no driving) -> wrongly damps the free-field.
    auto fC = [&](int s) { Eigen::VectorXd f = -ag(s * dt) * Mr; return f; };
    const auto RC = solve_newmark(M, Cb, K, fC, dt, ns, z, z);
    double scale = 0, dB = 0, dC = 0;
    for (size_t i = 0; i < RA.u.size(); ++i) {
        scale = std::max(scale, RA.u[i].cwiseAbs().maxCoeff());
        dB = std::max(dB, (RB.u[i] - RA.u[i]).cwiseAbs().maxCoeff());
        dC = std::max(dC, (RC.u[i] - RA.u[i]).cwiseAbs().maxCoeff());
    }
    std::printf("   free-field vs reference:  max|u_B - u_A|/scale = %.4f  (should be ~0)\n", dB / scale);
    std::printf("   absorbing-only vs ref:    max|u_C - u_A|/scale = %.4f  (should be large)\n", dC / scale);
    check(dB / scale < 0.02, "free-field lateral boundary reproduces the free-field (driving + dashpot)");
    check(dC / scale > 0.10, "absorbing-only sides wrongly damp the free-field (driving force is needed)");
}

// (i) 1D free-field column (solve_free_field_column): the driver's lateral free-field velocity v_ff.
// A uniform shear column (rigid base) shaken at the fundamental f_1 = Vs/(4H) has steady surface
// (relative) displacement |u_surf| = (4/pi) A/(w1^2 2 xi) (fundamental-mode SDOF, phi=sin(pi y/2H));
// off resonance it is far smaller. Validates the standalone site-response engine that feeds the
// free-field lateral boundary (D3b; docs/references/dynamic-seismic-formulation.md sec.8).
void test_freefield_column() {
    std::printf("-- (i) 1D free-field column (solve_free_field_column) = 1D site response --\n");
    const double H = 30.0, E = 208000.0, nu = 0.3, rho = 2.0;
    const double G = E / (2 * (1 + nu)), Vs = std::sqrt(G / rho);
    const double f1 = Vs / (4 * H), w1 = 2 * kPi * f1, xi = 0.05, A = 1.0;
    const double u_th = (4.0 / kPi) * A / (w1 * w1 * 2 * xi);
    // Uniform column: NE quadratic (3-node) segments, nodes 0..2*NE up the height (base = node 0).
    const int NE = 24, nn = 2 * NE + 1;
    std::vector<double> y(nn);
    for (int i = 0; i < nn; ++i) y[i] = H * i / (nn - 1);
    std::vector<katai::core::ShearColumnSegment> segs;
    for (int e = 0; e < NE; ++e) {
        katai::core::ShearColumnSegment s; s.nn = 3;
        s.node = {2 * e, 2 * e + 1, 2 * e + 2, 0, 0}; s.G = G; s.rho = rho;
        segs.push_back(s);
    }
    const auto ray = rayleigh_from_modes(f1, xi, 3 * f1, xi);
    auto peak_surf = [&](double freq, int periods, int steps_per_period) {
        const double dt = 1.0 / (freq * steps_per_period);
        const int steps = periods * steps_per_period;
        auto ag = [&](double t) { return A * std::sin(2 * kPi * freq * t); };
        const auto R = katai::core::solve_free_field_column(y, segs, ray.alpha, ray.beta, ag, dt, steps);
        double pk = 0.0;
        for (const auto& uu : R.u) pk = std::max(pk, std::fabs(uu[nn - 1]));   // surface = top node
        return pk;
    };
    const double res = peak_surf(f1, 30, 60);          // resonance, ~30 cycles -> near steady state
    const double off = peak_surf(f1 / 3.0, 30, 60);    // well below f_1
    std::printf("   Vs=%.1f  f_1=%.3f Hz  resonant |u_surf|=%.4f m  (theory %.4f, err %.1f%%)\n",
                Vs, f1, res, u_th, 100 * (res - u_th) / u_th);
    std::printf("   off-resonance |u_surf|=%.4f m  (ratio to resonance %.2f)\n", off, off / res);
    check(std::fabs(res - u_th) < 0.12 * u_th, "resonant surface disp = (4/pi) A/(w1^2 2xi)");
    check(off < 0.4 * res, "off-resonance response is much smaller (column amplifies at f_1)");
}

}  // namespace

int main() {
    std::printf("Dynamics core (Newmark + Rayleigh) verification\n\n");
    test_sdof();
    test_column_frequencies();
    test_column_transfer();
    test_mass_conservation();
    test_column_2d();
    test_absorbing_1d();
    test_absorbing_2d();
    test_freefield_lateral();
    test_freefield_column();
    if (g_failures == 0) {
        std::printf("\nOK: dynamics core = SDOF + 1D/2D site response + transfer + absorbing + free-field boundaries\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
