#pragma once
// Free-field 1D vertical shear column — the seismic lateral-boundary driver.
//
// The Lysmer-Kuhlemeyer free-field lateral boundary = the lateral dashpot C_b
// (assemble_boundary_dashpot, D3) + the driving force C_b·v_ff, where v_ff is the 1D site
// response of the boundary soil column (net boundary force C_b·(v_ff − v_2D);
// docs/references/dynamic-seismic-formulation.md §8, D3b). This header produces that v_ff:
// the motion of a vertical shear column RELATIVE to a RIGID base
//     m u'' + c u' + k u = −m·1·a_g(t),   c = α m + β k  (Rayleigh, same as 2D),
// with quadratic (tri6 edge = 3 nodes) / quartic (tri15 edge = 5 nodes) 1D Lagrange
// elements (k = ∫ G N' N' dy shear stiffness, m = ∫ ρ N N dy consistent mass). The same
// solve_newmark integrator. Verification: test_dynamics (i) — uniform column
// f_n=(2n-1)Vs/4H + surface transfer 2/(πξ) at resonance.

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <vector>

#include <Eigen/Dense>

#include <katai/analysis/dynamics.hpp>       // solve_newmark, NewmarkResult, DynamicsSolveFactory
#include <katai/fem/assembly/assembler.hpp>  // expand_to_full
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/math/sparse_matrix.hpp>

namespace katai::core {

// One 1D shear element of the column: nn (3 or 5) nodes node[0..nn) (ordered in y,
// corner→corner), shear modulus G and mass density ρ. Node indices refer to the y[] array
// passed to solve_free_field_column.
struct ShearColumnSegment {
    std::array<int, 5> node{};
    int nn = 3;
    double G = 0.0;
    double rho = 0.0;
};

// 1D Lagrange basis functions N and derivatives dN/dξ for nn nodes (equally spaced on [-1,1]).
inline void ff_lagrange_1d(int nn, double xi, double* N, double* dN) {
    double xk[5];
    for (int k = 0; k < nn; ++k) xk[k] = -1.0 + 2.0 * k / (nn - 1);
    for (int i = 0; i < nn; ++i) {
        double Ni = 1.0, dNi = 0.0;
        for (int m = 0; m < nn; ++m) {
            if (m == i) continue;
            // N_i = Π_{j≠i} (ξ−x_j)/(x_i−x_j);  dN_i/dξ = Σ_{m≠i} 1/(x_i−x_m) Π_{j≠i,m} (ξ−x_j)/(x_i−x_j)
            double term = 1.0 / (xk[i] - xk[m]);
            for (int j = 0; j < nn; ++j)
                if (j != i && j != m) term *= (xi - xk[j]) / (xk[i] - xk[j]);
            dNi += term;
            Ni *= (xi - xk[m]) / (xk[i] - xk[m]);
        }
        N[i] = Ni;
        dN[i] = dNi;
    }
}

struct FreeFieldColumnResult {
    std::vector<double> t;                 // time of each step (nsteps+1)
    std::vector<Eigen::VectorXd> u, v;     // full-nodal (size = y.size()); relative to the rigid base
    int base_node = -1;                    // the fixed (deepest) node
};

// 1D free-field site response. y = column node depths (any order; smallest y = base,
// fixed). segs = shear elements (with y[] indices). Relative motion under the horizontal
// base acceleration a_g(t): m u'' + c u' + k u = −m·1·a_g, c = α m + β k. Returns:
// full-nodal u/v at each step (base ~0).
inline FreeFieldColumnResult solve_free_field_column(
    const std::vector<double>& y, const std::vector<ShearColumnSegment>& segs,
    double alpha, double beta, const std::function<double(double)>& ag,
    double dt, int nsteps, const DynamicsSolveFactory& factory = {}) {
    FreeFieldColumnResult out;
    const int nc = static_cast<int>(y.size());
    if (nc < 2 || segs.empty()) return out;

    // Rigid base: fix the deepest (smallest y) node.
    int base = 0;
    for (int i = 1; i < nc; ++i) if (y[i] < y[base]) base = i;
    DofMap dofs(nc, 1);
    dofs.fix_node_component(base, 0);
    dofs.finalize();
    const int neq = dofs.equation_count();
    if (neq == 0) { out.base_node = base; return out; }

    // 5-point Gauss-Legendre (exact up to ~degree 8 for the quartic mass).
    static const double gp[5] = {-0.906179845938664, -0.538469310105683, 0.0,
                                 0.538469310105683, 0.906179845938664};
    static const double gw[5] = {0.236926885056189, 0.478628670499366, 0.568888888888889,
                                 0.478628670499366, 0.236926885056189};
    math::SparseMatrixBuilder bK(neq), bM(neq);
    for (const auto& s : segs) {
        const int nn = std::clamp(s.nn, 2, 5);
        double ke[5][5] = {}, me[5][5] = {};
        for (int q = 0; q < 5; ++q) {
            double N[5], dN[5];
            ff_lagrange_1d(nn, gp[q], N, dN);
            double J = 0.0;
            for (int i = 0; i < nn; ++i) J += dN[i] * y[s.node[i]];   // dy/dξ
            const double detJ = std::fabs(J);
            if (detJ < 1e-30) continue;
            const double wq = gw[q] * detJ;
            for (int i = 0; i < nn; ++i)
                for (int j = 0; j < nn; ++j) {
                    ke[i][j] += s.G * (dN[i] / J) * (dN[j] / J) * wq;  // ∫ G (dN/dy)² dy
                    me[i][j] += s.rho * N[i] * N[j] * wq;              // ∫ ρ N N dy
                }
        }
        for (int i = 0; i < nn; ++i) {
            const int ei = dofs.equation(dofs.global_dof(s.node[i], 0));
            if (ei < 0) continue;
            for (int j = 0; j < nn; ++j) {
                const int ej = dofs.equation(dofs.global_dof(s.node[j], 0));
                if (ej < 0) continue;
                bK.add_entry(ei, ej, ke[i][j]);
                bM.add_entry(ei, ej, me[i][j]);
            }
        }
    }
    const math::CsrMatrix K1 = bK.build(), M1 = bM.build();
    // Rayleigh C = α M + β K (same coefficients as the 2D solve).
    math::SparseMatrixBuilder bC(neq);
    auto add_scaled = [&](const math::CsrMatrix& A, double sc) {
        for (int r = 0; r < A.rows; ++r)
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
                bC.add_entry(r, A.col_indices[p], sc * A.values[p]);
    };
    add_scaled(M1, alpha);
    add_scaled(K1, beta);
    const math::CsrMatrix C1 = bC.build();

    // Horizontal base excitation: F = −a_g(t)·M·1 (influence r = 1, on the free DOFs).
    const Eigen::VectorXd Mr = M1 * Eigen::VectorXd::Ones(neq);
    auto force = [&](int step) -> Eigen::VectorXd { return -ag(step * dt) * Mr; };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
    const NewmarkResult R = solve_newmark(M1, C1, K1, force, dt, nsteps, z, z, 0.5, 0.25, factory);

    out.base_node = base;
    out.t = R.t;
    for (size_t s = 0; s < R.u.size(); ++s) {
        out.u.push_back(expand_to_full(dofs, R.u[s]));
        out.v.push_back(expand_to_full(dofs, R.v[s]));
    }
    return out;
}

// 1D free-field column with a COMPLIANT (absorbing) base — the lateral free-field driver of
// a 2D model with a compliant base (while the 2D base absorbs, the lateral columns must not
// impose the rigid-base free field; dynamic-seismic-formulation.md §11.2 item 3).
// Differences from the rigid version:
//   - The base node is NOT fixed; the base Lysmer dashpot c = ρV_s (unit-area column, point
//     term) is added to C and the input is applied through the SAME dashpot as the point
//     force 2·ρV_s·v_up(t) (Joyner-Chen factor 2; v_up = the ½·∫ integral of the within
//     a_g — SAME convention as 2D, both factors kept explicit).
//   - The solution is in TOTAL motion (NO −m·a_g inertial driver); u/v return TOTALS.
//     Starting from rest v_up(0)=0 → F(0)=0 → a(0)=0 in closed form (M is never factored).
// rhoVs_base = the ρV_s of the base half-space (SAME material as the 2D base dashpot = the
// deepest layer).
inline FreeFieldColumnResult solve_free_field_column_compliant(
    const std::vector<double>& y, const std::vector<ShearColumnSegment>& segs,
    double alpha, double beta, const std::function<double(double)>& ag,
    double dt, int nsteps, double rhoVs_base, const DynamicsSolveFactory& factory = {}) {
    FreeFieldColumnResult out;
    const int nc = static_cast<int>(y.size());
    if (nc < 2 || segs.empty() || rhoVs_base <= 0.0) return out;

    int base = 0;
    for (int i = 1; i < nc; ++i) if (y[i] < y[base]) base = i;
    DofMap dofs(nc, 1);       // no node is fixed (the base moves)
    dofs.finalize();
    const int neq = dofs.equation_count();

    static const double gp[5] = {-0.906179845938664, -0.538469310105683, 0.0,
                                 0.538469310105683, 0.906179845938664};
    static const double gw[5] = {0.236926885056189, 0.478628670499366, 0.568888888888889,
                                 0.478628670499366, 0.236926885056189};
    math::SparseMatrixBuilder bK(neq), bM(neq);
    for (const auto& s : segs) {
        const int nn = std::clamp(s.nn, 2, 5);
        double ke[5][5] = {}, me[5][5] = {};
        for (int q = 0; q < 5; ++q) {
            double N[5], dN[5];
            ff_lagrange_1d(nn, gp[q], N, dN);
            double J = 0.0;
            for (int i = 0; i < nn; ++i) J += dN[i] * y[s.node[i]];
            const double detJ = std::fabs(J);
            if (detJ < 1e-30) continue;
            const double wq = gw[q] * detJ;
            for (int i = 0; i < nn; ++i)
                for (int j = 0; j < nn; ++j) {
                    ke[i][j] += s.G * (dN[i] / J) * (dN[j] / J) * wq;
                    me[i][j] += s.rho * N[i] * N[j] * wq;
                }
        }
        for (int i = 0; i < nn; ++i) {
            const int ei = dofs.equation(dofs.global_dof(s.node[i], 0));
            for (int j = 0; j < nn; ++j) {
                const int ej = dofs.equation(dofs.global_dof(s.node[j], 0));
                bK.add_entry(ei, ej, ke[i][j]);
                bM.add_entry(ei, ej, me[i][j]);
            }
        }
    }
    const math::CsrMatrix K1 = bK.build(), M1 = bM.build();
    math::SparseMatrixBuilder bC(neq);
    auto add_scaled = [&](const math::CsrMatrix& A, double sc) {
        for (int r = 0; r < A.rows; ++r)
            for (int p = A.row_ptr[r]; p < A.row_ptr[r + 1]; ++p)
                bC.add_entry(r, A.col_indices[p], sc * A.values[p]);
    };
    add_scaled(M1, alpha);
    add_scaled(K1, beta);
    const int beq = dofs.equation(dofs.global_dof(base, 0));
    bC.add_entry(beq, beq, rhoVs_base);            // base absorber (unit-area point dashpot)
    const math::CsrMatrix C1 = bC.build();

    // v_up = ½·∫a_g (trapezoid; within→upward-travelling wave) — same as the 2D compliant driver.
    std::vector<double> vg(nsteps + 1, 0.0);
    for (int k = 1; k <= nsteps; ++k)
        vg[k] = vg[k - 1] + 0.5 * dt * (ag((k - 1) * dt) + ag(k * dt));
    auto force = [&](int step) -> Eigen::VectorXd {
        Eigen::VectorXd f = Eigen::VectorXd::Zero(neq);
        const double v_up = 0.5 * vg[std::min(step, nsteps)];
        f(beq) = rhoVs_base * (2.0 * v_up);        // Joyner-Chen factor 2 explicit
        return f;
    };
    const Eigen::VectorXd z = Eigen::VectorXd::Zero(neq);
    const NewmarkResult R = solve_newmark(M1, C1, K1, force, dt, nsteps, z, z, 0.5, 0.25, factory,
                                          {}, z /* a(0)=0: at rest + F(0)=0 */);

    out.base_node = base;
    out.t = R.t;
    for (size_t s = 0; s < R.u.size(); ++s) {
        out.u.push_back(expand_to_full(dofs, R.u[s]));
        out.v.push_back(expand_to_full(dofs, R.v[s]));
    }
    return out;
}

}  // namespace katai::core
