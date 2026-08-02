// Structural internal-force OUTPUTS (structural_forces.hpp) verification — the PLAXIS
// Output counterpart (N/Q/M diagram + anchor/geogrid axial force). All closed-form, no MKL
// needed (pure post-processing + a small dense solve).
//
//  (1) Cantilever plate beam, transverse tip load P → M(s)=P(L−s) LINEAR (PL at the root,
//      0 at the tip), |Q|=P constant, N=0 — the diagram recovers the FULL distribution
//      (not just the max). Timoshenko φ is quadratic = exact.
//  (2) Cantilever with axial tip N → N(s)=N constant, M=0.
//  (3) Anchor axial-force recovery mirrors the solver exactly: N=kk(U−Up), Fmax clamp.
//  (4) Geogrid diagram: uniform ε → N=EA·ε; compression → slack (N=0); N_p yield → N=N_p.
// (See docs/references/structural-plate-formulation.md §9.)
#include <katai/analysis/structural_forces.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include <Eigen/Dense>

#include <katai/fem/assembly/dof_map.hpp>
#include <katai/mesh/mesh.hpp>

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}

using namespace katai::core;
using katai::mesh::Mesh;

// Horizontal cantilever plate beam: ne 3-node elements along the x axis, node 0 clamped,
// (fx,fy,m) load at the tip. Global DOF layout: DofMap(nnodes,2) translations + an extra
// rotational DOF per node. Returns {mesh, dofs, plates, disp} (disp full global-DOF). For
// the diagram verification.
struct Cantilever {
    Mesh mesh;
    DofMap dofs;
    std::vector<PlateElement> plates;
    Eigen::VectorXd disp;
    Cantilever() : dofs(0, 2) {}
};

Cantilever build_cantilever(double L, int ne, const plate::PlateProps& p,
                            double fx, double fy, double m) {
    const int nnodes = 2 * ne + 1;
    const double h = L / (2 * ne);
    Cantilever c;
    c.mesh.node_count = nnodes;
    c.mesh.x.resize(nnodes);
    c.mesh.y.resize(nnodes);
    for (int i = 0; i < nnodes; ++i) { c.mesh.x[i] = i * h; c.mesh.y[i] = 0.0; }

    DofMap dofs(nnodes, 2);
    std::vector<int> rot(nnodes);
    for (int i = 0; i < nnodes; ++i) rot[i] = dofs.add_extra_dof();
    // Clamped: node 0 ux, uy, φ.
    dofs.fix_node_component(0, 0);
    dofs.fix_node_component(0, 1);
    dofs.fix(rot[0]);
    dofs.finalize();

    for (int e = 0; e < ne; ++e) {
        PlateElement pe;
        pe.nodes = {2 * e, 2 * e + 2, 2 * e + 1};  // A(ξ−1), B(ξ+1), mid(ξ0)
        pe.rot_dof = {rot[2 * e], rot[2 * e + 2], rot[2 * e + 1]};
        pe.props = p;  // trans_dof default −1 → mesh-node sharing (global_dof)
        c.plates.push_back(pe);
    }

    // Dense global K (total_dofs), solve the free subsystem.
    const int nt = dofs.total_dofs();
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nt, nt);
    for (const auto& pe : c.plates) {
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = c.mesh.x[pe.nodes[k]]; X(k, 1) = c.mesh.y[pe.nodes[k]]; }
        const auto Ke = plate::stiffness(X, p);
        int g[9];
        for (int k = 0; k < 3; ++k) {
            g[3 * k + 0] = dofs.global_dof(pe.nodes[k], 0);
            g[3 * k + 1] = dofs.global_dof(pe.nodes[k], 1);
            g[3 * k + 2] = pe.rot_dof[k];
        }
        for (int a = 0; a < 9; ++a)
            for (int b = 0; b < 9; ++b) K(g[a], g[b]) += Ke(a, b);
    }

    Eigen::VectorXd f = Eigen::VectorXd::Zero(nt);
    const int tip = nnodes - 1;
    f(dofs.global_dof(tip, 0)) = fx;
    f(dofs.global_dof(tip, 1)) = fy;
    f(rot[tip]) = m;

    std::vector<int> free;
    for (int d = 0; d < nt; ++d) if (!dofs.is_fixed(d)) free.push_back(d);
    const int nf = static_cast<int>(free.size());
    Eigen::MatrixXd Kff(nf, nf);
    Eigen::VectorXd ff(nf);
    for (int i = 0; i < nf; ++i) {
        ff(i) = f(free[i]);
        for (int j = 0; j < nf; ++j) Kff(i, j) = K(free[i], free[j]);
    }
    const Eigen::VectorXd uf = Kff.ldlt().solve(ff);
    c.disp = Eigen::VectorXd::Zero(nt);
    for (int i = 0; i < nf; ++i) c.disp(free[i]) = uf(i);
    c.dofs = std::move(dofs);
    return c;
}

void test_plate_bending_diagram() {
    // Stiff wall: E=3e7 kPa, d=0.5 m → EA=E d, EI=E d^3/12, nu=0.15.
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    plate::PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0, P = 100.0;

    Cantilever c = build_cantilever(L, 8, p, 0.0, -P, 0.0);  // transverse tip load −y
    const auto diag = plate_force_diagram(c.plates, c.mesh, c.dofs, c.disp, {}, 5);
    check(!diag.empty(), "diagram not empty");

    // M(s) = P(L−s) LINEAR (PL at the root, ~0 at the tip); |Q|=P constant; N≈0. s = arc length from the root.
    double max_M_err = 0.0, max_Q_err = 0.0, max_N = 0.0;
    for (const auto& st : diag) {
        const double M_exact = P * (L - st.s);
        max_M_err = std::max(max_M_err, std::fabs(std::fabs(st.M) - M_exact));
        max_Q_err = std::max(max_Q_err, std::fabs(std::fabs(st.Q) - P));
        max_N = std::max(max_N, std::fabs(st.N));
    }
    const auto env = force_envelope(diag);
    std::printf("  bending: max|M|=%.4f (PL=%.1f)  max_M_err=%.3e  max_Q_err=%.3e  max|N|=%.2e\n",
                env.max_abs_M, P * L, max_M_err, max_Q_err, max_N);
    check(close(env.max_abs_M, P * L, 1e-3), "root bending moment = PL");
    check(max_M_err < 1e-3 * P * L, "M(s) = P(L−s) LINEAR at all stations (distribution)");
    check(max_Q_err < 1e-3 * P, "Q(s) = P constant (distribution)");
    check(max_N < 1e-6 * P, "N ≈ 0 (pure bending)");
    // First station at the root (s≈0), last at the tip (s≈L): M≈0 at the tip.
    check(std::fabs(diag.back().M) < 1e-3 * P * L, "M ≈ 0 at the tip");
    check(close(diag.back().s, L, 1e-6), "last station arc length = L");
}

void test_plate_axial_diagram() {
    const double E = 3.0e7, d = 0.5, nu = 0.15;
    plate::PlateProps p; p.EA = E * d; p.EI = E * d * d * d / 12.0; p.nu = nu;
    const double L = 10.0, N0 = 500.0;

    Cantilever c = build_cantilever(L, 8, p, N0, 0.0, 0.0);  // axial tip load +x
    const auto diag = plate_force_diagram(c.plates, c.mesh, c.dofs, c.disp, {}, 5);
    double max_N_err = 0.0, max_M = 0.0;
    for (const auto& st : diag) {
        max_N_err = std::max(max_N_err, std::fabs(st.N - N0));
        max_M = std::max(max_M, std::fabs(st.M));
    }
    std::printf("  axial: N=%.4f (N0=%.1f) max_N_err=%.3e max|M|=%.2e\n",
                diag.front().N, N0, max_N_err, max_M);
    check(max_N_err < 1e-4 * N0, "N(s) = N0 constant (axial diagram)");
    check(max_M < 1e-6 * N0, "M ≈ 0 (pure axial)");
}

void test_anchor_force_recovery() {
    // fixed-end anchor: node_a free, far end a fixed point. EA/L · U, then Fmax clamp.
    Mesh mesh;
    mesh.node_count = 1;
    mesh.x = {0.0}; mesh.y = {0.0};
    DofMap dofs(1, 2);
    dofs.finalize();

    AnchorElement an;
    an.node_a = 0; an.node_b = -1;
    an.fixed_point = {3.0, 0.0};   // horizontal along +x
    an.EA = 2.0e5; an.L = 5.0;
    const double kk = an.EA / an.L;

    // A disp translating node_a by Ua along +x (the anchor elongates → tension). U along
    // dir·(Xb−Xa): dir=+x; g=[−dir, dir] → U = −Ua (node_a) ... solver sign: g_a=−dir →
    // U=−dir·u_a. If node_a moves +x the anchor SHORTENS → compression. For the tension
    // test move node_a to −x (so it elongates).
    const double Ua = -0.001;  // node_a to −x → distance grows → tension
    Eigen::VectorXd disp = Eigen::VectorXd::Zero(dofs.total_dofs());
    disp(dofs.global_dof(0, 0)) = Ua;
    // U = g_a·u = (−dir_x)·Ua = −1·(−0.001) = +0.001 ; N = kk·U > 0 (tension).
    const double U_expected = -1.0 * Ua;
    const auto r = anchor_force(an, mesh, dofs, disp);
    std::printf("  anchor elastic: N=%.4f (kk·U=%.4f) yielded=%d\n", r.N, kk * U_expected, (int)r.yielded);
    check(close(r.N, kk * U_expected, 1e-9), "anchor N = kk·(U) (elastic recovery)");
    check(!r.yielded, "elastic: not yielded");

    // Fmax_tens yield: clamp.
    an.Fmax_tens = 0.5 * kk * U_expected;  // smaller than expected → yields
    const auto r2 = anchor_force(an, mesh, dofs, disp);
    check(close(r2.N, an.Fmax_tens, 1e-12), "anchor N = Fmax_tens (yield clamp)");
    check(r2.yielded, "Fmax exceeded: yield flag");

    // committed U_p feedback: N = kk·(U − Up).
    an.Fmax_tens = -1.0;  // unbounded again
    const double Up = 0.0003;
    const auto r3 = anchor_force(an, mesh, dofs, disp, Up);
    check(close(r3.N, kk * (U_expected - Up), 1e-9), "anchor N = kk·(U−Up) (committed plastic)");
}

void test_geogrid_force_diagram() {
    // Horizontal geogrid (x axis), length Lg; end B moves +x by Δ → uniform ε=Δ/Lg → N=EA·ε.
    const double Lg = 4.0, EA = 1.0e5, Delta = 0.002;
    Mesh mesh;
    mesh.node_count = 3;
    mesh.x = {0.0, Lg, 0.5 * Lg};   // A, B, mid
    mesh.y = {0.0, 0.0, 0.0};
    DofMap dofs(3, 2);
    dofs.finalize();

    GeogridElement ge;
    ge.nodes = {0, 1, 2};
    ge.props.EA = EA; ge.props.Np = -1.0;  // pure tension-only

    Eigen::VectorXd disp = Eigen::VectorXd::Zero(dofs.total_dofs());
    disp(dofs.global_dof(1, 0)) = Delta;          // end B +x
    disp(dofs.global_dof(2, 0)) = 0.5 * Delta;    // mid at half (linear ux)
    const double eps = Delta / Lg, N_exact = EA * eps;

    const auto diag = geogrid_force_diagram(ge, mesh, dofs, disp);
    check(diag.size() == 2, "geogrid diagram has 2 Gauss stations");
    double max_err = 0.0;
    for (const auto& st : diag) max_err = std::max(max_err, std::fabs(st.N - N_exact));
    std::printf("  geogrid tension: N=%.4f (EA·ε=%.4f) max_err=%.3e\n",
                diag.front().N, N_exact, max_err);
    check(max_err < 1e-6 * N_exact, "geogrid N = EA·ε (tension diagram)");

    // Compression → slack (N=0): move end B to −x.
    Eigen::VectorXd dispc = Eigen::VectorXd::Zero(dofs.total_dofs());
    dispc(dofs.global_dof(1, 0)) = -Delta;
    dispc(dofs.global_dof(2, 0)) = -0.5 * Delta;
    const auto diagc = geogrid_force_diagram(ge, mesh, dofs, dispc);
    check(std::fabs(diagc.front().N) < 1e-12, "geogrid slack in compression (N=0)");

    // N_p yield: under the cap.
    ge.props.Np = 0.5 * N_exact;
    const auto diagp = geogrid_force_diagram(ge, mesh, dofs, disp);
    for (const auto& st : diagp) check(close(st.N, ge.props.Np, 1e-12), "geogrid N = N_p (yield)");
}

}  // namespace

int main() {
    test_plate_bending_diagram();
    test_plate_axial_diagram();
    test_anchor_force_recovery();
    test_geogrid_force_diagram();
    if (g_failures == 0) {
        std::printf("OK: structural force outputs verified "
                    "(plate M/Q/N diagram, anchor force, geogrid diagram)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
