// Plate M-N plastic hinge (the Mp/Np diamond, PLAXIS MMM sec 18.3;
// structural-plate-formulation.md sec 10) -- exact-class V&V:
//   (1) Return-map sweep: f <= 0, idempotence, and an INDEPENDENT oracle = energy-norm
//       nearest-point comparison against a coarse sampling of the diamond boundary
//       (closed form == CPP).
//   (2) FD tangent: inside a region the map is PIECEWISE LINEAR -> the consistent tangent
//       matches FD to machine precision.
//   (3) Element elastic-limit identity: caps unbounded -> quadrature f_int == K*u (round-off).
//   (4) BVP yield ONSET + saturation (the test_plate_soil harness): simple beam on soft
//       soil + central load. Closed form P_lim = 2*Mp/s_g, s_g = the critical GAUSS
//       point's distance from the support (PLAXIS: the check is AT THE STRESS POINT, not
//       the node). At 0.9*P_lim the plastic state is BIT-ZERO, at 1.1*P_lim it yields
//       (the onset bracket = verification of the closed form); at 1.5*P_lim the Gauss |M|
//       saturates EXACTLY at Mp + large softening. NOTE: on an elastic soil bed a single
//       hinge does NOT create a mechanism (excess load flows into the soil, the solve
//       converges) -- the Prandtl collapse-lf pattern does NOT apply here; measured, and
//       the test is pinned to the correct physics.
//   (5) M-N interaction BVP: axial pre-tension T = Np/2 -> the critical Gauss point sits
//       ON THE DIAMOND SURFACE (|N|/Np + |M|/Mp = 1) and the moment plateau
//       ~ Mp*(1 - T/Np); the Np-unbounded control reaches full Mp -> the drop truly comes
//       from the interaction.
//   (6) Moment plateau + diagram: the capped Gauss M expand to the stations by Lagrange.
//   (7) tri15 (5-node plate): the same onset + saturation pins.
#include <katai/analysis/nonlinear_solver.hpp>
#include <katai/analysis/structural_forces.hpp>
#include <katai/fem/assembly/dof_map.hpp>
#include <katai/fem/elements/plate.hpp>
#include <katai/materials/material_model.hpp>
#include <katai/geometry/rectangular_domain.hpp>
#include <katai/linsolve/direct_solver.hpp>
#include <katai/mesh/mesh.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using katai::core::DofMap;
using katai::core::MaterialModel;
using katai::core::MaterialType;
using katai::core::PlateElement;
using katai::core::PlateElement5;
using katai::geometry::RectangularDomain;
namespace linsolve = katai::linsolve;
using katai::mesh::Mesh;
namespace plate = katai::core::plate;

namespace {
int g_failures = 0;
void check(bool ok, const char* what) {
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what); ++g_failures; }
}
bool close(double a, double b, double tol) {
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(b));
}
Eigen::VectorXd solve_spd(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealSymmetricPositiveDefinite);
    s->factorize(k);
    return s->solve(r);
}
// The tangent can be singular/asymmetric in plasticity -> general (nonsymmetric) solve.
Eigen::VectorXd solve_gen(const katai::math::CsrMatrix& k, const Eigen::VectorXd& r) {
    auto s = linsolve::make_direct_solver(linsolve::MatrixType::RealNonsymmetric);
    s->factorize(k);
    return s->solve(r);
}

// ------------------------------------------------------------------------------------------
// (1) Return-map sweep + the independent CPP oracle.
// Oracle: coarsely sample the diamond BOUNDARY (4 edges x 4000 points), find the nearest
// boundary point in the energy norm E = ((N_tr-N)^2/EA + (M_tr-M)^2/EI)/2; the closed-form
// return cannot be farther than that (equality ~ sampling resolution). At an interior
// point the map is the identity.
// ------------------------------------------------------------------------------------------
void test_return_map_scan() {
    std::printf("-- (1) M-N return-map sweep + brute-force CPP oracle --\n");
    plate::PlateProps p;
    p.EA = 2.0e6; p.EI = 2.0e4; p.nu = 0.15; p.Mp = 150.0; p.Np = 3000.0;

    // Diamond boundary sampler: 4 edges corner to corner.
    const int nb = 4000;
    std::vector<double> bN(nb), bM(nb);
    for (int i = 0; i < nb; ++i) {
        const double t = 4.0 * i / nb;             // [0,4): 1 unit per edge
        const int edge = (int)t;
        const double s = t - edge;                 // in-edge parameter
        double N0, M0, N1, M1;
        switch (edge) {                            // (+Np,0)->(0,+Mp)->(-Np,0)->(0,-Mp)->start
            case 0:  N0 = p.Np; M0 = 0.0;  N1 = 0.0;   M1 = p.Mp;  break;
            case 1:  N0 = 0.0;  M0 = p.Mp; N1 = -p.Np; M1 = 0.0;   break;
            case 2:  N0 = -p.Np; M0 = 0.0; N1 = 0.0;   M1 = -p.Mp; break;
            default: N0 = 0.0;  M0 = -p.Mp; N1 = p.Np; M1 = 0.0;   break;
        }
        bN[i] = N0 + s * (N1 - N0);
        bM[i] = M0 + s * (M1 - M0);
    }
    auto energy = [&](double Ntr, double Mtr, double N, double M) {
        return 0.5 * ((Ntr - N) * (Ntr - N) / p.EA + (Mtr - M) * (Mtr - M) / p.EI);
    };

    bool all_f = true, all_idem = true, all_cpp = true, all_int = true;
    double worst_f = 0.0, worst_cpp = 0.0;
    const int ng = 81;
    for (int i = 0; i < ng; ++i)
        for (int j = 0; j < ng; ++j) {
            const double Ntr = -3.0 * p.Np + 6.0 * p.Np * i / (ng - 1);
            const double Mtr = -3.0 * p.Mp + 6.0 * p.Mp * j / (ng - 1);
            const double eps = Ntr / p.EA, kap = Mtr / p.EI;   // ep_c = kp_c = 0
            const auto r = plate::mn_return(p, eps, kap, 0.0, 0.0);
            const double f = std::fabs(r.N) / p.Np + std::fabs(r.M) / p.Mp - 1.0;
            if (f > 1e-12) { all_f = false; worst_f = std::max(worst_f, f); }
            const double ftr = std::fabs(Ntr) / p.Np + std::fabs(Mtr) / p.Mp - 1.0;
            if (ftr <= 0.0) {
                // Interior point: identity + the state unchanged (bit-for-bit).
                if (r.N != Ntr || r.M != Mtr || r.ep != 0.0 || r.kp != 0.0 || r.yielded)
                    all_int = false;
            } else {
                // Idempotence: same total strain with the returned plastic state -> same stress.
                const auto r2 = plate::mn_return(p, eps, kap, r.ep, r.kp);
                if (!close(r2.N, r.N, 1e-9) || !close(r2.M, r.M, 1e-9)) all_idem = false;
                // The independent CPP oracle.
                double ebest = 1e300;
                for (int b = 0; b < nb; ++b) ebest = std::min(ebest, energy(Ntr, Mtr, bN[b], bM[b]));
                const double emap = energy(Ntr, Mtr, r.N, r.M);
                // Sampling resolution edge length / nb -> energy gap O(h^2); leave slack.
                const double slack = 1e-4 * (1.0 + ebest);
                if (emap > ebest + slack) { all_cpp = false; worst_cpp = std::max(worst_cpp, emap - ebest); }
            }
        }
    std::printf("   %d x %d sweep: worst f = %.2e, worst CPP excess = %.2e\n",
                ng, ng, worst_f, worst_cpp);
    check(all_f, "return admissible everywhere (f <= 1e-12)");
    check(all_int, "interior-point identity bit-for-bit (state never written)");
    check(all_idem, "idempotence: re-evaluation with the returned state gives the same stress");
    check(all_cpp, "closed form == brute-force nearest point (energy norm)");

    // Single-cap cases: N unbounded -> N unchanged (bit), |M| <= Mp; M unbounded -> symmetric.
    plate::PlateProps pM = p; pM.Np = -1.0;
    plate::PlateProps pN = p; pN.Mp = -1.0;
    bool mcap_ok = true, ncap_ok = true;
    for (int i = 0; i < 41; ++i) {
        // Expected trials from the values the map actually SEES (the eps/kap round trip can
        // shift 1 ulp; the bit-identity claim is made against the map's own trial).
        const double eps = (-2.0 * p.Np + 4.0 * p.Np * i / 40.0) / p.EA;
        const double kap = (-2.0 * p.Mp + 4.0 * p.Mp * i / 40.0) / p.EI;
        const double Ntr = p.EA * eps, Mtr = p.EI * kap;
        const auto rm = plate::mn_return(pM, eps, kap, 0.0, 0.0);
        if (rm.N != Ntr || std::fabs(rm.M) > p.Mp * (1.0 + 1e-12)) mcap_ok = false;
        if (std::fabs(Mtr) > p.Mp && !close(std::fabs(rm.M), p.Mp, 1e-12)) mcap_ok = false;
        if (std::fabs(Mtr) > p.Mp && rm.ep != 0.0) mcap_ok = false;   // the axial state never drifts
        const auto rn = plate::mn_return(pN, eps, kap, 0.0, 0.0);
        if (rn.M != Mtr || std::fabs(rn.N) > p.Np * (1.0 + 1e-12)) ncap_ok = false;
        if (std::fabs(Ntr) > p.Np && !close(std::fabs(rn.N), p.Np, 1e-12)) ncap_ok = false;
        if (std::fabs(Ntr) > p.Np && rn.kp != 0.0) ncap_ok = false;
    }
    check(mcap_ok, "Mp-only: N bit-for-bit elastic, M exactly clamped, eps_p never drifts");
    check(ncap_ok, "Np-only: M bit-for-bit elastic, N exactly clamped, kap_p never drifts");

    // Caps unbounded -> pure elastic identity (bit-for-bit).
    plate::PlateProps pe = p; pe.Mp = -1.0; pe.Np = -1.0;
    const auto re = plate::mn_return(pe, 1.0e-3, 2.0e-3, 0.0, 0.0);
    check(re.N == p.EA * 1.0e-3 && re.M == p.EI * 2.0e-3 && !re.yielded,
          "uncapped: bit-for-bit elastic");
}

// ------------------------------------------------------------------------------------------
// (2) Consistent tangent vs central FD. Inside a region the map is PIECEWISE LINEAR -> FD
// matches exactly.
// ------------------------------------------------------------------------------------------
void test_tangent_fd() {
    std::printf("-- (2) consistent tangent vs central FD (region interiors) --\n");
    plate::PlateProps p;
    p.EA = 2.0e6; p.EI = 2.0e4; p.nu = 0.15; p.Mp = 150.0; p.Np = 3000.0;
    // Representative points (away from region BOUNDARIES): elastic interior, the 4
    // quadrant surfaces, the two corners.
    const struct { double Nt, Mt; const char* name; } pts[] = {
        {0.4 * p.Np, 0.3 * p.Mp, "elastic interior"},
        {0.9 * p.Np, 0.9 * p.Mp, "surface (+,+)"},
        {-0.9 * p.Np, 0.9 * p.Mp, "surface (-,+)"},
        {0.9 * p.Np, -0.9 * p.Mp, "surface (+,-)"},
        {-0.9 * p.Np, -0.9 * p.Mp, "surface (-,-)"},
        {0.05 * p.Np, 2.5 * p.Mp, "corner V_M"},
        {2.5 * p.Np, 0.05 * p.Mp, "corner V_N"},
    };
    const double he = 1e-8, hk = 1e-9;   // small steps that stay inside the region
    bool all_ok = true;
    for (const auto& q : pts) {
        const double eps = q.Nt / p.EA, kap = q.Mt / p.EI;
        const auto r0 = plate::mn_return(p, eps, kap, 0.0, 0.0);
        auto NM = [&](double e2, double k2, double& N, double& M) {
            const auto rr = plate::mn_return(p, e2, k2, 0.0, 0.0);
            N = rr.N; M = rr.M;
        };
        double Np1, Mp1, Nm1, Mm1;
        NM(eps + he, kap, Np1, Mp1); NM(eps - he, kap, Nm1, Mm1);
        const double dNde = (Np1 - Nm1) / (2 * he), dMde = (Mp1 - Mm1) / (2 * he);
        NM(eps, kap + hk, Np1, Mp1); NM(eps, kap - hk, Nm1, Mm1);
        const double dNdk = (Np1 - Nm1) / (2 * hk), dMdk = (Mp1 - Mm1) / (2 * hk);
        const double scale = p.EA;
        const bool ok = std::fabs(dNde - r0.Dnn) < 1e-4 * scale &&
                        std::fabs(dNdk - r0.Dnm) < 1e-4 * scale &&
                        std::fabs(dMde - r0.Dnm) < 1e-4 * scale &&
                        std::fabs(dMdk - r0.Dmm) < 1e-4 * scale;
        if (!ok) {
            std::printf("   %s: D=[%.4g %.4g; %.4g] FD=[%.4g %.4g %.4g %.4g]\n", q.name,
                        r0.Dnn, r0.Dnm, r0.Dmm, dNde, dNdk, dMde, dMdk);
            all_ok = false;
        }
    }
    check(all_ok, "D_ep == FD at every representative point (elastic/4 surfaces/2 corners)");
}

// ------------------------------------------------------------------------------------------
// (3) Element elastic-limit identity: caps unbounded -> quadrature f_int == K*u (round-off).
// ------------------------------------------------------------------------------------------
void test_element_elastic_identity() {
    std::printf("-- (3) element elastic-limit identity (quadrature vs K*u) --\n");
    plate::PlateProps p;
    p.EA = 5.0e6; p.EI = 8.5e3; p.nu = 0.2;   // caps unbounded by default
    plate::NodeCoords X;
    X << 0.0, 0.0, 2.0, 0.4, 1.0, 0.18;       // slightly curved 3-node
    plate::Dof u;
    u << 1e-3, -2e-3, 3e-4, -4e-4, 5e-4, -6e-4, 7e-5, 8e-4, -9e-5;
    plate::Dof f;
    plate::ElementMatrix K;
    std::array<double, 6> sc{}, st{};
    plate::internal_force_plastic(X, p, u, sc.data(), st.data(), f, &K);
    const plate::Dof f_el = plate::stiffness(X, p) * u;
    const plate::ElementMatrix K_el = plate::stiffness(X, p);
    check((f - f_el).norm() <= 1e-12 * (1.0 + f_el.norm()), "3-node: f quadrature == K*u (1e-12)");
    check((K - K_el).norm() <= 1e-12 * K_el.norm(), "3-node: K_T == K_elastic (1e-12)");
    bool state_zero = true;
    for (double v : st) if (v != 0.0) state_zero = false;
    check(state_zero, "3-node: an elastic call writes no state (bit 0)");

    plate::NodeCoords5 X5;
    X5 << 0.0, 0.0, 0.5, 0.05, 1.0, 0.08, 1.5, 0.05, 2.0, 0.0;
    plate::Dof5 u5;
    for (int i = 0; i < 15; ++i) u5(i) = std::sin(1.0 + i) * 1e-3;
    plate::Dof5 f5;
    plate::ElementMatrix5 K5;
    std::array<double, 10> sc5{}, st5{};
    plate::internal_force_plastic5(X5, p, u5, sc5.data(), st5.data(), f5, &K5);
    const plate::Dof5 f5_el = plate::stiffness5(X5, p) * u5;
    const plate::ElementMatrix5 K5_el = plate::stiffness5(X5, p);
    check((f5 - f5_el).norm() <= 1e-12 * (1.0 + f5_el.norm()), "5-node: f quadrature == K*u");
    check((K5 - K5_el).norm() <= 1e-12 * K5_el.norm(), "5-node: K_T == K_elastic");
}

// ------------------------------------------------------------------------------------------
// BVP harness (the test_plate_soil pattern): a W x H soft-soil block, a plate chain on top,
// ends supported (uy = 0), a point load at the center. free_right_ux: the top-right corner
// ux is left FREE (so an axial preload can be applied, the M-N interaction test).
// ------------------------------------------------------------------------------------------
constexpr double kW = 10.0, kH = 2.0;

struct Bvp {
    Mesh mesh;
    DofMap dofs{0, 2};
    std::vector<PlateElement> plates;
    std::vector<int> top;
    int mid_node = -1;
};
Bvp make_bvp(const plate::PlateProps& pr, bool free_right_ux) {
    Bvp b;
    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    b.mesh = katai::mesh::generate_structured_tri6(domain, 8, 2);
    b.dofs = DofMap(b.mesh.node_count, 2);
    std::vector<int> top(b.mesh.top_nodes.begin(), b.mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int c) { return b.mesh.x[a] < b.mesh.x[c]; });
    std::vector<int> rot(b.mesh.node_count, -1);
    for (int n : top) rot[n] = b.dofs.add_extra_dof();
    for (size_t e = 0; e * 2 + 2 < top.size(); ++e)
        b.plates.push_back(PlateElement{{top[2 * e], top[2 * e + 2], top[2 * e + 1]},
                                        {rot[top[2 * e]], rot[top[2 * e + 2]], rot[top[2 * e + 1]]},
                                        pr});
    b.top = top;
    b.mid_node = top[(top.size() - 1) / 2];
    for (int n : b.mesh.bottom_nodes) { b.dofs.fix_node_component(n, 0); b.dofs.fix_node_component(n, 1); }
    for (int n : b.mesh.left_nodes) b.dofs.fix_node_component(n, 0);
    for (int n : b.mesh.right_nodes)
        if (!(free_right_ux && n == top.back())) b.dofs.fix_node_component(n, 0);
    b.dofs.fix_node_component(top.front(), 1);
    b.dofs.fix_node_component(top.back(), 1);
    return b;
}
// The critical Gauss lever arm s_g = max_g min(x_g, W - x_g) (distance from the support).
double critical_gauss_lever(const Bvp& b) {
    double sg = 0.0;
    for (const auto& pe : b.plates) {
        plate::NodeCoords X;
        for (int k = 0; k < 3; ++k) { X(k, 0) = b.mesh.x[pe.nodes[k]]; X(k, 1) = b.mesh.y[pe.nodes[k]]; }
        for (double xi : plate::kBendGaussXi) {
            const Eigen::Vector3d N = plate::shape(xi);
            const double x = N(0) * X(0, 0) + N(1) * X(1, 0) + N(2) * X(2, 0);
            sg = std::max(sg, std::min(x, kW - x));
        }
    }
    return sg;
}

// Shared runner: central load P (+ optional right-end axial T), the result + mid deflection.
struct BvpRun {
    katai::core::NewtonResult nr;
    double u_mid = 0.0;
    double max_gauss_M = 0.0;      // max capped Gauss |M| (with the committed state)
    double crit_f = -1.0;          // |N|/Np + |M|/Mp - 1 at the most loaded Gauss (-1 if uncapped)
    double crit_N = 0.0;           // N at the same point
};
BvpRun run_bvp(const plate::PlateProps& props, const MaterialModel& soil, double P, double T,
               bool free_right_ux) {
    Bvp bb = make_bvp(props, free_right_ux);
    bb.dofs.finalize();
    Eigen::VectorXd f = Eigen::VectorXd::Zero(bb.dofs.equation_count());
    f(bb.dofs.equation(bb.dofs.global_dof(bb.mid_node, 1))) = -P;
    if (T != 0.0) f(bb.dofs.equation(bb.dofs.global_dof(bb.top.back(), 0))) = T;
    BvpRun out;
    out.nr = katai::core::solve_nonlinear(bb.mesh, bb.dofs, {soil}, f, solve_gen,
                                          {40, 60, 1e-9}, {}, {},
                                          katai::core::Structures{bb.plates, {}});
    out.u_mid = out.nr.displacement[bb.dofs.global_dof(bb.mid_node, 1)];
    double best = -1.0;
    for (size_t ei = 0; ei < bb.plates.size(); ++ei) {
        const auto& pe = bb.plates[ei];
        plate::NodeCoords X;
        plate::Dof u;
        for (int k = 0; k < 3; ++k) {
            X(k, 0) = bb.mesh.x[pe.nodes[k]]; X(k, 1) = bb.mesh.y[pe.nodes[k]];
            u(3 * k + 0) = out.nr.displacement[bb.dofs.global_dof(pe.nodes[k], 0)];
            u(3 * k + 1) = out.nr.displacement[bb.dofs.global_dof(pe.nodes[k], 1)];
            u(3 * k + 2) = out.nr.displacement[pe.rot_dof[k]];
        }
        std::array<double, 3> Ng, Mg;
        plate::gauss_forces_plastic(X, pe.props, u,
                                    out.nr.plate_plastic.data() + ei * plate::kPlasticStateSize,
                                    Ng, Mg);
        for (int q = 0; q < 3; ++q) {
            out.max_gauss_M = std::max(out.max_gauss_M, std::fabs(Mg[q]));
            if (props.Mp > 0.0 && props.Np > 0.0) {
                const double fq = std::fabs(Ng[q]) / props.Np + std::fabs(Mg[q]) / props.Mp - 1.0;
                if (fq > best) { best = fq; out.crit_f = fq; out.crit_N = Ng[q]; }
            }
        }
    }
    return out;
}

// ------------------------------------------------------------------------------------------
// (4) Yield onset + saturation + softening. Closed form P_lim = 2*Mp/s_g (the critical
// GAUSS lever). On an elastic soil bed a single hinge is not a mechanism (the load flows
// into the soil, the solve converges) -- pins: the onset bracket (0.9x bit-zero / 1.1x
// yields), full-Mp saturation, softening.
// ------------------------------------------------------------------------------------------
void test_bvp_limit_load() {
    std::printf("-- (4) BVP yield onset + Mp saturation + softening --\n");
    const double E = 3.0e7, d = 0.4, nu = 0.15;
    plate::PlateProps pr;
    pr.EA = E * d; pr.EI = E * d * d * d / 12.0; pr.nu = nu;
    pr.Mp = 100.0;    // plastic moment [kNm/m]
    const MaterialModel soil{MaterialType::LinearElastic, 10.0, 0.3};   // very soft soil

    Bvp b = make_bvp(pr, false);
    const double sg = critical_gauss_lever(b);
    const double P_lim = 2.0 * pr.Mp / sg;
    std::printf("   s_g = %.4f m -> P_lim = %.4f kN\n", sg, P_lim);

    // (a) Elastik-limit kimligi: 0.9*P_lim'de kapakli plate == kapaksiz plate.
    plate::PlateProps pr_el = pr; pr_el.Mp = -1.0;                     // kapaksiz ikiz
    const BvpRun r09 = run_bvp(pr, soil, 0.9 * P_lim, 0.0, false);
    const BvpRun r09e = run_bvp(pr_el, soil, 0.9 * P_lim, 0.0, false);
    check(r09.nr.converged && r09e.nr.converged, "0.9*P_lim: both runs converged");
    std::printf("   0.9*P_lim: u_cap = %.10e, u_el = %.10e (diff %.2e)\n", r09.u_mid, r09e.u_mid,
                std::fabs(r09.u_mid - r09e.u_mid));
    check(close(r09.u_mid, r09e.u_mid, 1e-9), "below yield capped == uncapped (1e-9)");
    bool zero09 = true;
    for (double v : r09.nr.plate_plastic) if (v != 0.0) zero09 = false;
    check(zero09, "0.9*P_lim: plastic state BIT-ZERO (not above the closed-form onset)");

    // (b) Onset bracket: yield starts at 1.1*P_lim -> P_lim = 2*Mp/s_g verified (+-10% bracket).
    const BvpRun r11 = run_bvp(pr, soil, 1.1 * P_lim, 0.0, false);
    bool any11 = false;
    for (double v : r11.nr.plate_plastic) if (v != 0.0) any11 = true;
    check(r11.nr.converged && any11, "1.1*P_lim: a hinge formed (onset bracket closed)");

    // (c) Saturation + softening: at 1.5*P_lim the Gauss |M| sits EXACTLY at Mp; the
    // deflection grows disproportionately.
    const BvpRun r15 = run_bvp(pr, soil, 1.5 * P_lim, 0.0, false);
    std::printf("   1.5*P_lim: conv=%d, max Gauss |M| = %.10f, u_mid = %.4e (0.9x: %.4e)\n",
                (int)r15.nr.converged, r15.max_gauss_M, r15.u_mid, r09.u_mid);
    check(r15.nr.converged, "1.5*P_lim converges (the soil bed carries -- not a mechanism)");
    check(close(r15.max_gauss_M, pr.Mp, 1e-9), "critical Gauss |M| saturated EXACTLY at Mp (clamp exact)");
    check(std::fabs(r15.u_mid) > 10.0 * std::fabs(r09.u_mid) * (1.5 / 0.9),
          "post-hinge softening: deflection >10x the proportional elastic");
    bool any_kp = false;
    for (size_t i = 1; i < r15.nr.plate_plastic.size(); i += 2)
        if (r15.nr.plate_plastic[i] != 0.0) any_kp = true;
    check(any_kp, "permanent kappa_p accumulated (the hinge is real)");
}

// ------------------------------------------------------------------------------------------
// (5) M-N interaction: axial pre-tension T = Np/2 -> the critical Gauss point sits ON THE
// DIAMOND SURFACE (|N|/Np + |M|/Mp = 1) and the moment plateau ~ Mp*(1 - N/Np); the
// Np-unbounded control reaches full Mp -> the drop truly comes from the interaction (not
// the moment cap).
// ------------------------------------------------------------------------------------------
void test_bvp_mn_interaction() {
    std::printf("-- (5) BVP M-N interaction (axial load lowers the moment plateau) --\n");
    const double E = 3.0e7, d = 0.4, nu = 0.15;
    plate::PlateProps pr;
    pr.EA = E * d; pr.EI = E * d * d * d / 12.0; pr.nu = nu;
    pr.Mp = 100.0; pr.Np = 2000.0;
    const MaterialModel soil{MaterialType::LinearElastic, 10.0, 0.3};
    const double T = 0.5 * pr.Np;   // axial pre-tension [kN/m]

    Bvp b = make_bvp(pr, true);
    const double sg = critical_gauss_lever(b);
    const double P_lim_full = 2.0 * pr.Mp / sg;   // full capacity without axial load

    const BvpRun r = run_bvp(pr, soil, 1.5 * P_lim_full, T, /*free_right_ux=*/true);
    std::printf("   conv=%d, critical f = %.3e, N = %.1f (T = %.1f), max Gauss |M| = %.4f "
                "(Mp*(1-T/Np) = %.1f)\n",
                (int)r.nr.converged, r.crit_f, r.crit_N, T, r.max_gauss_M,
                pr.Mp * (1.0 - T / pr.Np));
    check(r.nr.converged, "axial + bending run converged");
    check(std::fabs(r.crit_f) <= 1e-9, "critical Gauss ON THE DIAMOND SURFACE (|N|/Np + |M|/Mp = 1)");
    check(r.crit_N > 0.4 * pr.Np && r.crit_N < 0.6 * pr.Np,
          "N ~ T at the critical point (the axial preload carried into the plate)");
    check(r.max_gauss_M < 0.75 * pr.Mp,
          "moment plateau clearly below Mp (the diamond section; T/Np = 0.5)");

    // Control: Np unbounded -> the same loads reach full Mp (the drop comes from the M-N interaction).
    plate::PlateProps pr_noN = pr; pr_noN.Np = -1.0;
    const BvpRun rc = run_bvp(pr_noN, soil, 1.5 * P_lim_full, T, true);
    std::printf("   control (Np unbounded): conv=%d, max Gauss |M| = %.10f\n",
                (int)rc.nr.converged, rc.max_gauss_M);
    check(rc.nr.converged && close(rc.max_gauss_M, pr.Mp, 1e-9),
          "with Np unbounded the plateau is full Mp -> the drop truly comes from the interaction");
}

// ------------------------------------------------------------------------------------------
// (6) Diagram: the capped Gauss M expand to the stations by Lagrange (the PLAXIS
// stress-point extrapolation rule); an empty plastic state = an elastic report (D6b, the
// linear dynamic envelope rule).
// ------------------------------------------------------------------------------------------
void test_bvp_diagram() {
    std::printf("-- (6) diagram (Gauss->station Lagrange; empty state = elastic report) --\n");
    const double E = 3.0e7, d = 0.4, nu = 0.15;
    plate::PlateProps pr;
    pr.EA = E * d; pr.EI = E * d * d * d / 12.0; pr.nu = nu;
    pr.Mp = 100.0;
    const MaterialModel soil{MaterialType::LinearElastic, 10.0, 0.3};

    Bvp b = make_bvp(pr, false);
    const double sg = critical_gauss_lever(b);
    const double P_lim = 2.0 * pr.Mp / sg;
    const BvpRun r = run_bvp(pr, soil, 1.5 * P_lim, 0.0, false);
    check(r.nr.converged, "the 1.5*P_lim run converged");

    Bvp bb = make_bvp(pr, false);
    bb.dofs.finalize();
    const auto diag = katai::core::plate_force_diagram(bb.plates, bb.mesh, bb.dofs,
                                                       r.nr.displacement, r.nr.plate_plastic);
    double maxMst = 0.0;
    for (const auto& st : diag) maxMst = std::max(maxMst, std::fabs(st.M));
    std::printf("   capped diagram peak |M| = %.4f (Mp = %.1f)\n", maxMst, pr.Mp);
    check(maxMst > 0.90 * pr.Mp && maxMst < 1.10 * pr.Mp,
          "plastic diagram peak ~ Mp (Gauss->station Lagrange; 10% band)");

    // EMPTY plastic state -> ELASTIC report (no cap): the linear-dynamic envelope rule (D6b).
    const auto diag_el = katai::core::plate_force_diagram(bb.plates, bb.mesh, bb.dofs,
                                                          r.nr.displacement);
    double maxMel = 0.0;
    for (const auto& st : diag_el) maxMel = std::max(maxMel, std::fabs(st.M));
    std::printf("   elastic (empty-state) diagram peak |M| = %.1f\n", maxMel);
    check(maxMel > 2.0 * pr.Mp,
          "empty state = uncapped elastic report (M >> Mp appears at the hinged displacement)");
}

// ------------------------------------------------------------------------------------------
// (7) tri15 (5-node plate): the same limit-load closed form.
// ------------------------------------------------------------------------------------------
void test_bvp_tri15_limit() {
    std::printf("-- (7) tri15 / 5-node plate limit load --\n");
    const double E = 3.0e7, d = 0.4, nu = 0.15;
    plate::PlateProps pr;
    pr.EA = E * d; pr.EI = E * d * d * d / 12.0; pr.nu = nu;
    pr.Mp = 100.0;
    const MaterialModel soil{MaterialType::LinearElastic, 10.0, 0.3};

    const RectangularDomain domain{0.0, 0.0, kW, kH, 0};
    const Mesh mesh = katai::mesh::generate_structured_tri15(domain, 4, 1);
    DofMap dofs(mesh.node_count, 2);
    std::vector<int> top(mesh.top_nodes.begin(), mesh.top_nodes.end());
    std::sort(top.begin(), top.end(), [&](int a, int c) { return mesh.x[a] < mesh.x[c]; });
    std::vector<int> rot(mesh.node_count, -1);
    for (int n : top) rot[n] = dofs.add_extra_dof();
    std::vector<PlateElement5> plates5;
    for (size_t e = 0; e * 4 + 4 < top.size(); ++e) {
        PlateElement5 pe;
        for (int k = 0; k < 5; ++k) { pe.nodes[k] = top[4 * e + k]; pe.rot_dof[k] = rot[top[4 * e + k]]; }
        pe.props = pr;
        plates5.push_back(pe);
    }
    const int mid_node = top[(top.size() - 1) / 2];
    for (int n : mesh.bottom_nodes) { dofs.fix_node_component(n, 0); dofs.fix_node_component(n, 1); }
    for (int n : mesh.left_nodes) dofs.fix_node_component(n, 0);
    for (int n : mesh.right_nodes) dofs.fix_node_component(n, 0);
    dofs.fix_node_component(top.front(), 1);
    dofs.fix_node_component(top.back(), 1);
    dofs.finalize();

    double sg = 0.0;
    for (const auto& pe : plates5) {
        plate::NodeCoords5 X;
        for (int k = 0; k < 5; ++k) { X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]]; }
        for (double xi : plate::kBendGaussXi5) {
            const auto N = plate::detail::shape5(xi);
            double x = 0.0;
            for (int k = 0; k < 5; ++k) x += N(k) * X(k, 0);
            sg = std::max(sg, std::min(x, kW - x));
        }
    }
    const double P_lim = 2.0 * pr.Mp / sg;
    std::printf("   s_g = %.4f m -> P_lim = %.4f kN\n", sg, P_lim);

    auto run_at = [&](double P, katai::core::NewtonResult& out) {
        Eigen::VectorXd f = Eigen::VectorXd::Zero(dofs.equation_count());
        f(dofs.equation(dofs.global_dof(mid_node, 1))) = -P;
        out = katai::core::solve_nonlinear(mesh, dofs, {soil}, f, solve_gen,
                                           {40, 60, 1e-9}, {}, {},
                                           katai::core::Structures{{}, {}, {}, {}, plates5, {}});
    };
    katai::core::NewtonResult rA, rB;
    run_at(0.85 * P_lim, rA);
    run_at(1.5 * P_lim, rB);
    bool zeroA = true, anyB = false;
    for (double v : rA.plate5_plastic) if (v != 0.0) zeroA = false;
    for (double v : rB.plate5_plastic) if (v != 0.0) anyB = true;
    // Max capped Gauss |M| (with the committed state).
    double maxM = 0.0;
    for (size_t ei = 0; ei < plates5.size(); ++ei) {
        const auto& pe = plates5[ei];
        plate::NodeCoords5 X;
        plate::Dof5 u;
        for (int k = 0; k < 5; ++k) {
            X(k, 0) = mesh.x[pe.nodes[k]]; X(k, 1) = mesh.y[pe.nodes[k]];
            u(3 * k + 0) = rB.displacement[dofs.global_dof(pe.nodes[k], 0)];
            u(3 * k + 1) = rB.displacement[dofs.global_dof(pe.nodes[k], 1)];
            u(3 * k + 2) = rB.displacement[pe.rot_dof[k]];
        }
        std::array<double, 5> Ng, Mg;
        plate::gauss_forces_plastic5(X, pe.props, u,
                                     rB.plate5_plastic.data() + ei * plate::kPlasticStateSize5,
                                     Ng, Mg);
        for (double m : Mg) maxM = std::max(maxM, std::fabs(m));
    }
    std::printf("   0.85x: conv=%d state-zero=%d | 1.5x: conv=%d max Gauss |M| = %.10f\n",
                (int)rA.converged, (int)zeroA, (int)rB.converged, maxM);
    check(rA.converged && zeroA, "tri15: 0.85*P_lim elastic (state bit-zero)");
    check(rB.converged && anyB, "tri15: 1.5*P_lim formed a hinge");
    check(close(maxM, pr.Mp, 1e-9), "tri15: critical Gauss |M| saturated EXACTLY at Mp");
}

}  // namespace

int main() {
    test_return_map_scan();
    test_tangent_fd();
    test_element_elastic_identity();
    test_bvp_limit_load();
    test_bvp_mn_interaction();
    test_bvp_diagram();
    test_bvp_tri15_limit();
    if (g_failures == 0) {
        std::printf("\nOK: plate M-N hinge verified (return map CPP-exact + FD tangent + "
                    "elastic identity + onset 2Mp/s_g + full-Mp saturation + M-N diamond "
                    "surface + diagram + tri15)\n");
        return 0;
    }
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
}
