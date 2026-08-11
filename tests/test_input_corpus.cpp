// The input corpus (roadmap section 6.3): every case is a CHECKED-IN .k2d under tests/corpus/
// plus the same problem built programmatically, with three assertions per case:
//
//   (1) IDENTITY -- project_to_json(programmatic build) is byte-identical to the checked-in
//       file (the writer is deterministic), so the file IS the programmatic problem and a
//       schema change cannot silently strand the corpus behind the format;
//   (2) CONTRACT -- the file loads with zero reader notes and validates with zero errors
//       (warnings are printed, never hidden);
//   (3) ORACLE -- the solve runs FROM THE FILE-LOADED project only: the mesh comes from the
//       file's own settings (mesh_from_project(pr)) and the run kind from the file's
//       initial_procedure (initial_phase_from). The declared closed form must hold. Nothing
//       from the programmatic build reaches the solve, so "reproducible from the file" is
//       what is actually tested.
//
// Regenerating after a deliberate format change: run with KATAI_CORPUS_WRITE=1 in the
// environment; the test rewrites every case file and reports it loudly -- the git diff of
// tests/corpus/ is then the review artifact. The files are compact single-line JSON (the
// writer's canonical form); `python -m json.tool` pretty-prints them for reading.
//
// verify: KV-CON-002
//   oracle:   closed_form
//   source:   Terzaghi (1943) one-dimensional consolidation, degree-of-consolidation series solution
//   locator:  U(Tv) = 1 - sum_j (2/M^2) exp(-M^2 Tv), M = (2j+1) pi/2, Tv = cv t / H_dr^2, cv = k Eoed / gamma_w (stated in full)
//   quantity: degree of consolidation U(Tv) of a laterally confined column drained at the top, surcharge applied at t=0+, run from the checked-in tests/corpus/kv-con-002-terzaghi-column.k2d [-]
//   expected: the series above at Tv ~ 0.2 / 0.4 / 0.6 / 0.9 and at the final Tv; undrained excess pore ~ q at early time, ~ 0 at Tv ~ 2
//   band:     3% on U, as asserted below -- measured -1.2% .. -0.4% at the sampled Tv on the file's own 0.4 m tri6 mesh with 120 time steps
//
// verify: KV-FND-008
//   oracle:   closed_form
//   source:   classical elasticity: Boussinesq point solution integrated over a uniform strip on a half-plane
//   locator:  sigma_z = (q/pi) [alpha + sin(alpha) cos(theta1 + theta2)], theta_i = atan((x -+ a)/z) measured from the strip edges, alpha = theta1 - theta2 (stated in full)
//   quantity: vertical stress under a uniform strip load on a weightless elastic half-plane (independent of E and nu), run from the checked-in tests/corpus/kv-fnd-008-strip-load.k2d [kPa]
//   expected: the closed form above, evaluated at each sampled node's exact (x, z)
//   band:     3% at z = 2..6 m under the strip, as asserted below -- measured +0.9% .. +1.4%; a half-plane approximated by a 40 x 20 m finite domain plus discretisation
//
// verify: KV-FND-009
//   oracle:   closed_form
//   source:   Flamant (1892) line load on an elastic half-plane
//   locator:  sigma_z = (2 P / pi) z^3 / (x^2 + z^2)^2, x measured from the load line (stated in full)
//   quantity: vertical stress under a concentrated line load on a weightless elastic half-plane, tri15 elements, run from the checked-in tests/corpus/kv-fnd-009-flamant-line-load.k2d [kPa]
//   expected: the closed form above, evaluated at each sampled node's exact (x, z)
//   band:     3% at z = 2..6 m, as asserted below -- measured +0.0% .. +1.3%; the concentrated nodal load is singular at its point of application, so the comparison stays at depth
//
// verify: KV-NUM-003
//   oracle:   closed_form
//   source:   geostatic K0 equilibrium of submerged level ground; Terzaghi effective-stress principle
//   locator:  sigma'_v = -gamma' (H - y), gamma' = gamma_sat - gamma_w; sigma'_h = K0 sigma'_v with K0 = 1 - sin(phi) = 0.5; u_w = gamma_w (H - y); undisturbed ground does not displace (stated in full)
//   quantity: nodal effective stresses, hydrostatic pore pressure and residual displacement of a submerged block under the K0 procedure, run from the checked-in tests/corpus/kv-num-003-k0-geostatic-block.k2d [kPa; m]
//   expected: the closed forms above at every node
//   band:     1e-6 kPa absolute on the recovered stresses (measured ~5e-13: the K0 procedure sets this field by construction, the residual is round-off), 1e-9 kPa on pore pressure, 1e-6 m on displacement, as asserted below
//
// verify: KV-CST-001
//   oracle:   closed_form
//   source:   undrained (Skempton) 1D confined compression; docs/references/effective-stress-formulation.md
//   locator:  M_u = M' + K_w/n, K_w/n = 3 (nu_u - nu) / ((1 - 2 nu_u)(1 + nu)) K', nu_u = 0.495; u_y = -q H / M_u; sigma'_yy = -M' q / M_u (stated in full)
//   quantity: surface settlement and mid-column effective vertical stress of a weightless laterally confined column loaded undrained (A), run from the checked-in tests/corpus/kv-cst-001-undrained-column.k2d [m; kPa]
//   expected: the closed forms above
//   band:     2% on settlement and 3% on effective stress, as asserted below -- measured +0.0% on both (the 1D confined fields are exactly representable on this mesh)
//
// verify: KV-SLP-001
//   oracle:   published_benchmark
//   source:   Griffiths and Lane (1999), "Slope stability analysis by finite elements", Geotechnique 49(3); Rocscience Slide verification problem 1
//   locator:  homogeneous 1:2 slope on a foundation layer, gamma = 20.2 kN/m3, c' = 3 kPa, phi' = 19.6 deg, psi = 0; published FoS: Bishop 0.988, Spencer 0.987, Phase2 T6 0.997
//   quantity: slope factor of safety by phi-c reduction, run as the file's INITIAL procedure (initial_procedure = Safety) from the checked-in tests/corpus/kv-slp-001-griffiths-lane-slope.k2d [-]
//   expected: FoS ~ 0.99 (the multi-method consensus above)
//   band:     8% vs 0.99, as asserted below -- SRM-versus-LEM method scatter plus the file's own coarse mesh; measured FoS 1.010 (+2.1%); the mechanism must also displace
//
// verify: KV-EXC-001
//   oracle:   closed_form
//   source:   1D elastic unloading of a laterally confined column (oedometric heave) under the staged-construction SumMstage rule
//   locator:  heave u = +gamma_f h_exc H_rem / E_oed, E_oed = E (1 - nu) / ((1 + nu)(1 - 2 nu)); the base total stress sheds exactly the excavated weight (stated in full)
//   quantity: pit-floor heave and base total vertical stress after deactivating the upper layer in a staged phase, run from the checked-in tests/corpus/kv-exc-001-staged-excavation.k2d [m; kPa]
//   expected: u(y = 6) = +17 * 4 * 6 / E_oed; sigma_v(base) = -18 * 6; the initial full-geometry K0 phase does not displace
//   band:     2%, as asserted below -- measured +0.0% on the heave and -0.0% on the base stress on the file's own mesh
//
// verify: KV-DYN-002
//   oracle:   closed_form
//   source:   1D SH site response of a damped elastic shear column on a rigid base, fundamental-mode resonance (Kramer 1996, Geotechnical Earthquake Engineering, ch. 7; the same transfer-function amplification verified independently in test_dynamics D1)
//   locator:  f_1 = Vs/(4H); at resonance |a_surf| = 2 A/(pi xi) and |u_surf| = (4/pi) A/(w_1^2 2 xi), w_1 = 2 pi f_1 (stated in full)
//   quantity: peak surface displacement and acceleration of the column driven harmonically at f_1, and the off-resonance response at f_1/3, run from the checked-in tests/corpus/kv-dyn-002-resonant-column.k2d [m; m/s2]
//   expected: |u_surf| = (4/pi)/(w_1^2 2 xi) and |a_surf| = 2/(pi xi) (A = 1); the f_1/3 response is a small fraction of the resonant one; the initial K0 phase does not displace
//   band:     2% on |u_surf|, 1% on |a_surf|, off/res < 0.2 -- measured -0.6% / -0.2% / 0.140 on the file's own 0.4 m tri6 mesh
//
// verify: KV-FLW-001
//   oracle:   closed_form
//   source:   Charny (1951) exact-discharge theorem for the unconfined rectangular dam with a seepage face (Polubarinova-Kochina 1962, Theory of Ground Water Movement); formulation docs/references/seepage-formulation.md
//   locator:  q = k (h1^2 - h2^2) / (2 L), exact regardless of the free-surface shape; the crest above the phreatic surface stays unsaturated (stated in full)
//   quantity: free-surface seepage discharge through a rectangular dam between reservoir h1 and tailwater h2 with a seepage face, run from the checked-in tests/corpus/kv-flw-001-charny-unconfined-dam.k2d [m3/day per m]
//   expected: q = 0.5 (5^2 - 1^2) / (2 * 10) = 0.6; zero saturated crest nodes; global mass balance ~ 0
//   band:     2% on q (measured +1.02%), balance < 1e-9 (measured 8e-14), zero saturated crest nodes -- on the file's own 0.35 m tri6 mesh
//
// verify: KV-DYN-003
//   oracle:   published_benchmark
//   source:   El Centro 1940-05-18 NS (S00E) strong-motion record; provenance and download identity recorded in tests/data/elcentro-1940-ns.md
//   locator:  published characteristic PGA about 0.319 g at t of about 2 s; 5%-damped peak Sa within the widely published 2.0-3.5 x PGA band; the record travels IN the .k2d, so the file alone reproduces the run
//   quantity: the accelerogram stored in the checked-in tests/corpus/kv-dyn-003-el-centro-two-layer.k2d (PGA, peak timing, 5% spectral shape) and the compliant-base two-layer site response driven by it [g; m/s2]
//   expected: file-borne record: 1560 samples at 0.02 s, PGA 0.319 g at ~2 s, peak Sa/PGA in 2.0-3.5, Sa(3 s) < 0.15 g; the compliant-base run solves from the file alone and reproduces the verified path's surface response
//   band:     record identity near-exact (0.005 g; measured 0.31882 g at 2.02 s); spectrum against the published bands (measured 2.87x at 0.19 s, Sa(3 s) 0.118 g); surface response 5% reproduction band around the measured 5.149 m/s2; Sa(T_min = 0.05 s)/PGA measured 1.05, banded 0.9-1.3 (the physics direction rigid >= compliant is pinned in test_real_record)
//
// verify: KV-FND-010
//   oracle:   closed_form
//   source:   Prandtl (1921) bearing-capacity wedge solution for a weightless Tresca (phi = 0) half-space
//   locator:  q_ult = (2 + pi) c, i.e. N_c = 2 + pi ~= 5.14; the phase is loaded PAST collapse, so the honest non-convergence message carries the equilibrated fraction and load_factor * q_applied IS the incremental limit load (stated in full)
//   quantity: bearing capacity factor N_c of a flexible strip footing on a weightless c-only half-space, loaded past collapse in a staged phase, run from the checked-in tests/corpus/kv-fnd-010-prandtl-strip-footing.k2d [-]
//   expected: the collapse phase does NOT fully converge (that is the result); N_c = load_factor * 6c / c ~= 2 + pi; the weightless gravity initial does not displace
//   band:     2%, as asserted below -- measured +0.6% on the file's own 0.4 m tri15 mesh (the direct structured-mesh benchmark KV-FND-005 measures +1.1% at tri15)
//
// verify: KV-FND-011
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Gibson (1967)
//   locator:  Section 2.2, strip load on incompressible Gibson soil (E = 299 z, nu = 0.495); half-space closed form s = q / (2 dG/dz), uniform under the load
//   quantity: centreline surface settlement under a q = 10 kPa strip load on a 4 m Gibson layer, E(y) via the schema's E_inc/y_ref profile, run from the checked-in tests/corpus/kv-fnd-011-gibson-strip-load.k2d [m]
//   expected: 0.047 (PLAXIS, same finite layer); the half-space closed form gives 0.050 and the finite layer must sit BELOW it (shared bias)
//   band:     5% vs the PLAXIS finite-layer value, as asserted below -- measured -3.4% (0.0454) on the file's own 0.15 m tri15 mesh (the direct structured benchmark KV-FND-002 measures -4.0%)
//
// verify: KV-FND-012
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Giroud (1972)
//   locator:  Section 2.1, smooth rigid strip footing on elastic soil; F = 2 (1 + nu) G B s / rho, rho = 0.88
//   quantity: footing force at a prescribed settlement of 10 mm, via the schema's line prescribed displacement (v2) and the reaction output, run from the checked-in tests/corpus/kv-fnd-012-giroud-rigid-footing.k2d [kN/m]
//   expected: 15.15 (analytic); PLAXIS publishes 15.24; the weightless gravity initial does not displace; a footing node carries exactly the imposed u_y
//   band:     2% vs analytic and 3% vs PLAXIS, as asserted below -- measured +1.1% / +0.5% (15.32) on the file's own 0.5 m tri15 mesh (the direct structured benchmark KV-FND-001 measures +1.4% / +0.8%)
//
// verify: KV-FND-013
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); slip-line solution Cox (1962)
//   locator:  Section 3.1, bearing capacity of a smooth rigid circular footing (axisymmetric Mohr-Coulomb), run from the checked-in tests/corpus/kv-fnd-013-cox-circular-footing.k2d
//   quantity: limit pressure p_max from the axisymmetric reaction output at a prescribed settlement of 0.35 m, associated flow (psi = phi -- the slip-line solution is the associated limit load) [kPa]
//   expected: 225.6 (analytic, 141 c); PLAXIS publishes 220.0; the K0 initial phase does not displace; a footing node carries exactly the imposed u_y
//   band:     5% vs analytic, as asserted below -- measured +3.7% (233.9) on the file's own 0.25 m tri15 mesh (the direct structured benchmark KV-FND-003 measures +3.9%; a 0.5 m mesh over-predicts by ~9%, the coarse-mesh bearing bias)
//
// verify: KV-FND-014
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems); analytic solution Davis & Booker (1973)
//   locator:  Section 3.2, smooth strip footing on clay with strength increasing with depth, run from the checked-in tests/corpus/kv-fnd-014-davis-booker-strip-footing.k2d
//   quantity: limit pressure p_max [kPa] with c(z) = c0 + c_inc z and E(z) = E0 + E_inc z via the schema's c_inc / E_inc / y_ref profile, from the reaction output at a prescribed settlement of 30 mm
//   expected: 7.80 (analytic, rho [(2 + pi) c0 + B c_inc / 4]); PLAXIS publishes 7.86; the weightless gravity initial does not displace; a footing node carries exactly the imposed u_y
//   band:     5% vs analytic and 3% vs PLAXIS, as asserted below -- measured +1.4% / +0.6% (7.91) on the file's own 0.5 m tri15 mesh (the direct structured benchmark KV-FND-004 measures +2.8%)
//
// verify: KV-CST-002
//   oracle:   closed_form
//   source:   the Hardening Soil oedometric stiffness law as published in the PLAXIS Material Models Manual: E_oed = E_oed^ref ((c cos(phi) + sigma_1 sin(phi))/(c cos(phi) + p_ref sin(phi)))^m, integrated over one-dimensional primary loading; the same law is verified at the material point against the manual's own figures in test_hardening_soil and test_hs_berlin
//   locator:  with c = 0 the stiffness factor reduces to (sigma_1/p_ref)^m and d eps_1 = d sigma_1 / E_oed gives -eps_1 = (p_ref^m / E_oed^ref) [sigma_1^(1-m)]/(1-m) between the two stress levels (stated in full and evaluated in the test, not called from the material header)
//   quantity: settlement increment of a laterally confined weightless Hardening Soil column when the vertical stress steps from 50 to 200 kPa, run from the checked-in tests/corpus/kv-cst-002-hs-oedometer.k2d [m]
//   expected: the closed form above with E_oed^ref = 30 MPa, p_ref = 100 kPa, m = 0.5, H = 4 m
//   band:     3%, as asserted below -- measured +0.41% on the file's own 0.5 m tri6 mesh with the driver's 40 load steps. The first HS boundary-value case in the corpus: the model was already verified at the material point, this verifies the path from the file through the mesher, the cap return mapping and the load stepping. It also pins a phase convention: a phase reports displacement relative to its own start, so the loading phase's field IS the increment (the seating phase's 0 -> 50 kPa settlement of 0.0164 m is reported separately and is not comparable to the same integral, because the law's stiffness vanishes as sigma -> 0)

// verify: KV-STR-002
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems)
//   locator:  Section 3.3, sliding block for testing interfaces; the manual states the answer as failure force = width * c_w + weight * tan(phi_w) = 4 * 2.5 + 100 * 0.5 = 60 kN/m, and reports 60.4 kN/m for its own run. Block E = 30 GN/m2, nu = 0, gamma = 25 kN/m3 with K0 = 0; interface in a separate elastoplastic data set, E = 3 GN/m2, nu = 0.45, phi_w = 26.6 deg, c_w = 2.5 kN/m2; bottom fully fixed, u_x = 0.1 m prescribed on the left side with u_y free, everything else free. Width 4 m and a weight of 100 kN/m at gamma = 25 fix the block at 4 m x 1 m -- the manual prints the arithmetic rather than the height
//   quantity: the horizontal failure force, as the sum of the reactions on the pushed edge, run from the checked-in tests/corpus/kv-str-002-plaxis-sliding-block.k2d [kN/m]
//   expected: 60.076 -- the manual's own formula evaluated at the phi_w it specifies. Its printed 60 uses tan(phi_w) = 0.5, i.e. 26.565 deg, so the manual is 0.13% self-inconsistent and the closed form is quoted here at the INPUT it publishes; PLAXIS reports 60.4
//   band:     2% vs the closed form and 2% vs PLAXIS, as asserted below -- measured -0.59% (59.7202) on the file's own 0.25 m tri6 mesh, converging monotonically with refinement (-1.31% / -0.59% / -0.27% at 0.5 / 0.25 / 0.125 m). Four further checks make the number more than a coincidence: the block translates rigidly rather than shearing; the force is bit-identical when the imposed slip is doubled, so it is a plateau and not a stiffness reading; adhesion and friction each move the answer by exactly the closed form's amount, so the two terms are reproduced separately; and deleting the interface changes the answer by five orders of magnitude. That last one is a regression sentry: until 2026-08-10 the interface lay along a fixed boundary whose two split sides share coordinates, the boundary conditions fixed both, and the joint was welded shut in silence
//
// verify: KV-STR-003
//   oracle:   published_benchmark
//   source:   PLAXIS 2D Validation Manual, Version 8 (Bentley Systems)
//   locator:  Section 2.3, bending of beams. Two problems on a simply supported span of l = 2 m, with the characteristics of an HEB 200 steel beam, which in plane strain is a plate 1 m wide out of plane: EA = 1.64e6 kN, EI = 1200 kNm2, nu = 0.0, a single point load F = 100 kN at mid-span and a uniformly distributed load q = 100 kN/m. The manual publishes both extremes for both problems: point load M_max = 50.0 kNm and u_max = 13.96 mm, distributed load M_max = 50.0 kNm and u_max = 17.43 mm. Its build is reproduced as stated: the two beams are added to the bottom line of a block cluster with a spacing in between, point fixities at their end points, and the soil cluster deactivated so that only the beams remain on a very coarse mesh
//   quantity: mid-span deflection and peak bending moment of BOTH beams from one run, the deflections from the nodal field and the moments from each element's own force diagram, run from the checked-in tests/corpus/kv-str-003-plaxis-beam-bending.k2d [m; kNm/m]
//   expected: the Mindlin (Timoshenko) closed forms, stated in full and evaluated in the test rather than called from the plate header: w = F l^3/(48 EI) + F l/(4 kGA') and w = 5 q l^4/(384 EI) + q l^2/(8 kGA'), with the manual's own shear rigidity kGA' = k EA/(2(1+nu)), k = 5/6 (Material Models Manual Eq. 18-8). They evaluate to 13.96206 mm and 17.43428 mm, which is exactly where the published 13.96 and 17.43 come from -- a PLAXIS plate is shear-deformable, and Euler-Bernoulli alone would give 13.8896 and 17.3618. M_max = F l / 4 = q l^2 / 8 = 50 kNm
//   band:     1% vs the closed form and 2% vs PLAXIS on the deflections, 2% vs the published 50 kNm on the moments, as asserted below -- measured -0.000% on BOTH deflections (13.96206 / 17.43428 mm) and +0.00% / +1.04% on the moments (50.00000 / 50.52083 kNm) on the file's own 0.25 m tri6 mesh. Five further witnesses make the pair more than a coincidence: the moment DISTRIBUTION follows F s/2 and q s(l-s)/2 station by station (worst 0.0000% and 1.0417% of M_max, not just the peak); the distributed beam's peak overshoots by exactly q h^2/12 -- the parabola the element's linear curvature cannot hold inside one element -- reproduced to five figures at h = 0.5 / 0.25 / 0.125 m (52.08333 / 50.52083 / 50.13021), so the residual is a structural discretisation bias that vanishes with h, and "a very coarse mesh is sufficient" is true of the displacements but not of the peak moment; bending and shear are moved SEPARATELY (EI x 4 divides the bending term alone, EA x 100 drives the answer onto the Euler-Bernoulli limit 13.8896 mm) and the run follows the closed form to -0.000% in each; the two spans are bit-for-bit independent; and with the beams deleted the model has no free DOF at all and the run refuses with "Every DOF is fixed; nothing to solve". That last one is the regression sentry for the fault this case was built to find: until 2026-08-11 a deactivated soil cluster pinned every translation of the beams standing in it (fix_inactive_nodes did not ask whether a structure held the node), so the beams were welded to the outside world along their whole length -- the only free DOFs left were their rotations, which was enough for the solve to converge, report "ok" and hand back max|u| = 0.000000e+00 with no warning of any kind
//
// verify: KV-SLP-002
//   oracle:   published_benchmark
//   source:   Griffiths and Lane (1999), "Slope stability analysis by finite elements", Geotechnique 49(3), Example 1
//   locator:  homogeneous 2:1 slope with no foundation layer (D = 1): phi' = 20 deg, c'/gamma H = 0.05, psi = 0, nominal E' = 1e5 kPa and nu' = 0.3 (their stated values); published FoS: 1.4 by FE (non-convergence at the 1000-iteration ceiling), 1.380 by the Bishop and Morgenstern (1960) charts
//   quantity: slope factor of safety by phi-c reduction, run as the file's INITIAL procedure (initial_procedure = Safety) from the checked-in tests/corpus/kv-slp-002-griffiths-lane-example1.k2d, dimensionalised as H = 10 m, gamma = 20 kN/m3, c' = 10 kPa [-]
//   expected: FoS between the published pair 1.380 (Bishop-Morgenstern) and 1.4 (Griffiths-Lane FE)
//   band:     4% vs Bishop-Morgenstern 1.380, as asserted below -- measured FoS 1.384 (+0.3%) on the file's own 1.0 m tri6 mesh; the mechanism must also displace
#include <katai/analysis/response_spectrum.hpp>
#include <katai/jobs/mesh_builder.hpp>
#include <katai/jobs/driver.hpp>
#include <katai/jobs/flow_driver.hpp>
#include <katai/io/project_io.hpp>
#include <katai/io/validate.hpp>
#include <katai/model/project.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace m = katai::model;
namespace io = katai::io;

namespace {
constexpr double kPi = 3.14159265358979323846;
int g_failures = 0;
void check(bool ok, const char* what) {
    std::printf(ok ? "ok:   %s\n" : "FAIL: %s\n", what);
    if (!ok) ++g_failures;
}

int nearest_node(const katai::mesh::Mesh& mesh, double x, double y) {
    int best = 0;
    double bd = 1e300;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double dx = mesh.x[n] - x, dy = mesh.y[n] - y, d2 = dx * dx + dy * dy;
        if (d2 < bd) { bd = d2; best = n; }
    }
    return best;
}

// ---------------------------------------------------------------- the corpus harness --
struct CorpusCase {
    const char* file;                              // file name under tests/corpus/
    m::Project (*build)();                         // the programmatic path
    void (*oracle)(const m::Project&);             // solves from the FILE-LOADED project
};

void run_case(const CorpusCase& c) {
    std::printf("\n== %s ==\n", c.file);
    const std::string path = std::string(KATAI_CORPUS_DIR) + "/" + c.file;
    const m::Project built = c.build();
    const std::string want = m::project_to_json(built);

    if (std::getenv("KATAI_CORPUS_WRITE")) {
        std::string err;
        check(m::save_project(built, path, &err), "KATAI_CORPUS_WRITE: case file rewritten");
        std::printf("      REGENERATED %s -- review the git diff and commit it.\n", path.c_str());
    }

    // (1) identity: the checked-in bytes are exactly the programmatic build.
    std::ifstream f(path, std::ios::binary);
    check(f.good(), "checked-in case file exists");
    std::stringstream ss;
    ss << f.rdbuf();
    const std::string got = ss.str();
    check(got == want, "checked-in file is byte-identical to the programmatic build");

    // (2) contract: loads cleanly, validates with zero errors.
    m::Project loaded;
    std::string err;
    std::vector<io::Issue> notes;
    check(m::load_project(path, loaded, &err, &notes), "file loads");
    if (!err.empty()) std::printf("      (%s)\n", err.c_str());
    check(notes.empty(), "reader reports no forward-version notes");
    const io::ValidationReport rep = io::validate_project(loaded);
    for (const io::Issue& i : rep.issues)
        std::printf("      %s %s: %s\n",
                    i.severity == io::Severity::Error ? "[error]  " : "[warning]",
                    i.path.c_str(), i.message.c_str());
    check(rep.ok(), "validator reports no errors");
    check(m::project_to_json(loaded) == want, "file path == programmatic path (JSON-identical)");

    // (3) oracle, from the loaded project only.
    c.oracle(loaded);
}

// ------------------------------------------------- KV-CON-002: Terzaghi 1D column --
// Mirrors the physics of the integrated GUI-path Terzaghi test: nu = 0 so Eoed = E, top
// drainage only (H_dr = H), surcharge installed in the consolidation phase at t = 0+.
constexpr double kTzW = 1.0, kTzH = 12.0, kTzE = 1000.0, kTzK = 0.1, kTzQ = 10.0;

double terzaghi_U(double Tv) {
    double s = 0.0;
    for (int j = 0; j < 80; ++j) {
        const double M = (2 * j + 1) * kPi / 2.0;
        s += (2.0 / (M * M)) * std::exp(-M * M * Tv);
    }
    return 1.0 - s;
}

m::Project build_terzaghi() {
    m::Project pr;
    pr.name = "KV-CON-002 Terzaghi column";
    pr.x_min = 0.0; pr.x_max = kTzW;
    pr.y_min = 0.0; pr.y_max = kTzH;
    pr.has_water = false;   // pore = excess only (no hydrostatic background)
    pr.mesh.elem_size = 0.4;   // 0.5 * 0.4^2 = 0.08 m^2, the GUI-path test's target area
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Oedometer clay";
    s.model = m::SoilModel::LinearElastic;
    s.E = kTzE; s.nu = 0.0;
    s.gamma_unsat = 16.0; s.gamma_sat = 18.0; s.e_init = 0.5;
    s.kx = kTzK; s.ky = kTzK;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kTzW, kTzW, 0};
    P.y = {0, 0, kTzH, kTzH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed,
                   (int)m::FlowBCType::Head, (int)m::FlowBCType::Closed};   // top drains (p = 0)
    P.edge_head = {0.0, 0.0, kTzH, 0.0};
    pr.polygons.push_back(P);

    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Surcharge";
    L.x1 = 0; L.y1 = kTzH; L.x2 = kTzW; L.y2 = kTzH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kTzQ;
    pr.loads.push_back(L);

    pr.initial.load_active = {0};   // surcharge OFF in the initial K0 phase
    const double cv = kTzK * kTzE / katai::app::kGammaWater;   // nu = 0 -> Eoed = E
    m::Phase consol;
    consol.name = "Consolidation";
    consol.type = m::PhaseType::Consolidation;
    consol.duration = 2.0 * kTzH * kTzH / cv;   // Tv ~ 2: near-complete consolidation
    consol.time_steps = 120;
    consol.load_active = {1};                   // surcharge ON -> applied at t = 0+
    pr.phases.push_back(consol);
    return pr;
}

void oracle_terzaghi(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "column meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }

    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + consolidation phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    const auto& C = res[1];
    check(!C.consol_time.empty(), "consolidation produced a time series");
    if (C.consol_time.empty()) return;

    const double cv = kTzK * kTzE / katai::app::kGammaWater;
    const double s_inf = kTzQ * kTzH / kTzE;

    double pore_peak = 0.0;
    for (double p : C.consol_excess_pore) pore_peak = std::fmax(pore_peak, p);
    std::printf("      excess pore: peak = %.3f kPa (~q = %.1f), final = %.4f kPa\n",
                pore_peak, kTzQ, C.consol_excess_pore.back());
    check(pore_peak > 0.85 * kTzQ && pore_peak < 1.10 * kTzQ,
          "undrained excess pore generation ~ surcharge q");
    check(C.consol_excess_pore.back() < 0.05 * kTzQ, "excess pore dissipated by Tv ~ 2");

    std::printf("      Tv     U_FE     U_Terzaghi   err\n");
    int checked = 0;
    for (size_t i = 1; i < C.consol_time.size(); ++i) {
        const double Tv = cv * C.consol_time[i] / (kTzH * kTzH);
        const double Ufe = C.consol_settlement[i] / s_inf;
        const double Uth = terzaghi_U(Tv);
        if (std::fabs(Tv - 0.2) < 0.013 || std::fabs(Tv - 0.4) < 0.013 ||
            std::fabs(Tv - 0.6) < 0.013 || std::fabs(Tv - 0.9) < 0.013) {
            std::printf("      %.3f  %.4f   %.4f       %.1f%%\n", Tv, Ufe, Uth,
                        100.0 * (Ufe - Uth) / Uth);
            check(std::fabs(Ufe - Uth) < 0.03, "U_FE matches Terzaghi U(Tv) within 3%");
            ++checked;
        }
    }
    check(checked >= 3, "several Tv points sampled");

    const double Tvf = cv * C.consol_time.back() / (kTzH * kTzH);
    const double Uf = C.consol_settlement.back() / s_inf;
    std::printf("      final: Tv = %.2f  U_FE = %.4f  U_Terzaghi = %.4f\n", Tvf, Uf,
                terzaghi_U(Tvf));
    check(std::fabs(Uf - terzaghi_U(Tvf)) < 0.03, "final settlement matches Terzaghi U(Tv_final)");
}

// ------------------------------------- weightless elastic half-plane, shared geometry --
// 40 x 20 m domain, weightless LE soil (gamma = 0 is deliberately legal for verification
// materials), gravity-loading initial procedure: with zero self-weight the single-phase
// solve is exactly the elastic response to the applied load.
constexpr double kHpW = 40.0, kHpH = 20.0, kHpCx = 20.0, kHpQ = 100.0;

m::Project build_half_plane(const char* name) {
    m::Project pr;
    pr.name = name;
    pr.x_min = 0.0; pr.x_max = kHpW;
    pr.y_min = 0.0; pr.y_max = kHpH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = true;   // the load line/point becomes a refinement source

    m::Material s;
    s.name = "Weightless elastic";
    s.model = m::SoilModel::LinearElastic;
    s.E = 30000.0; s.nu = 0.3;         // sigma_z of both oracles is independent of E and nu
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Half-plane";
    P.material = 0;
    P.x = {0, kHpW, kHpW, 0};
    P.y = {0, 0, kHpH, kHpH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

// Solve the single-phase run exactly as a front end would run the file.
katai::app::SolveResult solve_single_phase(const m::Project& pr, katai::mesh::Mesh& mesh_out) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "meshed from the file's own settings");
    mesh_out = M.mesh;
    return katai::app::solve_gravity_le(pr, M.mesh,
                                        katai::app::initial_phase_from(pr.initial_procedure));
}

// ------------------------------------------------- KV-FND-008: uniform strip load --
constexpr double kStripA = 2.0;   // strip half-width; the strip spans x = 18..22 at the surface

m::Project build_strip() {
    m::Project pr = build_half_plane("KV-FND-008 strip load");
    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Strip q";
    L.x1 = kHpCx - kStripA; L.y1 = kHpH;
    L.x2 = kHpCx + kStripA; L.y2 = kHpH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kHpQ;
    pr.loads.push_back(L);
    return pr;
}

double strip_sigma_z(double x, double z) {   // x from the strip centre, z depth
    const double t1 = std::atan((x + kStripA) / z), t2 = std::atan((x - kStripA) / z);
    const double alpha = t1 - t2;
    return (kHpQ / kPi) * (alpha + std::sin(alpha) * std::cos(t1 + t2));
}

void oracle_strip(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "strip-load solve converged");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    std::printf("      z[m]   x[m]    sigma_z FE   closed form   err\n");
    for (double zt : {2.0, 3.0, 4.0, 6.0}) {
        const int n = nearest_node(mesh, kHpCx, kHpH - zt);
        const double x = mesh.x[n] - kHpCx, z = kHpH - mesh.y[n];
        const double fe = -R.stress.stress[n](1);          // compression positive here
        const double cf = strip_sigma_z(x, z);
        std::printf("      %.2f  %+.2f   %8.3f     %8.3f    %+.1f%%\n", z, x, fe, cf,
                    100.0 * (fe - cf) / cf);
        check(std::fabs(fe - cf) < 0.03 * cf, "sigma_z matches the strip closed form within 3%");
    }
}

// ------------------------------------------------ KV-FND-009: Flamant line load --
m::Project build_flamant() {
    m::Project pr = build_half_plane("KV-FND-009 Flamant line load");
    pr.mesh.order = 15;   // the corpus exercises tri15 through the file
    m::Load L;
    L.kind = m::LoadKind::Point;
    L.name = "Line load P";
    L.x1 = kHpCx; L.y1 = kHpH;
    L.qy1 = -kHpQ;
    pr.loads.push_back(L);
    return pr;
}

double flamant_sigma_z(double x, double z) {   // x from the load line, z depth
    const double r2 = x * x + z * z;
    return (2.0 * kHpQ / kPi) * z * z * z / (r2 * r2);
}

void oracle_flamant(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "Flamant solve converged");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    std::printf("      z[m]   x[m]    sigma_z FE   closed form   err\n");
    for (double zt : {2.0, 3.0, 4.0, 6.0}) {
        const int n = nearest_node(mesh, kHpCx, kHpH - zt);
        const double x = mesh.x[n] - kHpCx, z = kHpH - mesh.y[n];
        const double fe = -R.stress.stress[n](1);
        const double cf = flamant_sigma_z(x, z);
        std::printf("      %.2f  %+.2f   %8.3f     %8.3f    %+.1f%%\n", z, x, fe, cf,
                    100.0 * (fe - cf) / cf);
        check(std::fabs(fe - cf) < 0.03 * cf, "sigma_z matches the Flamant closed form within 3%");
    }
}

// ------------------------------------------ KV-NUM-003: submerged K0 geostatic block --
// 20 x 10 m block, water table at the surface: buoyant K0 effective stresses, hydrostatic
// pore pressure, and (the K0 procedure's defining property) no displacement.
constexpr double kK0H = 10.0, kK0W = 20.0;

m::Project build_k0_block() {
    m::Project pr;
    pr.name = "KV-NUM-003 K0 geostatic block";
    pr.x_min = 0.0; pr.x_max = kK0W;
    pr.y_min = 0.0; pr.y_max = kK0H;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Sand";
    s.model = m::SoilModel::LinearElastic;
    s.gamma_unsat = 17.0; s.gamma_sat = 20.0;
    s.E = 1.0e4; s.nu = 0.3;
    s.phi = 30.0; s.c = 1.0;   // k0_auto: K0 = 1 - sin(30) = 0.5
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Block";
    P.material = 0;
    P.x = {0, kK0W, kK0W, 0};
    P.y = {0, 0, kK0H, kK0H};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    pr.has_water = true;
    pr.wx = {0.0, kK0W};
    pr.wy = {kK0H, kK0H};   // water table at the surface
    return pr;
}

void oracle_k0_block(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "submerged K0 solve converged");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    std::printf("      max|u| = %.3e m\n", R.max_disp);
    check(R.max_disp < 1e-6, "undisturbed submerged ground does not displace");

    const double gamma_eff = 20.0 - katai::app::kGammaWater;
    double max_sv = 0.0, max_sh = 0.0, max_pore = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        const double d = kK0H - mesh.y[n];
        const double sv_ex = -gamma_eff * d;
        max_sv = std::fmax(max_sv, std::fabs(R.stress.stress[n](1) - sv_ex));
        max_sh = std::fmax(max_sh, std::fabs(R.stress.stress[n](0) - 0.5 * sv_ex));
        max_pore = std::fmax(max_pore, std::fabs(R.pore[n] - katai::app::kGammaWater * d));
    }
    std::printf("      max|sigma'_v err| = %.3e  max|sigma'_h err| = %.3e  max|pore err| = %.3e kPa\n",
                max_sv, max_sh, max_pore);
    check(max_sv < 1e-6 && max_sh < 1e-6,
          "buoyant K0 effective-stress field exact to round-off (1e-6 kPa)");
    check(max_pore < 1e-9, "hydrostatic pore pressure u_w = gamma_w (H - y)");
}

// ------------------------------------------- KV-CST-001: undrained confined column --
// Weightless 2 x 10 m column, Undrained (A), full-width surcharge: the Skempton 1D
// closed forms, with the engine's K_w/n built from nu_u = 0.495.
constexpr double kUcE = 1.0e4, kUcNu = 0.3, kUcNuU = 0.495, kUcH = 10.0, kUcW = 2.0,
                 kUcQ = 50.0;

m::Project build_undrained_column() {
    m::Project pr;
    pr.name = "KV-CST-001 undrained column";
    pr.x_min = 0.0; pr.x_max = kUcW;
    pr.y_min = 0.0; pr.y_max = kUcH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Weightless undrained clay";
    s.model = m::SoilModel::LinearElastic;
    s.drainage = m::Drainage::Undrained;
    s.E = kUcE; s.nu = kUcNu;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;   // only the surcharge acts
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kUcW, kUcW, 0};
    P.y = {0, 0, kUcH, kUcH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::NormallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::NormallyFixed};
    pr.polygons.push_back(P);

    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Surcharge";
    L.x1 = 0; L.y1 = kUcH; L.x2 = kUcW; L.y2 = kUcH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kUcQ;
    pr.loads.push_back(L);
    return pr;
}

void oracle_undrained_column(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "undrained column solve converged");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }

    const double Mp = kUcE * (1.0 - kUcNu) / ((1.0 + kUcNu) * (1.0 - 2.0 * kUcNu));
    const double Kp = kUcE / (3.0 * (1.0 - 2.0 * kUcNu));
    const double kwn = 3.0 * (kUcNuU - kUcNu) / ((1.0 - 2.0 * kUcNuU) * (1.0 + kUcNu)) * Kp;
    const double Mu = Mp + kwn;
    const double uy_ex = -kUcQ * kUcH / Mu;
    const double sig_ex = -Mp * kUcQ / Mu;

    double u_top = 0.0, sig_mid = 0.0;
    for (int n = 0; n < mesh.node_count; ++n) {
        if (mesh.y[n] > kUcH - 1e-6) u_top = std::fmin(u_top, R.disp[n * 2 + 1]);
        if (mesh.y[n] > 0.2 * kUcH && mesh.y[n] < 0.8 * kUcH)
            sig_mid = std::fmin(sig_mid, R.stress.stress[n](1));
    }
    std::printf("      u_top = %.5e m (exact %.5e, %+.1f%%)   sigma'_yy(mid) = %.3f kPa (exact %.3f, %+.1f%%)\n",
                u_top, uy_ex, 100.0 * (u_top - uy_ex) / uy_ex,
                sig_mid, sig_ex, 100.0 * (sig_mid - sig_ex) / sig_ex);
    check(std::fabs(u_top - uy_ex) < 0.02 * std::fabs(uy_ex),
          "undrained settlement = -qH/M_u within 2%");
    check(std::fabs(sig_mid - sig_ex) < 0.03 * std::fabs(sig_ex),
          "effective stress = -M'q/M_u within 3% (the rest is excess pore)");
}

// ---------------------------------------- KV-SLP-001: Griffiths and Lane slope FoS --
// The file's OWN initial procedure is Safety: loading the .k2d and running it IS the
// phi-c reduction. Geometry and parameters as published (1:2 slope on a foundation).
m::Project build_gl_slope() {
    m::Project pr;
    pr.name = "KV-SLP-001 Griffiths-Lane slope";
    pr.x_min = 20.0; pr.x_max = 70.0;
    pr.y_min = 20.0; pr.y_max = 35.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::Safety;
    pr.mesh.elem_size = 3.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Griffiths-Lane soil";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3;
    s.gamma_unsat = 20.2; s.c = 3.0; s.phi = 19.6; s.psi = 0.0;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Slope";
    P.material = 0;
    // CCW: base, right, top, slope face, foundation top, back.
    P.x = {20, 70, 70, 50, 30, 20};
    P.y = {20, 20, 35, 35, 25, 25};
    P.edge_bc = {(int)m::BCType::FullyFixed,        (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free,              (int)m::BCType::Free,
                 (int)m::BCType::Free,              (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

void oracle_gl_slope(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "Safety (phi-c reduction) ran from the file");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    const double ref = 0.99;
    std::printf("      FoS = %.3f  (published ~%.2f: Bishop 0.988 / Spencer 0.987 / T6 0.997)  err = %+.1f%%   mechanism max|u| = %.3e\n",
                R.fos, ref, 100.0 * (R.fos - ref) / ref, R.max_disp);
    check(std::fabs(R.fos - ref) < 0.08 * ref, "factor of safety within 8% of the benchmark");
    check(R.max_disp > 1e-6, "the failure mechanism displaces (a genuine slip surface)");
}

// -------------------------------------------- KV-EXC-001: staged excavation heave --
// Two stacked layers (lower 0..6 gamma 18, upper 6..10 gamma 17); the staged phase
// deactivates the upper layer and the pit floor heaves by the oedometric unloading.
constexpr double kExE = 1.0e4, kExNu = 0.3;

m::Project build_excavation() {
    m::Project pr;
    pr.name = "KV-EXC-001 staged excavation";
    pr.x_min = 0.0; pr.x_max = 20.0;
    pr.y_min = 0.0; pr.y_max = 10.0;
    pr.has_water = false;
    pr.mesh.elem_size = 1.5;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Lower stratum";
    s.model = m::SoilModel::LinearElastic;
    s.E = kExE; s.nu = kExNu;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0;
    pr.materials.push_back(s);
    m::Material f = s;
    f.name = "Upper stratum (excavated)";
    f.gamma_unsat = 17.0;
    pr.materials.push_back(f);

    m::SoilPolygon L;
    L.name = "Lower";
    L.material = 0;
    L.x = {0, 20, 20, 0}; L.y = {0, 0, 6, 6};
    L.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(L);
    m::SoilPolygon U;
    U.name = "Upper";
    U.material = 1;
    U.x = {0, 20, 20, 0}; U.y = {6, 6, 10, 10};
    U.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(U);

    m::Phase exc;
    exc.name = "Excavate";
    exc.poly_active = {1, 0};
    pr.phases.push_back(exc);
    return pr;
}

void oracle_excavation(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "two-layer model meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + excavation phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    check(res[0].max_disp < 1e-9, "initial full-geometry K0 phase does not displace");

    const double Eoed = kExE * (1.0 - kExNu) / ((1.0 + kExNu) * (1.0 - 2.0 * kExNu));
    const double uh_ex = +17.0 * 4.0 * 6.0 / Eoed;
    const auto& R = res[1];
    int best = 0;
    double bd = 1e300;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        const double d = std::hypot(R.mesh.x[n] - 10.0, R.mesh.y[n] - 6.0);
        if (d < bd) { bd = d; best = n; }
    }
    const double heave = R.disp[best * 2 + 1];
    double sv = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] < 1e-6) sv = std::fmin(sv, R.stress.stress[n](1));
    const double sv_ex = -18.0 * 6.0;
    std::printf("      heave u(y=6) = %.5e m (exact %.5e, %+.1f%%)   base sigma_v = %.2f kPa (exact %.2f, %+.1f%%)\n",
                heave, uh_ex, 100.0 * (heave - uh_ex) / uh_ex,
                sv, sv_ex, 100.0 * (sv - sv_ex) / std::fabs(sv_ex));
    check(std::fabs(heave - uh_ex) < 0.02 * uh_ex, "pit-floor heave = unloading/E_oed within 2%");
    check(std::fabs(sv - sv_ex) < 0.02 * std::fabs(sv_ex),
          "base stress sheds exactly the excavated weight within 2%");
}

// ------------------------------------------- KV-DYN-002: resonant shear column --
// The verified GUI-path site-response case (test_dynamic_gui) remodelled onto the
// schema: an elastic column in pure horizontal shear (VerticallyFixed sides), driven
// harmonically at its own fundamental frequency and, in a second phase, well below it.
// gamma/g = 2.0 exactly, so Vs = 200 m/s and f_1 = 2.5 Hz on round numbers.
constexpr double kDyW = 2.0, kDyH = 20.0, kDyE = 208000.0, kDyNu = 0.3, kDyGamma = 19.62;
constexpr double kDyA = 1.0, kDyXi = 0.05;

// Mirrored EXPRESSION FOR EXPRESSION by the Python DSL build (dsl_corpus.py): the
// file stores these derived numbers, and byte-identity across authors needs the
// same IEEE operations in the same order.
double dyn_f1() {
    const double G = kDyE / (2.0 * (1.0 + kDyNu));
    const double rho = kDyGamma / 9.81;
    return std::sqrt(G / rho) / (4.0 * kDyH);
}

m::Project build_resonant_column() {
    m::Project pr;
    pr.name = "KV-DYN-002 resonant shear column";
    pr.x_min = 0.0; pr.x_max = kDyW;
    pr.y_min = 0.0; pr.y_max = kDyH;
    pr.has_water = false;
    pr.mesh.elem_size = 0.4;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Shear column soil";
    s.model = m::SoilModel::LinearElastic;
    s.E = kDyE; s.nu = kDyNu;
    s.gamma_unsat = kDyGamma; s.gamma_sat = kDyGamma; s.e_init = 0.5;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kDyW, kDyW, 0};
    P.y = {0, 0, kDyH, kDyH};
    // Rigid base; VerticallyFixed sides (u_y = 0, u_x free) suppress cantilever
    // bending so the thin column deforms in pure shear -- the 1D SH column.
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(P);

    const double f1 = dyn_f1();
    auto dyn = [&](const char* name, double freq, double dur, int steps) {
        m::Phase p;
        p.name = name;
        p.type = m::PhaseType::Dynamic;
        p.seismic_wave = m::SeismicWave::Harmonic;
        p.seismic_amp = kDyA;
        p.seismic_freq = freq;
        p.damping_ratio = kDyXi;
        p.rayleigh_f1 = f1;
        p.rayleigh_f2 = 3.0 * f1;
        p.duration = dur;
        p.time_steps = steps;
        return p;
    };
    pr.phases.push_back(dyn("Resonance", f1, 8.0, 800));         // ~20 cycles -> steady state
    pr.phases.push_back(dyn("Off resonance", f1 / 3.0, 9.0, 720));
    return pr;
}

void oracle_resonant_column(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "column meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 3 && res[0].ok && res[1].ok && res[2].ok,
          "initial + two dynamic phases converged");
    if (res.size() != 3 || !res[1].ok || !res[2].ok) return;
    check(res[0].max_disp < 1e-9, "initial K0 phase does not displace");

    const double f1 = dyn_f1(), w1 = 2.0 * kPi * f1;
    const double u_ex = (4.0 / kPi) * kDyA / (w1 * w1 * 2.0 * kDyXi);
    const double a_ex = 2.0 / kPi / kDyXi * kDyA;
    const double u_res = res[1].max_disp, u_off = res[2].max_disp;
    const double a_res = res[1].dyn_peak_surface_a;
    std::printf("      f_1 = %.3f Hz   |u_surf| = %.4f m (exact %.4f, %+.1f%%)   "
                "|a_surf| = %.3f m/s^2 (exact %.3f, %+.1f%%)   off/res = %.3f\n",
                f1, u_res, u_ex, 100.0 * (u_res - u_ex) / u_ex,
                a_res, a_ex, 100.0 * (a_res - a_ex) / a_ex, u_off / u_res);
    check(std::fabs(u_res - u_ex) < 0.02 * u_ex,
          "resonant surface displacement = (4/pi) A/(w_1^2 2 xi) within 2%");
    check(std::fabs(a_res - a_ex) < 0.01 * a_ex,
          "peak surface acceleration = 2 A/(pi xi) within 1%");
    check(u_off < 0.2 * u_res, "the f_1/3 response is far below resonance (site selectivity)");
}

// ------------------------------------------- KV-FLW-001: Charny unconfined dam --
// The canonical free-surface benchmark (test_seepage_gui) remodelled onto the schema:
// a rectangular dam between reservoir h1 and tailwater h2, with a SEEPAGE FACE above
// the tailwater. Extra polygon vertices at (L, h2) and (0, h1) split the vertical
// faces the way a user would draw them, so each edge takes one flow condition.
// Charny's theorem fixes the discharge exactly, whatever the free surface does.
constexpr double kFwL = 10.0, kFwD = 6.0, kFwH1 = 5.0, kFwH2 = 1.0, kFwK = 0.5;

m::Project build_charny_dam() {
    m::Project pr;
    pr.name = "KV-FLW-001 Charny unconfined dam";
    pr.x_min = 0.0; pr.x_max = kFwL;
    pr.y_min = 0.0; pr.y_max = kFwD;
    pr.has_water = false;   // saturation comes from the flow solution, not a polyline
    pr.mesh.elem_size = 0.35;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Dam fill";
    s.model = m::SoilModel::LinearElastic;
    s.E = 1.0e4; s.nu = 0.3;
    s.gamma_unsat = 18.0; s.gamma_sat = 20.0; s.e_init = 0.5;
    s.kx = kFwK; s.ky = kFwK;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Dam";
    P.material = 0;
    P.x = {0, kFwL, kFwL, kFwL, 0, 0};
    P.y = {0, 0, kFwH2, kFwD, kFwD, kFwH1};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::Free, (int)m::BCType::Free};
    //             bottom                     right<h2               right>h2
    P.edge_flow = {(int)m::FlowBCType::Closed, (int)m::FlowBCType::Head, (int)m::FlowBCType::Seepage,
                   (int)m::FlowBCType::Closed, (int)m::FlowBCType::Closed, (int)m::FlowBCType::Head};
    //             top                         left>h1                    left<h1
    P.edge_head = {0.0, kFwH2, 0.0, 0.0, 0.0, kFwH1};
    pr.polygons.push_back(P);
    return pr;
}

void oracle_charny_dam(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "dam meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto F = katai::app::solve_groundwater_flow(pr, M.mesh);
    check(F.ok, "unconfined flow with a seepage face solved from the file");
    if (!F.ok) { std::printf("      (%s)\n", F.message.c_str()); return; }

    const double q_ex = kFwK * (kFwH1 * kFwH1 - kFwH2 * kFwH2) / (2.0 * kFwL);
    int n_top_sat = 0;
    for (int n = 0; n < M.mesh.node_count; ++n)
        if (M.mesh.y[n] > kFwD - 1e-6 && F.pore[n] > 1.0) ++n_top_sat;
    std::printf("      q = %.5f m3/day/m (Charny %.5f, %+.2f%%)   balance = %.2e   "
                "saturated crest nodes = %d   iters = %d\n",
                F.discharge, q_ex, 100.0 * (F.discharge - q_ex) / q_ex,
                F.balance_err, n_top_sat, F.iterations);
    check(std::fabs(F.discharge - q_ex) < 0.02 * q_ex,
          "discharge = Charny k(h1^2-h2^2)/(2L) within 2%");
    check(n_top_sat == 0, "crest stays unsaturated (free surface inside the dam)");
    check(F.balance_err < 1e-9, "global mass balance ~ 0");
}

// --------------------------------------- KV-DYN-003: El Centro, record in the file --
// The real-record feature's corpus witness (verified GUI path: test_real_record): the
// El Centro 1940 NS accelerogram travels IN the .k2d, so the checked-in file alone
// reproduces the run -- no loose side files. Two-layer profile (soft over stiff, a real
// impedance contrast) on a COMPLIANT base, exactly the product configuration.
// Both authors (this builder and the DSL build) read tests/data/elcentro-1940-ns.dat
// with the same parse and the same unit conversion, so the file's seventeen-digit
// record bytes agree across authors.
bool load_elcentro(std::vector<double>& acc, double& dt) {
    const std::string path = std::string(KATAI_TEST_DATA_DIR) + "/elcentro-1940-ns.dat";
    std::ifstream f(path);
    if (!f) return false;
    std::vector<double> t;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ss(line);
        double tv, av;
        if (ss >> tv >> av) { t.push_back(tv); acc.push_back(av * 9.81); }
    }
    if (t.size() < 2) return false;
    dt = t[1] - t[0];
    return true;
}

m::Project build_el_centro() {
    m::Project pr;
    pr.name = "KV-DYN-003 El Centro two-layer column";
    pr.x_min = 0.0; pr.x_max = 2.0;
    pr.y_min = 0.0; pr.y_max = 20.0;
    pr.has_water = false;
    pr.mesh.elem_size = 1.0;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material soft;
    soft.name = "Soft upper layer";
    soft.model = m::SoilModel::LinearElastic;
    soft.E = 2.0 * 1.3 * 25920.0; soft.nu = 0.3;
    soft.gamma_unsat = 1.8 * 9.81; soft.gamma_sat = soft.gamma_unsat; soft.e_init = 0.5;
    m::Material stiff = soft;
    stiff.name = "Stiff lower layer";
    stiff.E = 2.0 * 1.3 * 189000.0;
    stiff.gamma_unsat = 2.1 * 9.81; stiff.gamma_sat = stiff.gamma_unsat;
    pr.materials.push_back(soft);
    pr.materials.push_back(stiff);

    m::SoilPolygon Pt;
    Pt.name = "Upper";
    Pt.material = 0;
    Pt.x = {0, 2, 2, 0}; Pt.y = {12, 12, 20, 20};
    Pt.edge_bc = {(int)m::BCType::Free, (int)m::BCType::VerticallyFixed,
                  (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    m::SoilPolygon Pb;
    Pb.name = "Lower";
    Pb.material = 1;
    Pb.x = {0, 2, 2, 0}; Pb.y = {0, 0, 12, 12};
    Pb.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::VerticallyFixed,
                  (int)m::BCType::Free, (int)m::BCType::VerticallyFixed};
    pr.polygons.push_back(Pt);
    pr.polygons.push_back(Pb);

    std::vector<double> acc; double dt = 0.0;
    if (!load_elcentro(acc, dt)) return pr;   // byte-identity will fail loudly
    m::Phase p;
    p.name = "ElCentro";
    p.type = m::PhaseType::Dynamic;
    p.seismic_wave = m::SeismicWave::Record;
    p.accel_record = acc;
    p.record_dt = dt;
    p.seismic_amp = 1.0;
    p.damping_ratio = 0.05;
    p.rayleigh_f1 = 3.0;
    p.rayleigh_f2 = 9.0;
    p.duration = (double)(acc.size() - 1) * dt;
    p.time_steps = std::min((int)acc.size() - 1, 20000);
    p.seismic_compliant_base = true;
    pr.phases.push_back(p);
    return pr;
}

void oracle_el_centro(const m::Project& pr) {
    check(pr.phases.size() == 1 && !pr.phases[0].accel_record.empty(),
          "the accelerogram travelled IN the file");
    if (pr.phases.empty() || pr.phases[0].accel_record.empty()) return;
    const auto& rec = pr.phases[0].accel_record;
    const double dt = pr.phases[0].record_dt;

    // (1) The file-borne record IS El Centro: published PGA and timing.
    check(rec.size() == 1560 && std::fabs(dt - 0.02) < 1e-12,
          "1560 samples at dt = 0.02 s, from the file");
    double pga = 0.0, tpk = 0.0;
    for (size_t i = 0; i < rec.size(); ++i)
        if (std::fabs(rec[i]) > pga) { pga = std::fabs(rec[i]); tpk = i * dt; }
    std::printf("      file record: PGA = %.5f g @ t = %.2f s (published 0.319 g @ ~2 s)\n",
                pga / 9.81, tpk);
    check(std::fabs(pga / 9.81 - 0.319) < 0.005, "file-borne PGA = published 0.319 g (< 0.005 g)");
    check(tpk > 1.5 && tpk < 2.5, "peak at ~2 s, where El Centro's is");

    // (2) The 5%-damped spectrum of the file's record has the published shape.
    std::vector<double> periods;
    for (int i = 0; i <= 60; ++i) periods.push_back(0.02 + (3.5 - 0.02) * i / 60.0);
    const auto Sa = katai::core::response_spectrum(rec, dt, periods, 0.05);
    double sa_pk = 0.0, t_pk = 0.0, sa3 = 0.0;
    for (size_t i = 0; i < Sa.size(); ++i) {
        if (Sa[i] > sa_pk) { sa_pk = Sa[i]; t_pk = periods[i]; }
        if (std::fabs(periods[i] - 3.0) < 0.04) sa3 = Sa[i];
    }
    std::printf("      file spectrum: peak Sa = %.2fx PGA at T = %.2f s;  Sa(3 s) = %.3f g\n",
                sa_pk / pga, t_pk, sa3 / 9.81);
    check(sa_pk / pga > 2.0 && sa_pk / pga < 3.5, "peak amplification in the published 2.0-3.5x band");
    check(t_pk > 0.1 && t_pk < 1.0, "spectral peak in the short-to-mid period range");
    check(sa3 / 9.81 < 0.15, "long-period ordinate small (Sa(3 s) < 0.15 g)");

    // (3) The compliant-base product run, from the file alone.
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "two-layer profile meshed from the file's own settings");
    if (!M.ok) return;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + 1559-step real-record phase converged");
    if (res.size() != 2 || !res[1].ok) return;
    const auto& R = res[1];
    const double a_surf = R.dyn_peak_surface_a;
    check(!R.dyn_response_sa.empty() && R.dyn_response_sa.size() == R.dyn_period.size(),
          "surface response spectrum produced");
    if (R.dyn_response_sa.empty()) return;
    std::printf("      run: peak surface accel = %.3f m/s^2 (%.2f g)   Sa(T_min = %.3f s) = %.3f m/s^2\n",
                a_surf, a_surf / 9.81, R.dyn_period.front(), R.dyn_response_sa.front());
    // Reproduction pin, declared as such: the number the verified GUI path produced
    // (test_real_record (d)), now required from the file alone; the physics direction
    // rigid >= compliant is pinned there, not re-argued here.
    check(std::fabs(a_surf - 5.149) < 0.05 * 5.149,
          "surface response reproduces the verified path's number (5% reproduction band)");
    // The shortest tabulated period is 0.05 s, not T -> 0: for a broadband surface
    // motion Sa there sits near but not exactly at PGA (measured 1.05x).
    const double ratio = R.dyn_response_sa.front() / a_surf;
    check(ratio > 0.9 && ratio < 1.3, "short-period spectral ordinate is PGA-scale (0.9-1.3x)");
}

// ------------------------------------ KV-FND-010: Prandtl strip footing, from a file --
// The first GLOBAL LIMIT LOAD corpus case, remodelled from the direct-FE benchmark
// (test_prandtl, KV-FND-005): a flexible strip footing on a weightless, cohesive,
// frictionless (Tresca) half-space, LOADED PAST COLLAPSE in one staged phase. The
// phase's honest non-convergence is the oracle: the driver reports the equilibrated
// fraction, and load_factor * q_applied is the incremental limit load.
constexpr double kPrW = 6.0, kPrH = 4.0, kPrX0 = 2.4, kPrX1 = 3.6;   // footing B = 1.2
constexpr double kPrC = 10.0, kPrQ = 6.0 * kPrC;                     // past N_c ~ 5.14

m::Project build_prandtl_footing() {
    m::Project pr;
    pr.name = "KV-FND-010 Prandtl strip footing";
    pr.x_min = 0.0; pr.x_max = kPrW;
    pr.y_min = 0.0; pr.y_max = kPrH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;   // weightless -> nil
    pr.mesh.elem_size = 0.4;
    pr.mesh.order = 15;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Weightless Tresca clay";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 10000.0; s.nu = 0.3;
    s.c = kPrC; s.phi = 0.0; s.psi = 0.0;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = 0.5;
    s.tension_cutoff = false;   // the wedge solution has no cut-off (mirrors KV-FND-005)
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Half-space";
    P.material = 0;
    P.x = {0, kPrW, kPrW, 0};
    P.y = {0, 0, kPrH, kPrH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Footing pressure";
    L.x1 = kPrX0; L.y1 = kPrH; L.x2 = kPrX1; L.y2 = kPrH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kPrQ;
    pr.loads.push_back(L);

    pr.initial.load_active = {0};   // footing OFF in the initial phase
    m::Phase collapse;
    collapse.name = "Load to collapse";
    collapse.load_active = {1};
    pr.phases.push_back(collapse);
    return pr;
}

void oracle_prandtl_footing(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "half-space meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok, "initial phase converged");
    if (res.size() != 2) return;
    check(res[0].max_disp < 1e-9, "weightless gravity initial does not displace");

    // The RESULT is the honest non-convergence: the load exceeds the capacity, the
    // driver says so, and the equilibrated fraction is the limit load.
    check(!res[1].ok, "the collapse phase honestly does NOT fully converge");
    check(res[1].message.find("collapse") != std::string::npos,
          "the refusal message names the collapse mechanism");
    const double nc = res[1].load_factor * kPrQ / kPrC;
    const double nc_ex = 2.0 + kPi;
    std::printf("      N_c = %.3f (Prandtl 2+pi = %.3f, %+.1f%%)\n",
                nc, nc_ex, 100.0 * (nc - nc_ex) / nc_ex);
    check(std::fabs(nc - nc_ex) < 0.02 * nc_ex, "N_c = 2 + pi within 2%");
}

// --------------------------------- KV-FND-011: Gibson strip load, PLAXIS Validation 2.2 --
// The first PLAXIS-Validation remodel onto the schema (the direct-FE original is
// KV-FND-002): a strip load on an incompressible Gibson soil whose stiffness grows
// linearly from ~zero at the surface -- E(y) expressed through the schema's own
// E_inc / y_ref profile fields, the same modelling decision as PLAXIS's Advanced
// E-increment. Half-model: the left edge is the symmetry axis (x fixed), exactly
// the rectangle defaults. The three sibling cases of the quartet (2.1, 3.1, 3.2)
// are DISPLACEMENT-controlled rigid footings; the input contract carries no
// prescribed-displacement boundary yet, so they stay with the direct benchmark
// until that schema feature lands -- a remodel that changed the problem to fit
// the file would no longer be the published case.
constexpr double kGbQ = 10.0, kGbB = 1.0, kGbH = 4.0;   // q [kPa], loaded half-width, layer

m::Project build_gibson() {
    m::Project pr;
    pr.name = "KV-FND-011 Gibson strip load";
    pr.x_min = 0.0; pr.x_max = 7.0;
    pr.y_min = 0.0; pr.y_max = kGbH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;   // weightless -> nil
    pr.mesh.elem_size = 0.15;
    pr.mesh.order = 15;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Gibson soil";
    s.model = m::SoilModel::LinearElastic;
    s.E = 0.01; s.nu = 0.495;          // ~zero surface stiffness, near-incompressible
    s.E_inc = 299.0; s.y_ref = kGbH;   // E(y) = E + E_inc (y_ref - y) = 299 z
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = 0.5;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Layer";
    P.material = 0;
    P.x = {0, 7, 7, 0};
    P.y = {0, 0, kGbH, kGbH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::Load L;
    L.kind = m::LoadKind::Distributed;
    L.name = "Strip q";
    L.x1 = 0.0; L.y1 = kGbH; L.x2 = kGbB; L.y2 = kGbH;
    L.qx1 = L.qx2 = 0; L.qy1 = L.qy2 = -kGbQ;
    pr.loads.push_back(L);

    pr.initial.load_active = {0};
    m::Phase ph;
    ph.name = "Strip load";
    ph.load_active = {1};
    pr.phases.push_back(ph);
    return pr;
}

void oracle_gibson(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "Gibson layer meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + strip-load phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    check(res[0].max_disp < 1e-9, "weightless gravity initial does not displace");

    const auto& R = res[1];
    const int nc = nearest_node(R.mesh, 0.0, kGbH);   // footing centreline, surface
    const double s_katai = -R.disp[nc * 2 + 1];
    const double s_plaxis = 0.047, s_exact = 0.050;
    std::printf("      settlement: exact (half-space) %.4f | PLAXIS (4 m layer) %.4f | "
                "file run %.4f (%+.1f%% vs PLAXIS)\n",
                s_exact, s_plaxis, s_katai, 100.0 * (s_katai - s_plaxis) / s_plaxis);
    check(std::fabs(s_katai - s_plaxis) < 0.05 * s_plaxis,
          "settlement within 5% of the published PLAXIS finite-layer value");
    check(s_katai < s_exact, "the finite layer settles less than the half-space (shared bias)");
}

// ------------------------------ KV-FND-012: Giroud rigid footing, PLAXIS Validation 2.1 --
// The first DISPLACEMENT-CONTROLLED corpus case: schema v2's line prescribed
// displacement imposes the smooth rigid footing (u_y = -10 mm, u_x free) and the new
// reaction output reads the footing force back. Half-model; the left edge is the
// symmetry axis, exactly the rectangle defaults. E is derived from G expression for
// expression in the DSL build (byte identity across authors).
constexpr double kGrG = 500.0, kGrNu = 1.0 / 3.0, kGrS = 0.010, kGrB = 1.0;

m::Project build_giroud() {
    m::Project pr;
    pr.name = "KV-FND-012 Giroud rigid footing";
    pr.x_min = 0.0; pr.x_max = 7.0;
    pr.y_min = 0.0; pr.y_max = 4.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;   // weightless -> nil
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 15;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Elastic soil";
    s.model = m::SoilModel::LinearElastic;
    s.E = 2.0 * kGrG * (1.0 + kGrNu); s.nu = kGrNu;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0; s.e_init = 0.5;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Half-space";
    P.material = 0;
    P.x = {0, 7, 7, 0};
    P.y = {0, 0, 4, 4};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::PrescribedDisp D;
    D.name = "Footing";
    D.x1 = 0.0; D.y1 = 4.0; D.x2 = kGrB; D.y2 = 4.0;
    D.set_ux = false; D.ux = 0.0;      // smooth: u_x stays free under the footing
    D.set_uy = true;  D.uy = -kGrS;
    pr.disps.push_back(D);

    pr.initial.disp_active = {0};      // footing not yet pushed in the initial phase
    m::Phase indent;
    indent.name = "Indent";
    indent.disp_active = {1};
    pr.phases.push_back(indent);
    return pr;
}

void oracle_giroud(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "half-space meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + indentation phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    check(res[0].max_disp < 1e-9, "weightless gravity initial does not displace");

    // The footing force from the REACTION output: sum R_y over the nodes on the line,
    // doubled for the half-model (the direct benchmark integrates B^T sigma; the
    // reaction field is that same discrete internal force at the fixed dofs).
    const auto& R = res[1];
    check(R.reaction.size() == 2 * R.mesh.node_count, "the static phase reports reactions");
    if (R.reaction.size() != 2 * R.mesh.node_count) return;
    double ry = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-6 && R.mesh.x[n] <= kGrB + 1e-6) ry += R.reaction[2 * n + 1];
    const double F_katai = 2.0 * std::fabs(ry);
    const double F_exact = 15.15, F_plaxis = 15.24;
    std::printf("      F: analytic %.2f | PLAXIS %.2f | file run %.2f (%+.1f%% vs analytic, "
                "%+.1f%% vs PLAXIS)\n", F_exact, F_plaxis, F_katai,
                100.0 * (F_katai - F_exact) / F_exact, 100.0 * (F_katai - F_plaxis) / F_plaxis);
    check(std::fabs(F_katai - F_exact) < 0.02 * F_exact,
          "footing force within 2% of the Giroud analytic value");
    check(std::fabs(F_katai - F_plaxis) < 0.03 * F_plaxis,
          "footing force within 3% of the published PLAXIS number");
    // The imposed settlement really happened, bit-for-bit: a node ON the line carries
    // exactly u_y = -10 mm (the ramp completes at load factor 1, so the prescribed value
    // is written verbatim). max_disp is NOT the right pin -- u_x is free under a smooth
    // footing, so edge nodes move laterally too and the vector norm exceeds 10 mm.
    int on_line = -1;
    for (int n = 0; n < R.mesh.node_count && on_line < 0; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-9 && std::fabs(R.mesh.x[n] - 0.5) < 0.26) on_line = n;
    check(on_line >= 0 && std::fabs(R.disp[2 * on_line + 1] + kGrS) < 1e-12,
          "a footing node carries exactly the imposed u_y = -10 mm");
}

// ------------------------------ KV-FND-013: Cox circular footing, PLAXIS Validation 3.1 --
// The first AXISYMMETRIC corpus case: the same displacement-controlled machine as
// KV-FND-012 (line prescribed displacement + reaction output) in r-z kinematics. The
// left edge is the symmetry axis (r = 0). ASSOCIATED flow (psi = phi): Cox (1962) is a
// slip-line solution, i.e. the associated limit load -- comparing it against a
// non-associated run would mix a modelling difference into a verification number. The
// K0 initial (gamma = 16, K0 = 1 - sin phi = 0.5) seeds the geostatic state; the Indent
// phase pushes the smooth rigid footing to 0.35 m -- the soft soil (E = 2400) needs a
// large indentation to reach the collapse plateau.
constexpr double kCoxC = 1.6, kCoxR = 1.0, kCoxSettle = 0.35;

m::Project build_cox() {
    m::Project pr;
    pr.name = "KV-FND-013 Cox circular footing";
    pr.axisymmetric = true;
    pr.x_min = 0.0; pr.x_max = 5.0;
    pr.y_min = 0.0; pr.y_max = 4.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 0.25;   // 0.5 m puts two elements across the radius: +9% (recorded in KV-FND-003)
    pr.mesh.order = 15;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Cox soil";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 2400.0; s.nu = 0.20;
    s.c = kCoxC; s.phi = 30.0; s.psi = 30.0;   // associated (the slip-line assumption)
    s.gamma_unsat = 16.0; s.gamma_sat = 16.0;
    s.tension_cutoff = false;                  // plain Mohr-Coulomb, as in the slip-line solution
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Soil cylinder";
    P.material = 0;
    P.x = {0, 5, 5, 0};
    P.y = {0, 0, 4, 4};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::PrescribedDisp D;
    D.name = "Footing";
    D.x1 = 0.0; D.y1 = 4.0; D.x2 = kCoxR; D.y2 = 4.0;
    D.set_ux = false; D.ux = 0.0;      // smooth: u_x stays free under the footing
    D.set_uy = true;  D.uy = -kCoxSettle;
    pr.disps.push_back(D);

    pr.initial.disp_active = {0};      // footing not yet pushed in the K0 phase
    m::Phase indent;
    indent.name = "Indent";
    indent.disp_active = {1};
    pr.phases.push_back(indent);
    return pr;
}

void oracle_cox(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "axisymmetric soil cylinder meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "K0 initial + indentation phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    check(res[0].max_disp < 1e-9, "K0 geostatic initial does not displace");

    // The limit pressure from the REACTION output. Axisymmetric nodal forces are
    // per radian (the r-weighted assembly), so the footing force is 2 pi |Ry| and
    // p = 2 pi |Ry| / (pi R^2) = 2 |Ry| / R^2 -- the same conversion as the direct
    // benchmark KV-FND-003.
    const auto& R = res[1];
    check(R.reaction.size() == 2 * R.mesh.node_count, "the static phase reports reactions");
    if (R.reaction.size() != 2 * R.mesh.node_count) return;
    double ry = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-6 && R.mesh.x[n] <= kCoxR + 1e-6) ry += R.reaction[2 * n + 1];
    const double p_katai = 2.0 * std::fabs(ry) / (kCoxR * kCoxR);
    const double p_exact = 225.6, p_plaxis = 220.0;
    std::printf("      p_max: Cox %.1f | PLAXIS %.1f | file run %.1f (%+.1f%% vs Cox, "
                "%+.1f%% vs PLAXIS)\n", p_exact, p_plaxis, p_katai,
                100.0 * (p_katai - p_exact) / p_exact, 100.0 * (p_katai - p_plaxis) / p_plaxis);
    check(std::fabs(p_katai - p_exact) < 0.05 * p_exact,
          "limit pressure within 5% of the Cox exact collapse pressure");
    int on_line = -1;
    for (int n = 0; n < R.mesh.node_count && on_line < 0; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-9 && std::fabs(R.mesh.x[n] - 0.5) < 0.13) on_line = n;
    check(on_line >= 0 && std::fabs(R.disp[2 * on_line + 1] + kCoxSettle) < 1e-12,
          "a footing node carries exactly the imposed u_y = -0.35 m");
}

// ------------------ KV-FND-014: Davis & Booker c(z) strip footing, PLAXIS Validation 3.2 --
// Tresca (phi = 0) with c = 1 + 2z and E = 299 + 498z through the schema's
// c_inc / E_inc / y_ref profile -- the corpus twin of the direct benchmark KV-FND-004,
// driven by the same displacement-controlled machine as KV-FND-012. Weightless soil:
// the gravity initial is an exact nil.
constexpr double kDbBhalf = 1.0, kDbSettle = 0.03;

m::Project build_davis_booker() {
    m::Project pr;
    pr.name = "KV-FND-014 Davis-Booker strip footing";
    pr.x_min = 0.0; pr.x_max = 5.0;
    pr.y_min = 0.0; pr.y_max = 4.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;   // weightless -> nil
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 15;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Davis-Booker clay";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 299.0; s.nu = 0.3;
    s.c = 1.0; s.phi = 0.0; s.psi = 0.0;       // Tresca
    s.E_inc = 498.0; s.c_inc = 2.0; s.y_ref = 4.0;   // + per metre BELOW y_ref (the surface)
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;
    s.tension_cutoff = false;                  // plain Tresca, as in the analytic solution
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Clay layer";
    P.material = 0;
    P.x = {0, 5, 5, 0};
    P.y = {0, 0, 4, 4};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::PrescribedDisp D;
    D.name = "Footing";
    D.x1 = 0.0; D.y1 = 4.0; D.x2 = kDbBhalf; D.y2 = 4.0;
    D.set_ux = false; D.ux = 0.0;      // smooth: u_x stays free under the footing
    D.set_uy = true;  D.uy = -kDbSettle;
    pr.disps.push_back(D);

    pr.initial.disp_active = {0};      // footing not yet pushed in the initial phase
    m::Phase indent;
    indent.name = "Indent";
    indent.disp_active = {1};
    pr.phases.push_back(indent);
    return pr;
}

void oracle_davis_booker(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "clay layer meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + indentation phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    check(res[0].max_disp < 1e-9, "weightless gravity initial does not displace");

    // Average pressure under the footing from the REACTION output: p = |sum Ry| / (B/2).
    const auto& R = res[1];
    check(R.reaction.size() == 2 * R.mesh.node_count, "the static phase reports reactions");
    if (R.reaction.size() != 2 * R.mesh.node_count) return;
    double ry = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-6 && R.mesh.x[n] <= kDbBhalf + 1e-6) ry += R.reaction[2 * n + 1];
    const double p_katai = std::fabs(ry) / kDbBhalf;
    const double p_exact = 7.80, p_plaxis = 7.86;
    std::printf("      p_max: Davis-Booker %.2f | PLAXIS %.2f | file run %.2f (%+.1f%% vs "
                "analytic, %+.1f%% vs PLAXIS)\n", p_exact, p_plaxis, p_katai,
                100.0 * (p_katai - p_exact) / p_exact, 100.0 * (p_katai - p_plaxis) / p_plaxis);
    check(std::fabs(p_katai - p_exact) < 0.05 * p_exact,
          "limit pressure within 5% of the Davis-Booker exact collapse pressure");
    check(std::fabs(p_katai - p_plaxis) < 0.03 * p_plaxis,
          "limit pressure within 3% of the published PLAXIS number");
    int on_line = -1;
    for (int n = 0; n < R.mesh.node_count && on_line < 0; ++n)
        if (R.mesh.y[n] > 4.0 - 1e-9 && std::fabs(R.mesh.x[n] - 0.5) < 0.26) on_line = n;
    check(on_line >= 0 && std::fabs(R.disp[2 * on_line + 1] + kDbSettle) < 1e-12,
          "a footing node carries exactly the imposed u_y = -30 mm");
}

// ------------------------------ KV-SLP-002: Griffiths and Lane (1999) Example 1 --
// The paper's own first example (D = 1, no foundation layer), dimensionalised as
// H = 10 m, gamma = 20 kN/m3, c' = 10 kPa so that c'/gamma H = 0.05 exactly. Geometry
// from their Fig. 1: a 1.2H crest plateau, a 2H slope run, vertical rollers on the
// left boundary, full fixity at the base. psi = 0 and the nominal E' = 1e5 / nu' = 0.3
// are the paper's stated values; no tension crack modelling (plain Mohr-Coulomb), as
// in the paper.
m::Project build_gl_example1() {
    m::Project pr;
    pr.name = "KV-SLP-002 Griffiths-Lane Example 1";
    pr.x_min = 0.0; pr.x_max = 32.0;
    pr.y_min = 0.0; pr.y_max = 10.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::Safety;
    pr.mesh.elem_size = 1.0;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Example 1 soil";
    s.model = m::SoilModel::MohrCoulomb;
    s.E = 1.0e5; s.nu = 0.3;
    s.gamma_unsat = 20.0; s.gamma_sat = 20.0;
    s.c = 10.0; s.phi = 20.0; s.psi = 0.0;
    s.tension_cutoff = false;
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Slope";
    P.material = 0;
    // CCW: base, slope face (toe at x = 3.2H), crest plateau, left boundary.
    P.x = {0, 32, 12, 0};
    P.y = {0, 0, 10, 10};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free,       (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);
    return pr;
}

void oracle_gl_example1(const m::Project& pr) {
    katai::mesh::Mesh mesh;
    const auto R = solve_single_phase(pr, mesh);
    check(R.ok, "Safety (phi-c reduction) ran from the file");
    if (!R.ok) { std::printf("      (%s)\n", R.message.c_str()); return; }
    const double ref = 1.380;   // Bishop & Morgenstern (1960), the chart value the paper cites
    std::printf("      FoS = %.3f  (published: Bishop-Morgenstern 1.380, Griffiths-Lane FE 1.4)  "
                "err vs 1.380 = %+.1f%%   mechanism max|u| = %.3e\n",
                R.fos, 100.0 * (R.fos - ref) / ref, R.max_disp);
    check(std::fabs(R.fos - ref) < 0.04 * ref,
          "factor of safety within 4% of the Bishop-Morgenstern chart value");
    check(R.max_disp > 1e-6, "the failure mechanism displaces (a genuine slip surface)");
}


// ------------------------------------------- KV-CST-002: Hardening Soil oedometer (1D) --
// The first HARDENING SOIL boundary-value problem in the corpus. Until now the HS family was
// verified only at the material point (test_hardening_soil, test_hs_cap, test_hs_berlin against
// the Material Models Manual figures) -- the model was proven, the PATH from a .k2d file through
// the mesher, the assembler, the cap return mapping and the load stepping was not.
//
// A one-dimensional compression test is the right first case because it is the experiment that
// DEFINES the parameter under test: E_oed^ref is read off an oedometer, and the closed-form
// integral of the HS stiffness law is what a practitioner assumes when they type that number in.
// The column is weightless, so the vertical stress is the surcharge itself and uniform -- there
// is no depth integral to argue about, and what is compared is the model's own law against the
// finite element solution of it: an independent path through the whole machinery.
constexpr double kOedW = 1.0, kOedH = 4.0;        // column [m]
constexpr double kOedQ0 = 50.0, kOedQ1 = 200.0;   // vertical stress before / after [kPa]
constexpr double kOedEoedRef = 30000.0, kOedPref = 100.0, kOedM = 0.5;

m::Project build_hs_oedometer() {
    m::Project pr;
    pr.name = "KV-CST-002 HS oedometer";
    pr.x_min = 0.0; pr.x_max = kOedW;
    pr.y_min = 0.0; pr.y_max = kOedH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 0.5;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Hardening Soil sand";
    s.model = m::SoilModel::HardeningSoil;
    s.drainage = m::Drainage::Drained;
    s.gamma_unsat = 0.0; s.gamma_sat = 0.0;   // weightless: sigma_1 IS the surcharge
    s.E = 30000.0; s.nu = 0.2;                // fallback only; HS uses the reference moduli
    s.c = 0.0; s.phi = 35.0; s.psi = 5.0;
    s.E50ref = 30000.0; s.Eoedref = kOedEoedRef; s.Eurref = 90000.0;
    s.m = kOedM; s.p_ref = kOedPref; s.nu_ur = 0.2; s.Rf = 0.9;
    s.k0nc_auto = true;
    s.tension_cutoff = false;   // 1D compression never reaches it; keeps the run's report clean
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0, kOedW, kOedW, 0};
    P.y = {0, 0, kOedH, kOedH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    // Two surcharges on the same line: the seating stress, active from the start, and the
    // increment the second phase switches on. The schema activates loads per phase, it does not
    // rescale them, so a stress step is two loads rather than one changed number.
    m::Load L0;
    L0.kind = m::LoadKind::Distributed;
    L0.name = "Seating stress";
    L0.x1 = 0; L0.y1 = kOedH; L0.x2 = kOedW; L0.y2 = kOedH;
    L0.qx1 = L0.qx2 = 0; L0.qy1 = L0.qy2 = -kOedQ0;
    pr.loads.push_back(L0);
    m::Load L1 = L0;
    L1.name = "Load increment";
    L1.qy1 = L1.qy2 = -(kOedQ1 - kOedQ0);
    pr.loads.push_back(L1);

    pr.initial.load_active = {1, 0};   // seating only
    m::Phase step;
    step.name = "Load to 200 kPa";
    step.load_active = {1, 1};
    pr.phases.push_back(step);
    return pr;
}

void oracle_hs_oedometer(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "oedometer column meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "seating + loading phases converged");
    if (res.size() != 2 || !res[1].ok) {
        if (!res.empty()) std::printf("      (%s)\n", res.back().message.c_str());
        return;
    }

    // The closed form, written out here rather than called from the material header: with c = 0
    // the stiffness factor is (sigma/p_ref)^m (sin(phi) cancels), so
    //   -eps_1 = (p_ref^m / Eoed_ref) * [sigma^(1-m)]/(1-m), evaluated from sigma_0 to sigma_1.
    const double coef = std::pow(kOedPref, kOedM) / kOedEoedRef;
    const double eps = coef * (std::pow(kOedQ1, 1.0 - kOedM) - std::pow(kOedQ0, 1.0 - kOedM)) /
                       (1.0 - kOedM);
    const double want = eps * kOedH;   // settlement increment of the top [m]

    // A phase reports displacement relative to ITS OWN start, so the loading phase's field is
    // already the increment the closed form describes. Measured rather than assumed: the seating
    // phase settles 0.0164 m under 0 -> 50 kPa, and an accumulating phase would then report
    // about 0.038 m (the 0 -> 200 kPa integral) instead of the 0.019 m it does.
    const int top = nearest_node(res[1].mesh, 0.5 * kOedW, kOedH);
    const double got = -res[1].disp[top * 2 + 1];
    std::printf("      seating phase (0 -> 50 kPa) settles %.6f m; the loading phase reports its "
                "own increment\n", -res[0].disp[top * 2 + 1]);
    std::printf("      settlement increment %.6f m vs closed form %.6f m (%+.2f%%)\n", got, want,
                100.0 * (got - want) / want);
    check(std::fabs(got - want) / want < 0.03,
          "1D compression settlement matches the HS oedometric law within 3%");

    // The stress state has to be the one the law was integrated over, or the agreement above is
    // a coincidence: uniform vertical stress equal to the surcharge, top to bottom.
    double sv_min = 0.0, sv_max = -1e300;
    for (int n = 0; n < res[1].mesh.node_count; ++n) {
        sv_min = std::fmin(sv_min, res[1].stress.stress[n](1));
        sv_max = std::fmax(sv_max, res[1].stress.stress[n](1));
    }
    std::printf("      vertical stress range %.3f .. %.3f kPa (surcharge %.1f)\n", sv_min, sv_max,
                -kOedQ1);
    check(std::fabs(sv_min + kOedQ1) / kOedQ1 < 0.02 && std::fabs(sv_max + kOedQ1) / kOedQ1 < 0.02,
          "the column carries the surcharge as a uniform vertical stress");
}

}  // namespace

// ------------------------------ KV-STR-002: sliding block, PLAXIS Validation 3.3 ---------
// The first case in which an INTERFACE decides the answer, and it is the manual's own
// interface verification. A stiff block is pushed sideways until it slides on a Coulomb
// joint at its base; the failure force is the joint's capacity and nothing else.
//
// Width 4 m and "weight 100 kN/m" at gamma = 25 fix the block at 4 m x 1 m -- the manual
// prints the arithmetic, not the height. The interface lives in a SEPARATE data set, as the
// manual specifies, reached through the schema's iface_material.
//
// Why the case is here at all: it did not run until 2026-08-10. The interface lies ALONG a
// fixed boundary, whose two split sides sit at identical coordinates, and the boundary
// conditions -- which match by coordinate -- fixed both of them. The joint was assembled and
// then welded shut, in silence, and this file read 5.4e6 kN/m instead of 60. The last check
// below is the sentry for that: with no interface the answer is not slightly different, it
// is wrong by five orders of magnitude.
constexpr double kSbW = 4.0, kSbH = 1.0, kSbGamma = 25.0;
constexpr double kSbCw = 2.5, kSbPhiw = 26.6, kSbPush = 0.1;

// The manual's own formula, evaluated for whatever (c_w, gamma) it is handed.
double sliding_capacity(double cw, double gamma) {
    return kSbW * cw + (kSbW * kSbH * gamma) * std::tan(kSbPhiw * kPi / 180.0);
}

m::Project build_sliding_block_at(double cw, double gamma, double push) {
    m::Project pr;
    pr.name = "KV-STR-002 PLAXIS 3.3 sliding block";
    pr.x_min = 0.0; pr.x_max = kSbW;
    pr.y_min = 0.0; pr.y_max = kSbH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 0.25;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material block;                       // "The block is modelled as a stiff linear elastic material"
    block.name = "Concrete block";
    block.model = m::SoilModel::LinearElastic;
    block.E = 3.0e7; block.nu = 0.0;         // 30 GN/m2
    block.gamma_unsat = gamma; block.gamma_sat = gamma;
    block.k0_auto = false; block.k0 = 0.0;   // "self weight stresses are switched on using K0 = 0"
    pr.materials.push_back(block);

    m::Material joint;                       // "stored in a separate elastoplastic data set"
    joint.name = "Interface";
    joint.model = m::SoilModel::MohrCoulomb;
    joint.E = 3.0e6; joint.nu = 0.45;        // 3 GN/m2
    joint.c = cw; joint.phi = kSbPhiw; joint.psi = 0.0;
    joint.gamma_unsat = gamma; joint.gamma_sat = gamma;
    joint.k0_auto = false; joint.k0 = 0.0;
    joint.rinter_rigid = true;               // R_inter = 1: the data set IS the interface strength
    pr.materials.push_back(joint);

    m::SoilPolygon P;
    P.name = "Block";
    P.material = 0;
    P.x = {0, kSbW, kSbW, 0};
    P.y = {0, 0, kSbH, kSbH};
    // "Bottom nodes are fully fixed. All other nodes are entirely free."
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::Free,
                 (int)m::BCType::Free, (int)m::BCType::Free};
    pr.polygons.push_back(P);

    m::StructElement iface;
    iface.name = "Base interface";
    iface.kind = m::StructKind::Interface;
    iface.x1 = 0.0; iface.y1 = 0.0; iface.x2 = kSbW; iface.y2 = 0.0;
    iface.iface_material = 1;
    pr.structs.push_back(iface);

    m::PrescribedDisp D;                     // "a prescribed horizontal displacement of 0.1 m ...
    D.name = "Push";                         //  but the nodes at this side are free to move vertically"
    D.x1 = 0.0; D.y1 = 0.0; D.x2 = 0.0; D.y2 = kSbH;
    D.set_ux = true;  D.ux = push;
    D.set_uy = false; D.uy = 0.0;
    pr.disps.push_back(D);

    pr.initial.struct_active = {1};
    pr.initial.disp_active = {0};
    m::Phase shove;
    shove.name = "Push";
    shove.struct_active = {1};
    shove.disp_active = {1};
    pr.phases.push_back(shove);
    return pr;
}

m::Project build_sliding_block() { return build_sliding_block_at(kSbCw, kSbGamma, kSbPush); }

// The horizontal force the prescribed displacement applies: sum R_x over the pushed edge.
// Returns a negative number if the run did not get there, which every caller checks.
double sliding_force(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return -1.0;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[1].ok) return -1.0;
    const auto& R = res[1];
    if (R.reaction.size() != 2u * (size_t)R.mesh.node_count) return -1.0;
    double rx = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.x[n]) < 1e-9) rx += R.reaction[2 * n];
    return std::fabs(rx);
}

void oracle_sliding_block(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    check(M.ok, "block meshed from the file's own settings");
    if (!M.ok) { std::printf("      (%s)\n", M.message.c_str()); return; }
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    check(res.size() == 2 && res[0].ok && res[1].ok, "initial + push phases converged");
    if (res.size() != 2 || !res[1].ok) return;
    const auto& R = res[1];

    // (a) The block SLIDES: it translates by the imposed amount instead of distorting. Every
    // node of the block carries the same u_x, which is what "it hardly deforms" means.
    check(R.reaction.size() == 2u * (size_t)R.mesh.node_count, "the static phase reports reactions");
    if (R.reaction.size() != 2u * (size_t)R.mesh.node_count) return;
    double ux_lo = 1e300, ux_hi = -1e300;
    for (int n = 0; n < R.mesh.node_count; ++n) {
        if (R.mesh.y[n] < 1e-9) continue;                 // the seam's outer side is the fixed world
        ux_lo = std::fmin(ux_lo, R.disp[2 * n]);
        ux_hi = std::fmax(ux_hi, R.disp[2 * n]);
    }
    std::printf("      block u_x spans %.9f .. %.9f m (imposed %.3f)\n", ux_lo, ux_hi, kSbPush);
    check(std::fabs(ux_lo - kSbPush) < 1e-4 && std::fabs(ux_hi - kSbPush) < 1e-4,
          "the block translates rigidly on the joint rather than shearing against a welded base");

    // (b) The published comparison.
    double rx = 0.0;
    for (int n = 0; n < R.mesh.node_count; ++n)
        if (std::fabs(R.mesh.x[n]) < 1e-9) rx += R.reaction[2 * n];
    const double F = std::fabs(rx);
    const double F_exact = sliding_capacity(kSbCw, kSbGamma);   // 60.076 at phi_w = 26.6 deg
    const double F_manual = 60.0, F_plaxis = 60.4;
    std::printf("      F: manual's arithmetic %.1f | exact for phi_w = %.1f deg %.4f | PLAXIS %.1f"
                " | file run %.4f (%+.2f%% vs exact, %+.2f%% vs PLAXIS)\n",
                F_manual, kSbPhiw, F_exact, F_plaxis, F,
                100.0 * (F - F_exact) / F_exact, 100.0 * (F - F_plaxis) / F_plaxis);
    check(std::fabs(F - F_exact) < 0.02 * F_exact, "failure force within 2% of the closed form");
    check(std::fabs(F - F_plaxis) < 0.02 * F_plaxis,
          "failure force within 2% of the published PLAXIS number");

    // (c) It is a LIMIT load, not a stiffness reading: pushing twice as far must not push
    // twice as hard. On a plateau the two runs agree to the last bit.
    const double F_far = sliding_force(build_sliding_block_at(kSbCw, kSbGamma, 2.0 * kSbPush));
    std::printf("      pushed 2x as far: %.9f vs %.9f kN/m\n", F_far, F);
    check(F_far > 0.0 && std::fabs(F_far - F) <= 1e-9 * F,
          "the force is a plateau: twice the imposed slip gives the same failure force");

    // (d) The two TERMS of the manual's formula, moved one at a time. Matching one number can
    // be a coincidence of two compensating errors; reproducing adhesion and friction
    // separately cannot. Doubling gamma leaves c_w alone and vice versa, and the relative
    // deviation must not move -- the discretisation bias is structural, not per-case.
    const double bias = (F - F_exact) / F_exact;
    struct Term { const char* what; double cw, gamma; };
    for (const Term& t : {Term{"gamma doubled ", kSbCw, 2.0 * kSbGamma},
                          Term{"adhesion off  ", 0.0, kSbGamma},
                          Term{"friction alone", 0.0, 2.0 * kSbGamma}}) {
        const double f = sliding_force(build_sliding_block_at(t.cw, t.gamma, kSbPush));
        const double e = sliding_capacity(t.cw, t.gamma);
        std::printf("      %s: %10.6f vs closed form %10.6f (%+.3f%%)\n", t.what, f, e,
                    100.0 * (f - e) / e);
        check(f > 0.0 && std::fabs((f - e) / e - bias) < 1e-6,
              "the term moves by exactly the closed form's amount (same relative bias)");
    }

    // (e) The sentry for the fault this case was built to catch. Delete the interface and the
    // block is welded to its base: the answer is not slightly stiffer, it is off the scale.
    m::Project welded = build_sliding_block();
    welded.structs.clear();
    welded.initial.struct_active.clear();
    welded.phases[0].struct_active.clear();
    const double F_welded = sliding_force(welded);
    std::printf("      with the interface removed: %.6g kN/m\n", F_welded);
    check(F_welded > 100.0 * F,
          "removing the interface changes the answer by orders of magnitude (it is not inert)");
}

// ------------------------------ KV-STR-003: bending of beams, PLAXIS Validation 2.3 ------
// The manual's own PLATE verification, and the first case in the corpus where a plate -- not
// the soil -- carries the entire load. Both of the manual's problems live in one model,
// built as it builds them: a single point load on one beam and a uniformly distributed load
// on another, "added to the bottom line with a spacing in between", with the soil cluster
// deactivated so that only the beams remain, supported at their end points.
//
// Why the case exists: the parity register counted Plate as verified on the strength of
// KV-STR-001, which is an ANCHOR prestress case -- the plate had no case of its own. Building
// this one found the reason it had never run. With the soil deactivated every node of the
// beam touches only passive elements, and fix_inactive_nodes pinned all of them: the beam was
// welded to the outside world along its whole length, and the file read max|u| = 0.000000e+00
// while reporting "ok" without one word of warning. Check (a) below is the sentry for that --
// if the exemption is ever taken back, the deflection does not drift, it collapses to zero.
constexpr double kBmEA = 1.64e6;     // HEB 200 in plane strain: a plate 1 m wide out of plane
constexpr double kBmEI = 1200.0;
constexpr double kBmNu = 0.0;
constexpr double kBmL = 2.0;         // span
constexpr double kBmGap = 1.0;       // "with a spacing in between"
constexpr double kBmF = 100.0;       // point load at mid-span [kN]
constexpr double kBmQ = 100.0;       // distributed load [kN/m]

// The closed forms, written out here rather than called from the plate header: the comparison
// must be the FE answer against the law, not the law against itself. A PLAXIS plate is a
// Mindlin (Timoshenko) beam, so mid-span deflection is bending PLUS shear, with the manual's
// own shear rigidity kGA' = k EA / (2(1+nu)), k = 5/6 (MMM Eq. 18-8).
double beam_kGA(double EA, double nu) { return (5.0 / 6.0) * EA / (2.0 * (1.0 + nu)); }
double beam_defl_point(double EA, double EI, double nu) {
    return kBmF * kBmL * kBmL * kBmL / (48.0 * EI) + kBmF * kBmL / (4.0 * beam_kGA(EA, nu));
}
double beam_defl_udl(double EA, double EI, double nu) {
    return 5.0 * kBmQ * kBmL * kBmL * kBmL * kBmL / (384.0 * EI) +
           kBmQ * kBmL * kBmL / (8.0 * beam_kGA(EA, nu));
}
// Bending moment along the span. Both problems peak at 50 kNm: F l / 4 = q l^2 / 8.
double beam_moment_point(double s) {
    return s <= 0.5 * kBmL ? kBmF * s / 2.0 : kBmF * (kBmL - s) / 2.0;
}
double beam_moment_udl(double s) { return kBmQ * s * (kBmL - s) / 2.0; }

m::Project build_beams_at(double EA, double EI, bool q_load, double elem_size) {
    m::Project pr;
    pr.name = "KV-STR-003 PLAXIS 2.3 bending of beams";
    pr.x_min = 0.0; pr.x_max = 2.0 * kBmL + kBmGap;
    pr.y_min = 0.0; pr.y_max = 1.0;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = elem_size;   // the manual: "A very coarse mesh is sufficient"
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material soil;             // never activated: the clusters exist to carry the beams
    soil.name = "Soil (deactivated)";
    soil.model = m::SoilModel::LinearElastic;
    soil.E = 1.0e4; soil.nu = 0.3;
    soil.k0_auto = false; soil.k0 = 0.5;
    pr.materials.push_back(soil);

    m::PlateMaterial heb;
    heb.name = "HEB 200";
    heb.EA = EA; heb.EI = EI; heb.nu = kBmNu; heb.w = 0.0;
    pr.plates.push_back(heb);

    for (int b = 0; b < 2; ++b) {
        const double x0 = b == 0 ? 0.0 : kBmL + kBmGap;
        m::SoilPolygon P;
        P.name = b == 0 ? "Block F" : "Block q";
        P.material = 0;
        P.x = {x0, x0 + kBmL, x0 + kBmL, x0};
        P.y = {0.0, 0.0, 1.0, 1.0};
        // The manual's "point fixities on the end points of the beam": left end pinned, right
        // end on a roller. The cluster is inactive, so these two edges hold nothing except the
        // beam's own end nodes -- which is exactly what a point fixity is here.
        P.edge_bc = {(int)m::BCType::Free, (int)m::BCType::VerticallyFixed,
                     (int)m::BCType::Free, (int)m::BCType::FullyFixed};
        pr.polygons.push_back(P);

        m::StructElement S;
        S.kind = m::StructKind::Plate;
        S.name = b == 0 ? "Beam F" : "Beam q";
        S.x1 = x0; S.y1 = 0.0; S.x2 = x0 + kBmL; S.y2 = 0.0;
        S.material = 0;
        pr.structs.push_back(S);
    }

    m::Load F;                    // "a single point load ... on a beam"
    F.kind = m::LoadKind::Point;
    F.name = "F";
    F.x1 = 0.5 * kBmL; F.y1 = 0.0; F.x2 = 0.5 * kBmL; F.y2 = 0.0;
    F.qx1 = F.qx2 = 0.0; F.qy1 = -kBmF; F.qy2 = 0.0;
    pr.loads.push_back(F);

    m::Load Q;                    // "a uniformly distributed load on a beam"
    Q.kind = m::LoadKind::Distributed;
    Q.name = "q";
    Q.x1 = kBmL + kBmGap; Q.y1 = 0.0; Q.x2 = 2.0 * kBmL + kBmGap; Q.y2 = 0.0;
    Q.qx1 = Q.qx2 = 0.0; Q.qy1 = Q.qy2 = -kBmQ;
    pr.loads.push_back(Q);

    pr.initial.poly_active = {0, 0};
    pr.initial.struct_active = {1, 1};
    pr.initial.load_active = {0, 0};
    m::Phase load;
    load.name = "Load";
    load.poly_active = {0, 0};
    load.struct_active = {1, 1};
    load.load_active = {1, (char)(q_load ? 1 : 0)};
    pr.phases.push_back(load);
    return pr;
}

// 0.25 m: eight elements per span. Still the coarse mesh the manual asks for, and coarse
// enough that the peak moment carries a visible discretisation bias -- see (c).
constexpr double kBmH = 0.25;
m::Project build_beams() { return build_beams_at(kBmEA, kBmEI, true, kBmH); }

// One run of the model: each beam's mid-span deflection (downward positive) and the peak |M|
// of its own force diagram.
struct BeamRead {
    bool ok = false;
    double w_F = 0.0, w_q = 0.0, M_F = 0.0, M_q = 0.0, max_u = 0.0;
    std::vector<katai::core::ForceStation> diag_F, diag_q;
};
BeamRead read_beams(const m::Project& pr) {
    BeamRead r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return r;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[0].ok || !res[1].ok) return r;
    const auto& R = res[1];
    if (R.disp.size() != 2u * (size_t)R.mesh.node_count) return r;
    r.max_u = R.max_disp;
    const auto mid_deflection = [&](double x) {
        for (int n = 0; n < R.mesh.node_count; ++n)
            if (std::fabs(R.mesh.x[n] - x) < 1e-9 && std::fabs(R.mesh.y[n]) < 1e-9)
                return -R.disp[2 * n + 1];
        return -1.0;                       // no node there: every caller checks for it
    };
    r.w_F = mid_deflection(0.5 * kBmL);
    r.w_q = mid_deflection(kBmL + kBmGap + 0.5 * kBmL);
    for (const auto& sf : R.struct_forces) {
        if (sf.name == "Beam F") { r.M_F = sf.max_M; r.diag_F = sf.stations; }
        if (sf.name == "Beam q") { r.M_q = sf.max_M; r.diag_q = sf.stations; }
    }
    r.ok = true;
    return r;
}

void oracle_beams(const m::Project& pr) {
    const BeamRead R = read_beams(pr);
    check(R.ok, "both phases converged and the beams' force diagrams were produced");
    if (!R.ok) return;
    check(R.w_F > 0.0 && R.w_q > 0.0, "a mesh node sits at each beam's mid-span");
    if (R.w_F <= 0.0 || R.w_q <= 0.0) return;

    // (a) The two published deflections. The closed form is Timoshenko because a PLAXIS plate
    // is: the shear part is only 0.5% of the answer, but it is the whole difference between
    // 13.889 mm (Euler-Bernoulli) and the 13.96 mm the manual prints.
    const double w_F_cf = beam_defl_point(kBmEA, kBmEI, kBmNu);
    const double w_q_cf = beam_defl_udl(kBmEA, kBmEI, kBmNu);
    std::printf("      u_max point load: closed form %.5f mm | PLAXIS 13.96 | file run %.5f mm"
                " (%+.3f%% vs closed form)\n",
                1e3 * w_F_cf, 1e3 * R.w_F, 100.0 * (R.w_F - w_F_cf) / w_F_cf);
    std::printf("      u_max distributed: closed form %.5f mm | PLAXIS 17.43 | file run %.5f mm"
                " (%+.3f%% vs closed form)\n",
                1e3 * w_q_cf, 1e3 * R.w_q, 100.0 * (R.w_q - w_q_cf) / w_q_cf);
    check(std::fabs(R.w_F - w_F_cf) < 0.01 * w_F_cf, "point-load deflection within 1% of the closed form");
    check(std::fabs(R.w_q - w_q_cf) < 0.01 * w_q_cf, "distributed deflection within 1% of the closed form");
    check(std::fabs(R.w_F - 13.96e-3) < 0.02 * 13.96e-3, "point-load deflection within 2% of PLAXIS");
    check(std::fabs(R.w_q - 17.43e-3) < 0.02 * 17.43e-3, "distributed deflection within 2% of PLAXIS");

    // (b) The manual's OTHER published pair: M_max = 50.0 kNm in both problems. A deflection
    // can be right for the wrong reason (a compensating support condition); the moment is read
    // from the element's own force diagram, through a different code path, and pins the same
    // answer.
    std::printf("      M_max: manual 50.0 / 50.0 kNm | file run %.5f / %.5f kNm (%+.2f%% / %+.2f%%)\n",
                R.M_F, R.M_q, 100.0 * (R.M_F - 50.0) / 50.0, 100.0 * (R.M_q - 50.0) / 50.0);
    check(std::fabs(R.M_F - 50.0) < 0.02 * 50.0, "point-load peak moment within 2% of the manual's 50 kNm");
    check(std::fabs(R.M_q - 50.0) < 0.02 * 50.0, "distributed peak moment within 2% of the manual's 50 kNm");

    // (c) Not just the peak -- the whole DIAGRAM. One number can be a coincidence; a moment
    // field that follows F s / 2 and q s (l - s) / 2 station by station cannot. The point-load
    // beam's field is piecewise linear and the element reproduces it to round-off; the
    // distributed beam's is a parabola sampled by an element whose curvature is linear, so it
    // is banded on the scale of that interpolation.
    struct Shape { const char* what; const std::vector<katai::core::ForceStation>* d;
                   double (*law)(double); double x0, band; };
    for (const Shape& sh : {Shape{"point load ", &R.diag_F, beam_moment_point, 0.0, 0.02},
                            Shape{"distributed", &R.diag_q, beam_moment_udl, kBmL + kBmGap, 0.08}}) {
        double worst = 0.0, worst_s = 0.0;
        for (const auto& st : *sh.d) {
            const double s = st.x - sh.x0;
            if (s < 0.15 * kBmL || s > 0.85 * kBmL) continue;   // the ends carry M ~ 0
            const double e = std::fabs(st.M - sh.law(s)) / 50.0;
            if (e > worst) { worst = e; worst_s = s; }
        }
        std::printf("      M(s) %s: worst deviation from the closed form %.4f%% of M_max (at s = %.3f m)\n",
                    sh.what, 100.0 * worst, worst_s);
        check(worst < sh.band, "the moment DISTRIBUTION follows the closed form, not just its peak");
    }

    // (d) What the peak moment of the DISTRIBUTED beam does with mesh size, because it is the
    // one published number this file does not reproduce exactly. Its moment field is a
    // parabola and the element's curvature is linear, so the peak overshoots -- by exactly
    // q h^2 / 12, the parabola the element cannot represent inside one span of length h. That
    // is a structural bias, not scatter: it is reproduced to five figures at three mesh sizes
    // and it vanishes as h -> 0. The deflection does not do this (it is exact at every mesh),
    // so "a very coarse mesh is sufficient" is true of the manual's displacements and not of
    // its peak moment -- worth knowing before reading a wall's M off a coarse run.
    std::printf("      M_max(distributed) vs mesh: ");
    for (double h : {0.5, kBmH, 0.125}) {
        const BeamRead V = read_beams(build_beams_at(kBmEA, kBmEI, true, h));
        const double rule = 50.0 + kBmQ * h * h / 12.0;
        std::printf("h=%.3f: %.5f (rule %.5f) ", h, V.M_q, rule);
        check(V.ok && std::fabs(V.M_q - rule) < 1e-4 * rule,
              "the peak-moment overshoot is exactly q h^2 / 12 at this mesh size");
    }
    std::printf("\n");

    // (e) The two TERMS of the deflection, separated. Bending and shear are added by the same
    // formula, so matching the total once proves neither. Stiffening EI by 4 divides the
    // bending term by 4 and leaves the shear term alone; stiffening EA by 100 does the
    // opposite and drives the answer onto the Euler-Bernoulli limit (13.8896 mm, the number
    // the manual would have published if a PLAXIS plate were not shear-deformable). The closed
    // form predicts both, and the FE run must follow it in each.
    struct Term { const char* what; double EA, EI; };
    for (const Term& t : {Term{"EI x 4  ", kBmEA, 4.0 * kBmEI},
                          Term{"EA x 100", 100.0 * kBmEA, kBmEI}}) {
        const BeamRead V = read_beams(build_beams_at(t.EA, t.EI, true, kBmH));
        const double cf_F = beam_defl_point(t.EA, t.EI, kBmNu), cf_q = beam_defl_udl(t.EA, t.EI, kBmNu);
        std::printf("      %s: point %.6f mm (cf %.6f, %+.3f%%) | udl %.6f mm (cf %.6f, %+.3f%%)\n",
                    t.what, 1e3 * V.w_F, 1e3 * cf_F, 100.0 * (V.w_F - cf_F) / cf_F,
                    1e3 * V.w_q, 1e3 * cf_q, 100.0 * (V.w_q - cf_q) / cf_q);
        check(V.ok && std::fabs(V.w_F - cf_F) < 0.01 * cf_F && std::fabs(V.w_q - cf_q) < 0.01 * cf_q,
              "bending and shear each move the answer by exactly the closed form's amount");
    }

    // (f) The beams do not feel each other: switching the distributed load off leaves the
    // point-loaded beam bit-for-bit where it was. Two spans in one file, one answer each.
    const BeamRead A = read_beams(build_beams_at(kBmEA, kBmEI, false, kBmH));
    std::printf("      with q switched off, the F beam reads %.12f mm (was %.12f)\n",
                1e3 * A.w_F, 1e3 * R.w_F);
    check(A.ok && std::fabs(A.w_F - R.w_F) <= 1e-12 * R.w_F, "the two spans are independent");

    // (g) The sentry, and the clearest statement of what was wrong. Delete the beams and the
    // model has no degrees of freedom left at all: the deactivated soil is entirely fixed, and
    // the run says so and stops. That is the honest answer for this geometry. Before the fix,
    // the beams were in exactly that state -- every translation pinned -- and only their
    // rotation DOFs remained free, which was enough for the solve to succeed, report "ok" and
    // hand back max|u| = 0. The difference between refusal and silence was three extra DOFs.
    m::Project bare = build_beams();
    bare.structs.clear();
    bare.initial.struct_active.clear();
    bare.phases[0].struct_active.clear();
    const auto MB = katai::app::mesh_from_project(bare);
    bool bare_refused = false;
    std::string bare_msg;
    if (MB.ok) {
        const auto rb = katai::app::solve_phases(bare, MB.mesh,
                                                 katai::app::initial_phase_from(bare.initial_procedure));
        for (const auto& p : rb)
            if (!p.ok) { bare_refused = true; bare_msg = p.message; break; }
    }
    std::printf("      with the beams removed: %s\n",
                bare_refused ? bare_msg.c_str() : "(the run still solved something)");
    check(bare_refused, "without the beams there is nothing left to solve, and the run says so");
}

int main() {
    std::printf("Input corpus: checked-in .k2d == programmatic build, validated, solved from the file\n");
    const CorpusCase cases[] = {
        {"kv-con-002-terzaghi-column.k2d", build_terzaghi, oracle_terzaghi},
        {"kv-fnd-008-strip-load.k2d", build_strip, oracle_strip},
        {"kv-fnd-009-flamant-line-load.k2d", build_flamant, oracle_flamant},
        {"kv-num-003-k0-geostatic-block.k2d", build_k0_block, oracle_k0_block},
        {"kv-cst-001-undrained-column.k2d", build_undrained_column, oracle_undrained_column},
        {"kv-slp-001-griffiths-lane-slope.k2d", build_gl_slope, oracle_gl_slope},
        {"kv-exc-001-staged-excavation.k2d", build_excavation, oracle_excavation},
        {"kv-dyn-002-resonant-column.k2d", build_resonant_column, oracle_resonant_column},
        {"kv-flw-001-charny-unconfined-dam.k2d", build_charny_dam, oracle_charny_dam},
        {"kv-dyn-003-el-centro-two-layer.k2d", build_el_centro, oracle_el_centro},
        {"kv-fnd-010-prandtl-strip-footing.k2d", build_prandtl_footing, oracle_prandtl_footing},
        {"kv-fnd-011-gibson-strip-load.k2d", build_gibson, oracle_gibson},
        {"kv-fnd-012-giroud-rigid-footing.k2d", build_giroud, oracle_giroud},
        {"kv-fnd-013-cox-circular-footing.k2d", build_cox, oracle_cox},
        {"kv-fnd-014-davis-booker-strip-footing.k2d", build_davis_booker, oracle_davis_booker},
        {"kv-slp-002-griffiths-lane-example1.k2d", build_gl_example1, oracle_gl_example1},
        {"kv-cst-002-hs-oedometer.k2d", build_hs_oedometer, oracle_hs_oedometer},
        {"kv-str-002-plaxis-sliding-block.k2d", build_sliding_block, oracle_sliding_block},
        {"kv-str-003-plaxis-beam-bending.k2d", build_beams, oracle_beams},
    };
    for (const CorpusCase& c : cases) run_case(c);

    if (g_failures == 0) {
        std::printf("\nOK: every corpus case is reproducible from its checked-in .k2d\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
