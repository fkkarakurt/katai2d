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
// verify: KV-CST-009
//   oracle:   closed_form
//   source:   the Soft Soil logarithmic compression law as published in the PLAXIS Material Models Manual chapter 10: e_v = lambda* ln(p'/p'_0) in primary loading (Eq 10-5) and e_v^e = kappa* ln(p'/p'_0) in unloading/reloading (Eq 10-6, with the tangent bulk modulus K_ur = p'/kappa* of Eq 10-7); M is not an input but is derived from K0nc (Eq 10-13) so that primary one-dimensional compression reaches that K0nc. The same equations are verified at the material point and through a hand-built BVP in test_soft_soil
//   locator:  a laterally confined weightless column walks the SAME stress range three times -- 50 -> 200 kPa primary, back to 50, and up to 200 again -- so the two indices and the pre-consolidation memory are each read from one file. Primary loading holds K0 = K0nc, so there the mean stress is proportional to the vertical one and the vertical stress ratio stands in for p'/p'_0. UNLOADING does not: with nu_ur = 0.15 the lateral stress falls far less than the vertical one, the ratio of horizontal to vertical stress RISES (sec. 10.3.5, "a well-known phenomenon in overconsolidated materials"), and sec. 10.3.1 draws the consequence that kappa* has no exact relation to the one-dimensional swelling index. The swelling leg therefore has its own closed form, evaluated on the mean stresses the elastic one-dimensional path d(sigma_h) = nu_ur/(1-nu_ur) d(sigma_v) produces (stated in full and evaluated in the test, not called from the material header)
//   quantity: settlement of each of the three legs, and the lateral stress ratio reached in primary loading, run from the checked-in tests/corpus/kv-cst-009-soft-soil-oedometer.k2d [m; -]
//   expected: primary lambda* ln(4) H = 0.110904 m with lambda* = 0.02, H = 4 m; swelling and reloading kappa* ln(p'_0/p'_1) H = 0.010186 m with kappa* = 0.004; K0 -> K0nc = 1 - sin(25 deg) = 0.5774
//   band:     2% on primary loading, 3% on the two elastic legs and 3% on K0nc, as asserted below -- measured -0.01% / -1.01% / -1.63% and +2.11% on the file's own 0.5 m tri6 mesh. The naive oracle kappa* ln(sigma_v/sigma_v0) reads 0.0222 m against a measured 0.0101, a factor of 2.2: the ORACLE is wrong there, not the run, and the manual says so in two separate places -- writing the closed form for the leg the soil actually walks is part of the verification, not a detail of it. Three further witnesses: the primary/reload ratio (11.07) matches the ratio of the two closed forms (10.89), which is the model's memory for the pre-consolidation stress stated as a number -- same file, same load, same stress range, an order of magnitude less settlement the second time; K0nc is MEASURED rather than the M formula being checked against itself; and the two indices are moved one at a time (lambda* x2, kappa* x2), each moving only its own leg and by exactly the closed form's amount
//
// verify: KV-CST-010
//   oracle:   closed_form
//   source:   the Soft Soil Creep differential law as published in the PLAXIS Material Models Manual chapter 11 (Buisman 1936, Bjerrum 1967, Garlanger 1972, Vermeer & Neher 1999): the volumetric creep rate is (mu*/tau)(p_eq/p_p^eq)^beta with beta = (lambda*-kappa*)/mu* (Eq 11-23), and tau is ONE DAY because the standard oedometer's 24-hour stage is the definition of the normal-consolidation line (Eq 11-13/14). The same law is verified at the material point in test_soft_soil_creep
//   locator:  on the normal-consolidation line p_eq = p_p, the rate reduces to mu*/tau independently of the stress level, and the differential law integrates exactly to e_v^c(t) = mu* ln(1 + t/tau) (derivation in docs/references/soft-soil-creep-formulation.md sec. 5.1, stated in full and evaluated in the test). Ground under its own weight, seeded normally consolidated by the K0 procedure, is left to sit for 100 days: NO load changes in the measured phase, only time passes, and the strain is uniform over the column even though the stress is not
//   quantity: surface settlement after 100 days of creep under self-weight alone, run from the checked-in tests/corpus/kv-cst-010-soft-soil-creep-column.k2d [m]
//   expected: mu* ln(1 + 100/1) H = 0.018460 m with mu* = 0.001, H = 4 m
//   band:     3%, as asserted below -- measured +0.74% on the file's own 0.5 m tri6 mesh with 50 time steps. The fixture is what it is because the manual predicted two earlier attempts failing: a weightless column loaded from zero cannot be used, because the initial pre-consolidation stress sits at the model's minimum of one stress unit, the first load puts p_eq far above it, and with beta = 16 the rate (p_eq/p_p)^beta collapses the run -- sec. 11.11's warning about unrealistic initial creep rates at OCR = 1 arriving as an arithmetic fact; and a zero-duration phase is elastic, because this model has no instantaneous plastic component at all (all inelastic strain is time-dependent). Three further witnesses: the law is sampled across three decades of time (1 / 10 / 1000 days, -2.06% / +1.62% / +0.61%), which no linear-in-time creep law could match at once and which locates tau at one day; the settlement is linear in mu*; and the differential witness -- the SAME file with the same ground as plain Soft Soil, which has every feature of this model except the creep, moves EXACTLY 0.000e+00 m over the same hundred days, so what is measured is creep and not a slow numerical drift
//
// verify: KV-CST-008
//   oracle:   closed_form
//   source:   the Hardening Soil with small-strain stiffness degradation law as published in the PLAXIS Material Models Manual chapter 7 (modified Hardin-Drnevich after Santos & Correia 2001): secant G_s/G0 = 1/(1 + a |gamma|/gamma_ref) with a = 0.385 (Eq 7-3), tangent G_t = G0/(1 + a gamma/gamma_ref)^2 (Eq 7-8) cut off below at G_ur = E_ur/(2(1+nu_ur)) (Eq 7-9), and Masing's rule gamma_0.7,re-loading = 2 gamma_0.7,virgin-loading (Eq 7-11), which the manual applies as a constant factor "throughout loading" rather than at a detected reversal; the same equations are verified at the material point in test_hssmall
//   locator:  one-dimensional UNLOADING of a laterally confined column: excavating h_exc of a gamma = 20 kN/m3 column relieves d(sigma) = gamma h_exc uniformly over the remaining depth. With m = 0 the stiffness is stress-independent, so the strain is uniform and the heave is eps x H_rem; in one-dimensional strain gamma = eps (gamma = sqrt(3/2 e:e) with e_yy = 2eps/3, e_xx = e_zz = -eps/3). Inverting the secant law sigma = E_oed,0 eps/(1 + a eps/gamma_ref) gives eps = d(sigma)/(E_oed,0 - a d(sigma)/gamma_ref), with E_oed,0 from E_0 = 2(1+nu_ur) G_0 (stated in full and evaluated in the test, not called from the material header)
//   quantity: heave of the excavated floor after deactivating the upper 2 m of a 10 m HSsmall column, run from the checked-in tests/corpus/kv-cst-008-hssmall-unloading.k2d [m]
//   expected: the closed form above with G0^ref = 187.5 MPa (E_0 = 450 MPa = 5 E_ur), gamma_0.7 = 1.5e-4, E_ur^ref = 90 MPa, nu_ur = 0.2 -> 7.1322e-4 m
//   band:     2%, as asserted below -- measured -0.30% (7.1108e-4) on the file's own 1.0 m tri6 mesh. UNLOADING is the point: in the HS family a deviatoric LOADING path is plastic from the first increment (the K0 state is seeded onto the shear surface), so there is no elastic window in which to measure a stiffness, whereas excavating moves away from both surfaces and is purely quasi-elastic. That the window really is elastic is not assumed but shown: the plain-HS twin, the same file with the overlay switched off (G0 = 0), reproduces ITS closed form -- elastic unloading at E_ur -- to +0.000%, and it heaves 4.500x as much, so the overlay is not a small correction here. Two further witnesses: setting G0 = G_ur makes the overlay redundant (E_0 = E_ur) and the run falls back onto plain HS BIT-FOR-BIT, which pins the overlay to the elastic branch alone (a version that also touched the plastic moduli or the failure surface could not close that identity); and the law is sampled at three points of the degradation curve, h_exc = 1 / 2 / 4 m, holding to -0.14% / -0.30% / -0.72%. Those three also pin Masing's factor: riding the VIRGIN backbone (gamma_0.7 instead of 2 gamma_0.7) would read +5.7% / +12.9% / +34.8% high, a deviation that grows with strain -- the signature of the wrong threshold, not of discretisation. The build did ride the virgin backbone until 2026-08-11
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
//
// verify: KV-STR-005
//   oracle:   closed_form
//   source:   the geogrid as specified in the PLAXIS 2D Reference Manual sec 6.5: "geogrids are flexible elastic or elastoplastic elements that represent a grid or sheet of fabric. Geogrids can only sustain tensile forces, but not compressive forces", with the axial stiffness defined in Eq. 6-51 as "the ratio of the axial force F per unit width and the axial strain (eps = dl/l)", EA = F/eps
//   locator:  a geogrid's translational degrees of freedom ARE the soil's -- it is a conforming chain of mesh nodes with no rotation -- so making the soil strain uniformly makes the geogrid strain with it, exactly. A homogeneous weightless block held at u_x = 0 on one side and pulled to u_x = D on the other deforms affinely, u_x = D x/L, so every horizontal fibre carries eps = D/L; the reinforcement spans the full width, so its ends sit ON the two boundaries and its elongation is D whatever the soil does. Its tension is then constant along its length, which means it applies no body force to the soil and the affine field remains the exact solution -- the closed form is not an approximation of this problem, it is this problem
//   quantity: axial force in the reinforcement, and the tension cut-off it stops at, run from the checked-in tests/corpus/kv-str-005-geogrid-tension.k2d [kN/m]
//   expected: EA D/L = 5000 x 0.004/10 = 2.0 kN/m elastic, and N_p = 3.0 kN/m past the cut-off, stated in full and evaluated in the test rather than called from the geogrid header
//   band:     1%, as asserted below -- measured 2.000000 kN/m, +0.000%, at EVERY station of the diagram on the file's own 0.5 m tri6 mesh. The fixture proves itself first: with no reinforcement the edge reaction is 21.978022 kN/m against the plane-strain closed form E/(1-nu^2) eps H = 21.978022, +0.0000%, so the field really is affine and the fibre strain really is D/L. Five further witnesses: the force is linear in BOTH inputs moved separately (twice the stretch and twice the stiffness each double it); past N_p it stops at N_p and doubling the stretch again changes nothing, while the elastic twin at the same stretch carries what it was told to; compressing the block instead of stretching it leaves the reinforcement carrying NOTHING, which is the manual's other sentence about this element; and the soil reaction is bit-for-bit identical with and without the reinforcement, which is what a constant-tension member must do and is the affine argument itself, measured. Building this case found a silent-wrong that was not the geogrid's: EVERY structural element (plate, plate5, geogrid, anchor, interface, embedded beam) read a node held by a NON-ZERO prescribed displacement as if it had not moved, because the element loops built their displacement vector from FREE degrees of freedom only. The limit was declared in internal_forces.hpp behind a guard that had expired -- "prescribed u_bar is today a kernel/test path, not available as a deformation BC in the GUI" -- which schema v2 made false when it added `disps`. Before the fix this case read 7.6363 kN/m instead of 2.0000, and its converged field was not affine: the node beside the driven edge sat 3.79% below its own left neighbour
//
// verify: KV-STR-004
//   oracle:   closed_form
//   source:   the embedded beam (pile row) as specified in the PLAXIS 2D Reference Manual sec 5.6.3 (connection point) and sec 6.6.3 (skin and base resistance), with the element formulation in the Scientific Manual sec 7.5. The manual fixes the loading path for us twice over: the material data set carries "only the bearing capacity (skin resistance and base resistance)" and not the stiffness response, and "embedded beams are not meant to be used as laterally loaded piles and will therefore not show accurate failure loads when subjected to transverse forces" -- so the defining quantity of this element is its AXIAL capacity
//   locator:  a pile row is smeared over a metre of wall, so every per-pile quantity is divided by the out-of-plane spacing exactly as EA, EI and the pile weight are (the reason Eq 6-65 divides the interface stiffnesses by L_spacing). Its ultimate axial load is therefore the skin resistance over the embedded length plus the base resistance, all per metre of wall. The head is pushed far past that limit and the PILE's own axial force is read, not the applied load: the load is applied at a soil node that the hinged head shares, so the soil carries the remainder and the plateau lives in the pile's force diagram
//   quantity: axial force at the head of the pile row under a head load of 1500 kN/m, far above its capacity, run from the checked-in tests/corpus/kv-str-004-axial-pile-capacity.k2d [kN/m]
//   expected: (T_skin,max L + F_max,base)/L_spacing = (100 x 10 + 500)/2.5 = 600 kN/m, stated in full and evaluated in the test rather than called from the driver
//   band:     2%, as asserted below -- measured 600.0000 kN/m, -0.00%, on the file's own 1.0 m tri6 mesh. The fixture obeys two rules the manual states: the soil is MOHR-COULOMB and not Linear Elastic, because PLAXIS ignores the shaft resistance AND the spacing inside a linear elastic cluster (it counts that as structure rather than soil), so an LE fixture would measure the one case PLAXIS treats differently; and its cohesion is far above anything mobilised, so the plateau measured is the pile's declared capacity and not a soil bearing failure. Soil and pile are weightless, so the load carried is the load applied. Five further witnesses: doubling the head load leaves the pile force at 600.0000 while the head goes on settling (0.116 -> 0.479 m), which is a limit load and not a stiffness reading; the two capacity terms are moved one at a time and each moves the total by exactly its own share (base 500 -> 100 kN gives 440.0000 against 440.0000, skin 100 -> 50 kN/m gives 400.0000 against 400.0000); doubling the out-of-plane spacing halves the capacity to 300.0000, ratio 0.5000 exactly, which is the check that Eq 6-65's division by L_spacing reaches the capacities and not only the stiffnesses; halving the element size leaves the capacity where it was, so nothing here is discretisation (a prediction of this case's own draft, that the tied node's Newton-Cotes share of the skin could not mobilise, was refuted by that measurement and is recorded in the test); and the sentry -- with the connection FREE, which is what this engine did for every pile until 2026-08-13 with no way to ask for anything else, the pile carries 0.0000 kN/m at its head, because a point load is delivered to the nearest SOIL node and reaches a free pile top only through the skin springs. Building this case also found and fixed two constants: Eq 6-65's division by L_spacing was not applied at all, leaving every skin and foot spring 2.5x too stiff at the default spacing, and the foot used D/2 where Eq 6-67 defines R_eq = sqrt(12 EI/EA)/2 = 0.433 D for a solid circular pile. The stiffness function's only consumer is the driver -- the element test passes its springs by hand -- so that factor of 2.5 stood while all 150 tests were green
#include <katai/analysis/response_spectrum.hpp>
#include <katai/mesh/boundary_extraction.hpp>   // collect_chain: the chain a geogrid is built on
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

// ------------------------------ KV-CST-008: HSsmall unloading, MMM ch. 7 -----------------
// The first HSsmall boundary-value case: the small-strain stiffness decides the answer, and
// it is read from the file. The model was already verified at the material point
// (test_hssmall runs Eq 7-3/7-7/7-8/7-10 against the closed form) -- what was unwitnessed is
// the path a user walks, and that is where the fault was.
//
// Why UNLOADING: in the HS family a deviatoric LOADING path is plastic from the first
// increment (the K0 state is seeded onto the shear surface, so there is no elastic window to
// measure a stiffness in). Excavating is the opposite -- it moves away from both surfaces, so
// the response is purely quasi-elastic and it is exactly the stiffness this overlay sets.
// The K0 procedure produces no strain, so the phase starts at gamma_hist = 0, which is where
// E0 lives. The plain-HS twin (checked below) reproduces its own closed form EXACTLY, which
// is what proves the window is elastic rather than merely assumed to be.
constexpr double kHsG0 = 187500.0;      // G0^ref: E0 = 2(1+nu_ur) G0 = 450 MPa = 5 x Eur
constexpr double kHsG07 = 1.5e-4;       // gamma_0.7 (VIRGIN loading, as the manual defines it)
constexpr double kHsEur = 90000.0;
constexpr double kHsNu = 0.2;           // nu_ur
constexpr double kHsGamma = 20.0;       // unit weight -> d(sigma) = gamma h_exc
constexpr double kHsDepth = 10.0;       // model depth

// Written out here rather than called from the material header, so the comparison is the FE
// answer against the law and not the law against itself. One-dimensional unloading of a
// laterally confined column: the stress relief is uniform, m = 0 makes the stiffness
// stress-independent, so the strain is uniform and the heave is eps x H_rem. The tangent
// modulus rides the Hardin-Drnevich hyperbola (Eq 7-8) whose integral is the secant law
// (Eq 7-7), and gamma = eps in one-dimensional strain (gamma = sqrt(3/2 e:e) with
// e_yy = 2eps/3, e_xx = e_zz = -eps/3). Inverting sigma = E_oed,0 eps / (1 + a eps/gamma_ref)
// gives the closed form below. gamma_ref is the RELOADING threshold 2 gamma_0.7 (Eq 7-11).
double hss_oed(double E) { return E * (1.0 - kHsNu) / ((1.0 + kHsNu) * (1.0 - 2.0 * kHsNu)); }
double hss_heave(double dsigma, double H_rem, double masing) {
    const double E0 = 2.0 * (1.0 + kHsNu) * kHsG0;
    const double eps = dsigma / (hss_oed(E0) - 0.385 * dsigma / (masing * kHsG07));
    return eps * H_rem;
}
double hs_heave(double dsigma, double H_rem) { return dsigma * H_rem / hss_oed(kHsEur); }

// h_exc metres of the column are excavated in the staged phase. G0 = 0 builds the plain-HS
// twin (the schema's own switch: no G0, no small-strain overlay).
m::Project build_hss_at(double h_exc, double G0) {
    m::Project pr;
    pr.name = "KV-CST-008 HSsmall unloading";
    pr.x_min = 0.0; pr.x_max = 2.0; pr.y_min = 0.0; pr.y_max = kHsDepth;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 1.0;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "HSsmall";
    s.model = m::SoilModel::HSsmall;
    s.gamma_unsat = s.gamma_sat = kHsGamma;
    s.E = 1.0e4; s.nu = kHsNu;              // not read by this model (HS reads its own pair)
    s.c = 10.0; s.phi = 35.0; s.psi = 0.0;
    s.tension_cutoff = false;               // K2D-M001: HS does not read it; say so in the file
    s.E50ref = 3.0e4; s.Eoedref = 3.0e4; s.Eurref = kHsEur;
    s.m = 0.0;                              // stress-independent: keeps the closed form exact
    s.nu_ur = kHsNu; s.p_ref = 100.0; s.Rf = 0.9;
    s.G0ref = G0; s.gamma07 = kHsG07;
    pr.materials.push_back(s);

    const double y_cut = kHsDepth - h_exc;
    m::SoilPolygon keep;
    keep.name = "Remaining";
    keep.material = 0;
    keep.x = {0.0, 2.0, 2.0, 0.0};
    keep.y = {0.0, 0.0, y_cut, y_cut};
    keep.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                    (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(keep);
    m::SoilPolygon dig;
    dig.name = "Excavated";
    dig.material = 0;
    dig.x = {0.0, 2.0, 2.0, 0.0};
    dig.y = {y_cut, y_cut, kHsDepth, kHsDepth};
    dig.edge_bc = {(int)m::BCType::Free, (int)m::BCType::HorizontallyFixed,
                   (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(dig);

    pr.initial.poly_active = {1, 1};
    m::Phase cut;
    cut.name = "Excavate";
    cut.poly_active = {1, 0};
    pr.phases.push_back(cut);
    return pr;
}

constexpr double kHsExc = 2.0;
m::Project build_hss() { return build_hss_at(kHsExc, kHsG0); }

// The heave at the excavated floor. Negative = the run did not get there.
double hss_run(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return -1.0;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[0].ok || !res[1].ok) return -1.0;
    return res[1].max_disp;
}

void oracle_hss(const m::Project& pr) {
    const double H_rem = kHsDepth - kHsExc, dsig = kHsGamma * kHsExc;
    const double u = hss_run(pr);
    check(u > 0.0, "the K0 and excavation phases converged");
    if (u <= 0.0) return;

    // (a) The published law, as the manual's model rides it: the reloading curve, 2 gamma_0.7.
    const double cf = hss_heave(dsig, H_rem, 2.0);
    const double cf_virgin = hss_heave(dsig, H_rem, 1.0);
    std::printf("      heave: closed form %.6e m | file run %.6e m (%+.2f%%); riding the VIRGIN"
                " backbone instead would give %.6e (%+.1f%%)\n",
                cf, u, 100.0 * (u - cf) / cf, cf_virgin, 100.0 * (cf_virgin - cf) / cf);
    check(std::fabs(u - cf) < 0.02 * cf, "heave within 2% of the small-strain closed form");

    // (b) The plain-HS twin, from the same file with the overlay switched off (G0 = 0). It must
    // reproduce ITS closed form -- elastic unloading at E_ur -- and it does so exactly, which is
    // what establishes that this window is quasi-elastic and the comparison in (a) is a
    // stiffness measurement rather than a plasticity one. It is also the differential witness:
    // the overlay is not a small correction here, it is a factor of four.
    const double u_hs = hss_run(build_hss_at(kHsExc, 0.0));
    const double cf_hs = hs_heave(dsig, H_rem);
    std::printf("      plain HS twin: closed form %.6e m | run %.6e m (%+.3f%%), ratio to HSsmall %.3f\n",
                cf_hs, u_hs, 100.0 * (u_hs - cf_hs) / cf_hs, u_hs / u);
    check(u_hs > 0.0 && std::fabs(u_hs - cf_hs) < 1e-6 * cf_hs,
          "the plain-HS twin reproduces elastic unloading at E_ur exactly (the window IS elastic)");
    check(u_hs > 3.0 * u, "the small-strain overlay is not a small correction (>3x here)");

    // (c) Identity: with G0 set to G_ur the overlay has nothing to add (E0 = E_ur), and the run
    // must fall back onto plain HS BIT-FOR-BIT. This pins the overlay to the right place in the
    // model -- a version that also touched the plastic moduli, or the failure surface, would
    // not close this identity.
    const double u_id = hss_run(build_hss_at(kHsExc, kHsEur / (2.0 * (1.0 + kHsNu))));
    std::printf("      G0 = G_ur identity: %.12e vs plain HS %.12e\n", u_id, u_hs);
    check(u_id > 0.0 && std::fabs(u_id - u_hs) <= 1e-12 * u_hs,
          "G0 = G_ur reduces to plain HS bit-for-bit (the overlay sits on the elastic branch only)");

    // (d) Not one point of the curve -- three. The degradation law is a function of strain, so
    // three unloadings of different size sample three different points of it, and the closed
    // form must hold at each. This is what pins Masing's factor: on the virgin backbone the
    // deviations below would be +5% .. +21% instead of a fraction of a percent, growing with
    // the strain, which is the signature of the wrong threshold rather than of discretisation.
    for (double h : {1.0, 4.0}) {
        const double uu = hss_run(build_hss_at(h, kHsG0));
        const double c2 = hss_heave(kHsGamma * h, kHsDepth - h, 2.0);
        const double c1 = hss_heave(kHsGamma * h, kHsDepth - h, 1.0);
        std::printf("      h_exc = %.0f m: run %.6e | closed form %.6e (%+.2f%%) | virgin backbone %+.1f%%\n",
                    h, uu, c2, 100.0 * (uu - c2) / c2, 100.0 * (c1 - c2) / c2);
        check(uu > 0.0 && std::fabs(uu - c2) < 0.02 * c2,
              "the closed form holds at another point of the degradation curve");
    }
}

// ------------------------------ KV-CST-009: Soft Soil oedometer, MMM ch. 10 --------------
// The Soft Soil model's defining behaviour is logarithmic compression with a SEPARATE
// unloading line and a memory for the pre-consolidation stress, and this case walks the same
// stress range three times to read all three of those off one file: primary loading at
// lambda*, unloading at kappa*, and reloading at kappa* again because the cap remembers where
// it has been. The material point and a hand-built BVP were already verified (test_soft_soil);
// what was unwitnessed is the path from a .k2d through the mesher, the K0/gravity initial
// state, the cap seeding and the load stepping.
constexpr double kSsLam = 0.02;     // lambda*, modified compression index
constexpr double kSsKap = 0.004;    // kappa*, modified swelling index (lambda*/kappa* = 5)
constexpr double kSsPhi = 25.0;
constexpr double kSsH = 4.0;        // column height
constexpr double kSsSeat = 50.0;    // seating stress
constexpr double kSsFull = 200.0;   // loaded stress

constexpr double kSsNu = 0.15;      // nu_ur, the manual's default for this model

// The manual's own law (Eq 10-5/10-6), written out here rather than called from the material
// header: e_v = index * ln(p'/p'_0). In one-dimensional strain the column is laterally
// confined, so e_v IS the vertical strain, and weightless soil keeps the stress uniform, so
// the settlement is e_v * H exactly.
//
// PRIMARY loading stays on the K0nc line -- the model's M is derived precisely so that
// one-dimensional compression reaches K0nc and then holds it -- so there the mean stress is
// proportional to the vertical one and the vertical stress ratio can stand in for p'/p'_0.
double ss_settlement(double index) { return index * std::log(kSsFull / kSsSeat) * kSsH; }

// UNLOADING is not the mirror of that, and the manual says why twice. Sec. 10.3.5: with a small
// nu_ur the lateral stress falls far less than the vertical one in one-dimensional unloading,
// so the ratio of horizontal to vertical stress RISES -- "a well-known phenomenon in
// overconsolidated materials". Sec. 10.3.1 draws the consequence: there is no exact relation
// between kappa* and the one-dimensional swelling index Cs, "because the ratio of horizontal
// and vertical stresses changes during one-dimensional unloading". So the mean stress does NOT
// follow the vertical stress on this leg, and kappa* ln(sigma_v/sigma_v0) is the wrong closed
// form -- it over-predicts the swelling by a factor of 2.2 here. The elastic one-dimensional
// path gives d(sigma_h) = nu_ur/(1 - nu_ur) d(sigma_v), and the law is then evaluated on the
// mean stresses that path actually produces.
double ss_swell(double kap) {
    const double K0nc = 1.0 - std::sin(kSsPhi * kPi / 180.0);
    const double sh0 = K0nc * kSsFull, p0 = (kSsFull + 2.0 * sh0) / 3.0;
    const double sh1 = sh0 + (kSsNu / (1.0 - kSsNu)) * (kSsSeat - kSsFull);
    const double p1 = (kSsSeat + 2.0 * sh1) / 3.0;
    return kap * std::log(p0 / p1) * kSsH;
}

m::Project build_ss_at(double lam, double kap) {
    m::Project pr;
    pr.name = "KV-CST-009 Soft Soil oedometer";
    pr.x_min = 0.0; pr.x_max = 1.0; pr.y_min = 0.0; pr.y_max = kSsH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::GravityLoading;
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Soft Soil clay";
    s.model = m::SoilModel::SoftSoil;
    s.gamma_unsat = s.gamma_sat = 0.0;       // weightless: the applied stress is the whole story
    s.c = 0.0; s.phi = kSsPhi; s.psi = 0.0;  // the manual's default dilatancy for this model
    s.tension_cutoff = false;                // K2D-M001: the SS return does not read it
    s.lam_star = lam; s.kap_star = kap;
    s.nu_ur = 0.15;                          // the manual's default
    s.k0nc_auto = true;                      // M is derived from K0nc (Eq 10-13)
    s.k0_auto = true;                        // start ON the K0nc line, so it stays there
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0.0, 1.0, 1.0, 0.0};
    P.y = {0.0, 0.0, kSsH, kSsH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    m::Load seat;
    seat.kind = m::LoadKind::Distributed;
    seat.name = "Seating stress";
    seat.x1 = 0.0; seat.y1 = kSsH; seat.x2 = 1.0; seat.y2 = kSsH;
    seat.qx1 = seat.qx2 = 0.0; seat.qy1 = seat.qy2 = -kSsSeat;
    pr.loads.push_back(seat);
    m::Load inc;
    inc.kind = m::LoadKind::Distributed;
    inc.name = "Load increment";
    inc.x1 = 0.0; inc.y1 = kSsH; inc.x2 = 1.0; inc.y2 = kSsH;
    inc.qx1 = inc.qx2 = 0.0; inc.qy1 = inc.qy2 = -(kSsFull - kSsSeat);
    pr.loads.push_back(inc);

    pr.initial.load_active = {1, 0};              // the seating stress establishes the state
    m::Phase load;  load.name  = "Load to 200 kPa";   load.load_active  = {1, 1};
    m::Phase unl;   unl.name   = "Unload to 50 kPa";  unl.load_active   = {1, 0};
    m::Phase rel;   rel.name   = "Reload to 200 kPa"; rel.load_active   = {1, 1};
    pr.phases.push_back(load);
    pr.phases.push_back(unl);
    pr.phases.push_back(rel);
    return pr;
}

m::Project build_ss() { return build_ss_at(kSsLam, kSsKap); }

// The three phases' settlements (each reported relative to its own start, so each IS the
// increment of that leg) plus the lateral stress ratio reached in primary loading.
struct SsRead { bool ok = false; double load = 0, unload = 0, reload = 0, k0 = 0; };
SsRead read_ss(const m::Project& pr) {
    SsRead r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return r;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 4) return r;
    for (const auto& p : res) if (!p.ok) return r;
    r.load = res[1].max_disp; r.unload = res[2].max_disp; r.reload = res[3].max_disp;
    // K0 = sigma_h/sigma_v at the end of primary loading, averaged over the column's interior
    // (the loaded top and the fixed base are boundary-disturbed).
    double acc = 0.0; int n_acc = 0;
    for (int n = 0; n < res[1].mesh.node_count; ++n) {
        const double y = res[1].mesh.y[n];
        if (y < 0.5 || y > kSsH - 0.5) continue;
        const double sv = res[1].stress.stress[n](1), sh = res[1].stress.stress[n](0);
        if (sv < -1.0) { acc += sh / sv; ++n_acc; }
    }
    r.k0 = n_acc ? acc / n_acc : 0.0;
    r.ok = true;
    return r;
}

void oracle_ss(const m::Project& pr) {
    const SsRead R = read_ss(pr);
    check(R.ok, "the initial, loading, unloading and reloading phases all converged");
    if (!R.ok) return;

    // (a) Primary loading rides the lambda* line.
    const double cf_lam = ss_settlement(kSsLam), cf_kap = ss_swell(kSsKap);
    std::printf("      primary loading: closed form %.6f m | file run %.6f m (%+.2f%%)\n",
                cf_lam, R.load, 100.0 * (R.load - cf_lam) / cf_lam);
    check(std::fabs(R.load - cf_lam) < 0.02 * cf_lam, "primary settlement within 2% of lambda* ln(4)");

    // (b) The SAME stress range unloaded rides the kappa* line -- a different index, measured
    // separately on the same file. This is the model's "distinction between primary loading and
    // unloading/reloading" and it cannot be faked by a single stiffness. The comparison is
    // against the closed form of the leg the soil actually walks (see ss_swell): the naive
    // kappa* ln(sigma_v/sigma_v0) reads 0.0222 m against a measured 0.0101, and it is the
    // ORACLE that is wrong there, not the run.
    std::printf("      unloading:       closed form %.6f m | file run %.6f m (%+.2f%%)\n",
                cf_kap, R.unload, 100.0 * (R.unload - cf_kap) / cf_kap);
    check(std::fabs(R.unload - cf_kap) < 0.03 * cf_kap, "swelling within 3% of the kappa* leg");

    // (c) Reloading the same range stays on kappa* -- the cap REMEMBERS the pre-consolidation
    // stress it reached, so the ground does not compress a second time. Same file, same load,
    // same stress range, an order of magnitude less settlement: that ratio is the sharpest
    // statement of the model's memory, and both of its legs are pinned to closed forms above.
    std::printf("      reloading:       closed form %.6f m | file run %.6f m (%+.2f%%); "
                "primary/reload ratio %.3f (closed forms %.3f)\n",
                cf_kap, R.reload, 100.0 * (R.reload - cf_kap) / cf_kap,
                R.load / R.reload, cf_lam / cf_kap);
    check(std::fabs(R.reload - cf_kap) < 0.03 * cf_kap, "reloading within 3% of the kappa* leg");
    check(std::fabs(R.load / R.reload - cf_lam / cf_kap) < 0.05 * (cf_lam / cf_kap),
          "the pre-consolidation memory: primary/reload matches the two closed forms");

    // (d) M is not an input: the file gives K0nc and the manual derives M from it (Eq 10-13) so
    // that primary one-dimensional compression REACHES that K0nc. Measuring the lateral stress
    // ratio is what verifies the derivation -- pinning M by its formula would only check the
    // formula against itself.
    const double k0nc = 1.0 - std::sin(kSsPhi * kPi / 180.0);
    std::printf("      K0 in primary loading: %.6f (K0nc = 1 - sin(phi) = %.6f, %+.2f%%)\n",
                R.k0, k0nc, 100.0 * (R.k0 - k0nc) / k0nc);
    check(std::fabs(R.k0 - k0nc) < 0.03 * k0nc,
          "primary one-dimensional compression reaches K0nc (so M was derived correctly)");

    // (e) The law is linear in its index: doubling lambda* doubles the primary settlement and
    // leaves the unloading leg alone; doubling kappa* does the opposite. Two indices moved one
    // at a time -- the same separation of terms that KV-STR-002 used on adhesion and friction.
    const SsRead A = read_ss(build_ss_at(2.0 * kSsLam, kSsKap));
    const SsRead B = read_ss(build_ss_at(kSsLam, 2.0 * kSsKap));
    std::printf("      lambda* x2: primary %.6f (cf %.6f) | kappa* x2: unloading %.6f (cf %.6f)\n",
                A.load, ss_settlement(2.0 * kSsLam), B.unload, ss_swell(2.0 * kSsKap));
    check(A.ok && std::fabs(A.load - ss_settlement(2.0 * kSsLam)) < 0.02 * ss_settlement(2.0 * kSsLam),
          "doubling lambda* doubles the primary settlement, as the closed form says");
    check(B.ok && std::fabs(B.unload - ss_swell(2.0 * kSsKap)) < 0.03 * ss_swell(2.0 * kSsKap),
          "doubling kappa* doubles the swelling, as the closed form says");
}

// ------------------------------ KV-CST-010: Soft Soil Creep, MMM ch. 11 ------------------
// Secondary compression from a .k2d: the ground is not loaded at all in the measured phase,
// only TIME passes, and the settlement that appears is the whole point of this model.
//
// Getting the experiment right took two attempts, and the manual predicted both failures. A
// weightless column loaded from zero (the KV-CST-009 fixture) cannot be used: the initial
// pre-consolidation stress sits at the model's minimum of one stress unit, the first load puts
// p_eq far above it, and the creep rate goes as (p_eq/p_p)^beta with beta = (lambda*-kappa*)/mu*
// = 16 here -- the run collapses, which is sec. 11.11's warning about unrealistically high
// initial creep rates at OCR = 1 arriving as an arithmetic fact. A phase of zero duration is no
// use either: in this model there is no instantaneous plastic component at all (all inelastic
// strain is time-dependent), so a zero-duration phase is elastic. What works is what the model
// is FOR: ground under its own weight, seeded normally consolidated by the K0 procedure, left
// to sit. There p_eq = p_p everywhere, the rate is exactly mu*/tau regardless of depth, and the
// strain is uniform even though the stress is not.
constexpr double kScLam = 0.02, kScKap = 0.004, kScMu = 0.001;   // lambda*/mu* = 20 (sec. 11.8.1)
constexpr double kScH = 4.0, kScGamma = 15.0, kScDays = 100.0;

// Eq 11-13/14 as the manual reads them: tau is ONE DAY, because the standard oedometer's
// 24-hour stage is the definition of the normal-consolidation line. Under constant effective
// stress on that line the differential creep law integrates exactly (soft-soil-creep-
// formulation.md sec. 5.1): e_v^c(t) = mu* ln(1 + t/tau). Written out here, not called from
// the material header.
double ssc_creep(double mu, double days) { return mu * std::log(1.0 + days / 1.0) * kScH; }

m::Project build_ssc_at(double mu, double days, bool creep_model) {
    m::Project pr;
    pr.name = "KV-CST-010 Soft Soil Creep column";
    pr.x_min = 0.0; pr.x_max = 1.0; pr.y_min = 0.0; pr.y_max = kScH;
    pr.has_water = false;
    pr.initial_procedure = m::InitialProcedure::K0Procedure;
    pr.mesh.elem_size = 0.5;
    pr.mesh.order = 6;
    pr.mesh.auto_refine = false;

    m::Material s;
    s.name = "Soft Soil Creep clay";
    s.model = creep_model ? m::SoilModel::SoftSoilCreep : m::SoilModel::SoftSoil;
    s.gamma_unsat = s.gamma_sat = kScGamma;
    s.c = 0.0; s.phi = kSsPhi; s.psi = 0.0;
    s.tension_cutoff = false;
    s.lam_star = kScLam; s.kap_star = kScKap; s.mu_star = mu;
    s.nu_ur = kSsNu;
    s.k0nc_auto = true; s.k0_auto = true;      // seeded ON the K0nc line, normally consolidated
    pr.materials.push_back(s);

    m::SoilPolygon P;
    P.name = "Column";
    P.material = 0;
    P.x = {0.0, 1.0, 1.0, 0.0};
    P.y = {0.0, 0.0, kScH, kScH};
    P.edge_bc = {(int)m::BCType::FullyFixed, (int)m::BCType::HorizontallyFixed,
                 (int)m::BCType::Free, (int)m::BCType::HorizontallyFixed};
    pr.polygons.push_back(P);

    pr.initial.name = "K0";
    pr.initial.duration = 0.0;
    m::Phase wait;
    wait.name = "Creep";
    wait.duration = days;
    wait.time_steps = 50;
    pr.phases.push_back(wait);
    return pr;
}

m::Project build_ssc() { return build_ssc_at(kScMu, kScDays, true); }

double ssc_run(const m::Project& pr) {
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return -1.0;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[0].ok || !res[1].ok) return -1.0;
    return res[1].max_disp;
}

void oracle_ssc(const m::Project& pr) {
    const double u = ssc_run(pr);
    check(u > 0.0, "the K0 and creep phases converged");
    if (u <= 0.0) return;

    // (a) The published law at the file's own duration.
    const double cf = ssc_creep(kScMu, kScDays);
    std::printf("      creep at %.0f days: closed form %.6f m | file run %.6f m (%+.2f%%)\n",
                kScDays, cf, u, 100.0 * (u - cf) / cf);
    check(std::fabs(u - cf) < 0.03 * cf, "creep settlement within 3% of mu* ln(1 + t/tau)");

    // (b) The law is LOGARITHMIC, and one point cannot show that. Three more durations spanning
    // three decades do: a linear creep law fitted through any one of them would miss the others
    // by a factor of ten. The spread also locates tau -- at t = tau the settlement is mu* ln 2,
    // and it is the 24-hour oedometer stage that fixes tau at one day (Eq 11-13/14).
    for (double d : {1.0, 10.0, 1000.0}) {
        const double uu = ssc_run(build_ssc_at(kScMu, d, true));
        const double c = ssc_creep(kScMu, d);
        std::printf("      t = %6.0f d: run %.6e | closed form %.6e (%+.2f%%)\n", d, uu, c,
                    100.0 * (uu - c) / c);
        check(uu > 0.0 && std::fabs(uu - c) < 0.03 * c, "the creep law holds across three decades of time");
    }

    // (c) The settlement is linear in mu*, the one parameter this model adds.
    const double u2 = ssc_run(build_ssc_at(2.0 * kScMu, kScDays, true));
    std::printf("      mu* x2: %.6f m (closed form %.6f)\n", u2, ssc_creep(2.0 * kScMu, kScDays));
    check(u2 > 0.0 && std::fabs(u2 - ssc_creep(2.0 * kScMu, kScDays)) < 0.03 * ssc_creep(2.0 * kScMu, kScDays),
          "doubling mu* doubles the creep, as the closed form says");

    // (d) The differential witness: the same file with the SAME ground as plain Soft Soil --
    // which has every feature of this model except the creep -- must sit still for the same
    // hundred days. If it moved, the settlement above would be something other than creep.
    const double u_ss = ssc_run(build_ssc_at(kScMu, kScDays, false));
    std::printf("      plain Soft Soil for the same 100 days: %.3e m (creep run %.3e, ratio %.0f)\n",
                u_ss, u, u / std::fmax(u_ss, 1e-15));
    check(u_ss >= 0.0 && u_ss < 0.01 * u, "without the creep model, time alone moves nothing");
}

// ------------------------------ KV-STR-004: axial capacity of a pile row --------------------
// The embedded beam's first case from a `.k2d` file, and the manual chooses the loading path
// for us: "embedded beams are not meant to be used as laterally loaded piles and will therefore
// not show accurate failure loads when subjected to transverse forces" (Reference Manual sec
// 6.6.4), and the material data set carries "only the bearing capacity" -- skin and base. So the
// defining quantity is the AXIAL capacity, and it has an exact closed form.
//
// Why the case exists, and what building it found. The parity register had the embedded beam as
// "implemented, unverified", and everything the verification touched turned out to be wrong:
//   * Eq 6-65's division by L_spacing was not applied at all, so every skin and foot spring was
//     2.5x too stiff at the default spacing. The function's only consumer is the driver -- the
//     element test passes its stiffnesses by hand -- so the constant was wrong by a factor of
//     2.5 while all 150 tests were green.
//   * The foot used D/2 where Eq 6-67 defines R_eq = sqrt(12 EI/EA)/2, which for a solid
//     circular pile is 0.433 D, not 0.5 D.
//   * The connection point was always FREE. PLAXIS connects it HINGED to the soil node when no
//     structure shares it, and a point load is carried by the nearest SOIL node -- so a pile
//     row could not be loaded at its head at all. Check (e) is the sentry for that one.
//
// Fixture, and the two rules behind it. The soil is MOHR-COULOMB and not Linear Elastic because
// PLAXIS ignores the shaft resistance and the spacing inside a linear elastic cluster (it counts
// that as structure, not soil), so an LE fixture would measure the one case PLAXIS treats
// differently. Its cohesion is set far above anything the run mobilises, so the plateau measured
// is the PILE's declared capacity and not a soil bearing failure. Soil and pile are weightless,
// so the load carried is the load applied and nothing else.
constexpr double kPlL = 10.0;        // embedded length [m]
constexpr double kPlD = 0.4;         // pile diameter [m]
constexpr double kPlLs = 2.5;        // out-of-plane spacing [m]
constexpr double kPlE = 3.0e7;       // pile stiffness [kN/m2]
constexpr double kPlTskin = 100.0;   // skin resistance cap, per pile [kN/m]
constexpr double kPlFbase = 500.0;   // base resistance, per pile [kN]
constexpr double kPlLoad = 1500.0;   // head load [kN/m of wall] -- far above the capacity
constexpr double kPlH = 1.0;         // element size [m]
constexpr double kPlTop = 16.0, kPlW = 16.0, kPlX = 8.0;

// The closed form, written out here rather than read from the driver: the ultimate axial load of
// a pile row is its skin resistance over the embedded length plus its base resistance, and the
// row carries it per metre of WALL, so every per-pile quantity is divided by the spacing --
// exactly as EA, EI and the pile weight are (PLAXIS Reference sec 6.6.3, Eq 6-65's own reason).
double pile_capacity(double Tskin, double L, double Fbase, double Ls) {
    return (Tskin * L + Fbase) / Ls;
}

m::Project build_pile_at(double Tskin, double Fbase, double Ls, double load, int conn, double h) {
    m::Project pr;
    pr.name = "KV-STR-004 axial capacity of a pile row";
    pr.x_min = 0; pr.x_max = kPlW; pr.y_min = 0; pr.y_max = kPlTop;
    pr.has_water = false;                   // dry, and said so rather than left to a default
    pr.mesh.elem_size = h; pr.mesh.order = 6;

    m::Material s;
    s.name = "Sand";
    s.model = m::SoilModel::MohrCoulomb;
    s.gamma_unsat = s.gamma_sat = 0.0;      // weightless: the applied load is the only load
    s.E = 30000.0; s.nu = 0.3;
    s.c = 1000.0;                           // far above anything mobilised: the soil never fails
    s.phi = 30.0; s.psi = 0.0;
    pr.materials.push_back(s);

    m::EmbeddedBeamMaterial em;
    em.name = "Bored pile";
    em.E = kPlE; em.gamma = 0.0; em.diameter = kPlD; em.Lspacing = Ls;
    em.Tskin_max = Tskin; em.Fmax_base = Fbase;
    pr.embedded.push_back(em);

    m::SoilPolygon P;
    P.name = "Soil"; P.material = 0;
    P.x = {0, kPlW, kPlW, 0};
    P.y = {0, 0, kPlTop, kPlTop};
    P.edge_bc = {4, 2, 0, 2};               // base fixed, sides on rollers, surface free
    pr.polygons.push_back(P);

    m::StructElement S;
    S.kind = m::StructKind::EmbeddedBeam;
    S.name = "Pile row";
    S.x1 = kPlX; S.y1 = kPlTop; S.x2 = kPlX; S.y2 = kPlTop - kPlL;
    S.material = 0;
    S.conn = conn;                          // 0 hinged (PLAXIS's default), 1 free
    pr.structs.push_back(S);

    m::Load F;
    F.kind = m::LoadKind::Point;
    F.name = "Head load";
    F.x1 = kPlX; F.y1 = kPlTop; F.x2 = kPlX; F.y2 = kPlTop;
    F.qx1 = F.qx2 = 0.0; F.qy1 = -load; F.qy2 = 0.0;
    pr.loads.push_back(F);

    pr.initial.poly_active = {1};
    pr.initial.struct_active = {1};
    pr.initial.load_active = {0};
    m::Phase load_phase;
    load_phase.name = "Head load";
    load_phase.poly_active = {1};
    load_phase.struct_active = {1};
    load_phase.load_active = {1};
    pr.phases.push_back(load_phase);
    return pr;
}

m::Project build_pile() {
    return build_pile_at(kPlTskin, kPlFbase, kPlLs, kPlLoad, 0, kPlH);
}

// One run: the axial force the pile carries AT ITS HEAD (the station with the largest y), and
// the head settlement. The force diagram is per metre of wall, like the beam's own EA/EI.
struct PileRead {
    bool ok = false;
    double N_head = 0.0, u_head = 0.0, max_u = 0.0;
};
PileRead read_pile(const m::Project& pr) {
    PileRead r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return r;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[0].ok || !res[1].ok) return r;
    const auto& R = res[1];
    r.max_u = R.max_disp;
    // An embedded beam reports its diagram AS A PLATE would (kind 0: it produces N/Q/M), so the
    // element is found by name rather than by kind -- the model has exactly one structure.
    for (const auto& sf : R.struct_forces) {
        if (sf.name != "Pile row" || sf.stations.empty()) continue;
        const katai::core::ForceStation* head = &sf.stations.front();
        for (const auto& st : sf.stations)
            if (st.y > head->y) head = &st;
        r.N_head = std::fabs(head->N);
        r.u_head = -head->uy;
        r.ok = true;
    }
    return r;
}

void oracle_pile(const m::Project& pr) {
    const double cap = pile_capacity(kPlTskin, kPlL, kPlFbase, kPlLs);
    const PileRead r = read_pile(pr);
    check(r.ok, "pile row solves and reports a force diagram");
    if (!r.ok) return;
    std::printf("      N at the pile head: %.4f kN/m  (T_max L + F_max)/Ls = %.4f  (%+.2f%%)\n",
                r.N_head, cap, 100.0 * (r.N_head / cap - 1.0));
    std::printf("      head settlement %.6f m, max|u| %.6f m\n", r.u_head, r.max_u);
    // The head load is FAR above the capacity, so the pile is fully mobilised and what it carries
    // is its capacity, no matter how hard the head is pushed. The soil takes the remainder -- the
    // plateau is in the PILE's force, not in the total applied load, because the load is applied
    // at a soil node that the hinged head shares.
    check(std::fabs(r.N_head / cap - 1.0) < 0.02, "the pile carries exactly its declared capacity");

    // (a) PLATEAU. Doubling the load must not change what the pile carries: a limit load, not a
    // stiffness reading. This is the same instrument as the sliding block's plateau.
    const PileRead r2 = read_pile(build_pile_at(kPlTskin, kPlFbase, kPlLs, 2.0 * kPlLoad, 0, kPlH));
    std::printf("      load x2: N_head %.4f kN/m (settlement %.6f m)\n", r2.N_head, r2.u_head);
    check(r2.ok && std::fabs(r2.N_head - r.N_head) < 0.01 * cap,
          "load x2 leaves the pile force unchanged: a plateau, not a stiffness");
    check(r2.u_head > 1.2 * r.u_head, "and the head goes on settling, so the run is past the limit");

    // (b) THE TWO TERMS, SEPARATELY. Skin and base are independent inputs and each must move the
    // capacity by exactly its own share -- one number agreeing could be a coincidence of two
    // errors, two numbers moving independently cannot.
    // Note the schema's convention, which a first draft of this check walked straight into: ZERO
    // means UNLIMITED for both caps, not zero capacity, so the base is moved between two non-zero
    // values rather than switched off (with Fmax_base = 0 the pile read 724 kN/m -- an unbounded
    // base, exactly as documented).
    const PileRead rs = read_pile(build_pile_at(kPlTskin, 100.0, kPlLs, kPlLoad, 0, kPlH));
    const PileRead rb = read_pile(build_pile_at(50.0, kPlFbase, kPlLs, kPlLoad, 0, kPlH));
    const double cap_s = pile_capacity(kPlTskin, kPlL, 100.0, kPlLs);
    const double cap_b = pile_capacity(50.0, kPlL, kPlFbase, kPlLs);
    std::printf("      base 500->100: %.4f vs %.4f | skin halved: %.4f vs %.4f kN/m\n",
                rs.N_head, cap_s, rb.N_head, cap_b);
    check(rs.ok && std::fabs(rs.N_head / cap_s - 1.0) < 0.03, "the base term moves by its own share");
    check(rb.ok && std::fabs(rb.N_head / cap_b - 1.0) < 0.03, "the skin term moves by its own share");

    // (c) SPACING. Everything about a row is per metre of wall, so doubling the out-of-plane
    // spacing must halve the capacity EXACTLY. This is the check that the /L_spacing of Eq 6-65
    // reaches the capacities as well as the stiffnesses.
    const PileRead rl = read_pile(build_pile_at(kPlTskin, kPlFbase, 2.0 * kPlLs, kPlLoad, 0, kPlH));
    const double cap_l = pile_capacity(kPlTskin, kPlL, kPlFbase, 2.0 * kPlLs);
    std::printf("      spacing x2: %.4f vs %.4f kN/m (ratio to base case %.4f)\n",
                rl.N_head, cap_l, rl.N_head / r.N_head);
    check(rl.ok && std::fabs(rl.N_head / cap_l - 1.0) < 0.03, "twice the spacing, half the capacity");

    // (d) MESH INDEPENDENCE, and a prediction of this test's own that the measurement refuted.
    // The draft expected a residual: at a hinged connection the skin point ON the tied node
    // cannot slip -- both sides of that joint are the same degree of freedom -- so its
    // Newton-Cotes share of the skin capacity looked like it could never mobilise, which would
    // leave a deviation shrinking with the element size, the sibling of the beam case's q h^2/12.
    // It does not happen: the capacity is EXACT at both densities. A limit load carried by
    // caps rather than by stiffnesses has nothing left for the discretisation to bias, and the
    // load simply redistributes to the points that can still take it. Recorded because a
    // prediction that failed is worth as much as one that held.
    const PileRead rf = read_pile(build_pile_at(kPlTskin, kPlFbase, kPlLs, kPlLoad, 0, 0.5 * kPlH));
    std::printf("      h/2: N_head %.4f kN/m (%+.2f%%) vs h: %+.2f%%\n", rf.N_head,
                100.0 * (rf.N_head / cap - 1.0), 100.0 * (r.N_head / cap - 1.0));
    check(rf.ok && std::fabs(rf.N_head / cap - 1.0) < 0.02,
          "halving the element size leaves the capacity where it was: no mesh dependence");

    // (e) THE SENTRY. With the connection FREE -- which is what this engine did for every pile
    // until 2026-08-13, with no way to ask for anything else -- the head load goes into the soil
    // beside the pile and reaches it only through the skin springs. The pile then carries a small
    // fraction of its capacity and the head settles several times as far. If the hinged default
    // is ever taken back, this does not drift: it collapses.
    const PileRead rfree = read_pile(build_pile_at(kPlTskin, kPlFbase, kPlLs, kPlLoad, 1, kPlH));
    std::printf("      connection FREE: N_head %.4f kN/m (%.1f%% of capacity), settlement %.6f m "
                "vs hinged %.6f m\n",
                rfree.N_head, 100.0 * rfree.N_head / cap, rfree.u_head, r.u_head);
    check(rfree.ok && rfree.N_head < 0.5 * cap,
          "a FREE head cannot deliver the load to the pile: the fault this case was built to find");
}

// ------------------------------ KV-STR-005: a geogrid's axial force and tension cut-off ------
// The geogrid's first case from a `.k2d` file. The manual's definition IS the oracle: "the axial
// stiffness is the ratio of the axial force F per unit width and the axial strain
// (eps = dl/l)", EA = F/eps (Reference Manual Eq. 6-51), and "geogrids can only sustain tensile
// forces, but not compressive forces" (sec 6.5).
//
// The fixture turns that definition into something a file can ask for. A geogrid's translational
// degrees of freedom ARE the soil's -- the line is a conforming chain of mesh nodes with no
// rotation -- so if the soil is made to strain uniformly, the geogrid strains with it, exactly.
// A homogeneous block held at u_x = 0 on one side and pulled to u_x = D on the other deforms
// affinely: u_x = D x/L, and every horizontal fibre has eps = D/L. The geogrid spans the full
// width, so its ends sit ON the two boundaries and its elongation is D whatever the soil does;
// its tension is then constant along its length, which means it applies no body force to the
// soil and the affine field stays the exact solution. N = EA D/L is therefore not an
// approximation of this problem -- it is this problem.
constexpr double kGgL = 10.0;      // block width [m] = the geogrid's length
constexpr double kGgH = 5.0;       // block height [m]
constexpr double kGgY = 2.0;       // the geogrid's level [m]
constexpr double kGgEA = 5000.0;   // axial stiffness [kN/m]
constexpr double kGgNp = 3.0;      // tension cut-off [kN/m]
constexpr double kGgD = 0.004;     // imposed stretch [m] -> eps = 4e-4, N = 2.0 kN/m (elastic)
constexpr double kGgHm = 0.5;      // element size [m]

// Eq. 6-51 rearranged, with the manual's tension cut-off applied on top, written out here rather
// than called from the geogrid header.
double geogrid_force(double EA, double stretch, double L, double Np) {
    const double N = EA * stretch / L;
    if (N <= 0.0) return 0.0;                      // tension only: compression carries nothing
    return (Np > 0.0 && N > Np) ? Np : N;
}

m::Project build_geogrid_at(double EA, double Np, double stretch, bool plastic, bool with_grid,
                            double h) {
    m::Project pr;
    pr.name = "KV-STR-005 geogrid axial force and tension cut-off";
    pr.x_min = 0; pr.x_max = kGgL; pr.y_min = 0; pr.y_max = kGgH;
    pr.has_water = false;
    pr.mesh.elem_size = h; pr.mesh.order = 6;

    m::Material s;
    s.name = "Fill";
    s.model = m::SoilModel::LinearElastic;
    s.gamma_unsat = s.gamma_sat = 0.0;     // weightless: the imposed stretch is the only action
    s.E = 10000.0; s.nu = 0.3;
    pr.materials.push_back(s);

    m::GeogridMaterial gm;
    gm.name = "Reinforcement";
    gm.elastoplastic = plastic;
    gm.EA = EA; gm.Np = Np;
    pr.geogrids.push_back(gm);

    m::SoilPolygon P;
    P.name = "Fill"; P.material = 0;
    P.x = {0, kGgL, kGgL, 0};
    P.y = {0, 0, kGgH, kGgH};
    // Bottom held vertically and the left side horizontally: enough to remove the rigid body
    // motions and nothing more, so the block is free to contract laterally and the strain field
    // stays affine. The right side is driven by the prescribed displacement below.
    P.edge_bc = {3, 0, 0, 2};
    pr.polygons.push_back(P);

    if (with_grid) {
        m::StructElement S;
        S.kind = m::StructKind::Geogrid;
        S.name = "Reinforcement";
        S.x1 = 0.0; S.y1 = kGgY; S.x2 = kGgL; S.y2 = kGgY;
        S.material = 0;
        pr.structs.push_back(S);
    }

    m::PrescribedDisp D;
    D.name = "Stretch";
    D.x1 = kGgL; D.y1 = 0.0; D.x2 = kGgL; D.y2 = kGgH;
    D.set_ux = true;  D.ux = stretch;
    D.set_uy = false; D.uy = 0.0;
    pr.disps.push_back(D);

    const std::vector<char> on(with_grid ? 1u : 0u, 1);
    pr.initial.poly_active = {1};
    pr.initial.struct_active = on;
    pr.initial.disp_active = {0};
    m::Phase pull;
    pull.name = "Stretch";
    pull.poly_active = {1};
    pull.struct_active = on;
    pull.disp_active = {1};
    pr.phases.push_back(pull);
    return pr;
}

m::Project build_geogrid() {
    return build_geogrid_at(kGgEA, kGgNp, kGgD, false, true, kGgHm);
}

// One run: the geogrid's axial force (its stations are all at the same tension here, so the
// extreme is the value) and the horizontal reaction the prescribed edge has to supply.
struct GgRead {
    bool ok = false;
    double N = 0.0, Rx = 0.0;
};
GgRead read_geogrid(const m::Project& pr) {
    GgRead r;
    const auto M = katai::app::mesh_from_project(pr);
    if (!M.ok) return r;
    const auto res = katai::app::solve_phases(pr, M.mesh,
                                              katai::app::initial_phase_from(pr.initial_procedure));
    if (res.size() != 2 || !res[0].ok || !res[1].ok) return r;
    const auto& R = res[1];
    if (std::getenv("KATAI_GG_DUMP")) {
        std::vector<std::pair<double, double>> nu;
        for (int n = 0; n < R.mesh.node_count; ++n)
            if (std::fabs(R.mesh.y[n] - kGgY) < 1e-9 && R.mesh.x[n] > 9.4)
                nu.push_back({R.mesh.x[n], R.disp[2 * n]});
        std::sort(nu.begin(), nu.end());
        std::printf("      [%s] u_x:", pr.structs.empty() ? "bare " : "grid ");
        for (const auto& p : nu)
            std::printf(" [x=%.4f u=%.4e aff=%.4e]", p.first, p.second, kGgD * p.first / kGgL);
        std::printf("\n");
    }
    for (const auto& sf : R.struct_forces)
        if (sf.name == "Reinforcement" && !sf.stations.empty()) {
            if (std::getenv("KATAI_GG_DUMP")) {
                std::printf("      stations (%zu):", sf.stations.size());
                for (size_t i = 0; i < sf.stations.size(); ++i)
                    if (i < 4 || i + 4 >= sf.stations.size() || i == sf.stations.size() / 2)
                        std::printf(" [x=%.2f N=%.4f]", sf.stations[i].x, sf.stations[i].N);
                std::printf("\n");
                // The nodal field along the reinforcement near the driven end, against the affine
                // solution u_x = D x/L. This separates a wrong DISPLACEMENT from a wrong FORCE
                // RECOVERY: only one of the two can be the fault.
                std::vector<std::pair<double, double>> nu;
                for (int n = 0; n < R.mesh.node_count; ++n)
                    if (std::fabs(R.mesh.y[n] - kGgY) < 1e-9 && R.mesh.x[n] > 9.4)
                        nu.push_back({R.mesh.x[n], R.disp[2 * n]});
                std::sort(nu.begin(), nu.end());
                std::printf("      u_x along the fibre:");
                for (const auto& p : nu)
                    std::printf(" [x=%.4f u=%.3e aff=%.3e]", p.first, p.second,
                                kGgD * p.first / kGgL);
                std::printf("\n");
                // The DRIVEN EDGE itself: every node on x = L was told u_x = D. Any node there
                // that did not get it is a hole in the prescribed displacement, and the field
                // beside it has to bend around the hole.
                std::vector<std::pair<double, double>> ed;
                for (int n = 0; n < R.mesh.node_count; ++n)
                    if (std::fabs(R.mesh.x[n] - kGgL) < 1e-9) ed.push_back({R.mesh.y[n], R.disp[2 * n]});
                std::sort(ed.begin(), ed.end());
                int missed = 0;
                for (const auto& p : ed) if (std::fabs(p.second - kGgD) > 1e-12) ++missed;
                std::printf("      driven edge: %zu nodes, %d did NOT receive u_x = D:", ed.size(),
                            missed);
                for (const auto& p : ed)
                    if (std::fabs(p.second - kGgD) > 1e-12)
                        std::printf(" [y=%.4f u=%.3e]", p.first, p.second);
                std::printf("\n");
            }
            r.N = sf.stations[sf.stations.size() / 2].N;   // mid-span: away from the ends
        }
    if (R.reaction.size() == 2u * (size_t)R.mesh.node_count)
        for (int n = 0; n < R.mesh.node_count; ++n)
            if (std::fabs(R.mesh.x[n] - kGgL) < 1e-9) r.Rx += R.reaction[2 * n];
    r.ok = true;
    return r;
}

void oracle_geogrid(const m::Project& pr) {
    const double want = kGgEA * kGgD / kGgL;              // Eq. 6-51, elastic (below N_p)
    const GgRead r = read_geogrid(pr);
    const GgRead rn0 = read_geogrid(build_geogrid_at(kGgEA, kGgNp, kGgD, false, false, kGgHm));
    check(r.ok, "the reinforced block solves");
    if (!r.ok) return;
    // The soil's own answer first, because it is what makes the geogrid's readable: with no
    // reinforcement the edge reaction must be the plane-strain closed form E/(1-nu^2) eps H
    // exactly. If it is, the strain field IS affine and the geogrid's fibre strain is D/L.
    const double Rx_soil = 10000.0 / (1.0 - 0.3 * 0.3) * (kGgD / kGgL) * kGgH;
    std::printf("      unreinforced edge reaction %.6f kN/m vs E/(1-nu^2) eps H = %.6f (%+.4f%%)\n",
                std::fabs(rn0.Rx), Rx_soil, 100.0 * (std::fabs(rn0.Rx) / Rx_soil - 1.0));
    check(rn0.ok && std::fabs(std::fabs(rn0.Rx) / Rx_soil - 1.0) < 0.005,
          "the bare block reproduces the affine closed form, so the fibre strain is D/L");

    // The chain the driver turns into quadratic elements. It is built by collecting every mesh
    // node ON the line and sorting them along it, and the driver then reads them as
    // corner, mid, corner, mid, ... -- {chain[2e], chain[2e+2], chain[2e+1]}. Nothing checks
    // that assumption, so it is checked here: every element's middle node must be the midpoint
    // of its two ends, or the element is not the element the code thinks it is.
    {
        const auto MM = katai::app::mesh_from_project(pr);
        const auto chain = katai::mesh::collect_chain(MM.mesh, 0.0, kGgY, kGgL, kGgY);
        double worst = 0.0; int worst_e = -1;
        for (size_t e = 0; 2 * e + 2 < chain.size(); ++e) {
            const int a = chain[2 * e], b = chain[2 * e + 2], mid = chain[2 * e + 1];
            const double mx = 0.5 * (MM.mesh.x[a] + MM.mesh.x[b]);
            const double my = 0.5 * (MM.mesh.y[a] + MM.mesh.y[b]);
            const double d = std::hypot(MM.mesh.x[mid] - mx, MM.mesh.y[mid] - my);
            if (d > worst) { worst = d; worst_e = (int)e; }
        }
        std::printf("      chain: %zu nodes, %s; worst mid-node offset %.3e m at element %d\n",
                    chain.size(), chain.size() % 2 ? "odd (as required)" : "EVEN",
                    worst, worst_e);
        if (worst_e >= 0 && worst > 1e-9) {
            const int a = chain[2 * worst_e], b = chain[2 * worst_e + 2], mid = chain[2 * worst_e + 1];
            std::printf("        element %d: A x=%.4f  mid x=%.4f  B x=%.4f\n", worst_e,
                        MM.mesh.x[a], MM.mesh.x[mid], MM.mesh.x[b]);
        }
        check(chain.size() % 2 == 1 && worst < 1e-9,
              "every element on the chain has its middle node at its own midpoint");
    }
    std::printf("      geogrid N: %.6f kN/m   EA eps = %.6f   (%+.3f%%)\n",
                r.N, want, 100.0 * (r.N / want - 1.0));
    check(std::fabs(r.N / want - 1.0) < 0.01, "the geogrid carries EA times the imposed strain");

    // (a) LINEAR IN BOTH FACTORS, SEPARATELY. Eq. 6-51 has two inputs and each must move the
    // force by exactly its own factor; one agreeing number could be two errors cancelling.
    const GgRead r2 = read_geogrid(build_geogrid_at(kGgEA, kGgNp, 2.0 * kGgD, false, true, kGgHm));
    const GgRead re = read_geogrid(build_geogrid_at(2.0 * kGgEA, kGgNp, kGgD, false, true, kGgHm));
    std::printf("      stretch x2: %.6f (x%.4f) | EA x2: %.6f (x%.4f) kN/m\n",
                r2.N, r2.N / r.N, re.N, re.N / r.N);
    check(r2.ok && std::fabs(r2.N / r.N - 2.0) < 0.01, "twice the strain, twice the force");
    check(re.ok && std::fabs(re.N / r.N - 2.0) < 0.01, "twice the stiffness, twice the force");

    // (b) THE CUT-OFF IS A PLATEAU. Past N_p the force stops at N_p however hard the block is
    // pulled -- a limit, not a stiffness. The elastic twin at the same stretch shows what it
    // would have carried without the cut-off, so the cut-off is measured against something.
    const double big = 4.0 * kGgD;                       // EA eps = 8 kN/m against N_p = 3
    const GgRead rp = read_geogrid(build_geogrid_at(kGgEA, kGgNp, big, true, true, kGgHm));
    const GgRead rp2 = read_geogrid(build_geogrid_at(kGgEA, kGgNp, 2.0 * big, true, true, kGgHm));
    const GgRead rel = read_geogrid(build_geogrid_at(kGgEA, kGgNp, big, false, true, kGgHm));
    std::printf("      cut-off: N %.6f and %.6f kN/m at 1x and 2x the stretch (N_p = %.3f); "
                "elastic twin %.6f\n", rp.N, rp2.N, kGgNp, rel.N);
    check(rp.ok && std::fabs(rp.N / kGgNp - 1.0) < 0.01, "past N_p the geogrid carries N_p");
    check(rp2.ok && std::fabs(rp2.N - rp.N) < 0.01 * kGgNp, "and doubling the stretch changes nothing");
    check(rel.ok && rel.N > 2.0 * kGgNp, "while the elastic twin carries what it was told to");

    // (c) TENSION ONLY. The manual's other sentence about this element: it "can only sustain
    // tensile forces, but not compressive forces". Push the block instead of pulling it and the
    // reinforcement must carry nothing at all -- not a small number, nothing.
    const GgRead rc = read_geogrid(build_geogrid_at(kGgEA, kGgNp, -kGgD, false, true, kGgHm));
    std::printf("      block compressed: geogrid N = %.3e kN/m\n", rc.N);
    check(rc.ok && rc.N < 1e-6 * want, "in compression the geogrid carries nothing");

    // (d) THE REINFORCEMENT MUST NOT DISTURB THE SOIL AT ALL. A member carrying CONSTANT tension
    // has zero force divergence, so it loads nothing but its own two ends -- both of which are
    // held here. The soil's stress field, and with it the reaction its supports carry, must
    // therefore be bit-for-bit what the bare block gives. That is a sharper statement than the
    // draft of this check attempted: it tried to read the geogrid's share OUT of the reaction,
    // which `SolveResult.reaction` cannot supply, because v1 sums the SOIL contribution only and
    // says so ("a structural end force landing on a fixed node is not included in v1"). The
    // measurement then proves the affine argument the whole case rests on, rather than merely
    // agreeing with it.
    const double dR = std::fabs(r.Rx) - std::fabs(rn0.Rx);
    std::printf("      soil reaction with the grid %.6f / without %.6f -> difference %.3e kN/m\n",
                std::fabs(r.Rx), std::fabs(rn0.Rx), dR);
    check(std::fabs(dR) < 1e-9 * Rx_soil,
          "a constant-tension member leaves the soil field exactly where it found it");
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
        {"kv-cst-008-hssmall-unloading.k2d", build_hss, oracle_hss},
        {"kv-cst-009-soft-soil-oedometer.k2d", build_ss, oracle_ss},
        {"kv-cst-010-soft-soil-creep-column.k2d", build_ssc, oracle_ssc},
        {"kv-str-004-axial-pile-capacity.k2d", build_pile, oracle_pile},
        {"kv-str-005-geogrid-tension.k2d", build_geogrid, oracle_geogrid},
    };
    for (const CorpusCase& c : cases) run_case(c);

    if (g_failures == 0) {
        std::printf("\nOK: every corpus case is reproducible from its checked-in .k2d\n");
        return 0;
    }
    std::fprintf(stderr, "\n%d check(s) failed\n", g_failures);
    return 1;
}
