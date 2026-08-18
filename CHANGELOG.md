# Changelog

All notable changes to KATAI 2D. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
MAJOR.MINOR.PATCH.

## [0.8.1] - 2026-08-19

A correctness fix in what the results panel and the calculation report SHOW, a
typesetting pass over the notation they show it in, and the maintainer's rule
about naming other programs turned into a test.

### Fixed

- **A field nothing had computed was drawn, tabulated and printed as a result.**
  S, the degree of saturation, is produced by the unsaturated flow — a transient
  or a fully coupled phase — and by nothing else. Every other run leaves the
  array empty and the field reader answers 1.0 when asked anyway, and nothing
  filtered the ask: an ordinary static analysis offered S in the field selector,
  painted a uniformly saturated model, tabulated "min 1.000, max 1.000", and
  printed a saturation row in both reports — in the document that goes to a
  client, beside numbers the solver did produce. Above a phreatic surface that
  is not merely uncomputed, it is the one answer that is certainly wrong. A
  field is now offered only where the run produced it.
- **`tau_max` said more than it computed.** It is the radius of the Mohr circle
  of the IN-PLANE stresses; the out-of-plane sigma_zz is not carried in the
  nodal field, so wherever it falls outside the in-plane pair — which
  plasticity can do — the true maximum shear is larger than the number shown.
  The field is now named `τmax (in-plane)`.
- **The report's columns fanned out.** `printf` pads by bytes, and the typeset
  names are not one byte per column; the extreme-value table is now padded by
  characters, which also fixes a Turkish material or load name having done the
  same thing since long before this release.

### Changed

- **The notation is typeset.** The interface and the reports wrote symbols the
  way a source file has to write them — `sigma'_yy`, `phi'`, `lambda*`,
  `kN/m3`, `>=`, `[deg]`. They now read `σ′yy`, `φ′`, `λ*`, `kN/m³`, `≥`, `[°]`:
  225 pieces of user-visible text across the Studio and both report languages.
  The UI font loads Greek, the prime, super/subscripts and the mathematical
  operators to make it possible, measured against the font before it was used.
  Indices stay plain letters on purpose — Unicode has no subscript y and no
  subscript w, and subscripting most of a set while three of it cannot be is
  worse than subscripting none.
- **No other FE program is named in anything a user reads.** Twenty-two
  Studio tooltips and ten Python docstrings named one; every one of them was
  making a point that stands without the brand, and every one of them still
  makes it. Source comments and the validation record are deliberately out of
  scope: where a default came from is worth recording, and a comparison is
  worthless without naming what was compared against.

### Added

- `test_product_text` (and `test_studio_product_text`): one scanner over both
  trees, reading string literals and docstrings only — it lexes past comments
  rather than splitting lines on `//`, so it sees what a user sees. Checked
  against a fixture that must fail it.
- Four checks pinning the saturation rule from both sides, and one that
  measures the report's column alignment in characters.

## [0.8.0] - 2026-08-17

Verification release, and the first one in which the Python surface reaches the
whole engine. Six capabilities — the interface, the plate, the geogrid, the
embedded beam, HS small, and Soft Soil with its creep — were implemented and
tested at the element, but had never been run along the path a user actually
takes: from a file, through the mesher and the driver. Run that way for the
first time, **five of the six were wrong**, and one further cross-cutting fault
(the tension cut-off) came out of the same work. Every closure arrives with a
benchmark input checked in and a citation to the manual clause it implements:
the verification record grows from 45 declared cases over 17 input files to
**57 cases over 26 files**, and the suite to **152 tests**.

The project file moves from version 12 to version 14. This build reads every
older file unchanged; an older build refuses these, which is the point of the
guard — each bump marks an input an older build would have dropped in silence.

### Added

- **The Python surface reaches the structural elements.** `prj.structures`
  creates walls, anchors, geogrids, pile rows and bare interfaces, one call
  each, and each returns a handle that phases activate and deactivate exactly
  like a load or a region — which is the whole of staged construction. Before
  this, a script that needed a wall had to drop to the raw bindings and assemble
  the material list, the element list and the index between them by hand, so the
  pleasant layer covered the tutorial and gave up at the first real job.
  Interfaces come from the wall that needs them (`interfaces="both"`), a wall can
  be made a groundwater screen, and a pile states how its head attaches — hinged
  by default, as PLAXIS does.
- **What a structure carries is readable from a script.** `ForceStation` and
  `StructForce` are bound: the stations along an element with their arc length,
  position, N, Q, M and displacement, plus the element's name, kind, yield flag
  and extremes. A wall's bending moment and a pile's axial force existed in C++
  and in the GUI and nowhere else — which is the surface a parameter study, a
  benchmark or a thesis actually lives in.
- **Flow phases, design codes and phase water from Python.**
  `prj.phases.transient_flow(...)` and `.fully_coupled(...)`; `design=` selects
  the EC7 or TBDY 2018 partial-factor catalogue for a single phase, so a
  characteristic run and its design check live in one file; `water=` overrides
  the table or the phreatic line for one phase, which is staged dewatering.
- **`katai.summary()` and `katai.extremes()`.** A readable account of a scripted
  run: the phases named, the extremes tabulated *with the place each occurs*,
  the structural force envelopes, and every diagnostic the engine raised —
  because those change how the numbers should be read. A field that is uniform
  says so instead of inventing a location.
- **"Reset small strain" as a phase option** (`.k2d` v14 `phases[].resetsmall`,
  diagnostic `K2D-M005`), after Material Models Manual sec. 7.6. A surcharge
  placed and removed to leave an overconsolidation behind also leaves a strain
  history, and in the real soil ageing erased that long before the analysis
  began. Verified by KV-CST-012 against an oracle established before the option
  existed: with the reset the run lands on KV-CST-008's fresh-K₀ answer to
  +0.0012%; without it, on a run 4.50× softer. On plain Hardening Soil, which
  has no history to reset, the flag is bit-for-bit inert and the run says so.
- **Li & Dafalias dilatancy below the phase-transformation line** for HS small
  (Material Models Manual sec. 7.9.1, Eq. 7-19…7-23), measured against an oracle
  written from the five equations and sharing no code with the kernel: worst
  difference 0.00e+00 over the branch. Recorded honestly, because the manual
  disagrees with itself here — Fig. 7-10 plots this function at 1.29× the
  amplitude Eq. 7-19 gives, and a formula printed in the specification outranks
  a constant reverse-engineered from a raster plot.
- **Numerical-uncertainty cases on axes other than the mesh** (KV-NUM-009…011).
  One of them measured a prediction wrong: the Newmark scheme was argued to be
  second order and the sweep recovered three, so the rule that "the algebra
  decides the order" is now bounded rather than assumed.

### Fixed

Each of these produced a converged run, a green suite and a wrong number.

- **An interface drawn along a fixed boundary was welded shut.** Both sides of
  the split sat at identical coordinates and boundary conditions are applied by
  coordinate, so both were fully fixed: the block sheared elastically against
  its own base instead of sliding, and every check in the run reported success.
  PLAXIS's own sliding-block case (Validation Manual V8 §3.3) returned
  **5,401,612 kN/m where the manual's arithmetic gives 60**. Which side holds the
  support is now decided rather than copied. KV-STR-002 from the checked-in file:
  59.7202 kN/m, and deleting the interface still returns 5.4e6 — that check is
  the sentry that fails loudly if the rule is ever undone.
- **Deactivating the soil welded the beams standing in it.** The manual builds
  its beam-bending case by removing the soil cluster; every node of the
  remaining beam then touched no active element and was pinned in both
  translations. The solve converged, reported "ok" and handed back
  max|u| = 0.000000e+00 with no diagnostic. A node a plate runs through is now
  exempt; an axial-only element keeps the fixity. KV-STR-003 reproduces the
  manual's 13.96 mm and 17.43 mm.
- **HS small rode the virgin backbone.** Masing's rule (Eq. 7-11,
  γ₀.₇,reloading = 2 γ₀.₇,virgin) was not applied, degrading the stiffness twice
  as fast — measured **+5.7% / +12.9% / +34.8%** too much heave at three
  unloading sizes, a deviation that grows with strain, which is the signature of
  a wrong threshold rather than of discretisation. Sec. 7.5's ceiling on
  E₀/E_ur was missing too; it is capped now and `K2D-M004` states the G₀ the run
  actually used.
- **The tension cut-off did not reach the models it belongs to.** `K2D-M001` had
  declared since 2026-08-08 that only the Mohr-Coulomb return read it, while the
  schema switches it on by default as PLAXIS does — so every Hardening Soil, HS
  small, Soft Soil and Soft Soil Creep run in this engine allowed tension past
  σ_t, a systematic difference from the reference code in the unsafe direction.
  KV-CST-011 pins it, including the identity that a cut-off the tension never
  reaches is bit-identical to no cut-off at all.
- **Every embedded-beam interface spring was 2.5× too stiff.** Reference Manual
  Eq. 6-65 divides all three springs by the out-of-plane spacing, as EA, EI, the
  weight and both capacities are divided; only Eq. 6-66's dimensionless factors
  were implemented. The foot also used D/2 where Eq. 6-67 defines
  R_eq = √(12 EI/EA)/2. A pile row could not be loaded at its head, which is why
  no case had ever run.
- **Structural elements read a driven node as standing still**, so a geogrid
  under a prescribed-displacement fixture carried nothing.
- **A prescribed-displacement line set to zero raised a warning it had no
  business raising** (`K2D-A003`): nothing is understated when the imposed value
  is zero, and that is the only way this schema can support a plate at a point.
- **K₀ = 0 is accepted.** With free vertical sides it is the *only* initial state
  in equilibrium with them; negative K₀ is still refused, with a message that
  says why.
- **The installer compared a real hash against the number 55.** Without
  `-UseBasicParsing`, Windows PowerShell 5.1 returns an `application/octet-stream`
  body — which is what GitHub serves a `.sha256` attachment as — as a *byte
  array*, so splitting it on whitespace split an array of bytes and element zero
  was the first byte of the hash, as a number. The integrity check would have
  "passed" for any file whose hash began with a matching byte. All web reads now
  go through one helper that decodes explicitly, the checksum is extracted by
  matching 64 hex characters, and a checksum file that does not contain 64 hex
  characters is a hard failure — a checksum that cannot be read is not a checksum
  that passed. TLS 1.2 is forced, since 5.1 still negotiates 1.0 on some builds.
- **`CITATION.cff` carried no DOI**, so GitHub's "Cite this repository" button
  produced an entry with nothing persistent in it — the one thing a citation
  exists to provide. The archived record had a DOI all along; the file never
  named it.
- **A verification case was asserting a result that only one solver backend
  could produce.** KV-STR-004's mesh-independence probe halved the element size
  under the fixture's 1500 kN/m overload. On MKL that converges and the pile
  carries its 600.0000 kN/m capacity; on the vendored Eigen solver — same
  source, same file, same mesh — the run reports a collapse mechanism after
  equilibrating 10% of the load. Continuous integration builds the `portable`
  composition, so it had been red since 2026-08-13 while every local run was
  green. The probe now runs at 1300 kN/m, inside the range where both backends
  return the capacity exactly (measured: 1200, 1300 and 1400 all give
  600.0000 at both densities), and re-runs the coarse density at the same load
  so the comparison is like-for-like. **The backend disagreement itself is not
  closed** — it is characterised and declared in
  `docs/validation/numerical-uncertainty.md` §10, together with what was tried
  and did not work.

### Changed

- The Hardening Soil dilatancy rule had been written out four times, once per
  return mapping. It is written once now, because a rule that is right in three
  copies and stale in the fourth is a silently wrong answer on whichever stress
  path reaches the fourth.
- `docs/references/hssmall-formulation.md` is in English throughout.

### Known limits, stated rather than implied

- HS small accumulates a monotone strain-history scalar and detects no reversal
  *inside* a phase; the manual refers that transformation to Benz (2006), which
  is a source to obtain rather than to paraphrase. Between phases,
  `resetsmall` is the manual's own remedy and is implemented.

## [0.7.1] - 2026-08-09

Packaging release. No engine change: the same solver, the same results, the same
verification record.

### Fixed

- The Windows executable carries its identity. It shipped with an empty version
  resource -- no product name, no version, no publisher, no copyright -- which is
  poor practice on its own and, on an unsigned and newly published binary, is
  also what a reputation heuristic reads as "anonymous": Microsoft Defender
  flagged `katai2d-0.7.0-win64.zip` as `Trojan:Win32/Wacatac.B!ml`, a generic
  machine-learning label, and removed the file. The executable now declares its
  product name, version, publisher, copyright and origin, generated from
  `version.hpp` so it cannot drift from what `katai info` prints. The rebuilt
  binary scans clean.

  This is a mitigation, not a guarantee, and the honest statement of the
  situation is in the README: an unsigned binary from a young project has no
  reputation, the published SHA-256 is how you verify what you downloaded, and
  the Python wheel -- which was never flagged -- is the alternative route.

## [0.7.0] — 2026-08-09

Capability release. Ten inputs a geotechnical model needs, and could not
express, are now in the file format; four of them closed a silently wrong
answer rather than a missing convenience. Every closure arrives with a
verification case and a citation to the manual clause it implements, and the
verification record grows from 27 to 45 declared cases.

The project file moves from version 8 to version 12. This build reads every
older file unchanged; an older build refuses these, which is the point of the
guard — each bump marks an input that an older build would have dropped in
silence and solved a different problem for.

### Engine — inputs that decide the answer

- **Per-material undrained stiffness** (`materials[].und_mode`, `nu_u`,
  `skempton_B`; `.k2d` v9). The pore fluid's bulk stiffness Kw/n was derived
  from a fixed nu_u = 0.495 for every undrained material — PLAXIS's default
  applied to users who had entered something else. Either the equivalent
  undrained Poisson ratio or Skempton's B is now a per-material input, related
  by the Material Models Manual's own equations. The same change fixed a
  measured silent error: for the Hardening Soil family, K' was read from the
  `E`/`nu` boxes that model never reads, sizing the pore fluid by an untouched
  default — a factor of 7.6 on an ordinary data set. It now follows the
  model's own unload/reload pair, and states the remaining limitation
  (`K2D-M002`).
- **Undrained (C)** (`materials[].drainage` = 4; `.k2d` v10). A total-stress
  analysis: undrained stiffness and strength, no pore pressure generated or
  carried, K0 on total stress. Available for the Linear Elastic and
  Mohr-Coulomb models, as in the manual; refused by name elsewhere, and in
  consolidation and fully-coupled phases. EC7 factors its cohesion as an
  undrained strength (gamma_cu).
- **Wells and drains** (`hydros`, `phases[].hydro`; `.k2d` v11). Dewatering
  from inside the model: a well prescribes a discharge and stops at `h_min`, a
  drain holds a head (normal one-sided, or vacuum). Both switch per phase. In
  consolidation and fully-coupled phases a drain sets the excess pore pressure
  to zero and a well is reported as not applied (`K2D-A009`); in a transient
  flow phase both are refused by name.
- **Cross permeability of walls and interfaces** (`structs[].flow_barrier`,
  `hyd_res`; `.k2d` v12). A cut-off wall can block flow: the groundwater
  calculation splits its own mesh along the barrier so the two sides carry
  separate pore-pressure degrees of freedom, exactly as the Scientific Manual
  describes. Impermeable, semi-permeable (hydraulic resistance d/k) or fully
  permeable — the last being the default and the previous behaviour.
- **Prescribed boundary flux** (`polygons[].edge_flux`, `Flux` in `edge_flow`;
  `.k2d` v6). Rainfall, infiltration and recharge boundaries. The kernel had
  the edge integral; no path from a file reached it, and no seepage solver took
  an external right-hand side.
- **Dilatancy cut-off** (`materials[].dilatancy_cutoff`, `e_max`; `.k2d` v5).
  A dilating soil stops dilating at its critical void ratio. Without it a dense
  sand dilates without limit and its bearing capacity comes out too high — an
  unsafe number, produced quietly. Off by default, as in PLAXIS.
- **Per-phase water conditions** (`phases[].water_override`, `wx`, `wy`;
  `.k2d` v4) and **anchor prestress** (`anchors[].prestress`; `.k2d` v3).
  Between them they make an anchored, dewatered excavation expressible: the
  pit is dewatered before it is dug, and the anchors are locked off rather than
  slack. Every anchored excavation built before this was the model of a weaker
  structure.
- **Staged-construction target and "ignore undrained behaviour"**
  (`phases[].mstage`, `ignoreund`; `.k2d` v8), and **per-phase numerical
  controls** (`tol`, `loadsteps`, `maxiter`; `.k2d` v7) so that a published
  run carries the numerics it was computed with.

### Diagnostics — an input may be used differently than drawn, never in silence

- `SolveResult.diagnostics` with stable codes and a written contract
  (`docs/diagnostics.md`), surfaced by the CLI, the Python package and the
  GUI alike. Sixteen input paths that were silently discarded — a load drawn
  above the surface, a wall the mesh never sees, a prescribed displacement
  that catches no node — are now refused or reported.
- Codes added in this release: `K2D-M002`, `K2D-M003` (what a material's
  parameters are read as), `K2D-A005`–`K2D-A011` (limits a result must be
  read under, including the mesh dependence of a factor of safety and the
  one-sided hand-over of a flow field with a barrier in it).
- `scripts/check_silent_drop.py` keeps the driver's object loops honest: a
  `continue` is either a selection, or diagnosed, or declared.

### Verification and numerical uncertainty

- 45 declared cases (was 27), 17 checked-in `.k2d` benchmark inputs, 150
  automated tests.
- **A numerical uncertainty band for published numbers**: the Grid Convergence
  Index after Roache (1994) and Celik et al. (2008), in the ASME V&V 20 sense,
  with `mesh::refine_uniform` supplying the nested triplets it assumes. No
  geotechnical vendor manual defines a discretisation-error estimator; this one
  comes from the verification literature, and is itself verified against
  manufactured triplets whose answer is known.
- Measured with it: the Giroud rigid-footing benchmark converges to 15.244
  kN/m ± 0.21%, which is 0.03% from the published PLAXIS 2D value — the
  file's own mesh reads 0.5% higher, and the difference is the mesh, not the
  physics.
- `K2D-A005`: a factor of safety computed with a non-associated flow rule
  depends on the mesh and falls as the mesh is refined (−7.9% over a fourfold
  refinement on Griffiths & Lane). The run says so.
- A verification case that measured nothing was found and fixed: the
  strength-reduction search hard-coded its own trial tolerance, so a sweep over
  the tolerance had been setting a number nothing read. Re-measured with the
  control connected, the factor of safety is bit-identical at and below 1e-3
  but +2.0% at 1e-2 and +45.6% at 1e-1 — always on the unsafe side.

### Compatibility

- Reads `.k2d` versions 1–12. Files written by this build are refused by 0.6.x,
  by design.
- The `.res` results format is unchanged (version 5).
- No breaking change to the Python surface or the CLI contract; both gained the
  new inputs.

## [0.6.1] — 2026-08-07

Patch release: the published binaries catch up with the tree.

### Engine
- Axisymmetric staged (Plastic) phases: load and prescribed-displacement
  staging now runs in r–z kinematics from every front end — the
  displacement-controlled bearing-capacity workflow works axisymmetrically.
  Previously any axisymmetric project was refused by the phase driver, so
  the capability was unreachable from the CLI and Python. What is still not
  plumbed for axisymmetry is refused honestly, never integrated as plane
  strain: activation changes (excavation / fill) and consolidation, flow,
  dynamic and safety phases.

### Verification
- Three new published benchmarks in the corpus, each a checked-in `.k2d`
  solved from the file by the suite:
  - `KV-FND-013` — bearing capacity of a smooth rigid circular footing,
    axisymmetric Mohr–Coulomb (Cox 1962 slip-line solution; PLAXIS 2D
    Validation Manual section 3.1).
  - `KV-FND-014` — smooth strip footing on clay with strength increasing
    with depth (Davis & Booker 1973; PLAXIS 2D Validation Manual
    section 3.2).
  - `KV-SLP-002` — Griffiths & Lane (1999) Example 1, the homogeneous 2:1
    slope, factor of safety by phi–c reduction against the published pair
    (their FE 1.4, Bishop & Morgenstern charts 1.380).
- A benchmark walkthrough document: the three cases end to end from the
  command line, with the published values beside the computed ones
  (`docs/validation/three-published-benchmarks.md`).

## [0.6.0] — 2026-07-20

The first published version: the engine, the command line, the Python
surface and the verification record.

### Engine
- Staged construction (plastic phases with PLAXIS-style inherited
  activation), K0 procedure and gravity initial stress.
- Constitutive models: Linear Elastic, Mohr–Coulomb with a tension cut-off,
  Hardening Soil, HS-small, Soft Soil, Soft Soil Creep; drained,
  undrained (A/B) and non-porous drainage types.
- Groundwater: steady-state and transient flow with free surfaces and
  seepage faces, Biot consolidation (elastoplastic), fully coupled
  flow–deformation.
- Dynamics: seismic time histories (harmonic, Ricker, stored accelerogram
  records) with compliant-base and free-field boundaries, Rayleigh damping,
  optional full plasticity during shaking, response spectra.
- Safety: phi–c strength reduction with an honest lower-bound flag.
- Structural elements: plates with elastoplastic Mp/Np hinges, embedded
  beam rows, node-to-node anchors, geogrids, Coulomb interfaces.
- Elements: 6- and 15-node triangles; plane strain and axisymmetric.
- Design codes: EC7 (EN 1997-1) and TBDY 2018 partial-factor catalogues.

### Front ends
- `katai` command line: `solve` / `validate` / `info` with a documented,
  test-pinned exit-code contract (0/2/3/4/5).
- Python package: the engineer-facing `katai.Project` builder over the
  published facade, plus the same `katai` command as a pip console script.
- The `.k2d` input format (version 2) with a JSON Schema, and the `.res`
  results format (version 5); one contract shared by every front end,
  pinned byte-for-byte by the suite.

### Verification
- The full CTest suite runs on every composition of the build; the
  verification corpus (13 published benchmark cases as checked-in `.k2d`
  files) and the generated verification matrix with its bibliography ship
  in the repository.
- The `portable` preset reproduces the entire suite without any
  proprietary component.
