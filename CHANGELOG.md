# Changelog

All notable changes to KATAI 2D. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
MAJOR.MINOR.PATCH.

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
