# Changelog

All notable changes to KATAI 2D. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
MAJOR.MINOR.PATCH.

## [Unreleased]

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
