# KATAI 2D

A verification-first 2D finite element engine for geotechnical engineering:
soil, structural elements and groundwater in one staged-construction
workflow, driven from the command line or from Python.

**Status: pre-release.** The engine and its validation record are under
active development; interfaces and the file format may still change. Parts
of `docs/references/` and `docs/validation/` are in the process of being
translated to English — the numbers in them are current, the prose is
catching up.

## Capabilities

- **Analysis types** — staged construction (plastic), K0 procedure and
  gravity initial stress, Biot consolidation (elastoplastic), fully coupled
  flow–deformation, steady-state and transient groundwater flow, safety
  (phi–c strength reduction), and seismic dynamics: acceleration time
  histories with compliant-base and free-field boundaries, plus response
  spectra.
- **Constitutive models** — Linear Elastic, Mohr–Coulomb with a tension
  cut-off, Hardening Soil, HS-small, Soft Soil and Soft Soil Creep; drained,
  undrained (A/B) and non-porous drainage types.
- **Elements** — 6- and 15-node triangles, plane-strain and axisymmetric;
  plates with elastoplastic Mp/Np hinges, embedded beam rows, node-to-node
  anchors, geogrids and Coulomb interfaces.
- **Design codes** — EC7 (EN 1997-1) and TBDY 2018 partial-factor
  catalogues, applied through material factoring.

## Verification

The suite is the specification. Every capability is checked in CTest against
a closed-form solution, an independent computation path that shares no code
with the solver, or a published benchmark; self-consistency alone does not
count as verification. The record lives in [docs/validation/](docs/validation/),
the formulation notes with their primary sources in
[docs/references/](docs/references/), and the benchmark inputs are plain
`.k2d` files checked in with the tests.

The `portable` preset builds and runs the whole suite without any
proprietary component, so every published number can be reproduced with the
vendored Eigen solver alone.

## Building (Windows, MSVC)

Requirements: a recent MSVC toolchain, CMake ≥ 3.25, Ninja. Optional: Intel
oneMKL (PARDISO solver backend; the build falls back to Eigen without it),
Python ≥ 3.10 with `nanobind` for the bindings.

```powershell
.\scripts\build.ps1 -Configure -Test        # configure + build + full suite
```

`scripts/build.ps1` loads the MSVC environment, discovers oneMKL and a real
Python, and drives the CMake presets (`msvc-rwdi`, `msvc-debug`, `portable`,
`engine`). `scripts/check_composition.ps1` walks all of them.

## Using it

```
katai validate model.k2d      # schema + physics validation, no solve
katai solve model.k2d         # run the staged analysis, write results
katai info                    # version and build information
```

The Python package exposes the same facade — build a project in a script,
solve it, and read results as arrays — and writes the same `.k2d` files the
CLI consumes. The input format is documented in
[docs/k2d-format.md](docs/k2d-format.md) with a machine-readable schema in
[docs/k2d.schema.json](docs/k2d.schema.json).

## License

Apache License 2.0 — see [LICENSE](LICENSE). Vendored third-party software
is listed in [NOTICE](NOTICE) (Eigen, MPL-2.0). If you use KATAI 2D in
research, see [CITATION.cff](CITATION.cff).
