# KATAI 2D

[![CI](https://github.com/fkkarakurt/katai2d/actions/workflows/ci.yml/badge.svg)](https://github.com/fkkarakurt/katai2d/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-green.svg)](CITATION.cff)

A verification-first 2D finite element engine for geotechnical engineering:
soil, structural elements and groundwater in one staged-construction
workflow, driven from the command line or from Python.

**Status: pre-release.** The engine and its validation record are under
active development; interfaces and the file format may still change. The
longer-form formulation and validation notes are being prepared for the
project website; what ships in `docs/` — the format spec, the schema and
the verification matrix — is pinned to the code by the test suite.

## Install

**Command line** — one PowerShell line: downloads the latest release,
verifies its SHA-256 and puts `katai` on your PATH (no admin rights, no
dependencies — `katai.exe` is a single self-contained file):

```powershell
irm https://raw.githubusercontent.com/fkkarakurt/katai2d/main/install.ps1 | iex
```

or download `katai2d-<version>-win64.zip` from
[Releases](https://github.com/fkkarakurt/katai2d/releases) and unzip it
anywhere.

**Python** — one wheel serves every CPython ≥ 3.12 (Windows x64); download
it from [Releases](https://github.com/fkkarakurt/katai2d/releases):

```powershell
pip install katai2d-0.6.0-cp312-abi3-win_amd64.whl
```

The wheel ships the `katai` package *and* the same `katai` command line —
one `pip install` gives you both front ends, with the same commands and the
same exit codes as the native executable. Plain `pip install katai2d` from
PyPI is planned for the announcement.

Both artifacts are verified the way the source tree is: before they leave
the build they must reproduce the corpus numbers in a clean environment
(`scripts/package_cli.ps1`, `scripts/build_wheel.ps1`).

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
count as verification. The record is the generated
[verification matrix](docs/validation/verification-matrix.md) with its
[bibliography](docs/validation/references.bib) — both produced from the
reference declarations inside the tests themselves — and the benchmark
inputs are plain `.k2d` files checked in under [tests/corpus/](tests/corpus/).

For a guided tour, [Three published benchmarks, end to
end](docs/validation/three-published-benchmarks.md) reproduces Cox (1962),
Davis & Booker (1973) and Griffiths & Lane (1999) letter-for-letter from the
sources and solves each one with the `katai` command line, published value
beside computed value.

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
`engine`, `release`). `scripts/check_composition.ps1` walks the first four.

`scripts/package_cli.ps1` builds the end-user distribution: a self-contained
`katai.exe` (optimized, static CRT, Eigen backend — no runtime dependency at
all), zipped under `dist/` with its SHA-256. Unzip anywhere and run.

## Using it

### The command line

```
katai solve <file.k2d> [--out <file.res>]   # run the staged analysis (--out writes a .res)
katai validate <file.k2d>                   # schema + physics validation, no solve
katai info                                  # version and file-format information
```

A real session on a corpus case — the Griffiths & Lane (1999) slope, whose
published factor of safety the suite pins:

```text
> katai validate kv-slp-001-griffiths-lane-slope.k2d
warning: materials[0].gamma_sat: "Griffiths-Lane soil": the saturated unit weight
(20 kN/m3) is below the unsaturated one (20.2 kN/m3), which is physically unusual
OK: kv-slp-001-griffiths-lane-slope.k2d satisfies the input contract (1 warning(s))

> katai solve kv-slp-001-griffiths-lane-slope.k2d
phase 1/1: Initial phase
phase 1/1: ok  max|u| = 2.999373e-01 m  FoS = 1.010
solved 1 phase(s) in 8.16 s
```

Validation never solves, and a solve never proceeds past a refused input: a
readable file with a broken physics contract is an error, not a warning. The
exit codes are part of the contract, pinned by `test_cli`:

| code | meaning |
|---|---|
| 0 | success — `validate`: no errors; `solve`: every phase converged |
| 2 | usage error |
| 3 | a file could not be read, parsed or written |
| 4 | the input contract refused the project — a solve must not proceed |
| 5 | the solve failed: a phase did not converge, or the engine refused it |

> **Note for MKL builds:** MKL is linked dynamically by default, so `katai.exe`
> needs the MKL runtime directory on `PATH` (`<oneAPI>\mkl\latest\bin`). A
> console without it kills the process in the loader — it exits silently,
> printing nothing at all. `scripts/build.ps1` arranges `PATH` for its own
> session and CTest does so for the suite; for your own shell, prepend the
> MKL `bin` directory once. A `portable` or `release` build has no such
> dependency.

The verification corpus under `tests/corpus/` is a good first thing to run:
every `.k2d` there is a published benchmark the suite pins.

### Python

The `katai` package writes and runs the same contract the CLI reads: one
project description, and the suite pins the script-built project against the
checked-in `.k2d` byte for byte. Slope stability end to end:

```python
import katai

# Units everywhere: kN, m, day.
prj = katai.Project("Griffiths & Lane slope", mesh_size=3.0, auto_refine=False)

soil = prj.materials.mohr_coulomb("Clayey sand", E=1.0e5, nu=0.3,
                                  c=3.0, phi=19.6, gamma=20.2)

# Geometry drawn counter-clockwise, one fixity name per edge.
prj.geometry.polygon(
    [(20, 20), (70, 20), (70, 35), (50, 35), (30, 25), (20, 25)],
    material=soil,
    fix=["full", "horizontal", "free", "free", "free", "horizontal"])

prj.initial(procedure="safety")     # phi-c reduction of the gravity state

job = prj.run()
res = job.results()[-1]
print(f"FoS = {res.fos:.3f}")       # published ~0.99; 1.010 on this mesh

prj.save("slope.k2d")               # the same case, ready for `katai solve`
```

Staged construction deactivates regions phase by phase, by handle; a phase
lists only what changes and inherits the rest (excerpt from the
staged-excavation example):

```python
upper = prj.geometry.rectangle(0.0, 6.0, 20.0, 10.0, material=upper_m,
                               name="Upper", bottom="free")
prj.phases.plastic("Excavate", deactivate=[upper])
```

Complete, commented versions of these ship under
[python/examples/](python/examples/): slope stability, Terzaghi
consolidation (the settlement time series against U(Tv)) and staged
excavation (rebound against the closed-form solution). Each one runs in
CTest against its verification-corpus band, so an example that stops telling
the truth fails the build.

The input format is documented in
[docs/k2d-format.md](docs/k2d-format.md) with a machine-readable schema in
[docs/k2d.schema.json](docs/k2d.schema.json).

## License

Apache License 2.0 — see [LICENSE](LICENSE). Vendored third-party software
is listed in [NOTICE](NOTICE) (Eigen, MPL-2.0). If you use KATAI 2D in
research, see [CITATION.cff](CITATION.cff).
