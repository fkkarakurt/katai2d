# KATAI 2D

[![CI](https://github.com/fkkarakurt/katai2d/actions/workflows/ci.yml/badge.svg)](https://github.com/fkkarakurt/katai2d/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)
[![Cite](https://img.shields.io/badge/Cite-CITATION.cff-green.svg)](CITATION.cff)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20692282-blue.svg)](https://doi.org/10.5281/zenodo.20692282)

A verification-first 2D finite element engine for geotechnical engineering:
soil, structural elements and groundwater in one staged-construction workflow,
driven from the command line or from Python.

Every capability is pinned to a closed-form solution, an independent
computation path that shares no code with the solver, or a published benchmark
— **self-consistency does not count as verification** — and the whole record
ships with the source: **58 declared verification cases**, **26 benchmark
`.k2d` input files** anyone can rerun, **153 automated tests**.

**Status: pre-release.** The engine and its validation record are under active
development; interfaces and the file format may still change. What ships in
`docs/` — the format spec, the JSON Schema and the verification matrix — is
generated from the code and pinned to it by the test suite.

---

## Contents

[Install](#install) · [Quick start](#quick-start) · [Python API](#the-python-api)
· [Command line](#the-command-line) · [Capabilities](#capabilities)
· [Verification](#verification) · [Building](#building-from-source)
· [License](#license-and-citation)

---

## Install

Current release: **v0.8.1** (`.k2d` format v14).
Windows x64. Both artifacts are on the
[Releases](https://github.com/fkkarakurt/katai2d/releases) page.

### Command line

One PowerShell line — downloads the latest release, verifies its SHA-256 and
puts `katai` on your PATH. No admin rights, no dependencies: `katai.exe` is a
single self-contained file.

```powershell
irm https://raw.githubusercontent.com/fkkarakurt/katai2d/main/install.ps1 | iex
```

Or download `katai2d-0.8.1-win64.zip` and unzip it anywhere.

### Python

One abi3 wheel serves every CPython ≥ 3.12 on Windows x64:

```powershell
pip install katai2d-0.8.1-cp312-abi3-win_amd64.whl
```

The wheel ships the `katai` package **and** the same `katai` command line — one
`pip install` gives you both front ends, with the same commands and the same
exit codes as the native executable. Plain `pip install katai2d` from PyPI is
planned for the announcement.

Both artifacts are verified the way the source tree is: before they leave the
build they must reproduce the corpus numbers in a clean environment
(`scripts/package_cli.ps1`, `scripts/build_wheel.ps1`).

### Which version am I on?

```powershell
katai info                                   # or:
python -c "import katai; print(katai.__version__)"
```

Every front end reads one constant, so the number is the build's real identity.
If a name in this document is missing on your machine, check that number first
— the API grew in **0.8.0**, and the table says what arrived when:

| Needs at least | Surface |
|---|---|
| 0.8.0 | `prj.structures` (walls, anchors, geogrids, pile rows, interfaces) |
| 0.8.0 | `katai.summary()`, `katai.extremes()`, bound structural forces |
| 0.8.0 | `prj.phases.transient_flow(...)`, `.fully_coupled(...)`, `design=`, `water=` |
| 0.7.1 | materials, geometry, water, loads, displacements, dewatering, plastic / consolidation / safety / dynamic phases |

A file written by a newer build is refused by an older one rather than
misread — each `.k2d` version bump marks an input the older build would have
dropped in silence.

### If Windows warns about the download

The Windows executable is **not code-signed**, and it is published by a young
project, so it has no reputation with Microsoft SmartScreen or Defender. That
combination can produce a generic machine-learning verdict — v0.7.0 was flagged
as `Trojan:Win32/Wacatac.B!ml`, a label Defender applies to unknown unsigned
binaries rather than to any identified malware. It is a false positive, and it
is being addressed rather than argued with: since v0.7.1 the executable carries
its full version and publisher information, the detection is reported to
Microsoft, and code signing is on the roadmap.

You do not have to take that on trust. Verify what you downloaded against the
checksum published with the release:

```powershell
(Get-FileHash .\katai2d-0.8.1-win64.zip -Algorithm SHA256).Hash.ToLower()
# compare with katai2d-0.8.1-win64.zip.sha256 from the same release
```

and, if you want an independent opinion, upload the file to VirusTotal — it is
a public release asset, so nothing private is disclosed by doing so. The Python
wheel has never been flagged and is the alternative route.

---

## Quick start

A slope, end to end. Units are fixed everywhere and stated in every docstring:
**kN, m, day**.

```python
import katai

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
print(f"FoS = {job.results()[-1].fos:.3f}")   # 1.010; published ~0.99

prj.save("slope.k2d")               # the same case, ready for `katai solve`
```

The same analysis, from the command line, on the file the suite pins:

```text
> katai validate kv-slp-001-griffiths-lane-slope.k2d
warning: materials[0].gamma_sat: "Griffiths-Lane soil": the saturated unit weight
(20 kN/m3) is below the unsaturated one (20.2 kN/m3), which is physically unusual
OK: kv-slp-001-griffiths-lane-slope.k2d satisfies the input contract (1 warning(s))

> katai solve kv-slp-001-griffiths-lane-slope.k2d
phase 1/1: Initial phase
phase 1/1: ok  max|u| = 2.999373e-01 m  FoS = 1.010
  warning K2D-A005  Griffiths-Lane soil: Strength reduction with a non-associated
  flow rule (phi = 19.6 deg, psi = 0 deg): the factor of safety depends on the mesh
  and falls as the mesh is refined, because the shear band narrows with the
  elements. Quote it with the mesh it was computed on, and confirm it with a
  refinement study.
solved 1 phase(s) in 9.79 s
```

That warning is the house style: a number that must be qualified is qualified
where it is produced, not in a footnote somewhere else.

---

## The Python API

`katai.Project` is the engineer-facing builder. A script reads like the
engineering problem — materials, geometry, water, loads, structures, phases —
and every object is created in **one call** that hands back a **handle**.
Phases activate and deactivate those handles, which is the whole of staged
construction.

Everything built here becomes a plain schema project (`prj.build()`), so a
script and a saved `.k2d` describe the same run and `katai solve` executes it
unchanged: the suite pins a DSL-built corpus case **byte-identical** to its
checked-in file.

### The surface, at a glance

| Namespace | What it makes |
|---|---|
| `prj.materials` | `linear_elastic` · `mohr_coulomb` · `hardening_soil` · `hs_small` · `soft_soil` · `soft_soil_creep` |
| `prj.geometry` | `polygon(points, material=…, fix=[…])` · `rectangle(x0, y0, x1, y1, …)` |
| `prj.water` | `table(y)` · `phreatic_line(points)` |
| `prj.loads` | `line_load(a, b, qx=, qy=)` · `point_load(at, qx=, qy=)` |
| `prj.structures` | `plate` · `anchor` · `geogrid` · `pile` · `interface` |
| `prj.displacements` | `line(a, b, ux=, uy=)` — prescribed displacement |
| `prj.dewatering` | `well(a, b, q=, h_min=)` · `drain(a, b, head=)` |
| `prj.phases` | `plastic` · `consolidation` · `safety` · `transient_flow` · `fully_coupled` · `dynamic` |
| `prj` | `initial(procedure=…)` · `validate()` · `build()` · `save(path)` · `run()` |
| `katai` | `summary(job, prj)` · `extremes(result)` · `Refusal` · `run(project)` |

Stiffnesses are **per metre out of plane**, as everywhere else in plane strain.
Anchors and pile rows are the exception and say so: they are discrete members,
so they carry a `spacing` [m] and the engine divides by it.

### A real job: an anchored excavation

Ground, a diaphragm wall with interfaces on both faces, a prestressed anchor
row, and a staged excavation that installs the anchor and digs the pit — with
no `.k2d` file anywhere.

```python
import katai

prj = katai.Project("Anchored excavation", mesh_size=2.0, auto_refine=False)

sand = prj.materials.mohr_coulomb("Berlin sand", E=3.0e4, nu=0.3,
                                  c=5.0, phi=32.0, gamma=19.0)

# Ground, split so the pit can be dug out of it later.
prj.geometry.rectangle(0, 0, 30, 12, material=sand, name="Below formation")
pit = prj.geometry.rectangle(0, 12, 12, 20, material=sand, name="Pit")
prj.geometry.rectangle(12, 12, 30, 20, material=sand, name="Retained ground")

# EA and EI per metre out of plane; an interface on each face, so the soil can
# slip against the wall instead of being glued to it.
wall = prj.structures.plate((12, 20), (12, 6), EA=1.2e7, EI=1.0e5, w=5.0,
                            interfaces="both", name="Diaphragm wall")

# A discrete member: it carries its out-of-plane spacing and a lock-off force
# that holds the excavation before anything moves.
anchor = prj.structures.anchor((12, 18), (20, 14), EA=2.0e5, spacing=2.5,
                               prestress=300.0, name="Anchor row 1")

prj.initial(procedure="k0", exclude=[anchor])
prj.phases.plastic("Excavate to formation", activate=[anchor], deactivate=[pit])

job = prj.run()
print(katai.summary(job, prj))
```

`summary()` names the phases, tabulates the extremes **with the place each
occurs**, lists the structural force envelopes and repeats every diagnostic the
engine raised — because those change how the numbers should be read (rows
abbreviated here; the full print also carries `u_x`, `sig'_xx`, `sig_xy`,
`tau_max`, `p_w` and the shear column):

```text
=== Excavate to formation =======================================
  status        converged
  max |u|       4.329029e-02 m
  diagnostics
    NOTE  K2D-A004
          A structural element ends on a supported node. The reported reactions
          are the soil's contribution only; the element's own end force at that
          support is not included in this build.
  extremes                    min          at                 max          at
    |u|       [  m]             0  (0, 16)                 0.04329  (3, 12)
    u_y       [  m]   -0.00089196  (30, 18)               0.043216  (3, 12)
    sig'_yy   [kPa]       -376.15  (30, 0)                   5.888  (10, 12)
  structural forces (envelope |max|)
    Anchor row 1           anchor    N       195.4 kN/m
    Diaphragm wall         plate     N       325.7 kN/m   M       58.33 kNm/m
```

"The largest settlement is 43 mm" and "…and it is three metres from the wall"
are two different findings; a field that is uniform says so instead of
inventing a location.

Anything can also be read directly — the summary is a convenience, not a gate:

```python
last = job.results()[-1]
wall_force = next(s for s in last.struct_forces if s.name == "Diaphragm wall")
print(wall_force.max_M, wall_force.max_N, wall_force.yielded)
for st in wall_force.stations:          # N, Q, M and position along the element
    print(st.s, st.x, st.y, st.N, st.Q, st.M)
```

Results are arrays, not screenshots: displacements, stresses, pore pressures,
reactions and structural force envelopes come back as NumPy data, so a
parameter study is a `for` loop.

### Staged construction

Activation is **inherited** from the previous phase, PLAXIS-style: a phase
lists only what changes.

```python
upper = prj.geometry.rectangle(0.0, 6.0, 20.0, 10.0, material=m,
                               name="Upper", bottom="free")
prj.phases.plastic("Excavate", deactivate=[upper])
prj.phases.consolidation("Dissipate", duration=100.0, steps=50)
prj.phases.safety("FoS after excavation")
```

Phases also take the controls that make a run reproducible and a design check
possible — all written into the `.k2d`:

```python
prj.phases.plastic("Half the stage", apply_fraction=0.5)      # PLAXIS Sum-Mstage
prj.phases.plastic("Design check", design="ec7_da3")          # partial factors
prj.phases.plastic("Dewater", water=8.0)                      # phase water table
prj.phases.plastic("Fine", tolerance=1e-3, max_iterations=60)
```

### Refusals

A refusal is an exception carrying the engine's own words and the stable
`K2D-*` codes, so a script branches on the reason instead of parsing prose:

```python
try:
    job = prj.run()
except katai.Refusal as e:
    if any(d.code == "K2D-G003" for d in e.diagnostics):
        ...            # a distributed load was drawn off the soil
    print(e.report)    # field-path issues, when the refusal came from validation
```

### Examples

Complete, commented programs ship under [python/examples/](python/examples/):

| Example | What it checks itself against |
|---|---|
| [`slope_stability.py`](python/examples/slope_stability.py) | Griffiths & Lane (1999); FoS 1.010 |
| [`consolidation_terzaghi.py`](python/examples/consolidation_terzaghi.py) | the Terzaghi U(T_v) series; \|ΔU\| ≤ 0.0061 at the sampled times |
| [`staged_excavation.py`](python/examples/staged_excavation.py) | 1D elastic unloading by hand; +0.00% |
| [`anchored_excavation.py`](python/examples/anchored_excavation.py) | its own invariants (convergence, anchor ≥ lock-off, wall carries moment) |

Each one runs in CTest against its verification band, so an example that stops
telling the truth fails the build rather than the first reader who tries it.

---

## The command line

```
katai solve <file.k2d> [--out <file.res>]   # run the staged analysis
katai validate <file.k2d>                   # schema + physics validation, no solve
katai info                                  # version and file-format information
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

The verification corpus under [tests/corpus/](tests/corpus/) is a good first
thing to run: every `.k2d` there is a published benchmark the suite pins.

> **Note for MKL builds:** MKL is linked dynamically by default, so `katai.exe`
> needs the MKL runtime directory on `PATH` (`<oneAPI>\mkl\latest\bin`). A
> console without it kills the process in the loader — it exits silently,
> printing nothing at all. `scripts/build.ps1` arranges `PATH` for its own
> session and CTest does so for the suite; for your own shell, prepend the MKL
> `bin` directory once. The **released** artifacts and the `portable` and
> `release` presets have no such dependency.

---

## Capabilities

- **Analysis types** — staged construction (plastic), K₀ procedure and gravity
  initial stress, Biot consolidation (elastoplastic), fully coupled
  flow–deformation, steady-state and transient groundwater flow, safety
  (φ–c strength reduction), and seismic dynamics: acceleration time histories
  with compliant-base and free-field boundaries, plus response spectra.
- **Constitutive models** — Linear Elastic, Mohr–Coulomb with a tension
  cut-off, Hardening Soil, HS-small (with Masing's rule and the Li & Dafalias
  dilatancy below the phase-transformation line), Soft Soil and Soft Soil
  Creep; drained, undrained (A/B/C) and non-porous drainage types.
- **Elements** — 6- and 15-node triangles, plane-strain and axisymmetric;
  plates with elastoplastic Mp/Np hinges, embedded beam rows, node-to-node
  anchors, geogrids and Coulomb interfaces.
- **Design codes** — EC7 (EN 1997-1) and TBDY 2018 partial-factor catalogues,
  applied per phase through material factoring.

**Not there yet, stated rather than implied:** three-dimensional analysis (this
is 2D by design), a graphical interface (a separate, in-development track),
unsaturated soil mechanics and liquefaction models, operating systems other
than Windows x64, and `pip install katai2d` from PyPI.

---

## Verification

The suite is the specification. Every capability is checked in CTest against a
closed-form solution, an independent computation path that shares no code with
the solver, or a published benchmark.

The record is the generated
[verification matrix](docs/validation/verification-matrix.md) with its
[bibliography](docs/validation/references.bib) — both produced from the
reference declarations **inside the tests themselves**, so the table cannot
drift from the suite; a gate fails the build when it does. The benchmark inputs
are plain `.k2d` files checked in under [tests/corpus/](tests/corpus/).

A selection:

| Benchmark | Reference | Result |
|---|---|---|
| Prandtl strip footing N_c (φ = 0) | 2 + π (Prandtl 1921) | +0.6% |
| Slope factor of safety (φ–c reduction) | Griffiths & Lane (1999), published ≈ 0.99 | 1.010 |
| Rigid strip footing on elastic soil | Giroud (1972) 15.15 | +1.1% |
| Unconfined dam discharge with a seepage face | Charny (1951) exact theorem | +1.02% |
| Terzaghi 1D consolidation U(T_v) | Terzaghi series | −1.2% … −0.4% |
| Resonant column at f₁ — \|u\|, \|a\| | damped SH closed form (Kramer 1996) | −0.6% / −0.2% |
| El Centro 1940 NS record identity | published PGA ≈ 0.319 g | 0.31882 g at 2.02 s |
| Sparse solve vs dense LU, independent path | Eigen FullPivLU, no shared code | 4.4×10⁻¹⁶ |

For a guided tour, [Three published benchmarks, end to
end](docs/validation/three-published-benchmarks.md) reproduces Cox (1962),
Davis & Booker (1973) and Griffiths & Lane (1999) letter-for-letter from the
sources and solves each one with the `katai` command line, published value
beside computed value.

Where a number is uncertain, the uncertainty is measured rather than asserted:
[numerical-uncertainty.md](docs/validation/numerical-uncertainty.md) records
grid-convergence bands, and states plainly where the mesh turned out **not** to
be the axis the error was on.

The `portable` preset builds and runs the whole suite without any proprietary
component, so every published number can be reproduced with the vendored Eigen
solver alone.

---

## Building from source

Windows, MSVC. Requirements: a recent MSVC toolchain, CMake ≥ 3.25, Ninja.
Optional: Intel oneMKL (PARDISO backend; the build falls back to Eigen without
it), Python ≥ 3.10 with `nanobind` for the bindings — ≥ 3.12 to build the
single abi3 wheel that ships.

```powershell
.\scripts\build.ps1 -Configure -Test        # configure + build + full suite
```

`scripts/build.ps1` loads the MSVC environment, discovers oneMKL and a real
Python, and drives the CMake presets (`msvc-rwdi`, `msvc-debug`, `portable`,
`engine`, `release`). `scripts/check_composition.ps1` walks the first four.

`scripts/package_cli.ps1` builds the end-user distribution: a self-contained
`katai.exe` (optimized, static CRT, Eigen backend — no runtime dependency at
all), zipped under `dist/` with its SHA-256. `scripts/build_wheel.ps1` builds
the wheel and verifies it in a clean virtual environment before it is allowed
to leave the build.

---

## Documentation

| | |
|---|---|
| [`docs/k2d-format.md`](docs/k2d-format.md) | the input format, field by field |
| [`docs/k2d.schema.json`](docs/k2d.schema.json) | machine-readable JSON Schema |
| [`docs/diagnostics.md`](docs/diagnostics.md) | every `K2D-*` code and what it means |
| [`docs/validation/verification-matrix.md`](docs/validation/verification-matrix.md) | the generated verification record, with its [bibliography](docs/validation/references.bib) |
| [`docs/validation/numerical-uncertainty.md`](docs/validation/numerical-uncertainty.md) | how far a computed number sits from the exact solution of its own equations |
| [`docs/validation/plaxis-2d-validation-comparison.md`](docs/validation/plaxis-2d-validation-comparison.md) | four documented validation cases, rebuilt from their problem statements |
| [`CHANGELOG.md`](CHANGELOG.md) | what changed, and which wrong answers it fixed |

---

## License and citation

Apache License 2.0 — see [LICENSE](LICENSE). Vendored third-party software is
listed in [NOTICE](NOTICE) (Eigen, MPL-2.0).

If you use KATAI 2D in research, cite it through [CITATION.cff](CITATION.cff)
or the archived record: [10.5281/zenodo.20692282](https://doi.org/10.5281/zenodo.20692282).

If a number looks wrong, the issue tracker has a dedicated *numerical result*
template — state the computed value and the reference with its source, and it
will be resolved against primary sources.
