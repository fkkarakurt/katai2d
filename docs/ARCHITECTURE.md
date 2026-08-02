# Architecture

This document describes the tree that exists. Where reality and the intended design differ, the
difference is stated together with the stage that closes it — a document that describes a plan as if
it were the code is a silent-wrong of its own kind.

## The layer contract

Five layers, each with one job, each depending only downward.

| # | Layer | Owns | May depend on |
|---|---|---|---|
| 1 | **engine** — `math`, `linsolve`, `geometry`, `mesh`, `materials`, `fem`, `analysis`, `model`, `io` | discretisation, constitutive integration, assembly, phase solvers, file formats | layer 1 |
| 2 | **jobs** — `kernel/jobs` | the execution seam: the analysis/mesh/flow drivers over the schema (phase sequencing lives here); the job object with progress, cancellation and batch runs is the rest of Stage D | 1 |
| 3 | **api** — `kernel/api` | the published surface: the C++ facade front ends program against (`<katai/api/katai.hpp>`); the C ABI for user materials and the Python binding target come later | 1, 2 |
| 4 | **front ends** — `studio/render`, `studio/app`, `cli` *(planned)* | presentation and argument parsing, nothing else | 3, and 2 for progress |
| 5 | **extension points** *(planned)* | third-party code entering through a versioned contract | 3 |

Enforced by `scripts/check_architecture.py`, run as the CTest `check_architecture`:

- **`dag`** — dependencies point downward only, and never cycle inside a layer.
- **`api`** — only `include/katai/**` is public; no module reaches into another module's `src/` or a
  private header.
- **`physics`** — no constitutive, assembly or analysis *definition* above the engine. Calling into
  the engine is what a front end is for, so the check distinguishes a definition from a call.

Accepted violations are listed inside that script together with the stage that removes each one. An
exception there is a debt with an owner, not a permission.

## Module tree

```
kernel/
  math/      include/katai/math/       dense and sparse containers, thread pool,
                                       solve-error vocabulary, PARDISO adapter
                                       (compiled only when MKL is present)
  linsolve/  include/katai/linsolve/   DirectSolver interface: Eigen backend always,
                                       PARDISO optional, chosen at link time; every
                                       solve verified by residual
  geometry/  include/katai/geometry/   2D primitives: rectangular domain,
                                       planar arrangement (segments -> closed faces)
  mesh/      include/katai/mesh/       robust predicates, incremental Delaunay,
                                       Ruppert refinement, tri6/tri15 promotion
  materials/ include/katai/materials/  constitutive models as pure functions, plus
                                       the registry that resolves a model by name
  fem/       include/katai/fem/        elements/ and assembly/: element kernels,
                                       DOF management, the global assembler
  analysis/  include/katai/analysis/   phase solvers and post/ (stress recovery)
  model/     include/katai/model/      the project schema
  io/        include/katai/io/         the .k2d project file (versioned JSON) and the
                                       .res results file (versioned binary), one module
                                       so the two version disciplines live side by side
  jobs/      include/katai/jobs/       layer 2, the execution seam: driver.hpp (schema ->
                                       neutral resolution -> engine phase strategies),
                                       mesh_builder.hpp (project -> FE mesh),
                                       flow_driver.hpp (groundwater-flow runs),
                                       job.hpp (a run as an object: progress, cooperative
                                       cancel, timings; validates before it solves)
  api/       include/katai/api/        layer 3, the published facade: katai.hpp, the one
                                       header a front end programs against (section 7.3:
                                       total over the input, closed over the implementation)
studio/
  render/    include/katai/render/     camera, field view, GL programs
  app/                                 the ImGui application
tests/                                 unit, regression and V&V suite; on-demand studies
scripts/                               build helper and the source-tree gates
docs/                                  formulation notes, V&V records, this document
LaTeX/                                 Scientific Manual, Verification Manual
```

### Include convention

One root for the whole tree, angle-bracketed, as with `Eigen/` or `boost/`:

```cpp
#include <katai/analysis/staged_construction.hpp>
#include <katai/mesh/mesh.hpp>
```

A module's public headers live under `<module>/include/katai/<module>/`. Anything under
`<module>/src/` is private to that module.

### Namespace convention

`katai::<module>`, with implementation details in `katai::<module>::detail`. Domain sub-namespaces
(`katai::core::plate`, `katai::core::tri`, `katai::core::iface`, ...) group element and material
families. No unscoped and no suffixed detail namespaces.

One deliberate exception: the code inside `materials`, `fem` and `analysis` still lives in
`katai::core`. The Stage A2 split was scoped to moves and include roots; renaming the namespace is
about a thousand qualified uses plus three hundred using-declarations, and is decided as its own
step rather than smuggled into a mechanical batch. The `io` and `jobs` modules carry the same
debt: `project_io.hpp` still spells `katai::model`, while `results_io.hpp` and the three relocated
driver headers in `katai/jobs` spell `katai::app`, because every consumer names their APIs there —
the spelling joins that dedicated rename step (moving a file and renaming its namespace are kept
as separate changes in this repository).

## Inside the engine modules that came out of `core`

| Module | Contents |
|---|---|
| `materials/` | Linear Elastic, Mohr-Coulomb, Hardening Soil, HS-small, Soft Soil, Soft Soil Creep, water retention (van Genuchten / Mualem), the `MaterialModel` dispatch, and the constitutive **registry** — models resolved by canonical name, construction and honest drainage refusals living with the models, external registration through the same seam (the planned user-material ABI implements it from outside the tree). Depends on vendored Eigen only, nothing else in the tree |
| `fem/elements/` | tri6, tri15, axisymmetric variants, element traits, point location; structural elements: plate, interface, geogrid, embedded beam |
| `fem/assembly/` | DOF map, including extra structural DOFs, and the global assembler |
| `analysis/` | nonlinear solver, staged construction, initial stress and K0, strength reduction, seepage, transient and unsaturated flow, consolidation, coupled flow-deformation, dynamics (linear and nonlinear), structural dynamics, free field, response spectrum, design codes, embedded wall builder, internal and structural forces |
| `analysis/post/` | stress recovery to nodes |

Elements and assembly are one module on purpose: `embedded_beam` and `dof_map` reference each
other, which inside a module is a legitimate detail and between two modules would be a cycle.

## Where the tree does not yet tell the truth

| Reality | Consequence | Closed by |
|---|---|---|
| Layers 2 and 3 exist (`katai/jobs` with the Job object, `katai/api` with the facade), but the GUI still calls the layer-2 free functions directly rather than submitting Jobs through the facade | the tested surface and the GUI's surface are equal only by convention until the migration; `test_layer_includes` cannot be turned on yet | Stage D item 3 (GUI onto Job/api), then the gate's layer tests |
| Layer 5 (extension points) does not exist | third-party code has no versioned entry | the UMAT work, queued after the five priorities |
| `materials` + `fem` + `analysis` together are 39 headers to 11 sources | header-heavy modules re-parse in every consumer (the precompiled header pays most of that cost today) | Stage C — move non-template code into `.cpp` |
| The code in `materials`, `fem` and `analysis` is still in namespace `katai::core` | the namespace names a module that no longer exists | a dedicated rename step, decided separately (see Namespace convention above) |

Stage names refer to the engineering roadmap kept with the maintainer's notes.

## Dependency direction

```
   app ──┬──────────────┬────────────┬──────────────┬──────────────────────────┐
    │    │              │            │              │                          │
    ▼    ▼              ▼            ▼              ▼                          ▼
  render model      linsolve     analysis ──► fem ──► materials    io (model, analysis, mesh)
    │                   │            │         │
    ▼                   ▼            ▼         ▼
  mesh ──► geometry   math      (math, mesh) (math, mesh)
```

Measured, not intended: `materials` depends on nothing in the tree (vendored Eigen only); `fem` on
`materials`, `math` and `mesh`; `analysis` on all of those; `linsolve` on `math`; `mesh` on
`geometry`; `io` on `model`, `analysis` and `mesh` (the `.res` file stores the engine's result
types over a mesh); `math`, `geometry` and `model` have no internal dependencies; `render` depends
on `mesh` only. `analysis` does **not** depend on `linsolve` — the linear solver enters as a
callback built by the composition root. Nothing in the engine depends on `render` or `app`. Run
`python scripts/check_architecture.py --list` to print the graph resolved from the actual include
statements rather than from this diagram.

## Build and composition

CMake presets, Ninja, MSVC on Windows. `katai_math`, `katai_linsolve`, `katai_mesh`,
`katai_materials`, `katai_fem`, `katai_analysis` and `katai_render` are compiled libraries;
`katai_geometry`, `katai_model` and `katai_io` are header-only `INTERFACE` targets. A shared precompiled header
of third-party includes (Eigen and the standard library — never our own headers) is built once in
`katai_math` and reused by every target; measured 2026-07-28, Eigen was 93–99.9 % of the parse
volume of a typical translation unit.

What a program is made of is decided at the composition root and nowhere else: `KATAI_WITH_GUI`
gates the studio and the GL stack, `KATAI_WITH_MODEL` the project schema and the file formats,
`KATAI_WITH_MKL` the PARDISO adapter — and shared code never tests any of them (architecture gate, check `config`). The
`engine` preset turns schema and GUI off, and exactly the tests confined to the engine modules are
built: a test names its modules in `katai_add_test`, which both labels it (`ctest -L materials`
runs precisely that module's tests) and scopes it (a composition without a named module skips the
test and says so at configure time — a missing module is refused, never silently substituted).
`scripts/check_composition.ps1` walks every buildable row of the configuration matrix.

**The sparse solver is `katai::linsolve`, and Intel oneMKL is genuinely optional.** Every consumer
solves through the `DirectSolver` interface; which backend a program contains — the always-available
Eigen backend or the PARDISO adapter — is decided by which selection translation unit CMake compiles,
never by a preprocessor conditional. No test is gated on MKL except the two that unit-test the PARDISO
adapter itself. The `portable` preset (`KATAI_WITH_MKL=OFF`) builds and passes the whole suite on the
Eigen backend, which is what makes the validation record reproducible without a proprietary component.
Every solve is verified against the system it claims to solve (relative residual above 1e-6 is
refused), so solution quality does not depend on the backend either.

`scripts/build.ps1` loads the MSVC environment and configures plus builds a preset in one step, and
`ctest` runs the suite. See `docs/ENVIRONMENT.md`.
