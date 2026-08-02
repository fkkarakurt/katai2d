# Decisions

Architecture decision record. A decision is never deleted once taken: it is marked **superseded** and
the replacement is named, so the reason a thing is the way it is stays recoverable.

Commercial and licensing decisions are deliberately absent from this file. They are deferred and kept
with the maintainer's notes; nothing in the code depends on them.

## Locked

| # | Subject | Decision | Status |
|---|---|---|---|
| D2 | Development model | single maintainer | accepted |
| D3 | Platform | Windows first, then macOS and Linux; written to be portable throughout | accepted |
| D4 | Language | modern C++ | accepted |
| D5 | FEM core | written from scratch, with two exceptions (below) | accepted |
| D8 | Constitutive models | Linear Elastic, Mohr-Coulomb, Hardening Soil, HS-small, Soft Soil, Soft Soil Creep | accepted; extended beyond the original three |
| D9 | Structural elements | plate, anchor, geogrid, interface, embedded beam row | accepted |
| D10 | Element type | 6- and 15-node triangles, user's choice | accepted; 15-node matches the PLAXIS default |
| D11 | Dimensional models | plane strain and axisymmetric | accepted |
| D16 | GUI toolkit | Dear ImGui, docking branch (MIT) | accepted; Qt, wxWidgets and FLTK were considered and rejected |
| D17 | Render backend | OpenGL 4.1 core behind a thin `Renderer` abstraction | accepted; the application is solver-bound, not GPU-bound, so Vulkan or Metal remain open without urgency |

## Superseded

| # | Subject | Original decision | Superseded by |
|---|---|---|---|
| D6 | Linear solver | Intel MKL PARDISO | **Superseded by `katai::linsolve` (Stage A2, 2026-07-28), which is what D6 always intended:** one `DirectSolver` interface, an always-available Eigen backend, PARDISO as an optional adapter chosen at link time, and every solve verified by residual. No test but the two adapter unit tests is gated on MKL; the `portable` preset runs the whole suite on Eigen. See `ENVIRONMENT.md`. |
| D14 | Scripting | own scripting layer, leaning towards Lua | **Python**, in two layers: a binding module over the engine and an engineer-facing package above it. Lua was never implemented. |
| D15 | Build speed | Ninja + sccache + fast linker + PCH | Ninja, ccache and (since 2026-07-28) a shared precompiled header are in place — measured, Eigen was 93–99.9 % of the parse volume, and MKL's static link another 57 minutes, so D15's original instinct about the PCH and the linker was right and the "header-heavy core is the bottleneck" correction was itself wrong by an order of magnitude. Clean build measured 39.6 min. Moving non-template code into `.cpp` remains worthwhile, at its own stage. |

## D5 exception — what "from scratch" excludes

- **The linear solver.** Eigen for the default direct solve, PARDISO when available. Writing a
  competitive sparse direct solver is not the point of this project.
- **Nothing else.** Mesh generation *is* written from scratch — constrained Delaunay plus Ruppert
  refinement — because no permissive library of adequate quality was available. The algorithms follow
  Shewchuk; none of the code does.

## Dependency policy

Third-party components must be usable in a commercial product. Copyleft (GPL, LGPL) is excluded;
permissive licences (MIT, BSD, Apache, MPL, Boost, zlib) are acceptable. Dependencies are vendored
under `third_party/` rather than fetched at configure time, so a build works offline and every commit
compiles against exactly the sources it was tested with.

Actually vendored today:

| Purpose | Library | Licence |
|---|---|---|
| Dense linear algebra, small matrices, sparse direct solve | Eigen | MPL-2.0 |
| GUI | Dear ImGui | MIT |
| Window and GL context | GLFW | zlib/libpng |
| GL loader | glad | permissive |
| Sparse direct solver (optional backend, not vendored) | Intel oneMKL PARDISO | proprietary, redistributable |
| Icons | Font Awesome 6 Free | icons CC BY 4.0 (attribution required), fonts OFL-1.1, code MIT |

Considered and **not** adopted: Clipper2 and Boost.Geometry (polygon and geometry helpers — the
planar arrangement was written directly instead), sokol_gfx, Lua, OpenMP and oneTBB (the engine uses
its own thread pool).
