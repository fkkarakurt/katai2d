# Development environment

Recorded for the maintainer's machine; scanned 2026-06 and still current. Paths are specific to this
installation, the traps below are not.

## Toolchain

| Tool | Version / location | Note |
|---|---|---|
| MSVC | VS Community 2026 (18.6); toolset 14.51.36231 and 14.44 | primary compiler; `cl.exe` under `Hostx64\x64` |
| Windows SDK | 10.0.28000 and 10.0.26100 | |
| CMake | 3.29.2 (`D:\Strawberry\c\bin`) | |
| Ninja | 1.12.0 (`D:\Strawberry\c\bin`) | build engine |
| ccache | 4.9.1 (`D:\Strawberry\c\bin`) | supports MSVC, so sccache is not needed |
| git | `D:\Program Files\Git` | |
| Python | 3.12.10 | required by the source-tree gates |
| GCC | 13.2.0 MinGW-UCRT (Strawberry) | second compiler, useful for portability checks |

CMake, Ninja, ccache and GCC all come from `D:\Strawberry\c\bin`. MSVC builds must run in an
environment where `vcvars64.bat` has been loaded; `scripts/build.ps1` does that itself.

## Build

```
.\scripts\build.ps1                        # configure if needed, then build
.\scripts\build.ps1 -Configure             # force a reconfigure (after a CMake change)
.\scripts\build.ps1 -Test                  # build, then run the suite
.\scripts\build.ps1 -Preset msvc-rwdi      # select a preset
ctest --test-dir build/msvc-rwdi -j 6      # run the suite directly
```

The script loads the MSVC environment into the session once, derives the oneMKL location, and passes
it to CMake. Loading the compiler environment and building in one invocation matters: a build started
in a session without `vcvars64` will fail in a way that looks like a missing header.

`ccache` is wired in through `CMAKE_CXX_COMPILER_LAUNCHER`, with `CMP0141` selecting MSVC `/Z7`
(embedded debug information) so that object files stay cacheable.

**A clean build costs about 40 minutes** (measured 2026-07-28: 39.6 min for the 291-target
`portable` build), down from the one-to-three hours earlier full builds cost depending on cache
state. The reduction is two measured fixes, not one guess: a shared precompiled header — Eigen was
93–99.9 % of the parse volume of a typical translation unit, so the earlier claim that the
header-heavy `core` caused the cost was wrong by an order of magnitude (it is about 7 %) — and
linking MKL dynamically, which removed some 57 minutes of static link time across the suite's
executables. The header-heavy `core` remains a real defect with its own stage. Changes that
trigger a full rebuild are still batched rather than verified one at a time.

## Intel oneMKL

Provides PARDISO plus BLAS and LAPACK.

**It is genuinely optional since Stage A2** (this corrects the 2026-07-27 measurement below, which
was true when taken). Every consumer solves through `katai::linsolve::DirectSolver`; the Eigen
backend (`SimplicialLDLT` for SPD, `SparseLU` otherwise) is always compiled because Eigen is
vendored, and the PARDISO adapter is added only when MKL is found. The backend is chosen by which
selection translation unit CMake compiles, so no source file conditions on `KATAI_WITH_MKL`. Only
the two PARDISO adapter unit tests are gated on MKL; everything else builds and runs in the
`portable` preset (`KATAI_WITH_MKL=OFF`), which is what makes the validation record reproducible
without a proprietary component. History, for the record: until 2026-07-28 PARDISO was the only
sparse solver in the tree, a build without MKL could not solve, and 72 of 116 tests were gated off.

Two things are worth knowing before touching this:

- **Installed on `D:`, not `C:`** — `D:\Program Files (x86)\Intel\oneAPI\mkl\latest` (version
  2026.0.0). The directory of the same name on `C:` is installer residue.
- **The top-level `setvars.bat` is broken for installation paths containing a space.** With
  `Program Files (x86)` in the path it does not quote its component `vars.bat` calls and fails with
  `'vars.bat' is not recognized`. Nothing here depends on it: `scripts/build.ps1` derives `MKLROOT`
  itself and passes `-D MKL_DIR=...` to CMake, which `MKLConfig.cmake` resolves from the environment.
  The component's own `mkl\latest\env\vars.bat` does work when called directly.

CMake selects **dynamic + sequential + lp64**. Dynamic is a measurement, not a taste: statically,
`test_solve.exe` was 130.7 MB and took 31.4 s to link against 4.6 MB and 5.0 s dynamically
(measured 2026-07-28), repeated across the ~76 test executables that solve — linking was 57 of the
~175 minutes a full static build cost. The price is that the MKL runtime must be locatable at load
time; the test suite gets it through CTest's own environment (`KATAI_MKL_RUNTIME_DIR` on each
test's `PATH`), and `build.ps1` prepends it for anything run by hand. Set `MKL_LINK=static` for an
executable that must run where MKL is not installed. Parallel PARDISO would need
`MKL_THREADING=tbb_thread`.

Evidence that the link and the numerics are correct: `tests/test_pardiso_smoke.cpp` solves a 5×5
sparse SPD system through PARDISO with `maxerr = 0.0`.

## Dependencies

Vendored drop-in under `third_party/` — Eigen, GLFW, glad, Dear ImGui — with versions and licences
recorded in `third_party/README.md`. Nothing is fetched at configure time, so a build works offline
and a given commit always compiles against exactly the sources it was tested with.

## Optional, not installed

`clang-cl` as a Visual Studio component, for a second opinion on warnings and for portability
checking. Not required by anything.
