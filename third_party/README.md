# third_party — vendored dependencies

Everything vendored here is permissively licensed and commercial-safe.
Vendoring is drop-in: the sources are copied into the repository and their
inner `.git` directories removed, so the build needs no network access and no
submodule step — the root CMake build compiles everything.

| Directory | Library | Version (vendored 2026-06-02) | License | Used for |
|-----------|---------|-------------------------------|---------|----------|
| `eigen/`  | Eigen | dev/master (≈3.4.90, WORLD=3) | MPL-2.0 | dense/small matrix algebra (header-only) |

## Updating

Shallow-clone the upstream, delete its `.git`, replace the directory. The
licenses ship with the sources; NOTICE at the repository root records the
attribution.

> Note: Intel oneMKL is **not** vendored — when present on the system it is
> linked behind the `LinearSolver` abstraction (`KATAI_WITH_MKL`); without it
> every build and the whole test suite run on the vendored Eigen backend.
