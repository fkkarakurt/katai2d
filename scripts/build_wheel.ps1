# KATAI 2D -- build, pack and VERIFY the katai2d wheel.
#
# The wheel is what `pip install katai2d` delivers: the _core extension built
# from the 'wheel' preset (Release, Eigen solver -- no proprietary component,
# no external DLL), the pure-Python package, and the `katai` console script
# with the same commands and exit codes as the native executable.
#
# The extension targets the CPython STABLE ABI, so one abi3 wheel serves every
# Python >= 3.12; no per-version build matrix.
#
# The proof is the last step, not the packaging: the wheel installs into a
# FRESH venv with a scrubbed environment, and there
#   * `katai info` runs from the console script,
#   * a corpus case validates and solves with the pinned numbers,
#   * the deliberately-past-collapse corpus case exits 5 (the honest refusal),
#   * the slope example runs its corpus-band guard,
#   * the anchored excavation runs its own -- a wall, an anchor and interfaces
#     built through `prj.structures`, which is the half of the surface an
#     import check cannot reach and the half a real job starts at.
# "The wheel works" means "the numbers are still right", not "import succeeded".
#
# Usage:  .\scripts\build_wheel.ps1 [-SkipVerify]
[CmdletBinding()]
param(
    [switch]$SkipVerify
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$work = Join-Path $repo "build\wheel-pkg"

# ------------------------------------------------------------------- build --
& (Join-Path $PSScriptRoot "build.ps1") -Preset wheel -Configure -Target katai_python_core

$staged = Join-Path $repo "build\wheel\python\katai"
if (-not (Test-Path (Join-Path $staged "__init__.py"))) {
    throw "no staged package under $staged"
}
$pyd = Get-ChildItem $staged -Filter "_core*.pyd" | Select-Object -First 1
if (-not $pyd) { throw "no _core extension under $staged" }

$python = Join-Path $env:LOCALAPPDATA "Python\bin\python3.exe"
if (-not (Test-Path $python)) { throw "python3 shim not found: $python" }

# The ABI story decides the tag: a cpXY name pins one interpreter; a bare or
# abi3 name is the stable ABI and one wheel covers >= 3.12.
if ($pyd.Name -match "cp(\d+)") {
    $tag = "cp$($Matches[1])"
    $stable = $false
    $requires = ">=3.10"
} else {
    $tag = "cp312"
    $stable = $true
    $requires = ">=3.12"
}

# ------------------------------------- the python-labelled tests, this config --
# The bindings the wheel ships are exercised where they were built: the DSL
# byte-identity, the structures suite, the schema coverage, the stub pin and every example.
ctest --test-dir (Join-Path $repo "build\wheel") -R "python" --output-on-failure -j 4
if ($LASTEXITCODE -ne 0) { throw "the python test subset failed in the wheel configuration" }

# ----------------------------------------------------------------- staging --
if (Test-Path $work) { Remove-Item -Recurse -Force $work }
$root = New-Item -ItemType Directory "$work\src" | Select-Object -ExpandProperty FullName
New-Item -ItemType Directory "$root\katai" | Out-Null
Copy-Item "$staged\*" "$root\katai" -Recurse -Exclude "__pycache__"
Copy-Item (Join-Path $repo "LICENSE") $root
Copy-Item (Join-Path $repo "NOTICE") $root

$version = & $python -c "import sys; sys.path.insert(0, r'$root'); import katai; print(katai.__version__)"
if (-not $version) { throw "cannot read katai.__version__ from the staged package" }

# BOM-less on purpose: PowerShell 5.1's Out-File utf8 prepends a BOM and the
# TOML parser refuses the file at byte one.
$noBom = New-Object System.Text.UTF8Encoding($false)

# A PyPI-facing readme: short and self-contained (the repo README's relative
# links would dangle on the package page).
[System.IO.File]::WriteAllText("$root\README.md", @"
# KATAI 2D

A verification-first 2D finite element engine for geotechnical engineering:
staged construction, groundwater flow, consolidation, seismic dynamics and
safety analysis over elastoplastic and creep constitutive models. Every
capability is validated against closed-form solutions, published benchmarks
or independent computation paths, and the validation record ships with the
source.

``pip install katai2d`` provides the ``katai`` Python package (an
engineer-facing project builder over the published facade) and the ``katai``
command line (``solve`` / ``validate`` / ``info``, documented exit codes).
Units are fixed: kN, m, day.

Documentation, examples and the verification matrix:
https://github.com/fkkarakurt/katai2d
"@, $noBom)

[System.IO.File]::WriteAllText("$root\pyproject.toml", @"
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "katai2d"
version = "$version"
description = "Verification-first 2D geotechnical finite element analysis -- engine, Python API and CLI"
readme = "README.md"
requires-python = "$requires"
license = { text = "Apache-2.0" }
authors = [{ name = "Fatih Kucukkarakurt", email = "fatihkucukkarakurt@gmail.com" }]
dependencies = ["numpy"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: Apache Software License",
    "Operating System :: Microsoft :: Windows",
    "Programming Language :: Python :: 3",
    "Topic :: Scientific/Engineering",
]

[project.urls]
Repository = "https://github.com/fkkarakurt/katai2d"

[project.scripts]
katai = "katai.cli:main"

[tool.setuptools]
packages = ["katai"]
license-files = ["LICENSE", "NOTICE"]

[tool.setuptools.package-data]
katai = ["*.pyd", "*.pyi"]
"@, $noBom)

# has_ext_modules makes setuptools tag the wheel for this platform instead of
# calling a wheel that ships a .pyd "pure"; py_limited_api tags the stable ABI.
$bdist = if ($stable) { "options={'bdist_wheel': {'py_limited_api': '$tag'}}," } else { "" }
[System.IO.File]::WriteAllText("$root\setup.py", @"
from setuptools import setup
from setuptools.dist import Distribution


class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True


setup($bdist distclass=BinaryDistribution)
"@, $noBom)

# -------------------------------------------------------------------- pack --
& $python -m pip install --quiet --disable-pip-version-check build
if ($LASTEXITCODE -ne 0) { throw "pip install build failed" }
& $python -m build --wheel --outdir "$work\dist" $root
if ($LASTEXITCODE -ne 0) { throw "wheel build failed" }
$wheel = Get-ChildItem "$work\dist\*.whl" | Select-Object -First 1
$mb = [math]::Round($wheel.Length / 1MB, 1)

# ------------------------------------------------------------------ verify --
if (-not $SkipVerify) {
    $venv = "$work\verify-venv"
    & $python -m venv $venv
    $vpy = "$venv\Scripts\python.exe"
    $vkatai = "$venv\Scripts\katai.exe"
    & $vpy -m pip install --quiet --disable-pip-version-check $wheel.FullName
    if ($LASTEXITCODE -ne 0) { throw "pip install of the wheel failed" }

    # A clean room: no KATAI_DLL_PATH, no PYTHONPATH.
    $env:KATAI_DLL_PATH = ""
    $env:PYTHONPATH = ""
    $corpus = Join-Path $repo "tests\corpus"

    & $vkatai info
    if ($LASTEXITCODE -ne 0) { throw "console script: `katai info` failed" }

    & $vkatai validate (Join-Path $corpus "kv-slp-001-griffiths-lane-slope.k2d")
    if ($LASTEXITCODE -ne 0) { throw "console script: validate exit $LASTEXITCODE, expected 0" }

    $res = Join-Path $work "verify.res"
    & $vkatai solve (Join-Path $corpus "kv-num-003-k0-geostatic-block.k2d") --out $res
    if ($LASTEXITCODE -ne 0) { throw "console script: solve exit $LASTEXITCODE, expected 0" }
    if (-not (Test-Path $res)) { throw "console script: --out wrote no .res" }

    & $vkatai solve (Join-Path $corpus "kv-fnd-010-prandtl-strip-footing.k2d")
    if ($LASTEXITCODE -ne 5) { throw "console script: past-collapse case exit $LASTEXITCODE, expected the honest 5" }

    & $vpy (Join-Path $repo "python\examples\slope_stability.py")
    if ($LASTEXITCODE -ne 0) { throw "the installed wheel failed the slope example's corpus-band guard" }

    # The structural surface, from the INSTALLED wheel. `import katai` succeeding says
    # nothing about whether prj.structures is in the build someone downloaded -- and a
    # published wheel whose documented API was absent is exactly what this gate exists
    # to catch. The example checks its own invariants and exits non-zero if they fail.
    & $vpy (Join-Path $repo "python\examples\anchored_excavation.py")
    if ($LASTEXITCODE -ne 0) { throw "the installed wheel failed the anchored-excavation example (prj.structures)" }

    Write-Host "verified: fresh venv, scrubbed environment -- console script exit codes 0/0/0/5, the corpus-band guard and the structural surface all hold" -ForegroundColor Green
}

Write-Host "OK -> $($wheel.FullName) ($mb MB, tag $tag$(if ($stable) { '-abi3' }))" -ForegroundColor Green
