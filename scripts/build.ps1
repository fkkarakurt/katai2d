# KATAI 2D — dev build script (T-ENV-2).
# Loads the MSVC environment (vcvars64.bat), then configure + build via a CMake preset.
# Usage:  .\scripts\build.ps1 [-Preset msvc-rwdi] [-Configure] [-Test] [-Target <name>]
[CmdletBinding()]
param(
    [string]$Preset = "msvc-rwdi",
    [switch]$Configure,   # force the configure step (first setup / CMake change)
    [switch]$Test,        # run ctest after the build
    [int]$Jobs = 6,       # ctest parallelism
    [string]$Target = ""  # optional single build target (default: everything)
)

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot

$vcvars = "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path $vcvars)) { throw "vcvars64.bat bulunamadi: $vcvars" }

# ccache caches NOTHING here unless it is told to tolerate the precompiled header: without this,
# every PCH-using compilation is reported "Could not use precompiled header" and recompiled in
# full. Measured 2026-08-13: 66% of calls uncacheable, 5.3% hit rate; with it, a rebuilt
# translation unit went from ~50 s to 5 s. Exported per build rather than written into the
# developer's machine-wide ccache.conf, so a scripted build is correct on a fresh machine.
if (-not $env:CCACHE_SLOPPINESS) { $env:CCACHE_SLOPPINESS = "pch_defines,time_macros" }

# Import the vcvars64 environment into this PowerShell session (once).
# 2>nul: even without vswhere.exe on PATH, VS finds itself from its own location
# (cl.exe works); we suppress that harmless stderr noise.
if (-not $env:KATAI_VCVARS_LOADED) {
    cmd /c "`"$vcvars`" >nul 2>nul && set" | ForEach-Object {
        if ($_ -match '^([^=]+)=(.*)$') { Set-Item -Path "env:$($matches[1])" -Value $matches[2] }
    }
    $env:KATAI_VCVARS_LOADED = "1"
}

# --- Intel oneMKL yolu (T-ENV-1 / Karar D6) --------------------------------
# The top-level setvars.bat is broken on an install path with spaces
# ("Program Files (x86)") -- it does not quote its vars.bat calls. So we derive
# MKLROOT ourselves and pass MKL_DIR to CMake via -D. Candidate paths in priority order:
$mklDir = $null
$oneapiMkl = if ($env:ONEAPI_ROOT) { Join-Path $env:ONEAPI_ROOT "mkl\latest" } else { $null }
# Force an array with @(...): if a single candidate returns, this prevents
# Where-Object collapsing to a string and [0] yielding the first CHARACTER.
$mklCandidates = @(
    @(
        $env:MKLROOT,
        $oneapiMkl,
        "D:\Program Files (x86)\Intel\oneAPI\mkl\latest",
        "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
    ) | Where-Object { $_ -and (Test-Path (Join-Path $_ "lib\cmake\mkl\MKLConfig.cmake")) }
)
if ($mklCandidates) {
    $mklRoot = $mklCandidates[0]
    $mklDir  = Join-Path $mklRoot "lib\cmake\mkl"
    $env:MKLROOT = $mklRoot
    # Keep the MKL DLLs on PATH in case of a dynamic link (harmless when static).
    $env:PATH = (Join-Path $mklRoot "bin") + ";" + $env:PATH
    Write-Host "oneMKL: $mklRoot" -ForegroundColor Cyan
} else {
    Write-Host "oneMKL bulunamadi -> Eigen cozucuyle devam." -ForegroundColor Yellow
}

# --- Python for the source-tree gates --------------------------------------
# check_language / check_architecture need a real interpreter. From a bare shell
# CMake's find_package can land on the WindowsApps store alias, reject it, and
# silently drop both gates from the test set -- measured: the portable configure
# registered 115 tests instead of 117. Resolve a real python here and pass it
# explicitly, so the gate set does not depend on the caller's PATH.
# Prefer the VERSIONLESS python3 shim: the PyManager install auto-updates and
# DELETES versioned shims when it does (measured 2026-08-02: python3.12.exe
# vanished mid-session when the default moved to 3.14; a versioned pin dies with
# its version, the python3 shim survives updates).
$python = @(Get-Command python3.exe, python.exe -ErrorAction SilentlyContinue |
            Where-Object { $_.Source -notmatch 'WindowsApps' })
if (-not $python) {
    $localPy = "$env:LOCALAPPDATA\Python\bin\python3.exe"
    if (Test-Path $localPy) { $python = @(@{ Source = $localPy }) }
}

Push-Location $root
try {
    $buildDir = Join-Path $root "build\$Preset"
    if ($Configure -or -not (Test-Path (Join-Path $buildDir "CMakeCache.txt"))) {
        $cfgArgs = @("--preset", $Preset)
        if ($mklDir) { $cfgArgs += @("-D", "MKL_DIR=$mklDir") }
        # One FindPython module across the tree (nanobind requires the NEW module,
        # and mixing FindPython with FindPython3 corrupts each other's cache).
        if ($python) { $cfgArgs += @("-D", "Python_EXECUTABLE=$($python[0].Source)") }
        cmake @cfgArgs
        if ($LASTEXITCODE -ne 0) { throw "configure basarisiz" }
    }
    $buildArgs = @("--build", "--preset", $Preset)
    if ($Target) { $buildArgs += @("--target", $Target) }
    cmake @buildArgs
    if ($LASTEXITCODE -ne 0) { throw "build failed" }

    if ($Test) {
        # Parallel, matching how the suite is normally run. Serially the same 119
        # tests take about 2.5x as long (441 s vs 175 s measured), which is enough
        # of a difference to discourage running the full suite.
        ctest --test-dir $buildDir --output-on-failure -j $Jobs
        if ($LASTEXITCODE -ne 0) { throw "test basarisiz" }
    }
    # A preset may build no GUI or CLI target at all; name only the front-end
    # binaries that actually exist (never say "OK" to a path that is not there).
    $fronts = @("katai_app.exe", "katai.exe") |
        ForEach-Object { Join-Path "$buildDir\bin" $_ } |
        Where-Object { Test-Path $_ }
    if ($fronts) { Write-Host "OK -> $($fronts -join ', ')" -ForegroundColor Green }
    else { Write-Host "OK -> $buildDir (library and test targets)" -ForegroundColor Green }
}
finally { Pop-Location }
