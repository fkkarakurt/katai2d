# KATAI 2D -- package the self-contained CLI for end users.
# Builds the 'release' preset (optimized, static CRT, Eigen solver backend:
# no proprietary component and no runtime dependency at all -- every kernel
# module links statically, so the distribution is katai.exe plus the license
# texts) and produces dist/katai2d-<version>-win64.zip with its SHA-256.
#
# Usage:  .\scripts\package_cli.ps1
[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot

# Version from the one source of truth (<katai/api/version.hpp>).
$versionHpp = Get-Content (Join-Path $root "kernel\api\include\katai\api\version.hpp") -Raw
if ($versionHpp -notmatch 'kVersion\s*=\s*"([^"]+)"') { throw "kVersion not found in version.hpp" }
$version = $Matches[1]

& (Join-Path $PSScriptRoot "build.ps1") -Preset release -Configure -Target katai_cli

$exe = Join-Path $root "build\release\bin\katai.exe"
if (-not (Test-Path $exe)) { throw "katai.exe was not produced: $exe" }

$name  = "katai2d-$version-win64"
$dist  = Join-Path $root "dist"
$stage = Join-Path $dist $name
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
New-Item -ItemType Directory -Force $stage | Out-Null

Copy-Item $exe $stage
Copy-Item (Join-Path $root "LICENSE") $stage
Copy-Item (Join-Path $root "NOTICE") $stage
@"
KATAI 2D $version -- command-line front end (win64, self-contained)

Usage:
  katai validate model.k2d    schema + physics validation, no solve
  katai solve model.k2d       run the staged analysis, write results
  katai info                  version and build information

This build carries no runtime dependency: the solver backend is Eigen and
the C runtime is linked statically, so the executable runs as-is on a
machine with nothing installed. Exit codes and the .k2d input format are
documented in the source repository:

  https://github.com/fkkarakurt/katai2d

Licensed under the Apache License 2.0 (see LICENSE); third-party notices
in NOTICE.
"@ | Out-File -FilePath (Join-Path $stage "README.txt") -Encoding utf8

$zip = Join-Path $dist "$name.zip"
if (Test-Path $zip) { Remove-Item -Force $zip }
Compress-Archive -Path $stage -DestinationPath $zip
$sha = (Get-FileHash $zip -Algorithm SHA256).Hash.ToLower()
$sha + "  $name.zip" | Out-File -FilePath "$zip.sha256" -Encoding ascii
Write-Host "packaged: $zip" -ForegroundColor Green
Write-Host "sha256:   $sha" -ForegroundColor Green
