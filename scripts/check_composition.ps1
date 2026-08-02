# KATAI 2D -- walk the composition matrix (roadmap section 4.3, rule 2):
# every configuration must configure, build and pass its own test subset,
# because a configuration nobody built is a configuration nobody can trust.
#
# Rows buildable today: the default (full: engine, schema, IO, jobs, api,
# cli, python), portable (no proprietary component) and engine (physics
# without schema or IO). python/ skips cleanly on machines without a Python
# development environment -- the row still builds and its subset still runs.
#
# Usage:  .\scripts\check_composition.ps1 [-Rows msvc-rwdi,portable,engine] [-Jobs 6]
# Environment (vcvars, MKL discovery, real Python) comes from build.ps1 so the
# matrix and the daily build can never drift apart.
[CmdletBinding()]
param(
    [string[]]$Rows = @("msvc-rwdi", "portable", "engine"),
    [int]$Jobs = 6
)

$ErrorActionPreference = "Stop"
$build = Join-Path $PSScriptRoot "build.ps1"
$results = @()

foreach ($row in $Rows) {
    Write-Host "`n=== composition row: $row ===" -ForegroundColor Cyan
    $sw = [System.Diagnostics.Stopwatch]::StartNew()
    try {
        & $build -Preset $row -Test -Jobs $Jobs
        $results += [pscustomobject]@{ Row = $row; Result = "PASS"; Minutes = [math]::Round($sw.Elapsed.TotalMinutes, 1) }
    } catch {
        $results += [pscustomobject]@{ Row = $row; Result = "FAIL"; Minutes = [math]::Round($sw.Elapsed.TotalMinutes, 1) }
        Write-Host "row '$row' failed: $_" -ForegroundColor Red
    }
}

Write-Host "`n=== composition matrix ===" -ForegroundColor Cyan
$results | Format-Table -AutoSize | Out-String | Write-Host

if ($results | Where-Object { $_.Result -ne "PASS" }) {
    exit 1
}
Write-Host "every buildable row of the matrix configures, builds and passes." -ForegroundColor Green
