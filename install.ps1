# KATAI 2D -- install the self-contained CLI in one line:
#
#   irm https://raw.githubusercontent.com/fkkarakurt/katai2d/main/install.ps1 | iex
#
# Downloads the latest release zip, verifies its published SHA-256, unpacks it
# to %LOCALAPPDATA%\Programs\katai2d and puts that directory on the user PATH.
# No admin rights, no dependencies -- katai.exe is a single static file.
#
# Parameters (for testing or offline installs):
#   -ZipPath <file>   install from a local zip instead of downloading
#   -NoPath           skip the PATH update
[CmdletBinding()]
param(
    [string]$ZipPath = "",
    [switch]$NoPath
)

$ErrorActionPreference = "Stop"
$repo = "fkkarakurt/katai2d"
$dest = Join-Path $env:LOCALAPPDATA "Programs\katai2d"

if ($ZipPath) {
    $zip = $ZipPath
    if (-not (Test-Path $zip)) { throw "zip not found: $zip" }
} else {
    $release = Invoke-RestMethod "https://api.github.com/repos/$repo/releases/latest"
    $asset = $release.assets | Where-Object { $_.name -like "katai2d-*-win64.zip" } | Select-Object -First 1
    if (-not $asset) { throw "no katai2d-*-win64.zip asset in the latest release" }
    $zip = Join-Path $env:TEMP $asset.name
    Invoke-WebRequest $asset.browser_download_url -OutFile $zip

    # Verify against the published .sha256 when the release carries one.
    $shaAsset = $release.assets | Where-Object { $_.name -eq "$($asset.name).sha256" } | Select-Object -First 1
    if ($shaAsset) {
        $expected = ((Invoke-WebRequest $shaAsset.browser_download_url).Content -split '\s+')[0].Trim().ToLower()
        $actual = (Get-FileHash $zip -Algorithm SHA256).Hash.ToLower()
        if ($actual -ne $expected) { throw "SHA-256 mismatch: expected $expected, got $actual" }
        Write-Host "sha256 verified: $actual"
    }
}

# Unpack; the zip contains one katai2d-<version>-win64 folder -- flatten it.
$stage = Join-Path $env:TEMP "katai2d-install-$PID"
if (Test-Path $stage) { Remove-Item -Recurse -Force $stage }
Expand-Archive $zip $stage
$inner = Get-ChildItem $stage -Directory | Select-Object -First 1
$src = if ($inner) { $inner.FullName } else { $stage }

New-Item -ItemType Directory -Force $dest | Out-Null
Copy-Item "$src\*" $dest -Recurse -Force
Remove-Item -Recurse -Force $stage

if (-not $NoPath) {
    $userPath = [Environment]::GetEnvironmentVariable("Path", "User")
    if ($userPath -notlike "*$dest*") {
        [Environment]::SetEnvironmentVariable("Path", "$userPath;$dest", "User")
        $env:PATH = "$env:PATH;$dest"
        Write-Host "added to user PATH: $dest"
    }
}

& (Join-Path $dest "katai.exe") info
Write-Host "installed -> $dest" -ForegroundColor Green
Write-Host "open a NEW terminal and run: katai info" -ForegroundColor Green
