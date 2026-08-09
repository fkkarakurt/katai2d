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

# Windows PowerShell 5.1 still defaults to TLS 1.0 on some builds; GitHub serves 1.2+ only.
try { [Net.ServicePointManager]::SecurityProtocol =
          [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12 } catch {}

# Every web call goes through here, for two reasons learned the hard way.
#
# -UseBasicParsing: without it, Windows PowerShell 5.1 hands the response to the Internet
# Explorer engine to build a DOM, which stops with a "Script Execution Risk" prompt and then
# refuses. An installer that asks a question mid-run is an installer that fails in a script.
#
# The response body is decoded EXPLICITLY. For a text/plain asset PowerShell returns a string,
# but for application/octet-stream -- which is what GitHub serves a .sha256 attachment as -- it
# returns a BYTE ARRAY, and every string operation then silently works on bytes: splitting it
# yielded 55, the byte value of the character "7", which was compared against a real hash and
# reported as "expected 55". Decoding before parsing is the fix.
function Get-TextFromUrl([string]$url) {
    $r = Invoke-WebRequest $url -UseBasicParsing
    if ($r.Content -is [byte[]]) { return [Text.Encoding]::UTF8.GetString($r.Content) }
    return [string]$r.Content
}

if ($ZipPath) {
    $zip = $ZipPath
    if (-not (Test-Path $zip)) { throw "zip not found: $zip" }
} else {
    # /releases/latest ignores pre-releases; while every published version is
    # a pre-release, fall back to the newest release in the full list.
    $release = $null
    try { $release = Invoke-RestMethod "https://api.github.com/repos/$repo/releases/latest" } catch {}
    if (-not $release) {
        $release = @(Invoke-RestMethod "https://api.github.com/repos/$repo/releases") | Select-Object -First 1
    }
    if (-not $release) { throw "no published release in $repo" }
    $asset = $release.assets | Where-Object { $_.name -like "katai2d-*-win64.zip" } | Select-Object -First 1
    if (-not $asset) { throw "no katai2d-*-win64.zip asset in the latest release" }
    $zip = Join-Path $env:TEMP $asset.name
    Invoke-WebRequest $asset.browser_download_url -OutFile $zip -UseBasicParsing

    # Verify against the published .sha256 when the release carries one;
    # otherwise against the SHA-256 digest the release API reports for the
    # asset. Either way the download is checked before it is unpacked.
    $shaAsset = $release.assets | Where-Object { $_.name -eq "$($asset.name).sha256" } | Select-Object -First 1
    $expected = $null
    if ($shaAsset) {
        # "<64 hex>  <filename>" -- take the hex word, whatever whitespace or line ending
        # follows it. Refuse anything that is not 64 hex characters rather than comparing
        # against a fragment: a checksum that cannot be read is not a checksum that passed.
        $shaText = Get-TextFromUrl $shaAsset.browser_download_url
        if ($shaText -match '([0-9a-fA-F]{64})') {
            $expected = $Matches[1].ToLower()
        } else {
            throw "the published checksum file $($shaAsset.name) does not contain a SHA-256: '$shaText'"
        }
    } elseif ("$($asset.digest)" -match '^sha256:([0-9a-fA-F]{64})$') {
        $expected = $Matches[1].ToLower()
    }
    if ($expected) {
        $actual = (Get-FileHash $zip -Algorithm SHA256).Hash.ToLower()
        if ($actual -ne $expected) { throw "SHA-256 mismatch: expected $expected, got $actual" }
        Write-Host "sha256 verified: $actual"
    } else {
        Write-Warning "no published SHA-256 for $($asset.name); installing unverified"
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
