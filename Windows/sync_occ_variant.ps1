param(
  [Parameter(Mandatory=$true)][string]$RepoRoot,
  [Parameter(Mandatory=$true)][string]$Configuration
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Get-OccVariant([string]$Config) {
  switch ($Config.ToLowerInvariant()) {
    "debug" { return "debug" }
    "profile" { return "debug" }
    default { return "release" }
  }
}

$RepoRoot = (Resolve-Path $RepoRoot).Path
$variant = Get-OccVariant $Configuration
$sourceDir = Join-Path $RepoRoot ("Lib\occ_" + $variant)
$destDir = Join-Path $RepoRoot "Lib\occ"
$destStamp = Join-Path $destDir "occ_variant.stamp"
$sourceLib = Join-Path $sourceDir "lib\occ_main.lib"

if (-not (Test-Path $sourceLib)) {
  throw "Required OCC variant is missing: $sourceLib. Run Windows\\Build_Dependencies\\build_deps.ps1 first."
}

if ((Test-Path $destStamp) -and ((Get-Content $destStamp -Raw).Trim() -eq $variant)) {
  return
}

if (Test-Path $destDir) {
  Remove-Item -Recurse -Force $destDir
}

Copy-Item -Recurse -Force $sourceDir $destDir
Set-Content -Path $destStamp -Value $variant -Encoding ASCII
