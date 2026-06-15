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
$sourceLibDir = Join-Path $sourceDir "lib"

if (-not (Test-Path $sourceLib)) {
  throw "Required OCC variant is missing: $sourceLib. Run Windows\\Build_Dependencies\\build_deps.ps1 first."
}

if ((Test-Path $destStamp) -and ((Get-Content $destStamp -Raw).Trim() -eq $variant)) {
  $destStampInfo = Get-Item $destStamp
  $newestSourceLib = Get-ChildItem -Path $sourceLibDir -Filter *.lib -File |
    Sort-Object LastWriteTimeUtc -Descending |
    Select-Object -First 1
  if ($null -ne $newestSourceLib -and $destStampInfo.LastWriteTimeUtc -ge $newestSourceLib.LastWriteTimeUtc) {
    return
  }
}

if (Test-Path $destDir) {
  Remove-Item -Recurse -Force $destDir
}

Copy-Item -Recurse -Force $sourceDir $destDir
Set-Content -Path $destStamp -Value $variant -Encoding ASCII

# TBB debug builds name their artifacts tbb12_debug.{lib,dll}.
# Normalize to tbb12.{lib,dll} so all project configs can link uniformly.
if ($variant -eq "debug") {
  # TBB debug builds name their artifacts tbb12_debug.{lib,dll}.
  $debugLib = Join-Path $destDir "lib\tbb12_debug.lib"
  $canonLib = Join-Path $destDir "lib\tbb12.lib"
  if ((Test-Path $debugLib) -and -not (Test-Path $canonLib)) {
    Copy-Item $debugLib $canonLib
  }
  $debugDll = Join-Path $destDir "bin\tbb12_debug.dll"
  $canonDll = Join-Path $destDir "bin\tbb12.dll"
  if ((Test-Path $debugDll) -and -not (Test-Path $canonDll)) {
    Copy-Item $debugDll $canonDll
  }

  # fmt and spdlog debug builds use a 'd' suffix.
  foreach ($pair in @(@("fmtd.lib","fmt.lib"), @("spdlogd.lib","spdlog.lib"))) {
    $src = Join-Path $destDir ("lib\" + $pair[0])
    $dst = Join-Path $destDir ("lib\" + $pair[1])
    if ((Test-Path $src) -and -not (Test-Path $dst)) {
      Copy-Item $src $dst
    }
  }
}
