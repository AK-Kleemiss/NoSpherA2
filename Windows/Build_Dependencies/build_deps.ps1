param(
  [Parameter(Mandatory=$true)][string]$RepoRoot,
  [string]$Configuration = "Release",
  [string]$Platform = "x64",
  [switch]$InstallMKLIfMissing
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Info($m) { Write-Host "[deps] $m" }
function Fail($m) { Write-Error "[deps] $m"; exit 1 }

$RepoRoot = (Resolve-Path $RepoRoot).Path
$LibDir   = Join-Path $RepoRoot "Lib"

# -----------------------------
# 0) Tool checks
# -----------------------------
if (-not (Get-Command cmake -ErrorAction SilentlyContinue)) { Fail "cmake not found on PATH" }
if (-not (Get-Command msbuild -ErrorAction SilentlyContinue)) { Fail "msbuild not found on PATH (use VS Developer Prompt or ensure VS Build Tools installed)" }

if (-not (Get-Command rustc -ErrorAction SilentlyContinue)) { Fail "rustc not found. Install Rust: https://www.rust-lang.org/tools/install" }
Info ("Rust: " + (rustc --version))

# -----------------------------
# 1) Detect MKLROOT
# -----------------------------
function Find-MKLRoot {
  $candidates = @()
  if ($env:MKLROOT) { $candidates += $env:MKLROOT }
  $candidates += "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
  $candidates += "C:\Program Files\Intel\oneAPI\mkl\latest"
  $candidates += (Join-Path $env:USERPROFILE "Intel\oneAPI\mkl\latest")

  foreach ($p in $candidates) {
    if ($p -and (Test-Path (Join-Path $p "include\mkl.h"))) { return $p }
  }
  return $null
}

function Download-File($Url, $OutFile) {
  Write-Host "[deps] Downloading $Url"
  Write-Host "[deps] -> $OutFile"

  # Prefer curl.exe on Windows if available (matches your Makefile)
  $curl = (Get-Command curl.exe -ErrorAction SilentlyContinue)
  if ($curl) {
    $p = Start-Process -FilePath $curl.Source -ArgumentList @("-L", "-o", $OutFile, $Url) -Wait -PassThru
    if ($p.ExitCode -ne 0) { throw "curl failed with exit code $($p.ExitCode)" }
    return
  }

  # Fallback
  Invoke-WebRequest -Uri $Url -OutFile $OutFile -UseBasicParsing
}

function Install-MKL {
  param(
    [Parameter(Mandatory=$true)][string]$RepoRoot
  )
  $url = "https://registrationcenter-download.intel.com/akdlm/IRC_NAS/ae472ff5-aa01-4a72-a452-ce7b559ef041/intel-onemkl-2025.3.1.10_offline.exe"
  $exeName = "intel-onemkl-2025.3.1.10-windows.exe"
  $exePath = Join-Path $RepoRoot $exeName

  Write-Host "[deps] MKL not found, downloading/installing Intel oneMKL for Windows"
  if (-not (Test-Path $exePath)) {
    Download-File -Url $url -OutFile $exePath
  } else {
    Write-Host "[deps] Installer already present: $exePath"
  }
  
  Write-Host "[deps] Installing MKL, this will take some time! DO NOT CLOSE THE TERMINAL!"
  # Run silent install + accept EULA
  $args = @("-a", "-s", "--eula=accept")
  $p = Start-Process -FilePath $exePath -ArgumentList $args -Wait -PassThru
  if ($p.ExitCode -ne 0) {
    throw "oneMKL installer failed with exit code $($p.ExitCode)"
  }

  Write-Host "[deps] oneMKL installer finished successfully"
}

$found = Find-MKLRoot
if (-not $found) {
  if ($InstallMKLIfMissing) {
    Install-MKL -RepoRoot $RepoRoot
    $found = Find-MKLRoot
  }
}

if (-not $found) {
  throw "MKL not found. Set MKLROOT or install Intel oneMKL."
}

$genDir = Join-Path $RepoRoot "Windows"
$autoProps = Join-Path $genDir "MKL.auto.props"

$xml = @"
<Project xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
  <PropertyGroup>
    <MKLROOT>$found</MKLROOT>
  </PropertyGroup>
</Project>
"@

Set-Content -Path $autoProps -Value $xml -Encoding UTF8
Write-Host "[deps] Wrote $autoProps"

$env:MKLROOT = $found
Info "Found MKLROOT=$env:MKLROOT"

# -----------------------------
# 2) Build featomic (if needed)
# -----------------------------
$FeatomicOut = Join-Path $LibDir "featomic_install\lib\metatensor.lib"
if (Test-Path $FeatomicOut) {
  Info "featomic already built ($FeatomicOut)"
} else {
  Info "Building featomic..."
  $src = Join-Path $RepoRoot "featomic\featomic"
  $bld = Join-Path $src "build"

  if (Test-Path $bld) { Remove-Item -Recurse -Force $bld }
  New-Item -ItemType Directory $bld | Out-Null
  Push-Location $bld
  try {
    cmake -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="..\..\..\Lib\featomic_install" ..
    msbuild -nologo .\featomic.sln /p:Configuration=Release /p:Platform=$Platform
    msbuild -nologo .\INSTALL.vcxproj /p:Configuration=Release /p:Platform=$Platform
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $FeatomicOut)) { Fail "featomic build finished but output not found: $FeatomicOut" }
  Info "featomic OK"
}

# -----------------------------
# 3) Build LibCint (if needed)
# -----------------------------
$LibCintOut = Join-Path $LibDir "LibCint\lib\cint.lib"
if (Test-Path $LibCintOut) {
  Info "LibCint already built ($LibCintOut)"
} else {
  Info "Building LibCint..."
  $src = Join-Path $RepoRoot "libcint"
  $bld = Join-Path $src "build"
  if (-not (Test-Path $bld)) { New-Item -ItemType Directory $bld | Out-Null }

  Push-Location $bld
  try {
    cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX="..\..\Lib\LibCint" -DCMAKE_C_COMPILER=cl ..
    cmake --build . --config RELEASE
    cmake --install .
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $LibCintOut)) { Fail "LibCint build finished but output not found: $LibCintOut" }
  Info "LibCint OK"
}

$stamp = Join-Path $RepoRoot "Windows\Build_Dependencies\deps.stamp"
New-Item -ItemType Directory -Force (Split-Path $stamp) | Out-Null
Set-Content -Path $stamp -Value ("OK " + (Get-Date).ToString("s")) -Encoding ASCII

Info "All dependencies are ready."
exit 0