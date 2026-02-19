param(
  [Parameter(Mandatory=$true)][string]$RepoRoot,
  [string]$Configuration = "Release",
  [string]$Platform = "x64"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Info($m) { Write-Host "[deps] $m" }
function Fail($m) { Write-Error "[deps] $m"; exit 1 }

$RepoRoot = (Resolve-Path $RepoRoot).Path
$LibDir   = Join-Path $RepoRoot "Lib"

function Import-VsDevCmdEnvironment {
  param(
    [Parameter(Mandatory=$true)][string]$VsDevCmdBat,
    [string]$Arch = "x64"
  )

  if (-not (Test-Path $VsDevCmdBat)) {
    throw "VsDevCmd not found: $VsDevCmdBat"
  }

  # Call VsDevCmd.bat, then dump environment via `set` and import into this PS session
  $cmd = 'call "{0}" -arch={1} >nul && set' -f $VsDevCmdBat, $Arch

  cmd.exe /s /c $cmd | ForEach-Object {
    if ($_ -match '^(.*?)=(.*)$') {
      Set-Item -Path ("Env:\" + $matches[1]) -Value $matches[2]
    }
  }
}

function Add-ToPathFront([string]$p) {
  if (-not [string]::IsNullOrWhiteSpace($p) -and (Test-Path $p)) {
    $env:PATH = $p + ";" + $env:PATH
  }
}

function Ensure-RustToolchain {

  # ------------------------------------------------------------
  # Check if it is just loaded
  # ------------------------------------------------------------
  if ((Get-Command rustc -ErrorAction SilentlyContinue) -and
      (Get-Command cargo -ErrorAction SilentlyContinue)) {
    return
  }

  Info "Rust not found in PATH. Searching common installation locations..."

  # ------------------------------------------------------------
  # Check common directories
  # ------------------------------------------------------------
  $candidateBins = @(
    (Join-Path $env:USERPROFILE ".cargo\bin"),
    (Join-Path $env:LOCALAPPDATA "Programs\Rust\bin"),
    "C:\Program Files\Rust\bin",
    "C:\Program Files (x86)\Rust\bin"
  )

  foreach ($bin in $candidateBins) {
    if (Test-Path (Join-Path $bin "rustc.exe")) {
      Add-ToPathFront $bin
      if ((Get-Command rustc -ErrorAction SilentlyContinue) -and
          (Get-Command cargo -ErrorAction SilentlyContinue)) {
        Info "Rust toolchain found at: $bin"
        return
      }
    }
  }

  # ------------------------------------------------------------
  # Fail if not found
  # ------------------------------------------------------------
  Write-Host ""
  Write-Host "============================================================" -ForegroundColor Red
  Write-Host "  Rust toolchain not detected" -ForegroundColor Red
  Write-Host "============================================================" -ForegroundColor Red
  Write-Host ""
  Write-Host "This project requires Rust (rustc + cargo) to build a dependency."
  Write-Host ""
  Write-Host "  Please install Rust via rustup:"
  Write-Host "     https://www.rust-lang.org/tools/install"
  Write-Host "   After installation please restart the console"
  Write-Host ""
  Write-Host "Expected default location:"
  Write-Host "  %USERPROFILE%\.cargo\bin"
  Write-Host ""
  Write-Host "You can verify manually with:"
  Write-Host "  rustc --version"
  Write-Host "  cargo --version"
  Write-Host ""
  exit 1
}

function Try-Load-MKLFromAutoProps {
  param(
    [Parameter(Mandatory=$true)]
    [string]$RepoRoot
  )

  $propsPath = Join-Path $RepoRoot "Windows\MKL.auto.props"

  if (-not (Test-Path $propsPath)) {
    return
  }

  Info "Found MKL.auto.props - loading MKLROOT..."

  try {
    $content = Get-Content $propsPath -Raw

    # Simple regex extraction (no XML namespace headaches)
    if ($content -match '<MKLROOT>(.*?)</MKLROOT>') {
      $mklRoot = $matches[1].Trim()

      if ($mklRoot) {
        $env:MKLROOT = $mklRoot
        Info "MKLROOT loaded from MKL.auto.props: $mklRoot"
        return
      }
    }

    Info "MKL.auto.props exists but MKLROOT not found."
    return
  }
  catch {
    Info "Failed to read MKL.auto.props: $_"
    return
  }
}

function Ensure-MKL {
  # ------------------------------------------------------------
  # Check if it is just loaded
  # ------------------------------------------------------------
  Try-Load-MKLFromAutoProps -RepoRoot $RepoRoot
  if ($env:MKLROOT) {
    $inc = Join-Path $env:MKLROOT "include\mkl.h"
    if (Test-Path $inc) {
      Info "MKL detected via MKLROOT=$env:MKLROOT"
      return
    }
  }

  Info "MKL not detected. Searching common installation locations..."

  # ------------------------------------------------------------
  # Check common directories
  # ------------------------------------------------------------
  $candidateRoots = @(
    # Standard oneAPI "latest"
    "C:\Program Files (x86)\Intel\oneAPI\mkl\latest",
    "C:\Program Files\Intel\oneAPI\mkl\latest",
    (Join-Path $env:USERPROFILE "Intel\oneAPI\mkl\latest"),

    # Some installs use explicit version folders instead of "latest"
    "C:\Program Files (x86)\Intel\oneAPI\mkl",
    "C:\Program Files\Intel\oneAPI\mkl",
    (Join-Path $env:USERPROFILE "Intel\oneAPI\mkl")
  )

  # helper: validate a potential MKL root
  function Test-MKLRoot([string]$root) {
    if (-not $root) { return $false }

    # If this is the parent "...mkl", try to resolve a child version folder or "latest"
    if (-not (Test-Path (Join-Path $root "include\mkl.h"))) {
      # Try common subfolders: latest, and the newest-looking version folder
      $latest = Join-Path $root "latest"
      if (Test-Path (Join-Path $latest "include\mkl.h")) {
        $script:FoundMKL = $latest
        return $true
      }

      try {
        $ver = Get-ChildItem -Path $root -Directory -ErrorAction SilentlyContinue |
               Sort-Object Name -Descending |
               Select-Object -First 1
        if ($ver -and (Test-Path (Join-Path $ver.FullName "include\mkl.h"))) {
          $script:FoundMKL = $ver.FullName
          return $true
        }
      } catch {
        return $false
      }

      return $false
    }

    $script:FoundMKL = $root
    return $true
  }

  $script:FoundMKL = $null
  foreach ($p in $candidateRoots) {
    if (Test-MKLRoot $p) {
      $env:MKLROOT = $script:FoundMKL
      Info "MKL found at: $env:MKLROOT"
      return
    }
  }

  # ------------------------------------------------------------
  # Fail if not found
  # ------------------------------------------------------------
  Write-Host ""
  Write-Host "============================================================" -ForegroundColor Red
  Write-Host "  Intel oneMKL not detected" -ForegroundColor Red
  Write-Host "============================================================" -ForegroundColor Red
  Write-Host ""
  Write-Host "This project requires Intel oneMKL (headers + libraries)."
  Write-Host ""
  Write-Host "How MKL is detected:"
  Write-Host "  - Environment variable MKLROOT, AND"
  Write-Host "  - The file %MKLROOT%\include\mkl.h must exist"
  Write-Host ""
  Write-Host "Fix options:"
  Write-Host "  1) Please install Intel oneAPI oneMKL:"
  Write-Host "     https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html"
  Write-Host ""
  Write-Host "  2) If MKL is already installed somewhere else, set MKLROOT manually:"
  Write-Host '     setx MKLROOT "C:\Program Files (x86)\Intel\oneAPI\mkl\latest"'
  Write-Host ""
  Write-Host "After installing / setting MKLROOT:"
  Write-Host "  - Restart Visual Studio (important)"
  Write-Host "  - Or open a new terminal so the environment refreshes"
  Write-Host ""
  Write-Host "Common install locations we looked in:"
  Write-Host "  C:\Program Files (x86)\Intel\oneAPI\mkl\latest"
  Write-Host "  C:\Program Files\Intel\oneAPI\mkl\latest"
  Write-Host "  %USERPROFILE%\Intel\oneAPI\mkl\latest"
  Write-Host ""
  exit 1
}

# -----------------------------
# 0) Initialize Visual Studio build environment (so cl/cmake/msbuild exist)
# -----------------------------
$vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vswhere)) { Fail "vswhere.exe not found at: $vswhere" }

$vsPath = & $vswhere -latest -products * -requires Microsoft.Component.MSBuild -property installationPath
if ([string]::IsNullOrWhiteSpace($vsPath)) { Fail "Visual Studio installation not found (vswhere returned empty)." }

$vsdev = Join-Path $vsPath "Common7\Tools\VsDevCmd.bat"

$vsArch = $Platform
if ($vsArch -ieq "Win32") { $vsArch = "x86" }
Import-VsDevCmdEnvironment -VsDevCmdBat $vsdev -Arch $vsArch


# -----------------------------
# 0) Tool checks
# -----------------------------
if (-not (Get-Command cmake -ErrorAction SilentlyContinue)) { Fail "cmake not found on PATH" }
if (-not (Get-Command msbuild -ErrorAction SilentlyContinue)) { Fail "msbuild not found on PATH (use VS Developer Prompt or ensure VS Build Tools installed)" }

Ensure-RustToolchain
Info ("Rust: " + (rustc --version))

# -----------------------------
# 1) Detect MKLROOT
# -----------------------------
Ensure-MKL
$genDir = Join-Path $RepoRoot "Windows"
$autoProps = Join-Path $genDir "MKL.auto.props"

$xml = @"
<Project xmlns="http://schemas.microsoft.com/developer/msbuild/2003">
  <PropertyGroup>
    <MKLROOT>$env:MKLROOT</MKLROOT>
  </PropertyGroup>
</Project>
"@

Set-Content -Path $autoProps -Value $xml -Encoding UTF8
Write-Host "[deps] Wrote $autoProps"

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
    cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX="..\..\Lib\LibCint" -DCMAKE_C_COMPILER=cl -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded" ..
    cmake --build . --config RELEASE
    cmake --install .
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $LibCintOut)) { Fail "LibCint build finished but output not found: $LibCintOut" }
  Info "LibCint OK"
}
# -----------------------------
# 4) Build OCC (if needed)
# -----------------------------
$OccOut = Join-Path $LibDir "occ\lib\occ_main.lib"
if (Test-Path $OccOut) {
  Info "OCC already built ($OccOut)"
} else {
  Info "Building LibCint..."
  Push-Location $RepoRoot
  try {
    cmake --workflow --preset windows-clang-cl
    cmake --install .\build-windows-clang-clf
  } finally {
    Pop-Location
  }
  if (-not (Test-Path $OccOut)) { Fail "OCC build finished but output not found: $OccOut" }
  Info "LibCint OK"
}
$stamp = Join-Path $RepoRoot "Windows\Build_Dependencies\deps.stamp"
New-Item -ItemType Directory -Force (Split-Path $stamp) | Out-Null
Set-Content -Path $stamp -Value ("OK " + (Get-Date).ToString("s")) -Encoding ASCII

Info "All dependencies are ready."
exit 0
