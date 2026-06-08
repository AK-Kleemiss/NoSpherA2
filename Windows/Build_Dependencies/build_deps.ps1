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
    cmake -DCMAKE_BUILD_TYPE=Release -DFEATOMIC_FETCH_METATENSOR=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="..\..\..\Lib\featomic_install" ..
    cmake --build . --config $Configuration --parallel
    cmake --install . --config $Configuration
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
    cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DPYPZPX="ON" -DCMAKE_INSTALL_PREFIX="..\..\Lib\LibCint" -DCMAKE_C_COMPILER=cl -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded" ..
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
# $OccCrtStamp records which CRT the OCC libs were built with (/MT or /MD).
# If the stamp is absent or doesn't say "MD", force a full rebuild so that
# the tbbmalloc_proxy IAT-patching mechanism can work.
$OccOut      = Join-Path $LibDir "occ\lib\occ_main.lib"
$OccCrtStamp = Join-Path $LibDir "occ\lib\occ_crt.stamp"
$occNeedsBuild = (-not (Test-Path $OccOut)) -or
                 (-not (Test-Path $OccCrtStamp)) -or
                 ((Get-Content $OccCrtStamp -ErrorAction SilentlyContinue) -ne "MD")
if (-not $occNeedsBuild) {
  Info "OCC already built with /MD ($OccOut)"
} else {
  # Remove stale libs AND the cmake build cache.
  # Changing /MT → /MD requires a clean configure; an incremental cmake build
  # will leave stale /MT object files that cause loader failures in the DLL.
  if (Test-Path (Join-Path $LibDir "occ\lib")) {
    Remove-Item -Recurse -Force (Join-Path $LibDir "occ\lib")
  }
  $occBuildDir = Join-Path $RepoRoot "build-windows-clang-cl"
  if (Test-Path $occBuildDir) {
    Info "Removing stale OCC build cache ($occBuildDir)..."
    Remove-Item -Recurse -Force $occBuildDir
  }
  Info "Building OCC..."
  Push-Location $RepoRoot
  try {
    cmake --workflow --preset windows-clang-cl
    cmake --install .\build-windows-clang-cl
    $tbbLib = Get-ChildItem `
      -Path (Join-Path $RepoRoot "build-windows-clang-cl") `
      -Recurse `
      -Filter "tbb12.lib" `
      -ErrorAction SilentlyContinue |
      Select-Object -First 1 -ExpandProperty FullName
    Copy-Item $tbbLib -Destination (Join-Path $LibDir "occ\lib\tbb12.lib") -Force
  } finally {
    Pop-Location
  }
  if (-not (Test-Path $OccOut)) { Fail "OCC build finished but output not found: $OccOut" }
  Set-Content -Path $OccCrtStamp -Value "MD" -Encoding ASCII
  Info "OCC OK"
}

# -----------------------------
# 4a) Build tbbmalloc + tbbmalloc_proxy (shared build, needed for allocator proxy)
# -----------------------------
# tbbmalloc_proxy can only be built as a shared library (its CMakeLists.txt
# returns early when BUILD_SHARED_LIBS is OFF). We configure a separate
# build tree for the already-fetched oneTBB source, build only the two
# malloc targets, and copy the import libs + DLLs into Lib/occ/ and the
# Windows output directories.
$TbbMallocProxyOut = Join-Path $LibDir "occ\lib\tbbmalloc_proxy.lib"
if (Test-Path $TbbMallocProxyOut) {
  Info "tbbmalloc_proxy already built ($TbbMallocProxyOut)"
} else {
  $tbbSrc    = Join-Path $RepoRoot "build-windows-clang-cl\_deps\onetbb-src"
  $tbbBldDir = Join-Path $RepoRoot "build-tbb-malloc"

  if (-not (Test-Path $tbbSrc)) {
    Fail "oneTBB source not found at $tbbSrc - run the OCC build first so CMake fetches it."
  }

  Info "Building tbbmalloc + tbbmalloc_proxy (shared)..."
  if (Test-Path $tbbBldDir) { Remove-Item -Recurse -Force $tbbBldDir }
  New-Item -ItemType Directory $tbbBldDir | Out-Null

  Push-Location $tbbBldDir
  try {
    # Use the same compiler as the OCC build; MultiThreadedDLL (/MD) is required
    # for shared TBB on MSVC — the proxy DLL uses its own heap internally.
    cmake -G Ninja $tbbSrc `
          -DCMAKE_C_COMPILER=clang-cl.exe `
          -DCMAKE_CXX_COMPILER=clang-cl.exe `
          -DCMAKE_C_COMPILER_TARGET=x86_64-pc-windows-msvc `
          -DCMAKE_CXX_COMPILER_TARGET=x86_64-pc-windows-msvc `
          -DCMAKE_BUILD_TYPE=Release `
          -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreadedDLL `
          "-DCMAKE_CXX_FLAGS=/EHsc" `
          "-DCMAKE_C_FLAGS=/EHsc" `
          -DBUILD_SHARED_LIBS=ON `
          -DTBB_TEST=OFF `
          -DTBB_STRICT=OFF `
          -DTBB_DISABLE_HWLOC_AUTOMATIC_SEARCH=OFF

    cmake --build . --target tbbmalloc tbbmalloc_proxy
  } finally {
    Pop-Location
  }

  # Locate the built artifacts (filter out CMakeFiles stubs)
  $srcMallocLib = Get-ChildItem $tbbBldDir -Recurse -Filter "tbbmalloc.lib" `
      -ErrorAction SilentlyContinue |
      Where-Object { $_.FullName -notmatch "CMakeFiles" } |
      Select-Object -First 1 -ExpandProperty FullName

  $srcProxyLib = Get-ChildItem $tbbBldDir -Recurse -Filter "tbbmalloc_proxy.lib" `
      -ErrorAction SilentlyContinue |
      Where-Object { $_.FullName -notmatch "CMakeFiles" } |
      Select-Object -First 1 -ExpandProperty FullName

  $srcMallocDll = Get-ChildItem $tbbBldDir -Recurse -Filter "tbbmalloc.dll" `
      -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName

  $srcProxyDll = Get-ChildItem $tbbBldDir -Recurse -Filter "tbbmalloc_proxy.dll" `
      -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName

  if (-not $srcMallocLib)  { Fail "tbbmalloc.lib not found after build" }
  if (-not $srcProxyLib)   { Fail "tbbmalloc_proxy.lib not found after build" }
  if (-not $srcMallocDll)  { Fail "tbbmalloc.dll not found after build" }
  if (-not $srcProxyDll)   { Fail "tbbmalloc_proxy.dll not found after build" }

  # Import libs → Lib/occ/lib/ (used by the VS linker)
  Copy-Item $srcMallocLib  -Destination (Join-Path $LibDir "occ\lib\tbbmalloc.lib")       -Force
  Copy-Item $srcProxyLib   -Destination (Join-Path $LibDir "occ\lib\tbbmalloc_proxy.lib")  -Force

  # Runtime DLLs → all Windows output directories
  $outDirs = @(
    (Join-Path $RepoRoot "Windows\x64\Release"),
    (Join-Path $RepoRoot "Windows\x64\Debug"),
    (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Debug")
  )
  foreach ($d in $outDirs) {
    New-Item -ItemType Directory -Force $d | Out-Null
    Copy-Item $srcMallocDll  -Destination $d -Force
    Copy-Item $srcProxyDll   -Destination $d -Force
  }

  Info "tbbmalloc_proxy OK"
}

# -----------------------------
# 4b) Deploy MKL runtime DLLs (needed when OCC is built with /MD)
# -----------------------------
# With /MD, the MKL static threading layer emits a dynamic dependency on
# mkl_intel_thread.2.dll and mkl_core.2.dll.  Deploy them alongside the other
# runtime DLLs so the test runner can find them.
$mklBinDir = Join-Path $env:MKLROOT "bin"
$mklDllsNeeded = @("mkl_intel_thread.2.dll", "mkl_core.2.dll")
$outDirsMkl = @(
  (Join-Path $RepoRoot "Windows\x64\Release"),
  (Join-Path $RepoRoot "Windows\x64\Debug"),
  (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Debug")
)
foreach ($mklDll in $mklDllsNeeded) {
  $src = Join-Path $mklBinDir $mklDll
  if (-not (Test-Path $src)) {
    Write-Warning "MKL DLL not found: $src - skipping deployment (test runner may fail to load NoSpherA2_DLL.dll)"
    continue
  }
  foreach ($d in $outDirsMkl) {
    New-Item -ItemType Directory -Force $d | Out-Null
    Copy-Item $src -Destination $d -Force
  }
  Info "Deployed $mklDll"
}

$stamp = Join-Path $RepoRoot "Windows\Build_Dependencies\deps.stamp"
New-Item -ItemType Directory -Force (Split-Path $stamp) | Out-Null
Set-Content -Path $stamp -Value ("OK " + (Get-Date).ToString("s")) -Encoding ASCII

# -----------------------------
# 4) Build BasisSetConverter (if needed)
# -----------------------------
$converterExe = Join-Path $RepoRoot "Lib/BasisSetConverter.exe"
if (Test-Path $converterExe) {
  Info "BasisSetConverter already built ($converterExe)"
} else {
  $src = Join-Path $RepoRoot "BasisSetGenerator"
  $bld = Join-Path $src "build"
  if (-not (Test-Path $bld)) { New-Item -ItemType Directory $bld | Out-Null }

  Push-Location $src
  try {
    cmake -S . -B build
    cmake --build build --config Release
    cmake --install build --config Release --prefix "../Lib"
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $converterExe)) { Fail "BasisSetConverter build finished but executable not found: $converterExe" }
  Info "BasisSetConverter OK"
}

Info "All dependencies are ready."
exit 0
