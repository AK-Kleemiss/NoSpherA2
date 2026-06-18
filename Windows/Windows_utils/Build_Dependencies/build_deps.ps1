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
$WindowsUtilsDir = Join-Path $RepoRoot "Windows\Windows_utils"
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

  $propsPath = Join-Path $WindowsUtilsDir "MKL.auto.props"

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
  function Get-MKLIncludePath([string]$root) {
    if (-not $root) { return $null }

    $candidates = @(
      (Join-Path $root "include\mkl.h"),
      (Join-Path $root "intelmkl.static.win-x64\build\native\include\mkl.h")
    )

    foreach ($candidate in $candidates) {
      if (Test-Path $candidate) {
        return $candidate
      }
    }

    return $null
  }

  # ------------------------------------------------------------
  # Check if it is just loaded
  # ------------------------------------------------------------
  Try-Load-MKLFromAutoProps -RepoRoot $RepoRoot
  if ($env:MKLROOT) {
    $inc = Get-MKLIncludePath $env:MKLROOT
    if ($inc) {
      Info "MKL detected via MKLROOT=$env:MKLROOT ($inc)"
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

    if (Get-MKLIncludePath $root) {
      $script:FoundMKL = $root
      return $true
    }

    # If this is the parent "...mkl", try to resolve a child version folder or "latest"
    if (-not (Get-MKLIncludePath $root)) {
      # Try common subfolders: latest, and the newest-looking version folder
      $latest = Join-Path $root "latest"
      if (Get-MKLIncludePath $latest) {
        $script:FoundMKL = $latest
        return $true
      }

      try {
        $ver = Get-ChildItem -Path $root -Directory -ErrorAction SilentlyContinue |
               Sort-Object Name -Descending |
               Select-Object -First 1
        if ($ver -and (Get-MKLIncludePath $ver.FullName)) {
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
  Write-Host "  - Environment variable MKLROOT, AND either"
  Write-Host "  - The file %MKLROOT%\include\mkl.h must exist, OR"
  Write-Host "  - The file %MKLROOT%\intelmkl.static.win-x64\build\native\include\mkl.h must exist"
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

function Get-OccVariant([string]$Config) {
  switch ($Config.ToLowerInvariant()) {
    "debug" { return "debug" }
    "profile" { return "debug" }
    default { return "release" }
  }
}

function Get-OccPreset([string]$Variant) {
  if ($Variant -eq "debug") {
    return "windows-msvc-debug"
  }
  return "windows-clang-cl"
}

function Get-OccInstallDir([string]$Variant) {
  return Join-Path $LibDir ("occ_" + $Variant)
}

function Get-OccRuntimeStamp([string]$Variant) {
  if ($Variant -eq "debug") {
    return "MDd"
  }
  return "MD"
}

function Get-OccBuildDir([string]$Preset) {
  return Join-Path $RepoRoot ("build-" + $Preset)
}

function Get-OccFingerprintFile([string]$Variant) {
  return Join-Path (Get-OccInstallDir $Variant) "occ_inputs.sha256"
}

function Get-Sha256Hex([string]$Path) {
  if (Get-Command Get-FileHash -ErrorAction SilentlyContinue) {
    return (Get-FileHash $Path -Algorithm SHA256).Hash.ToLowerInvariant()
  }

  $stream = [System.IO.File]::OpenRead($Path)
  $sha = [System.Security.Cryptography.SHA256]::Create()
  try {
    return ([System.BitConverter]::ToString($sha.ComputeHash($stream))).Replace("-", "").ToLowerInvariant()
  } finally {
    $sha.Dispose()
    $stream.Dispose()
  }
}

function Get-OccRepoFingerprint([string]$Variant) {
  $hashInputs = New-Object System.Collections.Generic.List[string]
  $hashInputs.Add("variant=$Variant")
  $hashInputs.Add("preset=$(Get-OccPreset $Variant)")
  $hashInputs.Add("runtime=$(Get-OccRuntimeStamp $Variant)")
  $hashInputs.Add("platform=$Platform")
  $hashInputs.Add("nos_avx=$env:NOS_AVX")

  $inputPaths = @(
    (Join-Path $RepoRoot "occ"),
    (Join-Path $RepoRoot "CMakeLists.txt"),
    (Join-Path $RepoRoot "CMakePresets.json")
  )

  foreach ($path in $inputPaths) {
    if (-not (Test-Path $path)) {
      continue
    }

    $item = Get-Item $path
    if ($item.PSIsContainer) {
      $files = Get-ChildItem $path -Recurse -File | Sort-Object FullName
    } else {
      $files = @($item)
    }

    foreach ($file in $files) {
      $relativePath = $file.FullName.Substring($RepoRoot.Length).TrimStart('\')
      $fileHash = Get-Sha256Hex $file.FullName
      $hashInputs.Add("$relativePath=$fileHash")
    }
  }

  $combined = [string]::Join("`n", $hashInputs)
  $bytes = [System.Text.Encoding]::UTF8.GetBytes($combined)
  $sha = [System.Security.Cryptography.SHA256]::Create()
  try {
    return ([System.BitConverter]::ToString($sha.ComputeHash($bytes))).Replace("-", "").ToLowerInvariant()
  } finally {
    $sha.Dispose()
  }
}

function Get-WindowsRuntimeOutputDirs {
  return @(
    (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Debug"),
    (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Profile"),
    (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Release"),
    (Join-Path $RepoRoot "Windows\NoSpherA2\x64\Release + Copy")
  )
}

function Ensure-OccVariantBuilt([string]$Variant) {
  $preset = Get-OccPreset $Variant
  $installDir = Get-OccInstallDir $Variant
  $buildDir = Get-OccBuildDir $preset
  $occOut = Join-Path $installDir "lib\occ_main.lib"
  $crtStamp = Join-Path $installDir "lib\occ_crt.stamp"
  $fingerprintFile = Get-OccFingerprintFile $Variant
  $expectedRuntime = Get-OccRuntimeStamp $Variant
  $expectedFingerprint = Get-OccRepoFingerprint $Variant
  $currentFingerprint = if (Test-Path $fingerprintFile) {
    (Get-Content $fingerprintFile -ErrorAction SilentlyContinue | Out-String).Trim()
  } else {
    ""
  }
  $tbbImportLib = Join-Path $installDir "lib\tbb12.lib"
  $tbbRuntimeDlls = @(Get-ChildItem (Join-Path $installDir "bin") -Filter "tbb12*.dll" -ErrorAction SilentlyContinue)
  $occNeedsBuild = (-not (Test-Path $occOut)) -or
                   (-not (Test-Path $crtStamp)) -or
                   ((Get-Content $crtStamp -ErrorAction SilentlyContinue) -ne $expectedRuntime) -or
                   (-not (Test-Path $fingerprintFile)) -or
                   ($currentFingerprint -ne $expectedFingerprint) -or
                   (-not (Test-Path $tbbImportLib)) -or
                   ($tbbRuntimeDlls.Count -eq 0)

  if (-not $occNeedsBuild) {
    Info "OCC $Variant already built for current inputs ($occOut)"
    return
  }

  if (Test-Path $installDir) {
    Remove-Item -Recurse -Force $installDir
  }
  if (Test-Path $buildDir) {
    Info "Removing stale OCC build cache ($buildDir)..."
    Remove-Item -Recurse -Force $buildDir
  }

  Info "Building OCC ($Variant)..."
  Push-Location $RepoRoot
  try {
    cmake --workflow --preset $preset
    cmake --install $buildDir

    $tbbLib = Get-ChildItem $buildDir -Recurse -Include "tbb12.lib","tbb12_debug.lib" `
      -ErrorAction SilentlyContinue |
      Where-Object { $_.FullName -notmatch "CMakeFiles" } |
      Select-Object -First 1 -ExpandProperty FullName
    $tbbDll = Get-ChildItem $buildDir -Recurse -Include "tbb12.dll","tbb12_debug.dll" `
      -ErrorAction SilentlyContinue |
      Where-Object { $_.FullName -notmatch "CMakeFiles" } |
      Select-Object -First 1 -ExpandProperty FullName
    if (-not $tbbLib) { Fail "tbb12.lib not found after OCC $Variant build" }
    if (-not $tbbDll) { Fail "tbb12.dll not found after OCC $Variant build" }

    New-Item -ItemType Directory -Force (Join-Path $installDir "lib") | Out-Null
    New-Item -ItemType Directory -Force (Join-Path $installDir "bin") | Out-Null
    Copy-Item $tbbLib -Destination (Join-Path $installDir "lib\tbb12.lib") -Force
    Copy-Item $tbbDll -Destination (Join-Path $installDir ("bin\" + [System.IO.Path]::GetFileName($tbbDll))) -Force

    # Keep MSBuild project files stable by providing generic import-lib names
    # that point at the debug variants when OCC was built in debug mode.
    if ($Variant -eq "debug") {
      $debugAliases = @(
        @{ Source = "fmtd.lib"; Alias = "fmt.lib" },
        @{ Source = "spdlogd.lib"; Alias = "spdlog.lib" }
      )
      foreach ($entry in $debugAliases) {
        $srcLib = Join-Path $installDir ("lib\" + $entry.Source)
        $dstLib = Join-Path $installDir ("lib\" + $entry.Alias)
        if ((Test-Path $srcLib) -and (-not (Test-Path $dstLib))) {
          Copy-Item $srcLib -Destination $dstLib -Force
        }
      }
    }

    foreach ($d in (Get-WindowsRuntimeOutputDirs)) {
      New-Item -ItemType Directory -Force $d | Out-Null
      Copy-Item $tbbDll -Destination $d -Force
    }
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $occOut)) { Fail "OCC $Variant build finished but output not found: $occOut" }
  Set-Content -Path $crtStamp -Value $expectedRuntime -Encoding ASCII
  Set-Content -Path (Join-Path $installDir "occ_variant.stamp") -Value $Variant -Encoding ASCII
  Set-Content -Path $fingerprintFile -Value $expectedFingerprint -Encoding ASCII
  Info "OCC $Variant OK"
}

function Sync-ActiveOccVariant([string]$Variant) {
  $sourceDir = Get-OccInstallDir $Variant
  $destDir = Join-Path $LibDir "occ"
  $sourceLib = Join-Path $sourceDir "lib\occ_main.lib"

  if (-not (Test-Path $sourceLib)) {
    Fail "Cannot activate OCC ${Variant}: missing $sourceLib"
  }

  if (Test-Path $destDir) {
    Remove-Item -Recurse -Force $destDir
  }

  Copy-Item -Recurse -Force $sourceDir $destDir
  Set-Content -Path (Join-Path $destDir "occ_variant.stamp") -Value $Variant -Encoding ASCII
  Info "Activated OCC $Variant at $destDir"
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
$autoProps = Join-Path $WindowsUtilsDir "MKL.auto.props"

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
    cmake -DBUILD_SHARED_LIBS=0 -DCMAKE_BUILD_TYPE=RELEASE -DPYPZPX="ON" -DCMAKE_INSTALL_PREFIX="..\..\Lib\LibCint" -DCMAKE_C_COMPILER=cl -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreadedDLL" ..
    cmake --build . --config RELEASE
    cmake --install .
  } finally {
    Pop-Location
  }

  if (-not (Test-Path $LibCintOut)) { Fail "LibCint build finished but output not found: $LibCintOut" }
  Info "LibCint OK"
}
# -----------------------------
# 4) Build OCC variants and activate the requested one
# -----------------------------
Ensure-OccVariantBuilt "release"
Ensure-OccVariantBuilt "debug"
Sync-ActiveOccVariant (Get-OccVariant $Configuration)

# -----------------------------
# 4a) Deploy MKL runtime DLLs (needed when OCC is built with /MD)
# -----------------------------
# With /MD, the MKL threading/runtime layer can pull in additional MKL runtime
# DLLs at load time. Deploy the known transitive dependencies alongside the
# other runtime DLLs so the test runner can find them.
# $mklBinDir = Join-Path $env:MKLROOT "bin"
# $mklDllsNeeded = @("mkl_intel_thread.2.dll", "mkl_core.2.dll", "mkl_def.2.dll")
# $outDirsMkl = Get-WindowsRuntimeOutputDirs
# foreach ($mklDll in $mklDllsNeeded) {
#   $src = Join-Path $mklBinDir $mklDll
#   if (-not (Test-Path $src)) {
#     Write-Warning "MKL DLL not found: $src - skipping deployment (test runner may fail to load NoSpherA2_DLL.dll)"
#     continue
#   }
#   foreach ($d in $outDirsMkl) {
#     New-Item -ItemType Directory -Force $d | Out-Null
#     Copy-Item $src -Destination $d -Force
#   }
#   Info "Deployed $mklDll"
# }

$stamp = Join-Path $WindowsUtilsDir "Build_Dependencies\deps.stamp"
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
