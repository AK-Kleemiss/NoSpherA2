# CPU profile of one NoSpherA2 run from the command line, no IDE session needed.
#
#   powershell -File scripts\profile_cpu.ps1 -AppArgs "<NoSpherA2 arguments>" [-Exe <path>] [-WorkingDir <dir>]
#       [-Tool auto|vtune|vs] [-Output <dir or file.diagsession>] [-Top 30]
#
# Build the Profile|x64 configuration of Windows\NoSpherA2\NoSpherA2.sln first: release code generation with
# full symbols, so the samples resolve to the functions and lines of the binary users actually run.
#
# Tool vtune (default when Intel VTune is installed): user-mode sampling, prints the top functions by self time
# and leaves the result directory for "vtune -report hotspots -r <dir> -group-by source-line" and friends.
# Tool vs: the Visual Studio Standard Collector writes a .diagsession to open in Visual Studio. Its /loadConfig
# is broken on VS 18 (Newtonsoft.Json binding) and /loadAgent rejects agent options, so it is the 1 kHz CPU agent only.
# In both, func@0x... inside libiomp5md.dll is the OpenMP spin-wait of idle threads, not NoSpherA2 work.
param(
    [Parameter(Mandatory = $true)]
    [string]$AppArgs,
    [string]$Exe = "",
    [string]$WorkingDir = (Get-Location).Path,
    [ValidateSet("auto", "vtune", "vs")]
    [string]$Tool = "auto",
    [string]$Output = "",
    [int]$Top = 30
)
$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
if ($Exe -eq "") { $Exe = Join-Path $root "build\Profile_x64\NoSpherA2.exe" }
$Exe = (Resolve-Path $Exe).Path
$WorkingDir = (Resolve-Path $WorkingDir).Path
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$vtune = "${env:ProgramFiles(x86)}\Intel\oneAPI\vtune\latest\bin64\vtune.exe"
if ($Tool -eq "auto") { $Tool = if (Test-Path $vtune) { "vtune" } else { "vs" } }
$argv = [System.Management.Automation.PSParser]::Tokenize($AppArgs, [ref]$null) | Where-Object { $_.Type -ne "NewLine" } | ForEach-Object { $_.Content }
$sw = [System.Diagnostics.Stopwatch]::StartNew()
if ($Tool -eq "vtune") {
    if (-not (Test-Path $vtune)) { throw "vtune.exe not found at $vtune" }
    if ($Output -eq "") { $Output = Join-Path $WorkingDir "vtune_$stamp" }
    $Output = [System.IO.Path]::GetFullPath($Output)
    & $vtune -collect hotspots -knob sampling-mode=sw -result-dir $Output -app-working-dir $WorkingDir -quiet -- $Exe @argv
    if ($LASTEXITCODE -ne 0) { throw "vtune collect failed ($LASTEXITCODE)" }
    $sw.Stop()
    cmd /c "`"$vtune`" -report hotspots -r `"$Output`" -group-by function -column `"CPU Time:Self`" -format text -limit $Top -quiet 2>&1" | Where-Object { $_ -notmatch "^vtune: |Column filter" }
    cmd /c "`"$vtune`" -report hotspots -r `"$Output`" -group-by function -format csv -csv-delimiter comma -quiet -report-output `"$Output\hotspots.csv`" 2>&1" | Out-Null
    Write-Host ("Collect and finalize {0:N1} s, result in {1} (hotspots.csv inside; vtune-gui {1} for the source view)" -f $sw.Elapsed.TotalSeconds, $Output)
    exit 0
}
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$vsroot = & $vswhere -latest -products * -property installationPath
$collector = Join-Path $vsroot "Team Tools\DiagnosticsHub\Collector"
$vsdiag = Join-Path $collector "VSDiagnostics.exe"
if (-not (Test-Path $vsdiag)) { throw "VSDiagnostics.exe not found under $collector" }
$agent = "/loadAgent:4EA90761-2248-496C-B854-3C0399A591A4;DiagnosticsHub.CpuAgent.dll"
if ($Output -eq "") { $Output = Join-Path $WorkingDir "NoSpherA2_$stamp.diagsession" }
$Output = [System.IO.Path]::GetFullPath($Output)
$session = [guid]::NewGuid().ToString()
$exeName = [System.IO.Path]::GetFileNameWithoutExtension($Exe)
$before = @(Get-Process -Name $exeName -ErrorAction SilentlyContinue | ForEach-Object { $_.Id })
Push-Location $WorkingDir
try {
    $out = cmd /c "`"$vsdiag`" start $session `"/launch:$Exe`" `"/launchArgs:$AppArgs`" $agent 2>&1" | Out-String
    if ($LASTEXITCODE -ne 0 -or $out -notmatch "Running") { throw "VSDiagnostics start failed:`n$out" }
    Start-Sleep -Milliseconds 500
    do {
        Start-Sleep -Milliseconds 500
        $alive = @(Get-Process -Name $exeName -ErrorAction SilentlyContinue | Where-Object { $before -notcontains $_.Id })
    } while ($alive.Count -gt 0)
    $sw.Stop()
    $out = cmd /c "`"$vsdiag`" stop $session `"/output:$Output`" 2>&1" | Out-String
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $Output)) { throw "VSDiagnostics stop failed:`n$out" }
}
finally { Pop-Location }
Write-Host ("Run time {0:N1} s, profile written to {1} (open in Visual Studio)" -f $sw.Elapsed.TotalSeconds, $Output)
