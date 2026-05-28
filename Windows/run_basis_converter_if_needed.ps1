param(
    [Parameter(Mandatory = $true)]
    [string]$ConverterExe,
    [Parameter(Mandatory = $true)]
    [string]$SolutionDir,
    [Parameter(Mandatory = $true)]
    [string]$Configuration,
    [Parameter(Mandatory = $true)]
    [string]$Platform,
    [Parameter(Mandatory = $false)]
    [string]$PlatformToolset,
    [Parameter(Mandatory = $false)]
    [string]$MSBuildExe
)

$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$auxFile = Join-Path $repoRoot "Src\basis_data.cpp"
$basisDir = Join-Path $repoRoot "BasisSetGenerator\basis_sets"

$hasAuxFile = Test-Path $auxFile
$needsRegeneration = -not $hasAuxFile

if (-not (Test-Path $basisDir)) {
    if ($hasAuxFile) {
        Write-Host "Basis set source directory not found ($basisDir), but basis_data.cpp exists. Skipping BasisSetConverter."
        exit 0
    }
    throw "Basis set source directory not found: $basisDir and basis_data.cpp is missing at $auxFile"
}

if (-not $needsRegeneration) {
    $csvFiles = Get-ChildItem -Path $basisDir -Filter "*.csv" -File
    $auxTimestamp = (Get-Item $auxFile).LastWriteTimeUtc
    foreach ($csv in $csvFiles) {
        if ($csv.LastWriteTimeUtc -gt $auxTimestamp) {
            $needsRegeneration = $true
            break
        }
    }
}

if (-not $needsRegeneration) {
    Write-Host "Basis set data is up to date; skipping BasisSetConverter."
    exit 0
}

Write-Host "Generating basis_data.cpp via BasisSetConverter..."
Write-Host "Using BasisSetConverter executable: $ConverterExe"

$srcDir = Join-Path $repoRoot "Src\"

& $ConverterExe $basisDir $srcDir
if ($LASTEXITCODE -ne 0) {
    throw "BasisSetConverter execution failed with exit code $LASTEXITCODE"
}
