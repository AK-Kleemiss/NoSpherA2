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
    [string]$MSBuildExe
)

$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$auxFile = Join-Path $repoRoot "Src\auxiliary_basis.cpp"
$basisDir = Join-Path $repoRoot "Src\basis_set_helper\basis_sets"

$csvFiles = Get-ChildItem -Path $basisDir -Filter "*.csv" -File
$needsRegeneration = -not (Test-Path $auxFile)

if (-not $needsRegeneration) {
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

if (-not (Test-Path $ConverterExe)) {
    $converterProject = Join-Path $PSScriptRoot "BasisSetConverter\BasisSetConverter.vcxproj"
    if (-not [string]::IsNullOrWhiteSpace($MSBuildExe)) {
        if (-not (Test-Path $MSBuildExe)) {
            throw "Provided MSBuild executable was not found: $MSBuildExe"
        }
        $msbuildCmd = $MSBuildExe
    }
    else {
        $resolvedMsbuild = Get-Command msbuild.exe -ErrorAction SilentlyContinue
        if ($null -eq $resolvedMsbuild) {
            throw "Could not resolve msbuild.exe. Pass -MSBuildExe explicitly from MSBuild/CI."
        }
        $msbuildCmd = $resolvedMsbuild.Source
    }

    Write-Host "BasisSetConverter executable not found. Building $converterProject ($Configuration|$Platform)..."
    & $msbuildCmd $converterProject /m /nologo "/p:Configuration=$Configuration" "/p:Platform=$Platform"
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to build BasisSetConverter via msbuild."
    }
}

Write-Host "Generating auxiliary_basis.cpp via BasisSetConverter..."
& $ConverterExe $SolutionDir
if ($LASTEXITCODE -ne 0) {
    throw "BasisSetConverter execution failed with exit code $LASTEXITCODE"
}
