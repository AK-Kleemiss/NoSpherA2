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
$auxFile = Join-Path $repoRoot "Src\auxiliary_basis.cpp"
$basisDir = Join-Path $repoRoot "Src\basis_set_helper\basis_sets"

$hasAuxFile = Test-Path $auxFile
$needsRegeneration = -not $hasAuxFile

if (-not (Test-Path $basisDir)) {
    if ($hasAuxFile) {
        Write-Host "Basis set source directory not found ($basisDir), but auxiliary_basis.cpp exists. Skipping BasisSetConverter."
        exit 0
    }
    throw "Basis set source directory not found: $basisDir and auxiliary_basis.cpp is missing at $auxFile"
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

    if ([string]::IsNullOrWhiteSpace($PlatformToolset) -and -not [string]::IsNullOrWhiteSpace($env:PlatformToolset)) {
        $PlatformToolset = $env:PlatformToolset
    }

    Write-Host "BasisSetConverter executable not found. Building $converterProject ($Configuration|$Platform)..."
    $msbuildArgs = @(
        $converterProject,
        '/m',
        '/nologo',
        "/p:Configuration=$Configuration",
        "/p:Platform=$Platform"
    )
    if (-not [string]::IsNullOrWhiteSpace($PlatformToolset)) {
        $msbuildArgs += "/p:PlatformToolset=$PlatformToolset"
        Write-Host "Forwarding PlatformToolset=$PlatformToolset to BasisSetConverter build."
    }

    & $msbuildCmd @msbuildArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to build BasisSetConverter via msbuild."
    }

    # The vcxproj can emit the exe under a project-specific OutDir. Resolve a usable path.
    $candidateConverterPaths = @(
        $ConverterExe,
        (Join-Path $PSScriptRoot "BasisSetConverter\$Platform\BasisSetConverter\BasisSetConverter.exe"),
        (Join-Path $PSScriptRoot "BasisSetConverter\$Platform\$Configuration\BasisSetConverter.exe")
    )
    $resolvedConverter = $null
    foreach ($candidate in $candidateConverterPaths) {
        if (Test-Path $candidate) {
            $resolvedConverter = $candidate
            break
        }
    }
    if ($null -eq $resolvedConverter) {
        throw "BasisSetConverter was built but executable could not be found. Checked: $($candidateConverterPaths -join '; ')"
    }
    $ConverterExe = $resolvedConverter
}

Write-Host "Generating auxiliary_basis.cpp via BasisSetConverter..."
Write-Host "Using BasisSetConverter executable: $ConverterExe"
& $ConverterExe $SolutionDir
if ($LASTEXITCODE -ne 0) {
    throw "BasisSetConverter execution failed with exit code $LASTEXITCODE"
}
