# Windows Build Dependencies

This directory contains the build script for Windows dependencies, including Intel oneMKL.

## Security: MKL Installer Verification

To mitigate supply-chain security risks, the build script supports SHA256 checksum verification of the Intel oneMKL installer.

### Obtaining the Official Checksum

1. Visit the official Intel oneMKL download page:
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

2. Locate the SHA256 checksum for version `2025.3.1.10` Windows offline installer:
   - Look in the download details section or release notes
   - The checksum is typically displayed near the download link or in a separate verification file
   - Intel may provide a `.sha256` file or list checksums in the product documentation

3. Alternatively, if you already have a trusted copy of the installer from a previous verified download, compute its SHA256:
   ```powershell
   Get-FileHash -Path "intel-onemkl-2025.3.1.10-windows.exe" -Algorithm SHA256
   ```
   Note: This method only works if you already have a verified copy from a trusted source.

### Using Checksum Verification

When building dependencies with automatic MKL installation, provide the checksum:

```powershell
.\build_deps.ps1 -RepoRoot "C:\path\to\repo" -InstallMKLIfMissing -MKLInstallerSHA256 "YOUR_SHA256_CHECKSUM_HERE"
```

### Behavior

- **With checksum**: The script will verify the downloaded installer against the provided SHA256. If verification fails, the build aborts with an error.
  
- **Without checksum**: The script will display a security warning and continue after a 3-second delay. This is not recommended for production or CI environments.

### Example

```powershell
# Secure build with checksum verification
.\build_deps.ps1 -RepoRoot $PSScriptRoot\..\..\.. -InstallMKLIfMissing -MKLInstallerSHA256 "abc123def456..."
```

## Script Usage

```powershell
.\build_deps.ps1 -RepoRoot <path> [-Configuration <Release|Debug>] [-Platform <x64|Win32>] [-InstallMKLIfMissing] [-MKLInstallerSHA256 <checksum>]
```

### Parameters

- **RepoRoot** (required): Path to the repository root
- **Configuration**: Build configuration (default: Release)
- **Platform**: Target platform (default: x64)
- **InstallMKLIfMissing**: Automatically download and install Intel oneMKL if not found
- **MKLInstallerSHA256**: Expected SHA256 checksum for installer verification (recommended when using -InstallMKLIfMissing)
