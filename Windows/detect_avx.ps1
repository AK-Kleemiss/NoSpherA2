# Prints 'ON' if the build host CPU + OS support AVX, 'OFF' otherwise.
# Used by Directory.Build.targets to pick EnableEnhancedInstructionSet.
Add-Type -Namespace Nos -Name Cpu -MemberDefinition '[DllImport("kernel32.dll")] public static extern bool IsProcessorFeaturePresent(int feature);'
# 39 = PF_AVX_INSTRUCTIONS_AVAILABLE (returns false on Windows versions that
# do not know the constant, which safely falls back to SSE-only builds)
if ([Nos.Cpu]::IsProcessorFeaturePresent(39)) { Write-Output 'ON' } else { Write-Output 'OFF' }
