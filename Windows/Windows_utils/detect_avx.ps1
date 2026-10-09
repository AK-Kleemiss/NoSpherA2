# Prints the widest instruction set the build host CPU + OS support: AVX2, AVX or SSE4.2.
# Used by NoSpherA2_detect_AVX.props when NOS_AVX is not set.

# Add-Type compiles C# with csc, which fails on the LIB/INCLUDE paths set by
# VS developer environments. Clear them for this process only.
$env:LIB = ''
$env:INCLUDE = ''

Add-Type -Namespace Nos -Name Cpu -MemberDefinition '[DllImport("kernel32.dll")] public static extern bool IsProcessorFeaturePresent(int feature);'
# 39 = PF_AVX_INSTRUCTIONS_AVAILABLE, 40 = PF_AVX2_INSTRUCTIONS_AVAILABLE (false on Windows versions
# that do not know the constant, which safely falls back to the narrower build). There is no PF_ for
# FMA, which cmake/detect_avx.c checks as well: every AVX2 CPU shipped so far also has FMA3.
if (-not [Nos.Cpu]::IsProcessorFeaturePresent(39)) { Write-Output 'SSE4.2' }
elseif ([Nos.Cpu]::IsProcessorFeaturePresent(40)) { Write-Output 'AVX2' }
else { Write-Output 'AVX' }
