# Prints the CudaCompile CodeGeneration string for the GPUs in this machine,
# e.g. 'compute_75,sm_75', one entry per distinct compute capability, joined
# by ';'. Prints nothing when nvidia-smi is missing or fails, which makes the
# build fall back to the portable architecture list.
try {
    $caps = & nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>$null
    if ($LASTEXITCODE -ne 0 -or -not $caps) { exit 0 }
    $out = @()
    foreach ($c in $caps) {
        $v = ($c -replace '[^0-9.]', '')
        if ($v -match '^\d+\.\d+$') {
            $sm = $v.Replace('.', '')
            $out += "compute_$sm,sm_$sm"
        }
    }
    $out = $out | Select-Object -Unique
    if ($out) { Write-Output ($out -join ';') }
} catch { }
