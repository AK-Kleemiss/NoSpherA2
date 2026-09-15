# Stress the single-precision I tensor against the double CPU path.
#
# fp32 GEMM error grows with the accumulation length, which here is the number of grid
# points per block, so -acc is the direct lever: acc=2 is 102,920 points, acc=3 is
# 435,008. The gaussian-halt run is the long trajectory - more fitting steps means more
# chances for a small integral error to compound.
#
# Runs are sequential on purpose: they share tests/P1_test as a working directory and
# would overwrite each other's NoSpherA2.log and I_tensor_stream.bin.

$exe = "D:\Git\NoSpherA2\build\release-windows\bin\NoSpherA2.exe"
$dir = "D:\Git\NoSpherA2\tests\P1_test"
$out = "D:\Git\NoSpherA2\stress"
New-Item -ItemType Directory -Force -Path $out | Out-Null
Set-Location $dir

$cases = @(
    @{ name = "acc2";      args = @("-acc", "2") },
    @{ name = "acc3";      args = @("-acc", "3") },
    @{ name = "halt_acc2"; args = @("-acc", "2", "-xcw_gaussian_halt") }
)

foreach ($c in $cases) {
    foreach ($mode in @("cpu", "gpu")) {
        $a = @("-do_XCW", "-cif", "P1_test_NA2.cif", "-hkl", "P1_test.hkl",
               "-anom_disp", "anom_disp.txt", "-XCW_settings", "test_settings.txt") + $c.args
        if ($mode -eq "gpu") { $a += "-gpu_itensor" }
        $sw = [Diagnostics.Stopwatch]::StartNew()
        & $exe @a 2>&1 | Out-Null
        $sw.Stop()
        $log = Join-Path $out ("{0}_{1}.log" -f $c.name, $mode)
        Copy-Item "NoSpherA2.log" $log -Force
        $t = (Get-Content $log | Select-String "XCW integrals").Line
        Write-Output ("{0,-10} {1,-3} wall {2,7:N1} s   {3}" -f $c.name, $mode, $sw.Elapsed.TotalSeconds, $t)
    }
}
Write-Output "STRESS_DONE"
