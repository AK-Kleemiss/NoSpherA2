#!/usr/bin/env python3
"""Check that sf_gpu.cu only uses runtime names gpu_backend.h maps for both backends.

Without an AMD card in CI the HIP build is not exercised, so a missing mapping would
only surface on someone else's machine. This catches it from the source alone.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRCS = [os.path.join(ROOT, "Src", "core", n)
        for n in ("sf_gpu.cu", "itensor_gpu.cu", "salted_gpu.cu")]
HDR = os.path.join(ROOT, "Src", "core", "gpu_backend.h")


def main():
    src = "".join(io.open(f, encoding="utf-8").read() for f in SRCS if os.path.exists(f))
    hdr = io.open(HDR, encoding="utf-8").read()

    used = {n for n in re.findall(r"\bgpu[A-Za-z_][A-Za-z0-9_]*", src) if n != "gpu_backend"}
    hip = hdr.split("#ifdef NOSPHERA2_USE_HIP")[1].split("#else")[0]
    cuda = hdr.split("#else")[1].split("#endif")[0]
    defined = lambda t: set(re.findall(r"#define\s+(gpu[A-Za-z_][A-Za-z0-9_]*)", t))
    d_hip, d_cuda = defined(hip), defined(cuda)

    bad = []
    for name in sorted(used - d_hip):
        bad.append("%s used but not mapped in the HIP branch" % name)
    for name in sorted(used - d_cuda):
        bad.append("%s used but not mapped in the CUDA branch" % name)
    for name in sorted(d_hip ^ d_cuda):
        bad.append("%s is mapped in only one branch" % name)

    # The fp32:fp64 ratio is the one thing HIP has no equivalent for, so those three
    # names are expected. Any other cuda* symbol compiles under nvcc and breaks on AMD.
    allowed = ("cudaDeviceGetAttribute", "cudaDevAttrSingleToDoublePrecisionPerfRatio", "cudaSuccess")
    stray = [m for m in re.findall(r"\bcuda[A-Z][A-Za-z0-9_]*", src) if m not in allowed]
    for name in sorted(set(stray)):
        bad.append("%s is CUDA-only and is not inside the CUDA branch of sf_gpu_fp64_ratio" % name)

    if bad:
        for b in bad:
            sys.stderr.write("gpu backend check: %s\n" % b)
        return 1
    sys.stdout.write("gpu backend check: %d runtime names map in both backends\n" % len(used))

    if "--compile" in sys.argv:
        return compile_hip_branch()
    return 0


def compile_hip_branch():
    """Compile the HIP branch with nvcc against the shim in tests/hip_compile_check.

    HIP on an NVIDIA host is itself a header translation onto CUDA, so this checks the
    branch's names resolve and its device code is valid. It does not check AMD code
    generation, the AMD device math library, or that anything runs.
    """
    import glob
    import shutil
    import subprocess
    import tempfile

    # Take the compiler CMake configured, not the one on PATH. On this machine both PATH
    # and CUDA_PATH point at a toolkit too old for the host compiler, so guessing picks
    # a toolkit that cannot build the project at all.
    nvcc = None
    for cache in glob.glob(os.path.join(ROOT, "build", "*", "CMakeCache.txt")):
        for line in io.open(cache, encoding="utf-8", errors="replace"):
            if line.startswith("CMAKE_CUDA_COMPILER:"):
                candidate = line.split("=", 1)[1].strip()
                if os.path.exists(candidate):
                    nvcc = candidate
                break
        if nvcc:
            break
    if not nvcc:
        nvcc = shutil.which("nvcc")
    if not nvcc or not os.path.exists(nvcc):
        sys.stdout.write("gpu backend check: no nvcc configured, skipping the HIP compile\n")
        return 0
    shim = os.path.join(ROOT, "tests", "hip_compile_check")
    for src_path in SRCS:
      if not os.path.exists(src_path):
        continue
      with tempfile.TemporaryDirectory() as tmp:
        cmd = [nvcc, "-c", "-O3", "-DNOSPHERA2_USE_HIP",
               "-I", shim, "-I", os.path.join(ROOT, "Src", "core"),
               "-o", os.path.join(tmp, "hip_probe.obj"), src_path]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            sys.stderr.write("gpu backend check: %s fails under HIP\n%s\n"
                             % (os.path.basename(src_path), p.stderr))
            return 1
    sys.stdout.write("gpu backend check: %d kernels compile under HIP (shim, not real hipcc)\n"
                     % len(SRCS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
