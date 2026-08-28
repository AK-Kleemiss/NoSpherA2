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
SRC = os.path.join(ROOT, "Src", "core", "sf_gpu.cu")
HDR = os.path.join(ROOT, "Src", "core", "gpu_backend.h")


def main():
    src = io.open(SRC, encoding="utf-8").read()
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
