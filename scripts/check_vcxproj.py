"""The Visual Studio projects list their sources by hand while CMake globs them; a test file or a core
source added on the CMake side silently stays out of the .sln build. Exit 1 when the two disagree.
Run from anywhere: python scripts/check_vcxproj.py"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# (vcxproj, source directory globbed by CMake GLOB_RECURSE, files CMake filters out again)
PAIRS = [
    # the two nos_math halves are #included by nos_math.cpp, the GPU host files are per-backend
    ("Windows/NoSpherA2_LIB/NoSpherA2_LIB.vcxproj", "Src/core",
     {"vec_nos_math.cpp", "mat_nos_math.cpp", "hip_runtime_shim.cpp", "gpu_dispatch.cpp"}),
    ("Windows/Tests/Tests.vcxproj", "tests/src", set()),
]
EXTENSIONS = (".cpp", ".cc", ".cxx", ".cu")  # .cu files sit in the CudaCompile (CUDA) and CustomBuild (HIP) groups


def sources_on_disk(directory, extensions):
    found = set()
    for dirpath, _, files in os.walk(os.path.join(ROOT, directory)):
        for f in files:
            if f.endswith(extensions):
                rel = os.path.relpath(os.path.join(dirpath, f), os.path.join(ROOT, directory))
                found.add(rel.replace("\\", "/"))
    return found


def sources_in_vcxproj(vcxproj, directory):
    with open(os.path.join(ROOT, vcxproj), encoding="utf-8") as fh:
        text = fh.read()
    project_dir = os.path.dirname(os.path.join(ROOT, vcxproj))
    listed = set()
    for inc in re.findall(r'<(?:ClCompile|CudaCompile|CustomBuild) Include="([^"]+)"', text):
        full = os.path.normpath(os.path.join(project_dir, inc.replace("\\", "/")))
        rel = os.path.relpath(full, os.path.join(ROOT, directory)).replace("\\", "/")
        if not rel.startswith(".."):
            listed.add(rel)
    return listed


def main():
    bad = False
    for vcxproj, directory, excluded in PAIRS:
        disk, listed = sources_on_disk(directory, EXTENSIONS) - excluded, sources_in_vcxproj(vcxproj, directory) - excluded
        for label, files in (("missing from", disk - listed), ("stale in", listed - disk)):
            if files:
                bad = True
                print(f"{label} {vcxproj}: {', '.join(sorted(files))}")
    print("vcxproj sources match CMake" if not bad else "vcxproj sources differ from CMake")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
