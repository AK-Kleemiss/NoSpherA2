#!/usr/bin/env bash
# Line coverage of the gtest suite with gcc/gcov: configure with NOSPHERA2_COVERAGE (only the
# NoSpherA2 targets are instrumented, the dependencies stay plain), run ctest, then gcovr.
# usage: scripts/coverage.sh [build dir]   (gcovr: pip install gcovr, or uv tool install gcovr)
# Result: <build dir>/coverage/index.html and a per-file summary on stdout.
set -euo pipefail
src=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-$src/build/coverage-linux}
cmake --preset release-linux -S "$src" -B "$build" -DNOSPHERA2_COVERAGE=ON -DNOSPHERA2_BUILD_TESTS=ON -DNOSPHERA2_GPU_AUTO=OFF
cmake --build "$build" -j "${JOBS:-12}"
find "$build" -name '*.gcda' -delete
(cd "$build" && OCC_DATA_PATH="$src/occ/share" ctest -j "${CTEST_PARALLEL_LEVEL:-4}" --output-on-failure) || true
mkdir -p "$build/coverage"
gcovr -r "$src" --object-directory "$build" --gcov-ignore-parse-errors suspicious_hits.warn_once_per_file \
    --filter "$src/Src/" --filter "$src/app/" \
    --exclude '.*_gpu\.(cu|h)$' \
    --html-details "$build/coverage/index.html" --print-summary --sort uncovered-number
gcovr -r "$src" --object-directory "$build" --gcov-ignore-parse-errors suspicious_hits.warn_once_per_file --filter "$src/Src/" --filter "$src/app/" --exclude '.*_gpu\.(cu|h)$' | tail -n +4 | sort -k4 -n -r | head -40
