#!/bin/bash
# Task 5: the gtest suite on the patched tree, TWICE on the SAME binary - default (renat5) and
# NAO_LEGACY_CASCADE=1 (the pre-change cascade).  Two runs of one binary separate "the change broke
# a test" from "this test was already failing", without needing a second build.
#
# The suite locates its data by walking up for tests/tests.toml, so it must start in
# <repo>/tests/src or be told NOS_REPO_ROOT - run it anywhere else and a dozen tests fail on a
# missing file, which is the harness misconfigured and not a defect.  First attempt (594821) did
# exactly that: 4 failures, all std::filesystem::exists false, from the wrong cwd.
#SBATCH --job-name=renat5-tests
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=/work/akkleemiss/florian/nao_renat5/slurm-tests-%j.out

set -u
TREE=/work/akkleemiss/florian/nao_renat5/tree
T=$TREE/build/release-linux
MY=/work/akkleemiss/florian/nao_renat5
export NOS_REPO_ROOT=$TREE
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export LD_LIBRARY_PATH=$T/bin:${LD_LIBRARY_PATH:-}

echo "=== $(hostname), $(grep -c ^processor /proc/cpuinfo) cores, $OMP_NUM_THREADS threads"
md5sum "$T/bin/NoSpherA2_Tests"
ls -la "$TREE/tests/tests.toml"

for arm in renat5 legacy; do
    w=$MY/tests_$arm
    rm -rf "$w"; mkdir -p "$w"
    if [ "$arm" = legacy ]; then E=NAO_LEGACY_CASCADE=1; else E=NAO_LEGACY_CASCADE=; fi
    t0=$(date +%s)
    ( cd "$TREE/tests/src" && env "$E" "$T/bin/NoSpherA2_Tests" --gtest_output=xml:"$w/gtest.xml" \
        > "$w/gtest.log" 2>&1 )
    rc=$?
    echo "=== arm $arm rc=$rc  $(( $(date +%s) - t0 ))s"
    grep -E "^\[==========\] [0-9]+ tests? from|^\[  PASSED  \]|^\[  FAILED  \] [0-9]+ test" \
        "$w/gtest.log" | tail -4
    echo "--- distinct failing tests, arm $arm:"
    grep -E "^\[  FAILED  \] [A-Za-z][A-Za-z0-9_]*\." "$w/gtest.log" | sed 's/ (.*//' | sort -u
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
