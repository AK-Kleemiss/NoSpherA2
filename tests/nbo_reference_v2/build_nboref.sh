#!/bin/bash
# Build MY OWN NoSpherA2 at nbo_external_reference, so -nbo_parse carries the valency-spin fix.
# Nothing here touches the cluster share or any sibling session's tree.
#
# DO NOT sbatch THIS ON AKL: cmake and git live on the login node AKL007 only and no module
# provides them, so a compute node dies rc=127 at the configure step (job 594630 on AKL008,
# "cmake: command not found").  Run the build on the login node instead - an incremental
# build of this tree took 22.7 s there - and keep this file for the env and path record.
#SBATCH --job-name=nboref-build
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/build-%j.out

set -u
TREE=/work/akkleemiss/florian/nos_nboref
cd "$TREE" || exit 1
J=${SLURM_CPUS_PER_TASK:-8}
echo "=== build on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $J"
echo "=== $(git log --oneline -1)"
T0=$(date +%s)
cmake --preset release-linux > cmake_configure.log 2>&1
RC1=$?
echo "=== configure rc=$RC1 after $(( $(date +%s) - T0 ))s"
[ $RC1 -eq 0 ] || { tail -30 cmake_configure.log; exit 1; }
cmake --build --preset release-linux -j "$J" > cmake_build.log 2>&1
RC2=$?
echo "=== build rc=$RC2 after $(( $(date +%s) - T0 ))s total"
tail -5 cmake_build.log
ls -la build/release-linux/bin/NoSpherA2 2>/dev/null
exit $RC2
