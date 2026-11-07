#!/bin/bash
# Stage 2 for every accepted molecule, in one job.
#
#   sbatch nbo_all_stage.sh [mol ...]        # default: all 22
#
# One job rather than 22. Stage 2 measured 1 s per molecule on a Debug node, so 22 separate
# jobs spend three orders of magnitude more time waiting in the queue than computing - and a
# 22-job wave is exactly the situation where the "cancel before resubmitting the same name"
# trap bites. The loop calls nbo_stage.sh unchanged (its own #SBATCH lines are comments to a
# plain bash invocation), so there is one implementation of the four steps, not two.
#
# Threads come from this job's own allocation and each molecule's provenance_nbo.json records
# them along with the hostname, so a per-molecule timing stays reusable.
#SBATCH --job-name=nboref-all
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/all-%j.out

set -u
HERE=${NBOREF_BIN:-/work/akkleemiss/florian/nbo_ref_v2/bin}
source "$HERE/env.sh"

MOLS="$*"
if [ -z "$MOLS" ]; then
    MOLS="acetylene ammonia benzene ch3 ethane ethene formaldehyde formate hcn lif n2 ni_co_4 nitromethane no o2 ozone pf5 pyridine sf6 so2 ticl4 water"
fi

echo "=== stage 2 for all on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has ${SLURM_CPUS_PER_TASK:-1} threads"
FAILED=""
for mol in $MOLS; do
    T0=$(date +%s)
    bash "$HERE/nbo_stage.sh" "$mol"
    RC=$?
    echo "=== $mol rc=$RC wall=$(( $(date +%s) - T0 ))s"
    [ $RC -eq 0 ] || FAILED="$FAILED $mol"
done
echo "=== done on $(hostname), threads ${SLURM_CPUS_PER_TASK:-1}, failed:${FAILED:- none}"
[ -z "$FAILED" ]
