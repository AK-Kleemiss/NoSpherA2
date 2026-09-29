#!/bin/bash
# Both stages of ONE held-out molecule, in one small job.
#
#   sbatch holdout_stage.sh <molecule>              # baseline binary, holdout root
#   NOSPHERA2=/path/to/other sbatch holdout_stage.sh <molecule>
#
# The 22 ran stage 1 and stage 2 as two jobs because stage 1 needed 2 h and stage 2 needed 8
# threads.  Every held-out molecule has 2-5 atoms, so both fit inside one 40-minute job, and
# the partition is saturated by other people's work - one small job per molecule backfills
# where a pair of larger ones would queue.  The RECIPE is unchanged: this calls orca_opt.sh
# and nbo_stage.sh themselves, so there is no second copy of either command anywhere.
#
# The two deliberate differences from the 22, both recorded in every provenance JSON rather
# than smoothed over:
#   - 4 MPI ranks for ORCA means SLURM_CPUS_PER_TASK is 1, so -nbo_threads is 1 and not 8.
#     Thread count must not move a number; if it does, that is a finding, not a nuisance.
#   - NOSPHERA2 points at a NAMED baseline binary instead of the share, so the "which binary"
#     column of the report is answered by the file the job actually ran.
#SBATCH --job-name=holdout
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2_holdout/slurm/holdout-%j.out

set -u
MOL=$1
export ROOT=${ROOT:-/work/akkleemiss/florian/nbo_ref_v2_holdout}
export NBOREF_BIN=${NBOREF_BIN:-$ROOT/bin}
export NOSPHERA2=${NOSPHERA2:-$ROOT/bin/NoSpherA2_baseline}

echo "=== holdout $MOL on $(hostname), job ${SLURM_JOB_ID:-none}"
echo "=== root      $ROOT"
echo "=== binary    $NOSPHERA2  ($(cat "$(dirname "$NOSPHERA2")/DEPLOYED_COMMIT.txt" 2>/dev/null))"

bash "$NBOREF_BIN/orca_opt.sh" "$MOL"
RC1=$?
if [ $RC1 -ne 0 ]; then
    echo "=== stage 1 rc=$RC1, stage 2 not attempted"
    exit $RC1
fi
grep -q "THE OPTIMIZATION HAS CONVERGED" "$ROOT/$MOL/$MOL.out" || {
    echo "=== stage 1 finished but the optimisation did not converge, stage 2 not attempted"
    exit 4
}
bash "$NBOREF_BIN/nbo_stage.sh" "$MOL"
RC2=$?
echo "=== holdout $MOL stage1=$RC1 stage2=$RC2"
exit $RC2
