#!/bin/bash
# The one surviving original archive, analysed by both programs.
#
#   sbatch ch3_archive_stage.sh
#
# tests/nbo_reference/ch3_orca_reference.47 is the only .47 that survived the first
# collection, and it carries the original density, Fock and overlap matrices - the actual
# numbers the stored ch3.nbo.json was computed from. So this run needs no regeneration at
# all: gennbo 7 reads that archive, -nbo_native reads the same archive, and the comparison
# is native against an external program on the original wavefunction data.
#
# -nbo_native still wants a positional wavefunction (convenience.cpp reads it before looking
# at -nbo_47), so the regenerated ch3.gbw is passed for it. That is defensible only because
# ch3's regeneration reproduced the recorded energy to 1.3e-8 Ha - ch3 is the one molecule
# whose starting geometry WAS the original optimised geometry, read out of this same archive.
# If the native numbers here disagree with the ones from ch3's own archive in the normal
# stage, the positional wavefunction is being used for more than the archive claims, and
# that is itself the finding.
#SBATCH --job-name=nboref-ch3arch
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/ch3arch-%j.out

set -u
HERE=${NBOREF_BIN:-/work/akkleemiss/florian/nbo_ref_v2/bin}
source "$HERE/env.sh"
D="$ROOT/ch3_original"   # named so compare_all.py sees molecule "ch3_original"
mkdir -p "$D" || exit 1
cd "$D" || exit 1

THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
NATIVE_FLAGS=$(python3 "$HERE/config.py" --native-flags ch3)
echo "=== ch3 original archive on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS"
echo "=== native flags: $NATIVE_FLAGS -nbo_threads $THREADS"

cp -f "$ROOT/nbo_reference/ch3_orca_reference.47" ch3_original.47 || exit 2
cp -f "$ROOT/ch3/ch3.gbw" . || exit 2

T0=$(date +%s)
"$HERE/gennbo7" ch3_original
RCG=$?
T1=$(date +%s)
"$NOSPHERA2" -nbo_parse ch3_original.nbo -nbo_json ch3_original.gennbo.nbo.json \
    > parse_ch3_original.log 2>&1
RCP=$?
T2=$(date +%s)
# shellcheck disable=SC2086
"$NOSPHERA2" -nbo_native ch3.gbw -nbo_47 ch3_original.47 $NATIVE_FLAGS \
    -nbo_threads "$THREADS" -nbo_json ch3_original.native.nbo.json \
    > native_ch3_original.log 2>&1
RCN=$?
T3=$(date +%s)

provenance_json ch3_original_archive ch3 "$THREADS" "$((T3 - T0))" "$((RCG + RCP + RCN))" \
  "\"archive\": \"tests/nbo_reference/ch3_orca_reference.47 (the original, NBAS=49 OPEN)\",
  \"nospher_a2\": \"$NOSPHERA2\",
  \"nospher_a2_commit\": \"$(cat "$(dirname "$NOSPHERA2")/DEPLOYED_COMMIT.txt" 2>/dev/null | tr -d '\n')\",
  \"native_flags\": \"$NATIVE_FLAGS -nbo_threads $THREADS\",
  \"seconds\": {\"gennbo\": $((T1 - T0)), \"parse\": $((T2 - T1)), \"native\": $((T3 - T2))},
  \"exit_codes\": {\"gennbo\": $RCG, \"parse\": $RCP, \"native\": $RCN}" \
  > provenance_nbo.json

echo "=== rc gennbo=$RCG parse=$RCP native=$RCN  seconds $((T1 - T0))/$((T2 - T1))/$((T3 - T2))  threads=$THREADS host=$(hostname)"
exit $((RCG + RCP + RCN))
