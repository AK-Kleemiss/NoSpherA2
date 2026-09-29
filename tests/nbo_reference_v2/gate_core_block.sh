#!/bin/bash
# The 22-molecule acceptance gate for the step-4 core-block fix (nbo_external_reference 2114ee1).
#
#   sbatch gate_core_block.sh
#
# Three arms on the SAME wavefunctions and the SAME archives already staged in $ROOT, so nothing
# is regenerated and the gennbo side is untouched:
#
#   arm_core    the new binary, default -> the core gets its own step-4 block
#   arm_pooled  the new binary, NAO_CORE_POOLED=1 -> the old fully pooled block
#   arm_old     the PRESERVED 7a466b07 binary, no env var
#
# arm_pooled against arm_old is the control: if the env var really restores the old form and the
# rebuild moved nothing else, those two JSONs are identical.  Without it a difference between
# arm_core and arm_old could be any of the three commits in between rather than the partition.
#SBATCH --job-name=nboref-gate
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/gate-%j.out

set -u
BIN=/work/akkleemiss/florian/nbo_ref_v2/bin
source "$BIN/env.sh"
TREE=/work/akkleemiss/florian/nos_nboref
NEW=$TREE/build/release-linux/bin/NoSpherA2
OLD=$TREE/build/release-linux/bin/NoSpherA2.7a466b07_pooled
GATE=/work/akkleemiss/florian/nbo_ref_v2/gate_2114ee1
MOLS="acetylene ammonia benzene ch3 ethane ethene formaldehyde formate hcn lif n2 ni_co_4 nitromethane no o2 ozone pf5 pyridine sf6 so2 ticl4 water"
T=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=$T

echo "=== gate on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $T threads"
echo "=== new $NEW"
echo "=== old $OLD"
ls -la "$NEW" "$OLD"

run_arm() {
    # run_arm <arm name> <binary> [env assignment ...]
    arm=$1; exe=$2; shift 2
    dir=$GATE/$arm
    mkdir -p "$dir"
    fail=""
    for mol in $MOLS; do
        m=$dir/$mol
        mkdir -p "$m"
        src=$ROOT/$mol
        FLAGS=$(python3 "$BIN/config.py" --native-flags "$mol")
        t0=$(date +%s)
        # shellcheck disable=SC2086
        env "$@" "$exe" -nbo_native "$src/$mol.gbw" -nbo_47 "$src/${mol}_native.47" $FLAGS \
            -nbo_threads "$T" -nbo_json "$m/$mol.native.nbo.json" > "$m/native_$mol.log" 2>&1
        rc=$?
        ln -sf "$src/$mol.gennbo.nbo.json" "$m/$mol.gennbo.nbo.json"
        printf '{"stage":"gate_%s","molecule":"%s","threads":%d,"wall_seconds":%d,"exit_code":%d,"hostname":"%s","node_cores_total":%d,"binary":"%s","env":"%s","slurm_job_id":"%s","finished_utc":"%s"}\n' \
            "$arm" "$mol" "$T" "$(( $(date +%s) - t0 ))" "$rc" "$(hostname)" \
            "$(grep -c ^processor /proc/cpuinfo)" "$exe" "$*" "${SLURM_JOB_ID:-none}" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$m/provenance_gate.json"
        [ $rc -eq 0 ] && [ -s "$m/$mol.native.nbo.json" ] || fail="$fail $mol"
    done
    echo "=== $arm done, failed:${fail:- none}"
    [ -z "$fail" ]
}

run_arm arm_core   "$NEW" NAO_CORE_POOLED=
RC_CORE=$?
run_arm arm_pooled "$NEW" NAO_CORE_POOLED=1
RC_POOL=$?
run_arm arm_old    "$OLD" NAO_CORE_POOLED=
RC_OLD=$?

echo "=== CONTROL: arm_pooled vs arm_old, must be identical or the env var is not the old form"
same=0; diff_=0
for mol in $MOLS; do
    if cmp -s "$GATE/arm_pooled/$mol/$mol.native.nbo.json" "$GATE/arm_old/$mol/$mol.native.nbo.json"; then
        same=$((same + 1))
    else
        diff_=$((diff_ + 1)); echo "    DIFFERS: $mol"
    fi
done
echo "=== control: $same identical, $diff_ different of 22"

echo "=== DISCRIMINATION: arm_core vs arm_pooled, must differ somewhere or the arm is a no-op"
moved=0
for mol in $MOLS; do
    cmp -s "$GATE/arm_core/$mol/$mol.native.nbo.json" "$GATE/arm_pooled/$mol/$mol.native.nbo.json" \
        || moved=$((moved + 1))
done
echo "=== discrimination: $moved of 22 molecules changed"

for arm in arm_core arm_pooled arm_old; do
    echo
    echo "############################## compare_all $arm"
    python3 "$BIN/compare_all.py" "$GATE/$arm" --json "$GATE/$arm.summary.json" 2>&1 | tail -45
    echo "############################## nao_class_leak $arm"
    python3 "$BIN/nao_class_leak.py" "$GATE/$arm" 2>&1 | tail -35
done
echo "=== rc core=$RC_CORE pooled=$RC_POOL old=$RC_OLD"
exit $((RC_CORE + RC_POOL + RC_OLD))
