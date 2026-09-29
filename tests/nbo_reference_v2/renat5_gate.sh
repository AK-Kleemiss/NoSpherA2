#!/bin/bash
# The before/after for the renat5 port, on the SAME 22 wavefunctions and the SAME .47 archives
# already staged in $ROOT, so nothing is regenerated and the gennbo side is untouched.
#
#   arm_renat5   the new binary, default            -> spec step 5 on the Rydberg set
#   arm_legacy   the new binary, NAO_LEGACY_CASCADE=1 -> the pre-change cascade
#   arm_oldbin   the PRESERVED pre-change binary, no env var
#
# arm_legacy against arm_oldbin is the control: if the env var really restores the old cascade and
# the rebuild moved nothing else, those two JSONs are identical.  Without it a difference between
# arm_renat5 and arm_oldbin could be the rebuild rather than the port.
#
# Stage B dumps the AO -> NAO matrix (NAO_DUMP_C) for the 8 closed-loop molecules from BOTH arms,
# on the very .47 files the numpy side reads, so the port-fidelity comparison is on identical
# arrays and nothing is converted.
#SBATCH --job-name=renat5-gate
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=/work/akkleemiss/florian/nao_renat5/slurm-%j.out

set -u
MY=/work/akkleemiss/florian/nao_renat5
ROOT=/work/akkleemiss/florian/nbo_ref_v2
NEW=$MY/bin/NoSpherA2.renat5
OLD=$MY/bin/NoSpherA2.baseline
MOLS="acetylene ammonia benzene ch3 ethane ethene formaldehyde formate hcn lif n2 ni_co_4 nitromethane no o2 ozone pf5 pyridine sf6 so2 ticl4 water"
FID="ammonia benzene ethane lif pf5 sf6 so2 water"
T=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=$T
export LD_LIBRARY_PATH=$MY/bin:${LD_LIBRARY_PATH:-}

echo "=== renat5 gate on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $T threads"
ls -la "$NEW" "$OLD"
md5sum "$NEW" "$OLD"

run_arm() {
    arm=$1; exe=$2; shift 2
    dir=$MY/$arm
    mkdir -p "$dir"
    fail=""
    for mol in $MOLS; do
        m=$dir/$mol
        mkdir -p "$m"
        src=$ROOT/$mol
        FLAGS=$(python3 "$MY/bin/config.py" --native-flags "$mol")
        t0=$(date +%s)
        # shellcheck disable=SC2086
        ( cd "$m" && env "$@" "$exe" -nbo_native "$src/$mol.gbw" -nbo_47 "$src/${mol}_native.47" \
            $FLAGS -nbo_threads "$T" -nbo_json "$m/$mol.native.nbo.json" \
            > "$m/native_$mol.log" 2>&1 )
        rc=$?
        cp "$src/$mol.gennbo.nbo.json" "$m/$mol.gennbo.nbo.json"
        printf '{"stage":"renat5_%s","molecule":"%s","threads":%d,"wall_seconds":%d,"exit_code":%d,"hostname":"%s","node_cores_total":%d,"binary":"%s","env":"%s","slurm_job_id":"%s","finished_utc":"%s"}\n' \
            "$arm" "$mol" "$T" "$(( $(date +%s) - t0 ))" "$rc" "$(hostname)" \
            "$(grep -c ^processor /proc/cpuinfo)" "$exe" "$*" "${SLURM_JOB_ID:-none}" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$m/provenance_nbo.json"
        [ $rc -eq 0 ] && [ -s "$m/$mol.native.nbo.json" ] || fail="$fail $mol"
        printf '%-11s %-13s rc=%d %5ds\n' "$arm" "$mol" "$rc" "$(( $(date +%s) - t0 ))"
    done
    echo "=== $arm done, failed:${fail:- none}"
}

echo "################ stage A: the 22-molecule arms"
run_arm arm_renat5 "$NEW" NAO_LEGACY_CASCADE=
run_arm arm_legacy "$NEW" NAO_LEGACY_CASCADE=1
run_arm arm_oldbin "$OLD" NAO_LEGACY_CASCADE=

echo "=== CONTROL arm_legacy vs arm_oldbin, must be identical or the env var is not the old cascade"
same=0; diff_=0
for mol in $MOLS; do
    if cmp -s "$MY/arm_legacy/$mol/$mol.native.nbo.json" "$MY/arm_oldbin/$mol/$mol.native.nbo.json"; then
        same=$((same + 1))
    else
        diff_=$((diff_ + 1)); echo "    DIFFERS: $mol"
    fi
done
echo "=== control: $same identical, $diff_ different of 22"

echo "=== DISCRIMINATION arm_renat5 vs arm_legacy, must differ or the port is a no-op"
moved=0
for mol in $MOLS; do
    cmp -s "$MY/arm_renat5/$mol/$mol.native.nbo.json" "$MY/arm_legacy/$mol/$mol.native.nbo.json" \
        || moved=$((moved + 1))
done
echo "=== discrimination: $moved of 22 molecules changed"

echo "################ stage B: NAO_DUMP_C on the 8 closed-loop molecules, both arms"
for mol in $FID; do
    s=$MY/sides/$mol
    for arm in renat5 legacy; do
        w=$MY/fid/$arm/$mol
        mkdir -p "$w"
        if [ "$arm" = legacy ]; then E=NAO_LEGACY_CASCADE=1; else E=NAO_LEGACY_CASCADE=; fi
        ( cd "$w" && env "$E" NAO_DUMP_C=1 NAO_DUMP_CPRE=1 "$NEW" -nbo_native "$s/$mol.gbw" \
            -nbo_47 "$s/$mol.47" -nbo_threads "$T" -nbo_json "$w/$mol.native.nbo.json" \
            > "$w/native_$mol.log" 2>&1 )
        rc=$?
        grep '^NAOC ' "$w/NoSpherA2.log" > "$w/$mol.naoc.txt"
        grep '^NAOCPRE ' "$w/NoSpherA2.log" > "$w/$mol.naocpre.txt"
        printf '%-7s %-9s rc=%d NAOC=%s NAOCPRE=%s\n' "$arm" "$mol" "$rc" \
            "$(wc -l < "$w/$mol.naoc.txt")" "$(wc -l < "$w/$mol.naocpre.txt")"
    done
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname), $T threads"
