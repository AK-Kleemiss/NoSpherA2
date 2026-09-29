#!/bin/bash
# Re-run -nbo_native for the open-shell molecules with the binary that emits per-spin NAO tables,
# on the SAME stored .47 archive the first run read, and write the result into the fixed root only.
# The old root is not touched: the symlink reparse_nboref.sh left there is replaced by a real file
# inside $FIXED, so both runs stay on disk and comparable.
#
#   bash native_spin_rerun.sh [mol ...]      (default: every molecule with open_shell true)
set -u
ROOT=/work/akkleemiss/florian/nbo_ref_v2
FIXED=/work/akkleemiss/florian/nbo_ref_v2_fixed
NEW=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
THREADS=${THREADS:-8}

[ -x "$NEW" ] || { echo "no binary at $NEW"; exit 1; }
echo "=== native rerun on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, -nbo_threads $THREADS"
echo "=== $(cd /work/akkleemiss/florian/nos_nboref && git log --oneline -1) plus uncommitted per-spin NAO patch"

MOLS="$*"
if [ -z "$MOLS" ]; then
    MOLS=$(cd "$ROOT" && for d in */; do
        m=${d%/}
        grep -q '"open_shell": true' "$m/$m.gennbo.nbo.json" 2>/dev/null && echo "$m"
    done)
fi
echo "=== molecules: $MOLS"

T0=$(date +%s)
FAIL=0
for MOL in $MOLS; do
    D="$ROOT/$MOL"
    [ -f "$D/${MOL}_native.47" ] || { echo "SKIP $MOL: no stored ${MOL}_native.47"; continue; }
    mkdir -p "$FIXED/$MOL"
    rm -f "$FIXED/$MOL/$MOL.native.nbo.json"   #the symlink into the old root
    NATIVE_FLAGS=$(python3 "$ROOT/bin/config.py" --native-flags "$MOL" 2>/dev/null)
    T=$(date +%s)
    (cd "$D" && "$NEW" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" $NATIVE_FLAGS \
        -nbo_threads "$THREADS" -nbo_json "$FIXED/$MOL/$MOL.native.nbo.json" \
        > "$FIXED/$MOL/native_rerun_$MOL.log" 2>&1)
    RC=$?
    W=$(( $(date +%s) - T ))
    if [ $RC -ne 0 ] || [ ! -s "$FIXED/$MOL/$MOL.native.nbo.json" ]; then
        echo "FAIL $MOL rc=$RC after ${W}s"; FAIL=$((FAIL + 1))
    else
        echo "ok   $MOL in ${W}s, flags: ${NATIVE_FLAGS:-none}"
    fi
done
T1=$(date +%s)

cat > "$FIXED/provenance_native_spin_rerun.json" <<JSON
{
  "stage": "native_rerun_with_per_spin_nao_tables",
  "hostname": "$(hostname)",
  "node_cores_total": $(grep -c ^processor /proc/cpuinfo),
  "nbo_threads": $THREADS,
  "wall_seconds": $((T1 - T0)),
  "molecules": "$MOLS",
  "molecules_failed": $FAIL,
  "binary": "$NEW",
  "binary_commit_plus_uncommitted_patch": "$(cd /work/akkleemiss/florian/nos_nboref && git log --format=%H -1)",
  "input": "$ROOT/<mol>/<mol>_native.47 and <mol>.gbw, unchanged",
  "finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON
echo "=== $((T1 - T0))s total, $FAIL failed"
exit $FAIL
