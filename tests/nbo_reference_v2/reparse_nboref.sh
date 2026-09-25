#!/bin/bash
# Re-parse gennbo 7's 22 stored .nbo files with the FIXED parser, into a NEW root.
#
#   bash reparse_nboref.sh
#
# Nothing in $ROOT is overwritten: the fixed JSON is written to $FIXED/<mol>/, and the native
# side and provenance are symlinked in so the existing consumers (compare_all.py,
# nrt_spin_test.py, ratio_probe.py) run against it unchanged. Re-parsing needs only the stored
# .nbo text - no ORCA, no .gbw, no .47 - which is why that text was kept.
#
# The parse is text work on files of a few hundred KB, so this runs on the login node in
# seconds rather than through slurm; it still records its host and thread count.
set -u
ROOT=/work/akkleemiss/florian/nbo_ref_v2
FIXED=/work/akkleemiss/florian/nbo_ref_v2_fixed
NEW=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
THREADS=${OMP_NUM_THREADS:-1}

[ -x "$NEW" ] || { echo "no binary at $NEW - build first"; exit 1; }

echo "=== reparse on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, threads=$THREADS"
echo "=== binary $NEW"
echo "=== $(cd /work/akkleemiss/florian/nos_nboref && git log --oneline -1)"

mkdir -p "$FIXED"
FAIL=0
T0=$(date +%s)
for D in "$ROOT"/*/; do
    MOL=$(basename "$D")
    [ -f "$D/$MOL.nbo" ] || continue
    mkdir -p "$FIXED/$MOL"
    "$NEW" -nbo_parse "$D/$MOL.nbo" -nbo_json "$FIXED/$MOL/$MOL.gennbo.nbo.json" \
        > "$FIXED/$MOL/parse_$MOL.log" 2>&1
    RC=$?
    ln -sf "$D/$MOL.native.nbo.json" "$FIXED/$MOL/$MOL.native.nbo.json"
    ln -sf "$D/provenance_nbo.json" "$FIXED/$MOL/provenance_nbo.json"
    ln -sf "$D/$MOL.nbo" "$FIXED/$MOL/$MOL.nbo"
    if [ $RC -ne 0 ] || [ ! -s "$FIXED/$MOL/$MOL.gennbo.nbo.json" ]; then
        echo "FAIL $MOL rc=$RC"
        FAIL=$((FAIL + 1))
    fi
done
T1=$(date +%s)

cat > "$FIXED/provenance_reparse.json" << JSON
{
  "stage": "reparse_gennbo_json_with_fixed_parser",
  "hostname": "$(hostname)",
  "node_cores_total": $(grep -c ^processor /proc/cpuinfo),
  "threads": $THREADS,
  "wall_seconds": $((T1 - T0)),
  "molecules_failed": $FAIL,
  "binary": "$NEW",
  "binary_commit": "$(cd /work/akkleemiss/florian/nos_nboref && git log --format=%H -1)",
  "source_nbo_text": "$ROOT/<mol>/<mol>.nbo, unchanged",
  "native_side": "symlink to $ROOT/<mol>/<mol>.native.nbo.json, NOT re-run",
  "finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
JSON

echo "=== reparsed $(ls -d "$FIXED"/*/ | wc -l) molecules in $((T1 - T0))s, $FAIL failed"
cat "$FIXED/provenance_reparse.json"
exit $FAIL
