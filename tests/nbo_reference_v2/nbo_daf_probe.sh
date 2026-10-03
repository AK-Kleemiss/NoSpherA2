#!/bin/bash
# STEP B of the secondary task.  Pages 41-55 of the SHIPPED manual are the authority (`$NBO HELP
# $END` is not a keyword in NBO 7.0.9 - the binary answers "Unrecognized $NBO keyword: HELP"),
# and they leave exactly two candidates unmeasured on this install:
#
#   NBODAF   p54: "Request writing the NBO direct access file (DAF) to external disk file LFN n,
#            or, if 'n' is not present, to the default LFN 48."  The earlier probe used
#            NBODAF=60 and got nothing.  48 is ALSO the unit that AONAO=W48 is refused for as
#            "reserved", which is consistent with 48 being the DAF's own unit - so the earlier
#            null may have been the unit, not the keyword.
#   BOAO     p54, default LFN 49 - a documented keyword whose file MUST appear, used here as the
#            positive control, so a run that writes nothing cannot be read as "the keyword does
#            nothing".
#
# Every file each arm creates is listed with its size and its first line; the unit is READ off
# that listing, never assumed.  Arms are in separate directories so nothing can be inherited.
set -u
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
SRC=/work/akkleemiss/florian/nbo_ref_v2
OUT=${OUT:-/work/akkleemiss/florian/nbo_daf}
MOL=${MOL:-lif}
rm -rf "$OUT"; mkdir -p "$OUT"

echo "host=$(hostname) date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) mol=$MOL"
for f in /work/software/bin/NBO/bin/gennbo.i8.exe /work/software/bin/NBO/bin/nbo7.i8.exe; do
    printf 'BIN %-44s %10s  %s\n' "$f" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -c1-32)"
done

BASE="PAOPNAO=W29 AOPNAO=W32 AONAO=W33"

arm() {
    NAME=$1; KEYS=$2
    W=$OUT/$NAME
    rm -rf "$W"; mkdir -p "$W"
    cp "$SRC/$MOL/${MOL}_native.47" "$W/$MOL.47"
    python3 -c '
import sys
p, keys = sys.argv[1], sys.argv[2]
lines = open(p).read().splitlines(True)
for i, line in enumerate(lines):
    if line.strip().upper().startswith("$NBO"):
        lines[i] = " $NBO %s $END\n" % keys
        break
else:
    raise SystemExit("no $NBO line in " + p)
open(p, "w").writelines(lines)
' "$W/$MOL.47" "$KEYS" || return 1
    ( cd "$W" && bash "$V2/gennbo7" "$MOL" >gennbo.log 2>&1 ); RC=$?
    echo "--- arm=$NAME rc=$RC keylist: $KEYS"
    grep -n -i -E "unrecogni|unknown|illegal|invalid|ignor|not recogni|reserved|abort|error" \
        "$W/$MOL.nbo" "$W/$MOL.gennbo.err" 2>/dev/null | head -8 | sed 's/^/      complaint /'
    ( cd "$W" && for f in $(ls -A); do
        SZ=$(stat -c %s "$f")
        case "$f" in *.47|*.nbo|gennbo.log|*.err) H="";; *) H=$(head -c 120 "$f" | head -1 | tr -cd '\40-\176');; esac
        printf '      %-20s %10s  %s\n' "$f" "$SZ" "$H"
      done )
}

arm base    "$BASE"
arm boao    "$BASE BOAO=W50 SVEC"
arm daf_def "$BASE NBODAF"
arm daf48   "$BASE NBODAF=48"
arm daf50   "$BASE NBODAF=50"

echo "=== diff of the file SETS against base (what each extra keyword actually produced):"
cd "$OUT"
for a in boao daf_def daf48 daf50; do
    echo "  $a: $(comm -13 <(ls -A base | sort) <(ls -A $a | sort) | tr '\n' ' ')"
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
