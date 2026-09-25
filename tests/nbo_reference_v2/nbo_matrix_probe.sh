#!/bin/bash
# PROBE: which of the four never-used NBO 7 matrix dumps this install actually honours, on lif.
#
# Nothing physical is measured here.  The only question is which keywords the BINARY accepts,
# which unit each matrix lands on, and which files are non-empty - because this lane has twice
# mistaken a 0-byte file with rc=0 for an answer, and AONAO=W48 is REFUSED because 48 is
# reserved.  Per the manual p49-50 only PAOPNAO(29) AOPNAO(32) AONAO(33) AOPNHO(34) AONHO(35)
# AOPNBO(36) AONBO(37) AOPNLMO(38) AONLMO(39) AOMO(40) AONO(41) DMAO(42) NAONBO(30) have their
# own default LFN and EVERY OTHER MATRIX DEFAULTS TO 49 - so DMPNAO and SPNAO would collide on
# one file if their units were left implicit.  Whether the binary honours an explicit W<n> for
# an OPERATOR matrix (the manual documents W[n] generically, not per keyword) is exactly what
# arm C measures.  NEVER ASSUME A UNIT NUMBER: every file each arm creates is listed with its
# size and the unit is read off that listing.
#
# Arms, each in its own directory so no file can be inherited from the arm before it:
#   base   PAOPNAO=W29 AOPNAO=W32 AONAO=W33 NRTSYM=OFF
#   rpnao  base + RPNAO                (p54: "Revises PAO to PNAO transformation matrix by
#                                       post-multiplying by TRyd and Tred")
#   ops    base + DMPNAO=W51 SPNAO=W52 DMNAO=W53 DMAO=W54   (p48 lists all four)
#   daf    base + NBODAF=60            (p54: writes the NBO direct access file to LFN n)
set -u
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
SRC=/work/akkleemiss/florian/nbo_ref_v2
OUT=${OUT:-/work/akkleemiss/florian/nbo_mat}
MOL=${MOL:-lif}
mkdir -p "$OUT"

echo "host=$(hostname) date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) mol=$MOL"

arm() {
    NAME=$1; KEYS=$2
    W=$OUT/$MOL.$NAME
    rm -rf "$W"; mkdir -p "$W"
    cp "$SRC/$MOL/${MOL}_native.47" "$W/$MOL.47"
    python3 - "$W/$MOL.47" "$KEYS" <<'PY' || return 1
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
PY
    ( cd "$W" && bash "$V2/gennbo7" "$MOL" >gennbo.log 2>&1 ); RC=$?
    echo "--- arm=$NAME rc=$RC keylist: $KEYS"
    cat "$W/gennbo.log"
    # every file, with its size.  The unit number is read off THIS listing, never assumed.
    ( cd "$W" && for f in $(ls -A); do printf '    %-28s %10s\n' "$f" "$(stat -c %s "$f")"; done )
    # what the run itself says: an unrecognised keyword, and every matrix header it wrote.
    echo "  --- keylist echo / complaints:"
    grep -n -i -E "unrecogni|unknown|illegal|invalid|ignor|not recogni|SETFN|abort|error" "$W/$MOL.nbo" "$W/$MOL.gennbo.err" 2>/dev/null | head -20
    sed -n '1,40p' "$W/$MOL.nbo" | grep -n -E "\\\$NBO|NBO options|keylist" | head
    echo "  --- header line of every text unit written (first line of each non-.47/.nbo file):"
    ( cd "$W" && for f in $(ls -A); do
        case "$f" in *.47|*.nbo|gennbo.log|*.err) continue;; esac
        printf '    %-28s | %s\n' "$f" "$(head -c 200 "$f" | head -1 | tr -d '\r')"
      done )
}

BASE="PAOPNAO=W29 AOPNAO=W32 AONAO=W33 NRTSYM=OFF"
arm base  "$BASE"
arm rpnao "$BASE RPNAO"
arm ops   "$BASE DMPNAO=W51 SPNAO=W52 DMNAO=W53 DMAO=W54"
arm daf   "$BASE NBODAF=60"
echo "=== probe done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
