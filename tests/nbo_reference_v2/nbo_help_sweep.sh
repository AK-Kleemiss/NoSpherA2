#!/bin/bash
# STEP A of the secondary task: ask the BINARY itself what it can write.  `$NBO HELP $END` is
# the authority - not the manual, not memory of the manual.  Nothing is inferred here; the whole
# help text is echoed back and every keyword that mentions a matrix is listed with its own line,
# so the next step can probe by name rather than by recollection.
#
# The question this serves: is there ANY keyword that writes a matrix BETWEEN the Schmidt stage
# (step 3/4 boundary of the published cascade) and the final NAO re-diagonalisation?  Unit 32
# (AOPNAO) is before it and unit 33 (AONAO) is after it, and nothing in between has ever been
# read on this install.
set -u
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
SRC=/work/akkleemiss/florian/nbo_ref_v2
OUT=${OUT:-/work/akkleemiss/florian/nbo_help}
MOL=${MOL:-lif}
rm -rf "$OUT"; mkdir -p "$OUT"
cd "$OUT"

echo "host=$(hostname) date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) mol=$MOL"
for f in /work/software/bin/NBO/bin/gennbo.i8.exe /work/software/bin/NBO/bin/nbo7.i8.exe; do
    printf 'BIN %-44s %10s  %s\n' "$f" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -c1-32)"
done

cp "$SRC/$MOL/${MOL}_native.47" "$MOL.47"
# rewrite the $NBO line to HELP, with python so no quoting can mangle it
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
' "$MOL.47" "HELP"

bash "$V2/gennbo7" "$MOL"; echo "rc=$?"
echo "--- files written, with sizes (unit numbers are READ HERE, never assumed):"
for f in $(ls -A); do printf '    %-28s %10s\n' "$f" "$(stat -c %s "$f")"; done
echo "=== FULL HELP TEXT BEGINS ($(wc -l < $MOL.nbo) lines in $MOL.nbo) ==="
cat "$MOL.nbo"
echo "=== FULL HELP TEXT ENDS ==="
