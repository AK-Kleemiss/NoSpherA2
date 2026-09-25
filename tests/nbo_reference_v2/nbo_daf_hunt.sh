#!/bin/bash
# The daf48 arm prints "/NBODAF / : NBO direct access file written on lfn48" and leaves NO file
# called lif.48 in the working directory.  A printed confirmation is not a file - that is the
# same class of blind instrument this lane keeps collecting - so this hunt asks WHERE, if
# anywhere, that write landed: the working directory (any name), $TMPDIR, /tmp, /dev/shm and the
# user's home, all by mtime inside the window of the run itself.
set -u
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
SRC=/work/akkleemiss/florian/nbo_ref_v2
W=/work/akkleemiss/florian/nbo_daf_hunt
MOL=lif
rm -rf "$W"; mkdir -p "$W"; cd "$W"
cp "$SRC/$MOL/${MOL}_native.47" "$MOL.47"
python3 -c '
import sys
p, keys = sys.argv[1], sys.argv[2]
lines = open(p).read().splitlines(True)
for i, line in enumerate(lines):
    if line.strip().upper().startswith("$NBO"):
        lines[i] = " $NBO %s $END\n" % keys
        break
else:
    raise SystemExit("no $NBO line")
open(p, "w").writelines(lines)
' "$MOL.47" "AOPNAO=W32 AONAO=W33 NBODAF=48"

echo "host=$(hostname) TMPDIR=${TMPDIR:-unset}"
# The window is taken from MARKER FILES, not from date(1): the first attempt used `date +%s` and
# the filesystem stamped every output 23 s later - an NFS server clock ahead of the login node -
# so the find window missed the entire run and reported 0 hits everywhere.  A marker file is
# stamped by the same clock as the thing being searched for.
: > "$W/.t0"; sleep 1
bash "$V2/gennbo7" "$MOL" > gennbo.log 2>&1; echo "rc=$?"
sleep 1; : > "$W/.t1"
echo "window: $(stat -c %y "$W/.t0") .. $(stat -c %y "$W/.t1")   (marker files, filesystem clock)"
echo "sanity: date(1) says $(date +%Y-%m-%d\ %H:%M:%S) - if that differs from the markers, the"
echo "        node and the filesystem disagree and only the markers may be used."
echo "--- working dir, ALL names:"
ls -la --time-style=+%H:%M:%S
echo "--- anything touched in the window elsewhere (each location prints its own count):"
for D in "${TMPDIR:-/tmp}" /tmp /dev/shm "$HOME" /work/akkleemiss/florian; do
    [ -d "$D" ] || { echo "    $D  ABSENT"; continue; }
    N=$(find "$D" -maxdepth 2 -newer "$W/.t0" ! -newer "$W/.t1" 2>/dev/null | wc -l)
    echo "    $D  $N hit(s) (searched $(find "$D" -maxdepth 2 2>/dev/null | wc -l) entries)"
    find "$D" -maxdepth 2 -newer "$W/.t0" ! -newer "$W/.t1" 2>/dev/null | head -10 | sed 's/^/        /'
done
echo "    CONTROL: the working dir itself must be a hit, or the window is blind ->"
echo "        $(find "$W" -maxdepth 1 -newer "$W/.t0" ! -newer "$W/.t1" 2>/dev/null | wc -l) of $(ls -A "$W" | wc -l) entries in $W"
echo "--- the line the run itself printed about the DAF:"
grep -n -i "daf\|lfn" "$MOL.nbo" | head
