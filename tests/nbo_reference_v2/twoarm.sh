#!/bin/bash
# Both arms of ONE binary over the held-out 16, native side only.
#
# Why this script exists: the column committed in d259879e was taken with a file labelled
# "60055d4b" by a hand-written DEPLOYED_COMMIT.txt, and that file is byte-identical
# (md5 1e45bc35f8d5e9617ee42dc62ca7d180) to the renat5 TREATMENT build.  A label beside a
# binary is not a property of the binary, so every arm below prints the md5sum of the file it
# actually invoked and the value of NAO_LEGACY_CASCADE it actually saw, per molecule, into the
# results stream.  Nothing is carried forward from the earlier column.
#
# DISCRIMINATION CHECK first: if the two arms agree on the first molecule, the arms are not
# wired up and the whole wave is void, so the script stops instead of producing 32 files that
# look like a before/after and are not.
#SBATCH --job-name=twoarm
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
set -u
H=/work/akkleemiss/florian/nbo_ref_v2_holdout
SRC=${SRC:-$H}
OUT=${OUT:-$H/twoarm}
BIN=${BIN:-$H/bin/NoSpherA2_519cb79}
MOLS=${MOLS:-"bf4_minus bh3 c2h5 cf3 ch2 cl2 clf3 hi hs mg_h2 nacl nh4_plus si2h6 sih4 thiophene zn_cl2"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"
BINMD5=$(md5sum "$BIN" | awk '{print $1}')
echo "host=$(hostname) threads=$THREADS"
echo "binary=$BIN md5=$BINMD5 bytes=$(stat -c %s "$BIN")"
echo "arms: baseline=NAO_LEGACY_CASCADE=1  renat5=<unset>"

run_one () {   # $1=mol $2=arm
    local MOL=$1 ARM=$2 W D FLAGS RC T0
    D=$SRC/$MOL
    W=$OUT/$ARM/$MOL
    [ -f "$D/${MOL}_native.47" ] || { echo "SKIP $MOL/$ARM: no ${MOL}_native.47"; return 1; }
    rm -rf "$W"; mkdir -p "$W"
    cp "$D/${MOL}_native.47" "$D/$MOL.gbw" "$D/$MOL.gennbo.nbo.json" "$W/" \
        || { echo "SKIP $MOL/$ARM: copy failed"; return 1; }
    FLAGS=$(python3 "$H/bin/config.py" --native-flags "$MOL")
    T0=$(date +%s)
    (
        cd "$W" || exit 1
        if [ "$ARM" = baseline ]; then export NAO_LEGACY_CASCADE=1; fi
        # printed from the same subshell that runs the binary, so a missing knob is visible
        echo "  $MOL/$ARM md5=$(md5sum "$BIN" | awk '{print $1}')" \
             "NAO_LEGACY_CASCADE=${NAO_LEGACY_CASCADE:-<unset>}" \
             "NAO_env=$(env | grep -c '^NAO_')"
        # shellcheck disable=SC2086
        "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" $FLAGS \
            -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native_$MOL.log" 2>&1
    )
    RC=$?
    echo "  $MOL/$ARM rc=$RC $(( $(date +%s) - T0 ))s"
    return $RC
}

FIRST=$(echo $MOLS | awk '{print $1}')
echo "=== discrimination check on $FIRST"
run_one "$FIRST" baseline || { echo "VOID: baseline arm failed on $FIRST"; exit 2; }
run_one "$FIRST" renat5   || { echo "VOID: renat5 arm failed on $FIRST"; exit 2; }
A=$OUT/baseline/$FIRST/$FIRST.native.nbo.json
B=$OUT/renat5/$FIRST/$FIRST.native.nbo.json
if /work/akkleemiss/florian/fp_fdp/bin/python3 "$H/bin/armdisc.py" "$A" "$B"; then
    echo "discrimination: PASS (the arms produce different numbers)"
else
    echo "VOID: the two arms agree on $FIRST - NAO_LEGACY_CASCADE is not wired up, no column taken"
    exit 3
fi

for MOL in $MOLS; do
    [ "$MOL" = "$FIRST" ] && continue
    run_one "$MOL" baseline
    run_one "$MOL" renat5
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
