#!/bin/bash
# Re-run ONLY the native side of the already-stamped 22 with the NAMED baseline binary, so that
# "the shipped native's worst Rydberg error is 0.3161 e" can be checked against the binary the
# held-out set was measured with instead of against the share binary that produced the published
# number.  gennbo's side is not touched: its stored JSON is copied in and reused, so the only
# thing that changes between the published column and this one is which NoSpherA2 ran.
#SBATCH --job-name=xcheck
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00
set -u
H=/work/akkleemiss/florian/nbo_ref_v2_holdout
SRC=${SRC:-/work/akkleemiss/florian/nbo_ref_v2}
OUT=${OUT:-$H/xcheck}
BIN=${BIN:-$H/bin/NoSpherA2_baseline}
MOLS=${MOLS:-"lif water ammonia ethane benzene pf5 so2 sf6 pyridine ch3"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"
echo "host=$(hostname) threads=$THREADS binary=$BIN commit=$(cat "$(dirname "$BIN")/DEPLOYED_COMMIT.txt" 2>/dev/null)"
for MOL in $MOLS; do
    D=$SRC/$MOL
    W=$OUT/$MOL
    [ -f "$D/${MOL}_native.47" ] || { echo "SKIP $MOL: no ${MOL}_native.47"; continue; }
    rm -rf "$W"; mkdir -p "$W"
    cp "$D/${MOL}_native.47" "$D/$MOL.gbw" "$D/$MOL.gennbo.nbo.json" "$W/" || { echo "SKIP $MOL: copy failed"; continue; }
    FLAGS=$(python3 "$H/bin/config.py" --native-flags "$MOL")
    T0=$(date +%s)
    # shellcheck disable=SC2086
    ( cd "$W" && "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" $FLAGS \
        -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native_$MOL.log" 2>&1 )
    RC=$?
    echo "$MOL rc=$RC $(( $(date +%s) - T0 ))s flags: $FLAGS"
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
