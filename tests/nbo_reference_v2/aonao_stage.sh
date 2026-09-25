#!/bin/bash
#SBATCH --job-name=aonao-stage
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=/work/akkleemiss/florian/aonao_cmp/slurm-%j.out
#
# Produce the two sides of the first ARBITRATED INTERMEDIATE comparison: the AO -> NAO
# transformation matrix, from gennbo 7 and from -nbo_native, on the same .47.
#
#   gennbo side   $NBO AONAO=W $END  ->  <mol>.33, nine decimals (probe job 588007 verified the
#                 keyword, the default unit 33, and that AONAO=W48 is refused - 48 is reserved)
#   native side   NAO_DUMP_C=1       ->  NAOC lines in NoSpherA2.log, one per NAO with its
#                 (atom, l, m, shell, class, occupancy) and its AO coefficients
#
# Why this is a fair comparison at all: both sides read the SAME <mol>_native.47, so the AO
# order, the overlap S and the density P are literally the same arrays.  Nothing is converted.
#
# Why the gennbo keylist here is AONAO=W ALONE and not the reference keylist plus AONAO: NRT on
# benzene is minutes and the NAO stage is upstream of NBO, E2 and NRT, so the matrix cannot
# depend on them.  That is an assumption, so it is measured rather than asserted: this run's own
# .nbo is parsed and its NAO occupancy table must equal the stamped reference JSON's.  If it
# does not, the shortcut is wrong and the comparison is void - checked in aonao_compare.py.
#
# Molecules: the two ends of the leak plus a sanity case.  ethane is the worst Rydberg excess
# (6.2x gennbo), benzene the largest molecule with a bad one (3.7x), pf5/so2/sf6 the only three
# that are already right, water/ammonia cheap controls, lif small enough to read by eye.
set -u
OUT=/work/akkleemiss/florian/aonao_cmp
SRC=/work/akkleemiss/florian/nbo_ref_v2
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
MOLS=${MOLS:-"lif water ammonia ethane benzene pf5 so2 sf6"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"

echo "host=$(hostname)  threads=$THREADS  node_cores=$(grep -c ^processor /proc/cpuinfo)  job=${SLURM_JOB_ID:-none}  partition=${SLURM_JOB_PARTITION:-none}"
echo "binary $BIN"
[ -x "$BIN" ] || { echo "no binary - run build_nboref.sh first"; exit 1; }

for MOL in $MOLS; do
    D=$SRC/$MOL
    W=$OUT/$MOL
    [ -f "$D/${MOL}_native.47" ] || { echo "SKIP $MOL: no ${MOL}_native.47"; continue; }
    rm -rf "$W"; mkdir -p "$W"
    cp "$D/${MOL}_native.47" "$W/$MOL.47"
    cp "$D/$MOL.gbw" "$W/" 2>/dev/null

    # --- gennbo side: AONAO=W alone, so no NRT time is spent on a matrix it cannot affect.
    # sf6 needs one more keyword.  Its first run wrote a 0-byte lfn 33 with rc=0 because NBO stopped
    # at "SYMOPS: generated 48 symmetry operator(s) for Th but expected 24" after 0.02 CPU seconds -
    # before the NAO stage, so the abort has nothing to do with NRT even though NRTSYM=OFF is what
    # suppresses it.  It is kept out of the other seven's keylist deliberately: their numbers were
    # produced with AONAO=W alone, and a keylist change would make them a different run.
    EXTRA_KEYS=${EXTRA_KEYS:-}
    [ "$MOL" = sf6 ] && EXTRA_KEYS="NRTSYM=OFF"
    python3 - "$W/$MOL.47" "$EXTRA_KEYS" <<'PY' || exit 1
import sys
p = sys.argv[1]
extra = (" " + sys.argv[2].strip()) if len(sys.argv) > 2 and sys.argv[2].strip() else ""
lines = open(p).read().splitlines(True)
for i, line in enumerate(lines):
    if line.strip().upper().startswith("$NBO"):
        lines[i] = " $NBO AONAO=W%s $END\n" % extra
        break
else:
    raise SystemExit("no $NBO line in " + p)
open(p, "w").writelines(lines)
print("keylist now: " + lines[i].strip())
PY
    T0=$(date +%s)
    ( cd "$W" && bash "$V2/gennbo7" "$MOL" ) || echo "gennbo rc nonzero for $MOL"
    T1=$(date +%s)
    ( cd "$W" && "$BIN" -nbo_parse "$MOL.nbo" -nbo_json "$MOL.aonao.nbo.json" > "parse_$MOL.log" 2>&1 )
    RCP=$?

    # --- native side: the same .47, the dump on
    T2=$(date +%s)
    ( cd "$W" && NAO_DUMP_C=1 "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "$MOL.47" \
        -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native_$MOL.log" 2>&1 )
    RCN=$?
    T3=$(date +%s)
    ( cd "$W" && grep '^NAOC ' NoSpherA2.log > "$MOL.naoc.txt" || true )
    printf '%-9s gennbo %3ds parse rc=%d native %3ds rc=%d  lfn33=%s bytes  NAOC lines=%s\n' \
        "$MOL" "$((T1 - T0))" "$RCP" "$((T3 - T2))" "$RCN" \
        "$(stat -c %s "$W/$MOL.33" 2>/dev/null || echo MISSING)" \
        "$(wc -l < "$W/$MOL.naoc.txt" 2>/dev/null || echo 0)"
done

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname), $THREADS threads"
