#!/bin/bash
# Stage 2: the external reference and the in-house answer, on ONE accepted wavefunction.
#
#   sbatch nbo_stage.sh <molecule>
#
# Four steps, in this order:
#
#   1. NoSpherA2 -convert_to_47 writes <mol>.47 with the $NBO keylist config.py gives.
#   2. <mol>_native.47 is a copy taken BEFORE step 3, because for so2 the gennbo archive
#      gets a hand-written $NRTSTR appended and the native reader has no business seeing it.
#   3. gennbo 7 runs on <mol>.47 -> <mol>.nbo, and NoSpherA2 -nbo_parse turns that into
#      <mol>.gennbo.nbo.json with the FIXED parser (b4584267: the NRT RS column is a rank,
#      not a structure number, and an open-shell NAO table prints Spin where the old parser
#      read Energy).
#   4. NoSpherA2 -nbo_native runs the in-house analysis on the same wavefunction and the
#      same archive -> <mol>.native.nbo.json.
#
# Both JSONs then have independent provenance - one is gennbo 7's output, the other is our
# own code - which is the whole point: `compare_nbo.py --all .` inside tests/nbo_reference
# compares that dataset with itself and reports 22/22 at 0.00000, and that is not a gate.
#SBATCH --job-name=nboref-nbo
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/nbo-%j.out

set -u
MOL=$1
HERE=${NBOREF_BIN:-/work/akkleemiss/florian/nbo_ref_v2/bin}
source "$HERE/env.sh"
cd "$ROOT/$MOL" || exit 1

THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
KW=$(python3 "$HERE/config.py" --nbo-keywords "$MOL")
NATIVE_FLAGS=$(python3 "$HERE/config.py" --native-flags "$MOL")
echo "=== $MOL on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS"
echo "=== gennbo keylist : $KW"
echo "=== native flags   : $NATIVE_FLAGS -nbo_threads $THREADS"

rm -f "$MOL.47" "${MOL}_native.47" "$MOL.nbo" "$MOL.gennbo.nbo.json" "$MOL.native.nbo.json"

T0=$(date +%s)
"$NOSPHERA2" -convert_to_47 "$MOL.gbw" -nbo_keywords "$KW" > "convert47_$MOL.log" 2>&1
RC47=$?
T1=$(date +%s)
[ -f "$MOL.47" ] || { echo "no $MOL.47 written"; tail -20 "convert47_$MOL.log"; exit 3; }
cp "$MOL.47" "${MOL}_native.47"
python3 "$HERE/config.py" --nrtstr "$MOL" >> "$MOL.47"

"$HERE/gennbo7" "$MOL"
RCG=$?
T2=$(date +%s)
"$NOSPHERA2" -nbo_parse "$MOL.nbo" -nbo_json "$MOL.gennbo.nbo.json" > "parse_$MOL.log" 2>&1
RCP=$?
T3=$(date +%s)
# shellcheck disable=SC2086
"$NOSPHERA2" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" $NATIVE_FLAGS \
    -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native_$MOL.log" 2>&1
RCN=$?
T4=$(date +%s)

provenance_json nbo_stage "$MOL" "$THREADS" "$((T4 - T0))" "$((RC47 + RCG + RCP + RCN))" \
  "\"nospher_a2\": \"$NOSPHERA2\",
  \"nospher_a2_commit\": \"$(cat "$(dirname "$NOSPHERA2")/DEPLOYED_COMMIT.txt" 2>/dev/null | tr -d '\n')\",
  \"gennbo_bin\": \"$NBOBIN\",
  \"nbo_keywords\": \"$KW\",
  \"native_flags\": \"$NATIVE_FLAGS -nbo_threads $THREADS\",
  \"seconds\": {\"convert_47\": $((T1 - T0)), \"gennbo\": $((T2 - T1)), \"parse\": $((T3 - T2)), \"native\": $((T4 - T3))},
  \"exit_codes\": {\"convert_47\": $RC47, \"gennbo\": $RCG, \"parse\": $RCP, \"native\": $RCN}" \
  > provenance_nbo.json

echo "=== $MOL rc 47=$RC47 gennbo=$RCG parse=$RCP native=$RCN"
echo "=== seconds 47=$((T1 - T0)) gennbo=$((T2 - T1)) parse=$((T3 - T2)) native=$((T4 - T3))  threads=$THREADS host=$(hostname)"
exit $((RC47 + RCG + RCP + RCN))
