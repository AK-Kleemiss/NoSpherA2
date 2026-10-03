#!/bin/bash
#SBATCH --job-name=nbo-step3-stages
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/step3-%j.out
#
# Re-run the native side in place with NAO_DUMP_STEP3=1 so every molecule's directory carries the
# occupancies of stage 1 (pre-NAO) and stage 3 (after the Schmidt projection and OWSO) next to the
# final ones and gennbo's.  The dump only prints; the JSON it writes is the same one the reference
# run writes, with the same per-molecule flags from config.py, so this does not fork the data.
#
# The pre-registration and the acceptance gate for any later change to step 3 are in the docstring
# of step3_stages.py, which is also what reads these logs.
set -u
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
CFG=$V2/config.py
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
export NAO_DUMP_STEP3=1

echo "=== step3_stages on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS threads"
echo "=== binary $BIN"
echo "=== $(cd /work/akkleemiss/florian/nos_nboref && git log --oneline -1) plus the uncommitted pre_occ dump column"
[ -x "$BIN" ] || { echo "no binary"; exit 1; }

for D in /work/akkleemiss/florian/nbo_ref_v2/*/ /work/akkleemiss/florian/nbo_ref_radicals/*/; do
    MOL=$(basename "$D")
    [ -f "$D/$MOL.gbw" ] || continue
    cd "$D" || continue
    FLAGS=$(python3 "$CFG" --native-flags "$MOL" 2>/dev/null)
    [ -n "$FLAGS" ] || { echo "SKIP $MOL: no config entry"; continue; }
    $BIN -nbo_native "$MOL.gbw" $FLAGS -nbo_47 "${MOL}_native.47" -nbo_threads $THREADS \
         -nbo_json "$MOL.native.nbo.json" > step3_native_$MOL.log 2>&1
    echo "$MOL exit=$? step3_lines=$(grep -c '^STEP3 ' NoSpherA2.log)"
done

echo "=== analysis, 22 + 3 molecules in one table"
cd $V2 || exit 1
python3 step3_stages.py /work/akkleemiss/florian/nbo_ref_v2
python3 step3_stages.py /work/akkleemiss/florian/nbo_ref_radicals
