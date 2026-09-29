#!/bin/bash
#SBATCH --job-name=nboref-regen
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/regen-%j.out
#
# One truth per molecule.  nbo_ref_v2/ held gennbo JSONs written by the OLD parser (NRT "RS" read
# as a structure number, the composite alpha+beta valency table overwriting beta, the open-shell
# NAO Spin column read as Energy) while nbo_ref_v2_fixed/ held correct re-parses of the same .nbo
# text.  Two copies of one reference is how a measurement gets made against the wrong one, which
# already happened once here.  So: re-parse EVERY kept .nbo in place with the current binary, and
# re-run the native side in place too, so both files in a molecule's directory carry
# parser_version 2 and nothing in the tree predates the fix.
#
# This destroys no independence: -nbo_parse only re-reads gennbo 7's own output text, which is
# kept and unchanged.  The reference is still NBO 7's numbers; only our reading of them is new.
set -u
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
CFG=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2/config.py
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS

echo "=== nboref-regen on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS threads"
echo "=== binary $BIN"
echo "=== $(cd /work/akkleemiss/florian/nos_nboref && git log --oneline -1)"
[ -x "$BIN" ] || { echo "no binary"; exit 1; }

ok=0; bad=0
for D in /work/akkleemiss/florian/nbo_ref_v2/*/ /work/akkleemiss/florian/nbo_ref_radicals/*/; do
    MOL=$(basename "$D")
    [ -f "$D/$MOL.nbo" ] || continue
    cd "$D" || continue
    #The parse needs no flags - it re-reads gennbo's own text - so it must not be skipped when
    #config.py has no entry for the directory.  ch3_original, the surviving archive, is exactly
    #that case and it kept the old parse for one run because this check sat above the parse.
    $BIN -nbo_parse "$MOL.nbo" -nbo_json "$MOL.gennbo.nbo.json" > regen_parse_$MOL.log 2>&1
    p=$?
    FLAGS=$(python3 "$CFG" --native-flags "$MOL" 2>/dev/null)
    if [ -n "$FLAGS" ] && [ -f "$MOL.gbw" ]; then
        $BIN -nbo_native "$MOL.gbw" $FLAGS -nbo_47 "${MOL}_native.47" -nbo_threads $THREADS \
             -nbo_json "$MOL.native.nbo.json" > regen_native_$MOL.log 2>&1
        n=$?
    else
        n="no-native"  #no .gbw, or no config entry: a parse-only directory
    fi
    gv=$(python3 -c "import json;print(json.load(open('$MOL.gennbo.nbo.json')).get('parser_version'))" 2>/dev/null)
    nv=$(python3 -c "import json;print(json.load(open('$MOL.native.nbo.json')).get('parser_version'))" 2>/dev/null)
    echo "$MOL parse_exit=$p native_exit=$n gennbo_parser_version=$gv native_parser_version=$nv flags='$FLAGS'"
    if [ "$gv" = "2" ] && [ "$nv" = "2" ]; then ok=$((ok+1)); else bad=$((bad+1)); fi
done
echo "=== stamped both files: $ok molecules, incomplete: $bad"
