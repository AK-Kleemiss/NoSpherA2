#!/bin/bash
# CAVEAT 2026-09-25: the file this ran, bin/NoSpherA2_baseline, is byte-identical to the renat5
# TREATMENT build (md5 1e45bc35f8d5e9617ee42dc62ca7d180) - the name and the DEPLOYED_COMMIT.txt
# beside it were wrong.  No arm below sets NAO_LEGACY_CASCADE=1, so none of them is a baseline and
# the job cannot bracket the pre-change number.  What it DOES still prove is that environment
# variables reach the binary (the class_split control moves benzene) - nothing more.
# Four arms of ONE binary on benzene, with a POSITIVE CONTROL, because the first
# attempt at this (job 594818) set NAO_CORE_POOLED=1 and never proved the variable reached the
# binary - a silent knob reads exactly like a knob with no effect.  NAO_CLASS_SPLIT=1 is the
# control: it is known to move the number (the class-split arm was measured 5.6x worse), so if it
# moves here the plumbing works and a flat NAO_CORE_POOLED is a real exclusion; if it does NOT
# move, every env-gated exclusion in this tree is void and none may be quoted.
#SBATCH --job-name=naoarms
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
set -u
H=/work/akkleemiss/florian/nbo_ref_v2_holdout
SRC=$H/../nbo_ref_v2
OUT=$H/naoarms
BIN=$H/bin/NoSpherA2_baseline
MOLS=${MOLS:-"benzene ethane"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"
echo "host=$(hostname) threads=$THREADS binary=$BIN commit=$(cat "$H/bin/DEPLOYED_COMMIT.txt")"
for MOL in $MOLS; do
    for ARM in default core_pooled owso_off class_split; do
        W=$OUT/$MOL.$ARM
        rm -rf "$W"; mkdir -p "$W"
        cp "$SRC/$MOL/${MOL}_native.47" "$SRC/$MOL/$MOL.gbw" "$SRC/$MOL/$MOL.gennbo.nbo.json" "$W/" || continue
        FLAGS=$(python3 "$H/bin/config.py" --native-flags "$MOL")
        (
            cd "$W" || exit 1
            case $ARM in
                core_pooled) export NAO_CORE_POOLED=1 ;;
                owso_off) export NAO_OWSO_OFF=1 ;;
                class_split) export NAO_CLASS_SPLIT=1 ;;
            esac
            # the knob, as the process that runs the binary actually sees it - printed from the
            # same subshell, one line, so a missing arm is visible instead of assumed
            echo "  arm=$ARM env: $(env | grep -E '^NAO_' | tr '\n' ' ')none_beyond_this"
            # shellcheck disable=SC2086
            "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" $FLAGS \
                -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native.log" 2>&1
            echo "  arm=$ARM rc=$?"
        )
    done
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
