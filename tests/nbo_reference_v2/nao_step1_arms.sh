#!/bin/bash
#SBATCH --job-name=nao-step1-arms
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=/work/akkleemiss/florian/nao_step1_exp/slurm-%j.out
#
# Three arms of ONE experiment on the SAME 25 wavefunctions, all with NAO_DUMP_STEP3=1 so the
# class totals after step 3 are visible as well as the final ones:
#
#   baseline    what the code ships
#   owso_off    NAO_OWSO_OFF=1 - the occupancy weighting inside each class is dropped and the
#               plain Loewdin does the whole within-class orthogonalisation
#   prenao_net  NAO_PRENAO_NET=1 - step 1 builds the pre-NAOs from the atom's own block of P with
#               S^A as the metric (the net atomic population, which is how the 1985 paper reads)
#               instead of from (S P S)^A
#
# PRE-REGISTERED, before the run:
#
#   owso_off is NOT a candidate.  It throws away the point of the OWSO and its final numbers are
#   expected to be worse.  It is the test of one claim: after step 3 the Rydberg set is the
#   S-orthogonal complement of the span of the natural minimal pre-NAOs, so each class TOTAL at
#   that stage depends only on that span and on nothing inside step 3.  Prediction: every
#   molecule's s3_Ryd is unchanged to within 1e-6 e while fin_Ryd moves.  What refutes it: an
#   s3_Ryd that moves.  Then step 3's internals are a lever after all and the search cannot be
#   narrowed the way the next arm assumes.
#
#   prenao_net IS a candidate, and it gets no directional prediction.  The comment in nao.cpp
#   rejected it on epoxide's NPA charges (half an electron per carbon), and the class-split arm
#   showed NPA charges cannot see a class-partition error at all - so that rejection is not
#   evidence about the Rydberg metric either way.  It is judged by the gate already fixed in
#   step3_stages.py: pf5, so2 and sf6 keep a final Rydberg total within 0.90-1.05 of gennbo's,
#   no molecule gets worse by more than 0.01 e, and the no-valence-block excess falls.  Improving
#   the mean while flattening those three is a fudge factor, not a fix.
#
# Only the native side is re-run.  gennbo 7's JSONs are untouched and come in by symlink from
# nbo_ref_v2, the single stamped truth - NOT from nbo_ref_v2_fixed, which nao_split_arms.sh still
# points at and which config.load_nbo now refuses on purpose.
set -u
OUT=/work/akkleemiss/florian/nao_step1_exp
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
export NAO_DUMP_STEP3=1

echo "=== nao_step1_arms on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS threads"
echo "=== binary $BIN"
echo "=== commit $(cat $V2/../../.git/HEAD 2>/dev/null | head -c 60) plus the uncommitted nao.cpp knobs (git is absent on the compute nodes)"
[ -x "$BIN" ] || { echo "no binary"; exit 1; }

for ARM in baseline owso_off prenao_net; do
    for ROOT in /work/akkleemiss/florian/nbo_ref_v2 /work/akkleemiss/florian/nbo_ref_radicals; do
        for D in "$ROOT"/*/; do
            MOL=$(basename "$D")
            [ -f "$D/$MOL.gbw" ] || continue
            [ -f "$D/${MOL}_native.47" ] || continue
            W=$OUT/$ARM/$(basename "$ROOT")/$MOL
            mkdir -p "$W"
            cp "$D/$MOL.gbw" "$D/${MOL}_native.47" "$W/"
            ln -sf "$D/$MOL.gennbo.nbo.json" "$W/$MOL.gennbo.nbo.json"
            FLAGS=$(python3 "$V2/config.py" --native-flags "$MOL" 2>/dev/null)
            [ -n "$FLAGS" ] || { echo "SKIP $MOL: no config entry"; continue; }
            case $ARM in
                owso_off)   KNOB=NAO_OWSO_OFF=1 ;;
                prenao_net) KNOB=NAO_PRENAO_NET=1 ;;
                *)          KNOB=NAO_BASELINE=1 ;;
            esac
            T0=$(date +%s)
            ( cd "$W" && env $KNOB "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" \
                $FLAGS -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" \
                > "native_$MOL.log" 2>&1 )
            RC=$?
            printf '%-11s %-12s rc=%d %4ds step3_lines=%s\n' "$ARM" "$MOL" "$RC" \
                "$(( $(date +%s) - T0 ))" "$(grep -c '^STEP3 ' "$W/NoSpherA2.log" 2>/dev/null)"
        done
    done
done

cd "$V2" || exit 1
for ARM in baseline owso_off prenao_net; do
    for ROOT in nbo_ref_v2 nbo_ref_radicals; do
        echo "=== ARM $ARM / $ROOT"
        python3 step3_stages.py "$OUT/$ARM/$ROOT"
    done
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname), $THREADS threads"
