#!/bin/bash
#SBATCH --job-name=nao-split-arms
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --output=/work/akkleemiss/florian/nao_split_exp/slurm-%j.out
#
# The two arms of ONE experiment on the SAME 22 wavefunctions:
#
#   baseline   step 4 re-diagonalises valence and Rydberg together (what the code ships)
#   split      NAO_CLASS_SPLIT=1, each class re-diagonalised on its own - the variant the
#              comment in nao.cpp says was rejected for "inflating the Rydberg occupancies
#              tenfold"
#
# Pre-registered expectation, written before the first run: SPLIT WILL BE WORSE.  The top
# eigenvalue of the m-averaged (atom, l) block is an upper bound on any single shell's own
# diagonal, so the valence shell can only come out with MORE population from one block than
# from a split one; splitting cannot reduce a Rydberg excess, it has to increase it.  If split
# comes out BETTER than the current +96.2 % mean Rydberg inflation, the one-block form is not
# the fix it is documented as and this whole line of reasoning is wrong.
#
# Only the native side is re-run.  gennbo 7's JSONs are untouched - they are the independent
# side and re-running them would destroy the only thing that makes this a gate.
set -u
ROOT=/work/akkleemiss/florian/nbo_ref_v2
OUT=/work/akkleemiss/florian/nao_split_exp
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
HERE=$ROOT/bin
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS

echo "=== nao_split_arms on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS threads"
echo "=== binary $BIN"
echo "=== $(cd /work/akkleemiss/florian/nos_nboref && git log --oneline -1) plus the uncommitted nao.cpp diagnostic knobs"
[ -x "$BIN" ] || { echo "no binary"; exit 1; }

mkdir -p "$OUT"
for ARM in baseline split; do
    for D in "$ROOT"/*/; do
        MOL=$(basename "$D")
        [ -f "$D/$MOL.gbw" ] || continue
        [ -f "$D/${MOL}_native.47" ] || continue
        W=$OUT/$ARM/$MOL
        mkdir -p "$W"
        cp "$D/$MOL.gbw" "$D/${MOL}_native.47" "$W/"
        FLAGS=$(python3 "$HERE/config.py" --native-flags "$MOL")
        T0=$(date +%s)
        if [ "$ARM" = split ]; then
            ( cd "$W" && NAO_CLASS_SPLIT=1 "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" \
                $FLAGS -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" \
                > "native_$MOL.log" 2>&1 )
        else
            ( cd "$W" && "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "${MOL}_native.47" \
                $FLAGS -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" \
                > "native_$MOL.log" 2>&1 )
        fi
        RC=$?
        T1=$(date +%s)
        # the gennbo side comes in by symlink so nao_class_leak.py runs on this tree unchanged.
        # It used to point into nbo_ref_v2_fixed, which was the correct parse while $ROOT held the
        # old one.  $ROOT is the single stamped truth now and the mirror is retired unstamped, so a
        # rerun would have stopped in config.load_nbo - the guard doing its job, but the link may as
        # well be right.
        ln -sf "$D/$MOL.gennbo.nbo.json" "$W/$MOL.gennbo.nbo.json"
        printf '%-8s %-12s rc=%d %4ds\n' "$ARM" "$MOL" "$RC" "$((T1 - T0))"
    done
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname), $THREADS threads"
