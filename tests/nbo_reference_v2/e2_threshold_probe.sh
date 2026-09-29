#!/bin/bash
# Is the open-shell E2 factor of two real? Make both programs print the numbers they suppressed.
#
#   sbatch e2_threshold_probe.sh [threshold_kcal]        # default 0.02
#
# gennbo prints "None above threshold" for ch3's E2 table in both spin blocks at its open-shell
# default of 0.25 kcal/mol, while -nbo_native prints 12 entries at 0.26-0.40. An empty table
# against 12 entries is *consistent* with native doubling a per-spin quantity and proves
# nothing - 12 matched pairs at ratio 2.00 would. So this lowers the print threshold on BOTH
# sides, writes the two JSONs into a directory named like any other molecule, and leaves
# ratio_probe.py to do the pairing:
#
#   python bin/ratio_probe.py <root>/e2_probe_<threshold>/ch3_lowe2 ch3_lowe2 --expect 2.0
#
# The threshold is an argument because the answer moved with it. At 0.02 both sides print their
# beta BD->BD* rows and native prints NO alpha ones, and the NBO tables say why: native's alpha
# BD* orbitals hold 0.00002 e against gennbo's 0.00056, so the same off-diagonal Fock element
# that fills them puts the alpha E2 rows near 0.006 kcal. Run it at 0.001 to make native print
# them and check that number instead of arguing about it.
#
# NBO's own defaults carry the same factor - 0.5 kcal closed shell against 0.25 open shell -
# which is why this is a convention question and not a bug report until the pairs land. Nothing
# here changes a convention; it only measures one.
#SBATCH --job-name=nboref-e2probe
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/e2probe-%j.out

set -u
HERE=${NBOREF_BIN:-/work/akkleemiss/florian/nbo_ref_v2/bin}
source "$HERE/env.sh"
LOW=${1:-0.02}
D="$ROOT/e2_probe_$LOW/ch3_lowe2"
mkdir -p "$D" || exit 1
cd "$D" || exit 1
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
echo "=== e2 threshold probe on $(hostname), node has $(grep -c ^processor /proc/cpuinfo) cores, this job has $THREADS threads"
echo "=== both sides at E2 print threshold $LOW kcal/mol"

# NBO 7 takes the E2 print threshold as a keyword value and an unrecognised spelling is a hard
# error, so try the documented one first and fall back. Only the threshold is added: the rest of
# the keylist is whatever stage 2 already put in the archive.
OK=""
for spec in "E2PERT=$LOW" "E2PERT=$LOW NBOSUM"; do
    cp -f "$ROOT/ch3/ch3_native.47" ch3_lowe2.47 || exit 2
    python3 "$HERE/set_nbo_keylist.py" ch3_lowe2.47 "$spec" || exit 2
    "$HERE/gennbo7" ch3_lowe2 && OK="$spec" && break
    echo "--- keylist '$spec' did not run, tail of its output:"
    tail -5 ch3_lowe2.nbo 2>/dev/null
done
[ -n "$OK" ] || { echo "no accepted E2PERT spelling"; exit 3; }
echo "=== accepted keylist: $OK"
grep -n "Threshold for printing" ch3_lowe2.nbo | head -4
echo "=== 'None above threshold' occurrences now: $(grep -c 'None above threshold' ch3_lowe2.nbo)"

T0=$(date +%s)
"$NOSPHERA2" -nbo_parse ch3_lowe2.nbo -nbo_json ch3_lowe2.gennbo.nbo.json > parse.log 2>&1
RCP=$?
T1=$(date +%s)
# The native side at the SAME threshold: a shorter table is not a better one, and comparing a
# table cut at 0.25 against one cut at 0.02 would manufacture "missing" entries out of a
# printing difference.
"$NOSPHERA2" -nbo_native "$ROOT/ch3/ch3.gbw" -nbo_47 ch3_lowe2.47 -nbo_e2min $LOW \
    -nrt -nrt_e2 1.0 -nbo_threads "$THREADS" -nbo_json ch3_lowe2.native.nbo.json \
    > native.log 2>&1
RCN=$?
T2=$(date +%s)

printf '%s\n' \
  "{\"host\": \"$(hostname)\", \"threads\": $THREADS," \
  " \"node_cores_total\": $(grep -c ^processor /proc/cpuinfo)," \
  " \"slurm_job_id\": \"${SLURM_JOB_ID:-none}\"," \
  " \"nbo_keywords\": \"$OK\", \"e2_print_threshold_kcal\": $LOW," \
  " \"native_flags\": \"-nbo_e2min $LOW -nrt -nrt_e2 1.0 -nbo_threads $THREADS\"," \
  " \"nospher_a2\": \"$NOSPHERA2\"," \
  " \"seconds\": {\"parse\": $((T1 - T0)), \"native\": $((T2 - T1))}," \
  " \"exit_codes\": {\"parse\": $RCP, \"native\": $RCN}}" > provenance_nbo.json

echo "=== rc parse=$RCP native=$RCN"
python3 -c "
import json
for s in ('gennbo', 'native'):
    r = json.load(open('ch3_lowe2.%s.nbo.json' % s))
    print(s, 'e2 entries', len(r.get('e2', [])), 'thresholds', r.get('thresholds'))
"
python3 "$HERE/ratio_probe.py" . ch3_lowe2 --expect 2.0 --json ratio_ch3_lowe2.json
echo "=== done on $(hostname), threads $THREADS"
exit $((RCP + RCN))
