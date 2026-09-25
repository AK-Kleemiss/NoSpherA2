#!/bin/bash
#SBATCH --job-name=aopnao
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00
#SBATCH --output=/work/akkleemiss/florian/aopnao_cmp/slurm-%j.out
#
# BISECT the NAO chain: the AO -> pre-NAO matrix from both sides, on the same eight .47 files.
#
# Why a bisection and not a hypothesis.  The AO -> NAO comparison (jobs 589060/589634) licensed
# exactly one clause: the operator handed to the final m-averaged within-block diagonalisation
# already carries the wrong occupancy multiset.  That is a LOCATION.  Turning it into a named
# step is what went wrong once already - step 3's Schmidt projection was named from evidence that
# only bracketed it, and the claim died within the hour.  A pre-NAO is upstream of every
# orthogonalisation, so these two matrices bracket the whole cascade and the interval halves
# whichever way the answer falls.
#
# PRE-REGISTERED, before the numbers exist - what each outcome licenses, and what it does not:
#
#   A. pre-NAOs DISAGREE (spectra differ beyond the floor, or the within-block mixing is not a
#      rank-paired near-identity).  Licenses: the two sides differ already at step 1's atomic
#      (atom, l) eigenproblem, upstream of every orthogonalisation; therefore no downstream
#      arbitration - OWSO on/off, the Schmidt priority chain, the class partition - is
#      interpretable until that is fixed, and a candidate fix is screened at step 1.
#      Does NOT license naming WHICH ingredient of step 1 is wrong (the m-averaging, gross vs net,
#      the S^-1/2 metric, the shell grouping).  That is a further bisection with its own arms.
#
#   B. pre-NAOs AGREE (both gates pass and the mixing is a rank-paired near-identity).  Licenses:
#      the disagreement is created strictly between PNAO and the m-averaged NAO block, i.e. inside
#      the orthogonalisation cascade and the class input it takes - steps 2 to 4 - and nothing
#      earlier.  Does NOT by itself separate Schmidt from OWSO; but the interval is then small
#      enough that NAO_OWSO_OFF and the class partition are two ARMS over a bracketed interval
#      rather than a guess, and a named step becomes reportable.
#
#   C. NO resolution passes the two-sided self-check.  Then nothing is quoted.  A metric that
#      fails its own side's self-check is not admissible for having been announced in advance -
#      pre_occ summing to 1041.75 e against N = 664 is what that costs - so the result is "the
#      PNAO matrices on disk cannot arbitrate", plus which property failed on which side.
#
# The self-check itself, carried over from the AONAO comparison, and it has to be different here
# because NBO 7 prints no pre-NAO occupancy table for either side to be checked against.  What
# both sides DO have is the defining property of a pre-NAO, which each can be tested against
# alone: inside one (atom, l) block, m-averaged, the columns are S-orthonormal and diagonalise
# S P S.  Native by construction, gennbo by NBO's own definition.  Both pass -> the comparison is
# readable at that resolution; either fails -> outcome C.
#
# Keylist: AOPNAO=W AONAO=W NRTSYM=OFF, the SAME line for all eight.  Three deliberate choices.
#   - AONAO in the same run, so the NAO matrix comes from the identical input as the PNAO one and
#     the two ends of the bracket cannot differ by their keylist.  It also re-measures the AONAO
#     numbers under NRTSYM=OFF, which is a regression check on the keylist change for free.
#   - NRTSYM=OFF everywhere, not only on sf6.  Job 589969 showed it does not reach this stage, and
#     sf6 needed it (0-byte lfn 33 with rc=0 after "SYMOPS ... expected 24"), so the old runs left
#     one set member with a different input.  Same input for every member, or the set is not a set.
#   - NRT is NOT requested: the NAO stage is upstream of NBO/E2/NRT, so it cannot depend on them,
#     and benzene's NRT is minutes.  That is an assumption, so it is MEASURED: this run's .nbo is
#     parsed and its NAO occupancy table must equal the stamped reference JSON's.
#
# MEASURED, job 591483 on AKL010, all eight rc=0: the PNAO matrix comes back on unit 32 (header
# " PNAOs in the AO basis:"), 31 kB to 766 kB, and sf6's is non-empty this time - NRTSYM=OFF in the
# shared keylist is what fixed the 0-byte lfn 33 it used to produce.  OUTCOME B: the pre-NAOs AGREE,
# eigenspace-projected mixing 1.000000 over 379 shells, 0 rank mismatches, both sides' self-check at
# 5e-10 / 5.8e-09.  The same run's NAO table still matches the stamped reference at 0.0e+00, so the
# uniform keylist changed nothing downstream either.  One oddity recorded and not chased: so2 also
# wrote a 6.4 MB unit 48, which nothing here reads.
#
# The lfn number of the PNAO matrix is NOT assumed.  AONAO=W came back on unit 33 and AONAO=W48
# was refused because 48 is reserved, so the default for AOPNAO is whatever this install says it
# is: every file the run creates is listed with its size, and the comparison reads the number off
# that list.  Guessing a unit is how a 0-byte file with rc=0 gets mistaken for an answer.
set -u
OUT=${OUT:-/work/akkleemiss/florian/aopnao_cmp}
SRC=/work/akkleemiss/florian/nbo_ref_v2
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
BIN=/work/akkleemiss/florian/nos_nboref/build/release-linux/bin/NoSpherA2
MOLS=${MOLS:-"lif water ammonia ethane benzene pf5 so2 sf6"}
THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$THREADS
mkdir -p "$OUT"

echo "host=$(hostname)  threads=$THREADS  node_cores=$(grep -c ^processor /proc/cpuinfo)  job=${SLURM_JOB_ID:-none}  partition=${SLURM_JOB_PARTITION:-none}"
echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# The build is incremental on purpose: only Src/core/nao.cpp changed (NAO_DUMP_CPRE), so this is
# one translation unit plus a link, which is backfill-sized.  A failed build must not stage
# anything, or the old binary's dump would be reported as the new one's.
# ... and it is built on the LOGIN node before this job is submitted, because the compute nodes
# carry no cmake, git or c++.  SKIP_BUILD defaults to 1 here for that reason; what makes the
# omission safe is the binary check below, which refuses to stage a binary without the new dump.
if [ "${SKIP_BUILD:-1}" != 1 ]; then
    bash "$V2/build_nboref.sh" || { echo "build failed - not staging"; exit 1; }
fi
[ -x "$BIN" ] || { echo "no binary at $BIN"; exit 1; }
echo "binary $BIN  mtime $(stat -c %y "$BIN")"
strings "$BIN" | grep -q NAO_DUMP_CPRE || { echo "binary does not carry NAO_DUMP_CPRE - stale build, not staging"; exit 1; }

for MOL in $MOLS; do
    D=$SRC/$MOL
    W=$OUT/$MOL
    [ -f "$D/${MOL}_native.47" ] || { echo "SKIP $MOL: no ${MOL}_native.47"; continue; }
    rm -rf "$W"; mkdir -p "$W"
    cp "$D/${MOL}_native.47" "$W/$MOL.47"
    cp "$D/$MOL.gbw" "$W/" 2>/dev/null

    python3 - "$W/$MOL.47" <<'PY' || exit 1
import sys
p = sys.argv[1]
lines = open(p).read().splitlines(True)
for i, line in enumerate(lines):
    if line.strip().upper().startswith("$NBO"):
        lines[i] = " $NBO AOPNAO=W AONAO=W NRTSYM=OFF $END\n"
        break
else:
    raise SystemExit("no $NBO line in " + p)
open(p, "w").writelines(lines)
print("keylist now: " + lines[i].strip())
PY

    T0=$(date +%s)
    ( cd "$W" && bash "$V2/gennbo7" "$MOL" ) || echo "gennbo rc nonzero for $MOL"
    T1=$(date +%s)
    ( cd "$W" && "$BIN" -nbo_parse "$MOL.nbo" -nbo_json "$MOL.aopnao.nbo.json" > "parse_$MOL.log" 2>&1 )
    RCP=$?

    T2=$(date +%s)
    ( cd "$W" && NAO_DUMP_CPRE=1 NAO_DUMP_C=1 "$BIN" -nbo_native "$MOL.gbw" -nbo_47 "$MOL.47" \
        -nbo_threads "$THREADS" -nbo_json "$MOL.native.nbo.json" > "native_$MOL.log" 2>&1 )
    RCN=$?
    T3=$(date +%s)
    ( cd "$W" && grep '^NAOCPRE ' NoSpherA2.log > "$MOL.naocpre.txt" || true )
    ( cd "$W" && grep '^NAOC ' NoSpherA2.log > "$MOL.naoc.txt" || true )

    printf '%-9s gennbo %3ds parse rc=%d native %3ds rc=%d  NAOCPRE=%s NAOC=%s lines\n' \
        "$MOL" "$((T1 - T0))" "$RCP" "$((T3 - T2))" "$RCN" \
        "$(wc -l < "$W/$MOL.naocpre.txt" 2>/dev/null || echo 0)" \
        "$(wc -l < "$W/$MOL.naoc.txt" 2>/dev/null || echo 0)"
    # every file gennbo wrote, with its size: this is where the PNAO matrix's lfn number comes
    # from, and a 0-byte file here is the sf6 failure mode showing itself instead of hiding.
    echo "  files: $(cd "$W" && ls -A | grep -v -E "^($MOL\.(47|gbw|nbo)|NoSpherA2\.log|.*\.txt|.*\.json|.*\.log)$" \
        | while read -r f; do printf '%s(%s) ' "$f" "$(stat -c %s "$W/$f")"; done)"
done

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname), $THREADS threads"
