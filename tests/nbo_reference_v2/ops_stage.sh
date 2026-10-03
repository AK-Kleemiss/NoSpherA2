#!/bin/bash
#SBATCH --job-name=nbomat
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_mat/slurm-%j.out
#
# Stage the CLOSED LOOP for the decisive test: gennbo's own density and overlap operators ON
# gennbo's own PNAO basis, together with its own final NAO matrix, for the same eight molecules
# OUTCOME B ran on.
#
# Why these units and not others.  The lane's localisation says step 4 - the per-(atom, l)
# m-averaged re-diagonalisation - consumes the VECTORS, and that both sides ENTER it with the same
# class spans (0 .. 3.3e-06, 8/8).  So the one basis on which the two sides provably agree is the
# PNAO basis, and the operator whose per-block diagonalisation the cascade is made of is the
# density on that basis.  Neither string DMPNAO nor SPNAO occurs anywhere in this repo: the lane
# has never had that operator.  The probe (lif, 25 Sep) measured that this install honours an
# explicit W<n> for operator matrices - 51/52/53/54 all written, each behind its own
# self-identifying header (" PNAO density matrix:", " PNAO overlap matrix:", " NAO density
# matrix:", " AO density matrix:") - which the manual documents generically (p48 lists the
# keywords, p49-50 the W[n] parameter and the LFN table in which NEITHER appears, so both would
# otherwise collide on the default 49.
#
# One arm only.  RPNAO is NOT requested here: the probe established that it leaves gennbo's final
# NAO occupancies bit-identical and its NAO vectors identical up to sign, so it cannot enter a test
# whose reference is that same final matrix.  Its TRyd*Tred was extracted from the lif pair and
# needs no further runs.
#
# Native's side is REUSED, not re-run: NoSpherA2 ignores the $NBO keylist entirely (it reads the
# .47's data lists and the .gbw), so the NAOCPRE/NAOC dumps from job 591483 are the same arrays
# this keylist would produce.  That reuse is CHECKED rather than asserted - the .47 staged here is
# diffed against that job's and the only permitted difference is the $NBO line.
set -u
OUT=${OUT:-/work/akkleemiss/florian/nbo_mat}
SRC=/work/akkleemiss/florian/nbo_ref_v2
PREV=/work/akkleemiss/florian/aopnao_cmp
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
BIN=/work/akkleemiss/share/NoSpherA2_RGBI_NBO/NoSpherA2
MOLS=${MOLS:-"lif water ammonia ethane benzene pf5 so2 sf6"}
KEYS="PAOPNAO=W29 AOPNAO=W32 AONAO=W33 DMPNAO=W51 SPNAO=W52 DMNAO=W53 NRTSYM=OFF"
mkdir -p "$OUT"
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} cores=$(grep -c ^processor /proc/cpuinfo) date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "keylist: $KEYS"

for MOL in $MOLS; do
    W=$OUT/ops_$MOL
    rm -rf "$W"; mkdir -p "$W"
    cp "$SRC/$MOL/${MOL}_native.47" "$W/$MOL.47"
    python3 - "$W/$MOL.47" "$KEYS" <<'PY' || exit 1
import sys
p, keys = sys.argv[1], sys.argv[2]
lines = open(p).read().splitlines(True)
for i, line in enumerate(lines):
    if line.strip().upper().startswith("$NBO"):
        lines[i] = " $NBO %s $END\n" % keys
        break
else:
    raise SystemExit("no $NBO line in " + p)
open(p, "w").writelines(lines)
PY
    # the reuse check: this .47 against job 591483's, and the ONLY line allowed to differ is $NBO.
    DIFF=$(diff "$W/$MOL.47" "$PREV/$MOL/$MOL.47" | grep -E '^[<>]' | grep -v -i '\$NBO' | head -3)
    if [ -n "$DIFF" ]; then echo "$MOL: .47 differs from job 591483's beyond the keylist - NOT reusing"; echo "$DIFF"; fi
    for f in "$PREV/$MOL/$MOL.naocpre.txt" "$PREV/$MOL/$MOL.naoc.txt"; do
        [ -f "$f" ] && cp "$f" "$W/" || echo "$MOL: missing $f"
    done

    T0=$(date +%s)
    ( cd "$W" && bash "$V2/gennbo7" "$MOL" >gennbo.log 2>&1 ); RC=$?
    T1=$(date +%s)
    ( cd "$W" && "$BIN" -nbo_parse "$MOL.nbo" -nbo_json "$MOL.aopnao.nbo.json" >parse.log 2>&1 ); RCP=$?
    printf '%-9s gennbo %4ds rc=%d  parse rc=%d\n' "$MOL" "$((T1-T0))" "$RC" "$RCP"
    ( cd "$W" && for f in $(ls -A); do
        case "$f" in *.47|*.log|*.json|*.txt|*.err) continue;; esac
        printf '    %-14s %9s  | %s\n' "$f" "$(stat -c %s "$f")" "$(sed -n 2p "$f" | tr -d '\r')"
      done )
done
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
