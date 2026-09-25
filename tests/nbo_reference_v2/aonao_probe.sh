#!/bin/bash
#SBATCH --job-name=aonao-probe
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --output=/work/akkleemiss/florian/aonao_probe/slurm-%j.out
#
# Does NBO 7 on THIS install emit the AO -> NAO transformation, and with what spelling?
#
# Why it matters: every comparison so far has been of FINAL tables (NPA charges, NAO
# occupancies, class totals), and a final table cannot say which orbital has the wrong shape.
# An arbitrated AO -> NAO matrix can: paired by rank within (atom, l), a column-by-column
# comparison names the block.  Without it the remaining freedom in nao.cpp - which vectors
# step 3 hands to each (atom, l) block - has no external arbiter at all.
#
# `strings nbo7.i4.exe` shows the keyword exists in five modes (print / print N columns /
# write to lfn / read from lfn / checkpoint), so the keyword AONAO is CONFIRMED to exist.
# What is NOT confirmed is the spelling of the write form or its default unit number, so the
# arms below are hypotheses, one per candidate, and HELP is asked first because this version's
# own help text is the manual for this version.
#
# LiF, 45 basis functions, 2 atoms: the printed matrix is small enough to read by eye, and
# small enough that a wrong keyword fails in seconds instead of minutes.
set -u
OUT=/work/akkleemiss/florian/aonao_probe
SRC=/work/akkleemiss/florian/nbo_ref_v2/lif
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
mkdir -p "$OUT"
cd "$OUT"

echo "host=$(hostname)  cpus_given=${SLURM_CPUS_ON_NODE:-none}  node_cores=$(grep -c ^processor /proc/cpuinfo)  job=${SLURM_JOB_ID:-none}"
echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"

run_arm() {
    # run_arm <name> <keyword ...>
    local name=$1; shift
    local W="$OUT/$name"
    rm -rf "$W"; mkdir -p "$W"
    cp "$SRC/lif.47" "$W/lif.47"
    if [ $# -gt 0 ]; then
        python3 "$V2/set_nbo_keylist.py" "$W/lif.47" "$@" || return 1
    fi
    ( cd "$W" && ls -A > .before && \
      NBOMEM=4gb bash "$V2/gennbo7" lif; echo "rc=$?" ) 2>&1 | sed "s/^/[$name] /"
    echo "[$name] keylist: $(sed -n 2p "$W/lif.47")"
    echo "[$name] files created besides lif.47/.nbo/.err:"
    ( cd "$W" && ls -A | grep -v -x -e lif.47 -e lif.nbo -e lif.gennbo.err -e .before | sed "s/^/[$name]   /" )
    echo "[$name] .nbo size $(stat -c %s "$W/lif.nbo") bytes; lines matching a NAO matrix header:"
    grep -n -i 'AO to NAO\|NAO transformation\|AONAO' "$W/lif.nbo" | head -10 | sed "s/^/[$name]   /"
    echo "[$name] any error/unknown-keyword line:"
    grep -n -i 'unknown\|illegal\|invalid\|error\|not recognized' "$W/lif.nbo" | head -10 | sed "s/^/[$name]   /"
    echo
}

# Arm 0: this version's own help text.  It is the only manual available on the node and it
# states the exact syntax of every matrix keyword, so it is asked before anything is built.
run_arm help HELP

# Arm 1: the plain form.  The help string "(1x,'      /AONAO  / : Print the AO to NAO
# transformation')" says this prints it into the .nbo output, which is text I already parse.
run_arm print AONAO

# Arms 2 and 3: the write form, two spellings.  =W with no number tests whether the default
# unit is used; =W48 tests an explicit logical file number.
run_arm w AONAO=W
run_arm w48 AONAO=W48

echo "probe done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
