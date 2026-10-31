#!/bin/bash
# Stage 1: re-optimise one reference molecule at the level of theory index.json recorded.
#
#   sbatch orca_opt.sh <molecule>
#
# Writes <name>.out, the .gbw next to it, and provenance_orca.json carrying the hostname,
# the cores this job was given and the core count of the node it landed on.
#SBATCH --job-name=nboref-orca
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=/work/akkleemiss/florian/nbo_ref_v2/slurm/orca-%j.out

set -u
MOL=$1
# NOT BASH_SOURCE: slurm copies the batch script into /var/spool/slurmd/job<id>/, so the
# script's own directory at run time is not the one it was submitted from.
HERE=${NBOREF_BIN:-/work/akkleemiss/florian/nbo_ref_v2/bin}
source "$HERE/env.sh"
cd "$ROOT/$MOL" || exit 1

load_modules
NTASKS=${SLURM_NTASKS:-1}

START=$(date +%s)
"$ORCA" "$MOL.inp" > "$MOL.out" 2>&1
RC=$?
END=$(date +%s)

provenance_json orca_opt "$MOL" "$NTASKS" "$((END - START))" "$RC" \
  "\"orca_path\": \"$ORCA\",
  \"orca_module\": \"$ORCA_MODULE\",
  \"orca_mpi_procs\": $NTASKS" > provenance_orca.json

echo "ORCA $MOL rc=$RC wall=$((END - START))s host=$(hostname) mpi_procs=$NTASKS node_cores=$(grep -c ^processor /proc/cpuinfo)"
exit $RC
