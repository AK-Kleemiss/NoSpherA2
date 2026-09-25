#!/bin/bash
# Shared environment for the NBO external-reference regeneration, sourced by both stages.
#
# TRAP: /work/akkleemiss/share/gennbo is the NBO 7 wrapper everyone points at, and on THIS
# cluster it is broken - it sets NBOBIN=/opt/nbo7/bin, which does not exist on any node, and
# fails only once it has already consumed its arguments. The working wrapper is
# /work/software/bin/NBO/gennbo (NBOBIN=/work/software/bin/NBO/bin). It also sets
# NBOMEM=100gb and, worse, does `rm -i $JOB.nbo` on a stale output, which in a batch job
# waits forever on a prompt nobody answers. gennbo7 in this directory is a bash wrapper
# around the same two binaries without either problem.
#
# TRAP: the orca module is NOT on the default MODULEPATH, on the login node or on a compute
# node. It lives in the spack lmod tree below. `module load orca` alone reports "unknown".

export SPACK_LMOD=/work/software/spack/share/spack/lmod/linux-rocky9-x86_64/Core
export NBOBIN=/work/software/bin/NBO/bin
export NOSPHERA2=/work/akkleemiss/share/NoSpherA2_RGBI_NBO/NoSpherA2
export ROOT=${ROOT:-/work/akkleemiss/florian/nbo_ref_v2}
export ORCA_MODULE=orca/avx2-6.1.1-xqm3jnz

load_modules() {
    if ! type module >/dev/null 2>&1; then
        source /usr/share/lmod/lmod/init/bash 2>/dev/null
    fi
    module use "$SPACK_LMOD"
    module load "$ORCA_MODULE"
}

# hostname, the cores this job was GIVEN and the cores the node HAS - an audit on a sibling
# branch chased a 6.5x phantom regression because a 4-thread job landed on a 96-core node and
# nothing in any log said so.
provenance_json() {
    # provenance_json <stage> <molecule> <threads> <wall_seconds> <exit_code> [extra json]
    cat <<JSON
{
  "stage": "$1",
  "molecule": "$2",
  "threads": $3,
  "wall_seconds": $4,
  "exit_code": $5,
  "hostname": "$(hostname)",
  "node_cores_total": $(grep -c ^processor /proc/cpuinfo),
  "slurm_job_id": "${SLURM_JOB_ID:-none}",
  "slurm_partition": "${SLURM_JOB_PARTITION:-none}",
  "slurm_cpus_on_node": "${SLURM_CPUS_ON_NODE:-none}",
  "finished_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"${6:+,
  $6}
}
JSON
}
