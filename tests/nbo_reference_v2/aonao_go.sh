#!/bin/bash
#SBATCH --job-name=aonao-go
#SBATCH --partition=Active
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:40:00
#SBATCH --output=/work/akkleemiss/florian/aonao_cmp/go-%j.out
#
# Build + stage in ONE small allocation, because the two-job chain (588263/588264, 32 cpus /
# 4 h and 8 cpus / 2 h) was given an estimated start of 15 Oct: partition Active is saturated
# with HWerner's S_QM_CC arrays and my own 43 ELI_complex jobs.  The 1-second probe job 588007
# got through at 2 cpus, so the constraint is the size of the ask, not access - a 4-cpu,
# 40-minute job is backfill-sized.
#
# The build is INCREMENTAL on purpose: /work/.../nos_nboref/build/release-linux/bin/NoSpherA2
# is from 06:55 and the only source touched since is Src/core/nao.cpp (NAO_DUMP_C), so this is
# one translation unit plus a link, not the 32-cpu full build the separate job was sized for.
# If that assumption is wrong the build simply takes longer and the 40 minutes is what fails,
# which is a visible failure rather than a silent one.
set -u
V2=/work/akkleemiss/florian/nos_nboref/tests/nbo_reference_v2
echo "=== aonao_go on $(hostname), job ${SLURM_JOB_ID:-none}, cpus ${SLURM_CPUS_PER_TASK:-?}, node cores $(grep -c ^processor /proc/cpuinfo)"
bash "$V2/build_nboref.sh" || { echo "build failed - not staging"; exit 1; }
bash "$V2/aonao_stage.sh"
