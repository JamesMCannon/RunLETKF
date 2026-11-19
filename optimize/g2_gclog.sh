#!/bin/bash
#SBATCH --job-name=letkf_g2_gclog
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# GC pause distribution: the A40 configuration (40 members, 40 threads), one call, GC logging on.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

H_MEMBERS=40 H_CALLS=1 H_GCLOG=true H_OUT=g2/gclog jl gclog 40 optimize/harness.jl

echo "== GC pauses"
"$JULIA" --project="$PROJ" optimize/gcpauses.jl "$LOGDIR/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out"
