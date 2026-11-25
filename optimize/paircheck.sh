#!/bin/bash
#SBATCH --job-name=pc
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%A_%a.out

# One pair-check array task (optimize/paircheck.jl). Requeued tasks skip finished units.
# Set by optimize/submit_paircheck.sh: PC_STATE, PC_LABEL, PC_NTASKS.
set -uo pipefail
cd /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
source optimize/jobcommon.sh
echo "== pc $PC_LABEL task $SLURM_ARRAY_TASK_ID/$PC_NTASKS, restart ${SLURM_RESTART_COUNT:-0}, start $(date -Is)"
JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" optimize/paircheck.jl
STATUS=$?
echo "== pc end $(date -Is), status $STATUS"
exit $STATUS
