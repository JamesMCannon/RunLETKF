#!/bin/bash
#SBATCH --job-name=h2export
#SBATCH --qos=preemptable
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# Writes h2/$H2_LABEL/$H2_NAME.jls (optimize/h2_export.jl) for Mesh_Testing's lmp_truth_states.jl.
# Usage: H2_STATE=g1/state_t1 H2_LABEL=g1_state_t1 [H2_MEMBERS=1] [H2_NAME=smoke] sbatch optimize/h2_export.sh
set -uo pipefail
cd /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
source optimize/jobcommon.sh
echo "== h2 export $H2_LABEL ${H2_NAME:-inputs}, members ${H2_MEMBERS:-all}, start $(date -Is)"
JULIA_NUM_THREADS=1 srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" optimize/h2_export.jl
STATUS=$?
echo "== h2 export end $(date -Is), status $STATUS"
exit $STATUS
