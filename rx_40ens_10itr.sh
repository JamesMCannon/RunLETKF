#!/bin/bash
#SBATCH --job-name=rx_40e
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --nodelist=bhpc-c5-u7-22,bhpc-c5-u7-23
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --output=/scratch/alpine/jaca8590/Julia/State_Vec_Updates/rx_phase/%x_%j.out
#SBATCH --error=/scratch/alpine/jaca8590/Julia/State_Vec_Updates/rx_phase/%x_%j.err

# change if CURC uses a different module name

# Ensure JULIA threads equals cpus-per-task
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

export JULIA_DEPOT_PATH="/projects/jaca8590/Julia/Julia_depot/"

export ENS_SIZE=40
export ITRS=10

echo "JULIA_DEPOT_PATH in this job: $JULIA_DEPOT_PATH"

echo "Job $SLURM_JOB_ID starting on $(date)"
echo "Node list: $SLURM_NODELIST"
echo "Tasks total: $SLURM_NTASKS; tasks/node: $SLURM_NTASKS_PER_NODE; cpus/task: $SLURM_CPUS_PER_TASK"
echo "JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo "Running $ENS_SIZE member ensemble for $ITRS itrs"

JULIA="/projects/jaca8590/software/julia-1.10.10/bin/julia"

srun --export=ALL --mpi=none $JULIA /projects/jaca8590/Julia/State_Vector_Updates/RX_phase_testing/Main.jl

echo "Job $SLURM_JOB_ID finished at $(date)"
