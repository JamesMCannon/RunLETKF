#!/bin/bash
#SBATCH --job-name=JLBatchTest
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=120:00:00
#SBATCH --mem=0
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Julia_tests/%x_%j.out
#SBATCH --error=/scratch/alpine/jaca8590/Julia/Julia_tests/%x_%j.err

# change if CURC uses a different module name

# Ensure JULIA threads equals cpus-per-task
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

export JULIA_DEPOT_PATH="/projects/jaca8590/Julia/Julia_depot/"

echo "JULIA_DEPOT_PATH in this job: $JULIA_DEPOT_PATH"

echo "Job $SLURM_JOB_ID starting on $(date)"
echo "Node list: $SLURM_NODELIST"
echo "Tasks total: $SLURM_NTASKS; tasks/node: $SLURM_NTASKS_PER_NODE; cpus/task: $SLURM_CPUS_PER_TASK"
echo "JULIA_NUM_THREADS=$JULIA_NUM_THREADS"

JULIA="/projects/jaca8590/software/julia-1.10.10/bin/julia"

# Run Julia script (one Julia process per node; addprocs(SlurmManager()) inside script will spawn workers)
#srun --export=ALL --mpi=none $JULIA /projects/jaca8590/Julia/HPC_Tests/batch_submission_tests/DCDL_BleedSims_v1.1.jl
srun --export=ALL --mpi=none $JULIA /projects/jaca8590/Julia/HPC_Tests/Julia_tests/Main.jl

echo "Job $SLURM_JOB_ID finished at $(date)"
