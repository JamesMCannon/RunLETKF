#!/bin/bash
#SBATCH --job-name=xyamp_11rho_100e
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=7-00:00:00
#SBATCH --mem=0
#SBATCH --output=/scratch/alpine/jaca8590/Julia/State_Vec_Updates/tx_pwrs/%x_%j.out
#SBATCH --error=/scratch/alpine/jaca8590/Julia/State_Vec_Updates/tx_pwrs/%x_%j.err

# change if CURC uses a different module name

# Ensure JULIA threads equals cpus-per-task
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

export JULIA_DEPOT_PATH="/projects/jaca8590/Julia/Julia_depot/"

export ENS_SIZE=100
export ITRS=10
export RHO=1.1
export SHUFFLE_XY=false
export DO_XY=true
export DO_TX=false
export DO_PHASE=false
export DO_AMP=true

echo "JULIA_DEPOT_PATH in this job: $JULIA_DEPOT_PATH"

echo "Job $SLURM_JOB_ID starting on $(date)"
echo "Node list: $SLURM_NODELIST"
echo "Tasks total: $SLURM_NTASKS; tasks/node: $SLURM_NTASKS_PER_NODE; cpus/task: $SLURM_CPUS_PER_TASK"
echo "JULIA_NUM_THREADS=$JULIA_NUM_THREADS"
echo "Running $ENS_SIZE member ensemble for $ITRS itrs with initial TX std_dev of $STDDEV_TX_KW kW. Shuffling TX: $SHUFFLE_TX"

JULIA="/projects/jaca8590/software/julia-1.10.10/bin/julia"

srun --export=ALL --mpi=none $JULIA /projects/jaca8590/Julia/State_Vector_Updates/TX_pwrs_testing/Main.jl

echo "Job $SLURM_JOB_ID finished at $(date)"
