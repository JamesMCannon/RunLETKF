#!/bin/bash
#SBATCH --job-name=letkf_g2_m1
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=03:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# M1: load balance at production size, 100 members on 40 threads, one call.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

H_MEMBERS=100 H_CALLS=1 H_OUT=g2/M1 jl M1 40 optimize/harness.jl
