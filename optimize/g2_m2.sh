#!/bin/bash
#SBATCH --job-name=letkf_g2_m2
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=05:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# M2: thread scaling with one member per thread (members = threads), one call each.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

for n in 1 5 10 20 40; do
    H_MEMBERS=$n H_CALLS=1 H_OUT=g2/M2_t$n jl "M2 t$n" $n optimize/harness.jl
done
