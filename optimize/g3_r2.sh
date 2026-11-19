#!/bin/bash
#SBATCH --job-name=letkf_g3_r2
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=03:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# R2: --heap-size-hint; 40 members, 40 threads, one process, one call each.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh
export H_MEMBERS=40 H_CALLS=1

H_OUT=g3/r2/base jl "R2 base" 40 optimize/harness.jl
for h in 50 100 150; do
    JFLAGS="--heap-size-hint=${h}G" H_OUT=g3/r2/heap$h jl "R2 heap$h" 40 optimize/harness.jl
done

echo "== compare"
"$JULIA" --project="$PROJ" optimize/compare_g3.jl g1/A40 g3/r2/base g3/r2/heap50 g3/r2/heap100 g3/r2/heap150
