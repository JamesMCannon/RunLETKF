#!/bin/bash
#SBATCH --job-name=letkf_g3_r3b
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=04:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# R3 screening, part B: 8x5 and 40x1 layouts, and numtasks = 1 in one 40-thread process.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

layout g3/r3b/base 1 40 40 240 || echo "layout failed"
layout g3/r3b/L8x5 8 5 5 300 || echo "layout failed"
layout g3/r3b/L40x1 40 1 1 600 || echo "layout failed"
H_NUMTASKS=1 layout g3/r3b/L1x40_nt1 1 40 40 240 || echo "layout failed"

echo "== compare"
"$JULIA" --project="$PROJ" optimize/compare_g3.jl g1/A40 g3/r3b/base g3/r3b/L8x5 g3/r3b/L40x1 g3/r3b/L1x40_nt1
