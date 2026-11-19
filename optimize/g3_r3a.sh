#!/bin/bash
#SBATCH --job-name=letkf_g3_r3a
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# R3 screening, part A: process x thread layouts for members 1-40, concurrent processes, one call each.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

layout g3/r3a/base 1 40 40 240 || echo "layout failed"
layout g3/r3a/L2x20 2 20 20 300 || echo "layout failed"
layout g3/r3a/L4x10 4 10 10 300 || echo "layout failed"
H_NUMTASKS=1 layout g3/r3a/L4x10_nt1 4 10 10 300 || echo "layout failed"

echo "== compare"
"$JULIA" --project="$PROJ" optimize/compare_g3.jl g1/A40 g3/r3a/base g3/r3a/L2x20 g3/r3a/L4x10 g3/r3a/L4x10_nt1
