#!/bin/bash
#SBATCH --job-name=letkf_p1
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# P1: effect of Newton polishing (AVIDModes finalize) on SafeMesh roots. 100 members on 40 threads,
# compared with g2/M1 (same members, unpolished); then one member on one thread for cost.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

P_MEMBERS=100 P_OUT=p1/M100 jl "P1 100 members" 40 optimize/p1.jl
P_MEMBERS=1 P_OUT=p1/cost jl "P1 cost" 1 optimize/p1.jl

echo "== P1 report"
"$JULIA" --project="$PROJ" optimize/p1_report.jl p1/M100 g2/M1 p1/cost
