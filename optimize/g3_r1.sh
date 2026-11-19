#!/bin/bash
#SBATCH --job-name=letkf_g3_r1
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=04:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# R1: GC mark threads (default 20) and concurrent sweeping; 40 members, 40 threads, one process, one call each.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh
export H_MEMBERS=40 H_CALLS=1

H_OUT=g3/r1/base jl "R1 base" 40 optimize/harness.jl
for g in 10 30 40; do
    JFLAGS="--gcthreads=$g" H_OUT=g3/r1/gc$g jl "R1 gc$g" 40 optimize/harness.jl
done
JFLAGS="--gcthreads=20,1" H_OUT=g3/r1/gc20s1 jl "R1 gc20s1" 40 optimize/harness.jl

echo "== compare"
"$JULIA" --project="$PROJ" optimize/compare_g3.jl g1/A40 g3/r1/base g3/r1/gc10 g3/r1/gc30 g3/r1/gc40 g3/r1/gc20s1
