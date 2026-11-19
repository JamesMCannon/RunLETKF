#!/bin/bash
#SBATCH --job-name=letkf_g3b
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# R1/R2 replicate: all eight GC configurations, 40 members, 40 threads, one process, one call each,
# in this replicate's fixed random order. Usage: sbatch optimize/g3b_r1r2.sh <replicate 1-4>
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh
export H_MEMBERS=40 H_CALLS=1

REP=$1
declare -A FLAGS=(
    [base]=""
    [gc40]="--gcthreads=40"
    [gc20s1]="--gcthreads=20,1"
    [gc40s1]="--gcthreads=40,1"
    [h150]="--heap-size-hint=150G"
    [gc40_h150]="--gcthreads=40 --heap-size-hint=150G"
    [gc20s1_h150]="--gcthreads=20,1 --heap-size-hint=150G"
    [gc40s1_h150]="--gcthreads=40,1 --heap-size-hint=150G"
)
case $REP in
    1) ORDER="gc20s1 gc40s1 base gc40 gc20s1_h150 h150 gc40_h150 gc40s1_h150" ;;
    2) ORDER="gc40 gc40s1 h150 gc40_h150 gc20s1_h150 base gc20s1 gc40s1_h150" ;;
    3) ORDER="gc40_h150 gc20s1_h150 gc40 gc20s1 gc40s1 gc40s1_h150 base h150" ;;
    4) ORDER="h150 gc20s1 gc20s1_h150 gc40s1 gc40_h150 gc40 gc40s1_h150 base" ;;
    *) echo "replicate must be 1-4" >&2; exit 1 ;;
esac
echo "== replicate $REP, order: $ORDER"

for cfg in $ORDER; do
    JFLAGS="${FLAGS[$cfg]}" H_OUT=g3b/rep$REP/$cfg jl "rep$REP $cfg" 40 optimize/harness.jl || echo "$cfg failed"
done
