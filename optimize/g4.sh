#!/bin/bash
#SBATCH --job-name=letkf_g4
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# G4 replicate <rep 1-3>: MVIA's multi-process ensemble_model! (4x10, 8x5) against threads only
# (1x40), 100 members, polished SafeMesh, one node. ym must equal P1's polished ym (p1/M100).
set -euo pipefail
REP=${1:?usage: g4.sh <replicate 1-3>}
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

case $REP in
    1) ORDER="L4x10 L1x40 L8x5" ;;
    2) ORDER="L8x5 L4x10 L1x40" ;;
    3) ORDER="L1x40 L8x5 L4x10" ;;
    *) echo "replicate must be 1-3"; exit 1 ;;
esac
echo "== replicate $REP, order: $ORDER"

export H_MEMBERS=100 H_CALLS=1
for cfg in $ORDER; do
    out=g4/rep$REP/$cfg
    case $cfg in
        L1x40) H_OUT=$out jl "rep$REP $cfg" 40 optimize/harness.jl || echo "$cfg failed" ;;
        L4x10) H_OUT=$out H_PROCS=4 H_PROC_THREADS=10 jl "rep$REP $cfg" 1 optimize/harness.jl || echo "$cfg failed" ;;
        L8x5)  H_OUT=$out H_PROCS=8 H_PROC_THREADS=5 jl "rep$REP $cfg" 1 optimize/harness.jl || echo "$cfg failed" ;;
    esac
done

echo "== compare"
"$JULIA" --project="$PROJ" optimize/compare_g4.jl p1/M100 g4/rep$REP
