#!/bin/bash
#SBATCH --job-name=letkf_g2_m34
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=02:30:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# M3: time and allocation per forward-model stage; M4: allocation sites. One member, one thread.
set -euo pipefail
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

S_MEMBER=1 S_OUT=g2/M3 jl M3 1 optimize/stages.jl
S_MEMBER=1 A_PATHSET=warm A_RATE=0.001 A_OUT=g2/M4 jl M4 1 optimize/allocs.jl
