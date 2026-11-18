#!/bin/bash
#SBATCH --job-name=letkf_g0
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

set -euo pipefail

JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
cd "$PROJ"

export PATHS=minimum TOD=day ENS_SIZE=4 ITRS=1
export DO_AMP=true DO_PHASE=true DO_S2=true DO_S3=true

echo "== node $(hostname)"
lscpu | grep -E "Model name|^CPU\(s\)"

echo "== code"
echo "RunLETKF $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"
echo "uncommitted changes:"
git status --porcelain
git diff --stat
tr -d '\r' < Manifest.toml | awk '/^\[\[deps\.(LongwaveModePropagator|ModifiedVLFInversionAlgorithms|RootsAndPoles)\]\]/{p=$0} p && /git-tree-sha1/{print p, $0; p=""}'
"$JULIA" --project="$PROJ" -e 'using Pkg; Pkg.status(["LongwaveModePropagator", "ModifiedVLFInversionAlgorithms", "RootsAndPoles"]; mode=Pkg.PKGMODE_MANIFEST)'

echo "== precompile, start $(date -Is)"
JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK "$JULIA" --project="$PROJ" -e 'using Pkg; Pkg.precompile()'
echo "== precompile end $(date -Is)"

run_letkf() {  # name, threads
    export NEW_FOLDER="g0/$1" JULIA_NUM_THREADS=$2
    rm -rf "$PROJ/$NEW_FOLDER"
    echo "== run $1, $2 threads, start $(date -Is)"
    srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" Main.jl
    echo "== run $1 end $(date -Is)"
}

run_letkf A 40
run_letkf B 40
run_letkf C 8

"$JULIA" --project="$PROJ" optimize/g0_compare.jl g0/A g0/B g0/C