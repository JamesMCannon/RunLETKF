#!/bin/bash
#SBATCH --job-name=letkf_g1_bench
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=05:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

set -euo pipefail

JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
cd "$PROJ"

export PATHS=AVID TOD=day ENS_SIZE=100
export DO_AMP=true DO_PHASE=true DO_S2=true DO_S3=true
export H_STATE=g1/state_t1 H_PATHSET=bench4

echo "== node $(hostname)"
lscpu | grep -E "Model name|^CPU\(s\)"
free -g | head -2
echo "== code"
echo "RunLETKF $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"
echo "uncommitted changes:"
git status --porcelain
git diff --stat
tr -d '\r' < Manifest.toml | awk '/^\[\[deps\.(LongwaveModePropagator|ModifiedVLFInversionAlgorithms|RootsAndPoles)\]\]/{p=$0} p && /git-tree-sha1/{print p, $0; p=""}'

echo "== precompile, start $(date -Is)"
"$JULIA" --project="$PROJ" -e 'using Pkg; Pkg.precompile()'
echo "== precompile end $(date -Is)"

harness() {  # label, threads; H_* set by caller
    echo "== $1, $2 threads, start $(date -Is)"
    JULIA_NUM_THREADS=$2 srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" optimize/harness.jl
    echo "== $1 end $(date -Is)"
}

H_MODE=fidelity H_T=0 H_MEMBERS=8 H_CALLS=1 H_OUT=g1/fidelity harness fidelity 40
H_MODE=bench H_T=1 H_MEMBERS=40 H_CALLS=3 H_OUT=g1/A40 harness A40 40
H_MODE=bench H_T=1 H_MEMBERS=40 H_CALLS=1 H_OUT=g1/B20 harness B20 20

"$JULIA" --project="$PROJ" optimize/compare_ym.jl g1/A40 g1/B20
