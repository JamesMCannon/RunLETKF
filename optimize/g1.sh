#!/bin/bash
#SBATCH --job-name=letkf_g1_state
#SBATCH --qos=blanca-lair
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=14:00:00
#SBATCH --begin=2026-11-18T19:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

set -euo pipefail

JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
cd "$PROJ"

export PATHS=AVID TOD=day ENS_SIZE=100 ITRS=1
export DO_AMP=true DO_PHASE=true DO_S2=true DO_S3=true
export NEW_FOLDER=g1/state_t1
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK

if [ -e "$PROJ/$NEW_FOLDER" ]; then
    echo "$PROJ/$NEW_FOLDER exists; refusing to overwrite or resume" >&2
    exit 1
fi

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

echo "== run, $JULIA_NUM_THREADS threads, start $(date -Is)"
srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" -e '
    t = @timed include("Main.jl")
    g = t.gcstats
    println("== total ", round(t.time; digits=1), " s; GC ", round(t.gctime; digits=1), " s (",
        round(100t.gctime/t.time; digits=1), "%); allocated ", round(t.bytes/2^30; digits=1),
        " GiB; pauses ", g.pause, ", full sweeps ", g.full_sweep, "; maxrss ",
        round(Sys.maxrss()/2^30; digits=1), " GiB")'
echo "== run end $(date -Is)"