# Sourced by the optimize/*.sh job scripts: paths, environment, provenance, precompile.
JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
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

jl() {  # label, threads, script; H_*/S_*/A_* set by caller
    echo "== $1, $2 threads, start $(date -Is)"
    JULIA_NUM_THREADS=$2 srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" "$3"
    echo "== $1 end $(date -Is)"
}
