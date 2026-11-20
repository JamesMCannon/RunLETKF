#!/bin/bash
# G4 submission: a setup job instantiates and precompiles once, then the three replicates run as a
# job array, then one job compares all replicates. Usage: bash optimize/submit_g4.sh
set -euo pipefail
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
COMMON=(--qos=blanca-lair --exclude=bhpc-c5-u7-23 --nodes=1 --ntasks=1 --parsable)
cd "$PROJ"

if [ -n "$(squeue -h -u "$USER" -n letkf_g3b)" ] && [ "${FORCE:-0}" != 1 ]; then
    echo "letkf_g3b is still running: its later configurations read this directory when they start."
    echo "Wait until its last configuration has warmed up, or rerun with FORCE=1."
    exit 1
fi
[ -f p1/M100/p1.jld2 ] || { echo "missing reference p1/M100/p1.jld2"; exit 1; }
[ -e g4 ] && { echo "g4/ exists; move it aside first (the harness refuses to overwrite)"; exit 1; }
for f in harness.jl compare_g4.jl g4.sh; do
    [ -f optimize/$f ] || { echo "missing optimize/$f"; exit 1; }
done

SETUP=$(sbatch "${COMMON[@]}" --job-name=letkf_g4_setup --cpus-per-task=40 --mem=0 --time=01:00:00 \
    --output="$LOGDIR/%x_%j.out" <<'JOB'
#!/bin/bash
set -euo pipefail
JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
export JULIA_NUM_PRECOMPILE_TASKS=$SLURM_CPUS_PER_TASK
cd /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
echo "== node $(hostname), start $(date -Is)"
echo "RunLETKF $(git rev-parse HEAD) ($(git rev-parse --abbrev-ref HEAD))"
git status --porcelain
tr -d '\r' < Manifest.toml | awk '/^\[\[deps\.(LongwaveModePropagator|ModifiedVLFInversionAlgorithms|RootsAndPoles)\]\]/{p=$0} p && /git-tree-sha1/{print p, $0; p=""}'
"$JULIA" --project=. -e '
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
    using Distributed, ModifiedVLFInversionAlgorithms
    SAFEMESH isa PolishedMesh || error("MVIA in the Manifest has no polishing: update ModifiedVLFInversionAlgorithms")
    hasmethod(ensemble_model!, Tuple{Any,Any,NamedTuple}, (:pids,)) ||
        error("MVIA in the Manifest has no multi-process ensemble_model!")
    println("setup ok")'
echo "== end $(date -Is)"
JOB
)
ARRAY=$(sbatch --parsable --array=1-3 --dependency=afterok:$SETUP --kill-on-invalid-dep=yes \
    --output="$LOGDIR/%x_%A_%a.out" optimize/g4.sh)
ANALYSIS=$(sbatch "${COMMON[@]}" --job-name=letkf_g4_analysis --cpus-per-task=1 --mem=16G --time=00:30:00 \
    --dependency=afterok:$SETUP,afterany:$ARRAY --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%j.out" <<'JOB'
#!/bin/bash
set -euo pipefail
JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
cd /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
"$JULIA" --project=. optimize/compare_g4.jl p1/M100 g4/rep1 g4/rep2 g4/rep3 | tee g4/analysis.txt
JOB
)
echo "setup $SETUP -> replicates $ARRAY (array 1-3) -> analysis $ANALYSIS"
echo "logs: $LOGDIR/letkf_g4_setup_$SETUP.out, letkf_g4_${ARRAY}_{1,2,3}.out, letkf_g4_analysis_$ANALYSIS.out"
