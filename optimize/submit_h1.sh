#!/bin/bash
# H1 submission: setup (instantiate, precompile for JULIA_CPU_TARGET, check MVIA and LMP) -> minimum
# pathset, 1x40 against 4x10 over 3 iterations, which must match exactly -> full AVID, one iteration,
# 100 and 50 members, 1x40 and 4x10 -> comparison report with the G5 runs. Results in h1/.
# Usage: bash optimize/submit_h1.sh
set -euo pipefail
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
COMMON=(--qos=blanca-lair --exclude=bhpc-c5-u7-23 --nodes=1 --ntasks=1 --parsable)
cd "$PROJ"

[ -e h1 ] && { echo "h1/ exists; move it aside first (runs refuse to overwrite or resume)"; exit 1; }
for f in g5.sh g5run.jl compare_g5.jl jobcommon.sh; do
    [ -f optimize/$f ] || { echo "missing optimize/$f"; exit 1; }
done
grep -q 'G5_OUT' optimize/g5.sh || { echo "optimize/g5.sh does not take G5_OUT"; exit 1; }
for d in g5/avid100_1x40 g5/avid100_4x10 g5/avid50_1x40 g5/avid50_4x10; do
    [ -f $d/g5run.jld2 ] || { echo "missing $d/g5run.jld2 (G5 timings for the report)"; exit 1; }
done

ENVSETUP='JULIA=/projects/jaca8590/software/julia-1.10.10/bin/julia
export JULIA_DEPOT_PATH=/projects/jaca8590/Julia/Julia_depot/
export JULIA_CPU_TARGET="generic;skylake-avx512,clone_all;znver2,clone_all"
cd /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF'

SETUP=$(sbatch "${COMMON[@]}" --job-name=letkf_h1_setup --cpus-per-task=40 --mem=0 --time=02:00:00 \
    --output="$LOGDIR/%x_%j.out" <<JOB
#!/bin/bash
set -euo pipefail
$ENVSETUP
export JULIA_NUM_PRECOMPILE_TASKS=\$SLURM_CPUS_PER_TASK
echo "== node \$(hostname), start \$(date -Is), JULIA_CPU_TARGET=\$JULIA_CPU_TARGET"
echo "RunLETKF \$(git rev-parse HEAD) (\$(git rev-parse --abbrev-ref HEAD))"
git status --porcelain
tr -d '\r' < Manifest.toml | awk '/^\[\[deps\.(LongwaveModePropagator|ModifiedVLFInversionAlgorithms|RootsAndPoles)\]\]/{p=\$0} p && /git-tree-sha1/{print p, \$0; p=""}'
"\$JULIA" --project=. -e '
    using Pkg
    Pkg.instantiate()
    Pkg.precompile()
    using Distributed, ModifiedVLFInversionAlgorithms
    const MVIA = ModifiedVLFInversionAlgorithms
    AVIDMESH isa PolishedMesh{TiledMesh} || error("MVIA in the Manifest has no F1r mesh")
    MVIA.dmodal_finite() || error("LMP in the Manifest has the old dmodalequation")
    all(s -> isdefined(MVIA, s), (:segmentmids, :bearing)) || error("MVIA in the Manifest lacks midpoint sampling or local bearings")
    hasmethod(ensemble_model!, Tuple{Any,Any,NamedTuple}, (:pids,)) ||
        error("MVIA in the Manifest has no multi-process ensemble_model!")
    println("setup ok")'
echo "== end \$(date -Is)"
JOB
)
MIN=$(sbatch --parsable --array=1-2 --job-name=letkf_h1 --time=06:00:00 \
    --export=ALL,G5_OUT=h1,G5_CFGS="min100_1x40 min100_4x10" \
    --dependency=afterok:$SETUP --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%A_%a.out" optimize/g5.sh)
CHECK=$(sbatch "${COMMON[@]}" --job-name=letkf_h1_check --cpus-per-task=1 --mem=16G --time=00:30:00 \
    --dependency=afterok:$MIN --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%j.out" <<JOB
#!/bin/bash
set -euo pipefail
$ENVSETUP
"\$JULIA" --project=. optimize/compare_g5.jl h1/min100_4x10:h1/min100_1x40 | tee h1/check.txt
JOB
)
AVID=$(sbatch --parsable --array=1-4 --job-name=letkf_h1 --time=24:00:00 \
    --export=ALL,G5_OUT=h1,G5_CFGS="avid100_1x40 avid100_4x10 avid50_1x40 avid50_4x10" \
    --dependency=afterok:$CHECK --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%A_%a.out" optimize/g5.sh)
ANALYSIS=$(sbatch "${COMMON[@]}" --job-name=letkf_h1_analysis --cpus-per-task=1 --mem=32G --time=00:30:00 \
    --dependency=afterok:$CHECK,afterany:$AVID --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%j.out" <<JOB
#!/bin/bash
set -euo pipefail
$ENVSETUP
"\$JULIA" --project=. optimize/compare_g5.jl h1/avid100_4x10:h1/avid100_1x40 h1/avid50_4x10:h1/avid50_1x40 \
    h1/min100_4x10:h1/min100_1x40 g5/avid100_4x10:g5/avid100_1x40 g5/avid50_4x10:g5/avid50_1x40 | tee h1/analysis.txt
JOB
)
echo "setup $SETUP -> minimum pathset $MIN (array 1-2) -> check $CHECK -> full AVID $AVID (array 1-4) -> analysis $ANALYSIS"
echo "logs in $LOGDIR: letkf_h1_setup_$SETUP, letkf_h1_${MIN}_{1,2}, letkf_h1_check_$CHECK, letkf_h1_${AVID}_{1..4}, letkf_h1_analysis_$ANALYSIS"