#!/bin/bash
# Pair check (optimize/PairCheck.jl) of each LETKF state in PC_STATES on preemptable nodes. Per state:
# a load test on a preemptable node -> an array of restartable tasks, each running its share of the
# (member, path) units -> a summary in pc/<label>/. Rerunning resumes: finished units are skipped.
# Usage: bash optimize/submit_paircheck.sh            (both states, full arrays)
#        PC_PILOT=1 bash optimize/submit_paircheck.sh (task 1 of each array only, to time the units)
set -euo pipefail
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
STATES=${PC_STATES:-"g1/state_t1 h1/avid100_4x10"}
NTASKS=${PC_NTASKS:-60}
EXPECTED=${PC_EXPECTED:-2200}
PRE=(--qos=preemptable --exclude=bhpc-c5-u7-23 --nodes=1 --ntasks=1 --parsable)
cd "$PROJ"

for f in PairCheck.jl paircheck.jl paircheck.sh paircheck_summary.jl loadtest.jl jobcommon.sh; do
    [ -f optimize/$f ] || { echo "missing optimize/$f"; exit 1; }
done
[ -d /projects/jaca8590/Julia/Julia_depot_runletkf/compiled ] ||
    { echo "job depot not built; run the H1 setup job first"; exit 1; }
for s in $STATES; do
    [ -d "$s/results" ] || { echo "no $s/results"; exit 1; }
done
ARRAY=1-$NTASKS
[ "${PC_PILOT:-0}" = 1 ] && ARRAY=1

for STATE in $STATES; do
    LABEL=${STATE//\//_}
    LOAD=$(sbatch "${PRE[@]}" --job-name=pc_load_$LABEL --cpus-per-task=4 --mem=16G --time=00:30:00 \
        --output="$LOGDIR/%x_%j.out" <<JOB
#!/bin/bash
set -uo pipefail
cd $PROJ
source optimize/jobcommon.sh
LT=\$(mktemp)
"\$JULIA" --project="\$PROJ" -t 4 optimize/loadtest.jl 2>&1 | tee "\$LT" || exit 1
if grep -q -E "missing from the cache|does not support precompilation" "\$LT"; then
    echo "package caches are inconsistent on this node; not starting the array"; exit 1
fi
JOB
)
    ARR=$(sbatch --parsable --array=$ARRAY --job-name=pc_$LABEL \
        --export=ALL,PC_STATE=$STATE,PC_LABEL=$LABEL,PC_NTASKS=$NTASKS \
        --dependency=afterok:$LOAD --kill-on-invalid-dep=yes optimize/paircheck.sh)
    SUM=$(sbatch "${PRE[@]}" --job-name=pc_sum_$LABEL --cpus-per-task=1 --mem=16G --time=00:30:00 \
        --dependency=afterany:$ARR --kill-on-invalid-dep=yes --output="$LOGDIR/%x_%j.out" <<JOB
#!/bin/bash
set -uo pipefail
cd $PROJ
source optimize/jobcommon.sh
PC_LABEL=$LABEL PC_EXPECTED=$EXPECTED "\$JULIA" --project="\$PROJ" optimize/paircheck_summary.jl | tee pc/$LABEL/summary.txt
JOB
)
    echo "$STATE: load test $LOAD -> array $ARR ($ARRAY of $NTASKS) -> summary $SUM; results in pc/$LABEL/"
done
