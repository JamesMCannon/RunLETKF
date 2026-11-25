#!/bin/bash
# Pair-check progress for each state: units finished, task states of the latest array, completion
# rate and ETA, failed units, and units with roots absent from production so far.
# Live view, refreshed every 60 s:  watch -n 60 bash optimize/pcstatus.sh
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
STATES=${PC_STATES:-"g1/state_t1 h1/avid100_4x10"}
EXPECTED=${PC_EXPECTED:-2200}

echo "pair check at $(date '+%a %d %b %H:%M:%S')"
for STATE in $STATES; do
    LABEL=${STATE//\//_}
    UNITS=$PROJ/pc/$LABEL/units
    done=$(find "$UNITS" -name '*.jld2' 2>/dev/null | wc -l)
    recent=$(find "$UNITS" -name '*.jld2' -mmin -10 2>/dev/null | wc -l)
    echo
    printf "%s: %d of %d units (%d%%), %d in the last 10 min" "$LABEL" "$done" "$EXPECTED" $((100 * done / EXPECTED)) "$recent"
    if ((recent > 0 && done < EXPECTED)); then
        eta=$(( (EXPECTED - done) * 600 / recent ))
        printf ", ETA %d min (%s)" $((eta / 60)) "$(date -d "@$(( $(date +%s) + eta ))" '+%H:%M')"
    fi
    echo
    newest=$(ls -t "$LOGDIR"/pc_"$LABEL"_*_*.out 2>/dev/null | head -1)
    [ -z "$newest" ] && { echo "  no task logs yet"; continue; }
    ARR=$(basename "$newest" | sed -E 's/.*_([0-9]+)_[0-9]+\.out$/\1/')
    echo "  array $ARR tasks: $(sacct -X -n -P -j "$ARR" -o State 2>/dev/null | sed 's/ by .*//' | sort | uniq -c | awk '{printf "%s %s, ", $2, $1}' | sed 's/, $//')"
    logs=("$LOGDIR"/pc_"$LABEL"_"$ARR"_*.out)
    failed=$(cat "${logs[@]}" 2>/dev/null | grep -c "Error: unit failed")
    absent=$(cat "${logs[@]}" 2>/dev/null | grep -E "^member .* absent$" | awk '$(NF-1) > 0' )
    echo "  failed units in this array's logs: $failed"
    if [ -n "$absent" ]; then
        echo "  units with roots absent from production: $(echo "$absent" | wc -l)"
        echo "$absent" | sed 's/^/    /' | head -10
    else
        echo "  units with roots absent from production: 0"
    fi
    [ -f "$PROJ/pc/$LABEL/summary.txt" ] && { echo "  summary:"; sed 's/^/    /' "$PROJ/pc/$LABEL/summary.txt"; }
done
