#!/bin/bash
# One-shot status of the latest H1 submission: Slurm state of every job and, for each Main.jl run,
# the iteration, ensemble progress (1x40 runs print a progress meter; 4x10 runs do not), an ETA,
# node memory and counts of polishing notes and errors. Live view, refreshed every 60 s:
#   watch -n 60 bash optimize/h1status.sh
PROJ=/projects/jaca8590/Julia/Optimize_LETKF/RunLETKF
LOGDIR=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs
NAMES=letkf_h1_setup,letkf_h1,letkf_h1_check,letkf_h1_analysis

hms() { local s=${1%.*}; ((s < 0)) && s=0; printf "%d:%02d:%02d" $((s/3600)) $((s%3600/60)) $((s%60)); }

# G5 wall time (s) of the same configuration, from its log.
g5total() {
    local f; f=$(grep -l "^== $1: " "$LOGDIR"/letkf_g5_*.out 2>/dev/null | tail -1)
    [ -n "$f" ] && awk '/^== total /{print $3; exit}' "$f"
}

runinfo() {  # $1 = log file
    local log=$1 cfg itrs done starts first cur now el frac prog eta mean g5 mem notes errs
    cfg=$(awk '/^== [a-z]+[0-9]+_[0-9x]+: PATHS=/{sub(/:$/, "", $2); print $2; exit}' "$log")
    [ -z "$cfg" ] && { echo "      (Main.jl not started yet)"; return; }
    itrs=$(grep -m1 -o 'ITRS=[0-9]*' "$log" | cut -d= -f2)
    done=$(grep -c 'Info: Elapsed' "$log")
    mapfile -t starts < <(grep -a -o 'start_time = [0-9T:.-]*' "$log" | awk '{print $3}')
    now=$(date +%s)
    prog=$(tr '\r' '\n' < "$log" | grep -a 'Progress:' | tail -1 | sed 's/\x1b\[[0-9;]*[A-Za-z]//g')
    frac=""
    [[ $prog =~ ([0-9]+)% ]] && frac=${BASH_REMATCH[1]}
    if grep -q "^== $cfg end" "$log"; then
        echo "      $cfg finished: $(grep -m1 "^== total" "$log")"
        return
    fi
    echo -n "      $cfg: iteration $((done + 1))/$itrs"
    if ((${#starts[@]} > done)); then
        first=$(date -d "${starts[0]}" +%s); cur=$(date -d "${starts[done]}" +%s)
        el=$((now - cur))
        echo -n ", this iteration $(hms $el)"
        [ -n "$frac" ] && echo -n ", members ${frac}%"
        eta=""
        if ((done > 0)); then
            mean=$(( (cur - first) / done ))
            if [ -n "$frac" ] && ((frac > 0)); then
                eta=$(( (itrs - done - 1) * mean + el * (100 - frac) / frac ))
            else
                eta=$(( (itrs - done) * mean - el ))
            fi
        elif [ -n "$frac" ] && ((frac > 0)); then
            eta=$(( itrs * el * 100 / frac - el ))
        fi
        if [ -n "$eta" ]; then
            echo -n ", ETA $(hms $eta) ($(date -d "@$((now + eta))" '+%a %H:%M'))"
        else
            echo -n ", ETA unknown until an iteration completes"
        fi
        echo
    else
        echo " (setting up)"
    fi
    [ -n "$prog" ] && echo "      progress meter: $(echo "$prog" | tr -s ' ' | cut -c1-110)"
    g5=$(g5total "$cfg")
    [ -n "$g5" ] && echo "      G5 wall time for $cfg: $(hms "$g5")"
    mem="$PROJ/h1/$cfg/peakmem.txt"
    [ -f "$mem" ] && read -r b p < "$mem" && echo "      node memory: start $((b / 1024)) GiB, peak $((p / 1024)) GiB"
    notes=$(grep -a -c -E 'root dropped|kept unpolished|guard90 of 90' "$log")
    errs=$(grep -a -c -E 'ERROR|Error' "$log")
    echo "      polishing notes $notes, error lines $errs"
}

echo "H1 status at $(date '+%a %d %b %H:%M:%S')"
SETUP=$(sacct -u "$USER" -S now-14days -X -n -P --name=letkf_h1_setup -o JobID 2>/dev/null | sort -n | tail -1)
[ -z "$SETUP" ] && { echo "no letkf_h1_setup job found in the last 14 days"; exit 0; }
sacct -u "$USER" -S now-14days -X -n -P --name="$NAMES" -o JobID,JobName,State,Elapsed,NodeList |
while IFS='|' read -r id name state elapsed node; do
    (( ${id%%_*} < SETUP )) && continue
    printf "%-16s %-18s %-11s %10s  %s\n" "$id" "$name" "$state" "$elapsed" "$node"
    if [ "$name" = letkf_h1 ] && [[ $id =~ ^[0-9]+_[0-9]+$ ]] && [ "$state" != PENDING ]; then
        log="$LOGDIR/letkf_h1_${id}.out"
        [ -f "$log" ] && runinfo "$log"
    fi
done
for f in check analysis; do
    if [ -f "$PROJ/h1/$f.txt" ]; then echo; echo "== h1/$f.txt"; tail -n 8 "$PROJ/h1/$f.txt"; fi
done
