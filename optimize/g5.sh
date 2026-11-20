#!/bin/bash
#SBATCH --job-name=letkf_g5
#SBATCH --qos=blanca-lair
#SBATCH --exclude=bhpc-c5-u7-23
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --time=16:00:00
#SBATCH --output=/scratch/alpine/jaca8590/Julia/Optimize_LETKF/logs/%x_%j.out

# G5 run <config>: Main.jl end to end, results in g5/<config>. A config is <paths><members>_<layout>:
# paths min (2 paths, 3 iterations) or avid (22 paths, 1 iteration); layout 1x40 or 4x10.
# As an array task, the config is entry SLURM_ARRAY_TASK_ID of G5_CFGS.
set -euo pipefail
CFG=${1:-}
if [ -z "$CFG" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CFGS=(${G5_CFGS:-})
    CFG=${CFGS[$((SLURM_ARRAY_TASK_ID - 1))]:-}
fi
[[ $CFG =~ ^(min|avid)([0-9]+)_(1x40|4x10)$ ]] || { echo "usage: g5.sh <min|avid><members>_<1x40|4x10>"; exit 1; }
SET=${BASH_REMATCH[1]} MEMBERS=${BASH_REMATCH[2]} LAYOUT=${BASH_REMATCH[3]}
source /projects/jaca8590/Julia/Optimize_LETKF/RunLETKF/optimize/jobcommon.sh

case $SET in
    min)  export PATHS=minimum ITRS=3 ;;
    avid) export PATHS=AVID ITRS=1 ;;
esac
export ENS_SIZE=$MEMBERS NEW_FOLDER=g5/$CFG
case $LAYOUT in
    1x40) export DISTRIBUTED=false ;;
    4x10) export DISTRIBUTED=true NPROCS=4 PROC_THREADS=10 ;;
esac

OUT=$PROJ/$NEW_FOLDER
[ -e "$OUT" ] && { echo "$OUT exists; refusing to overwrite or resume"; exit 1; }
mkdir -p "$OUT"

# Node memory in use (MiB): at start and the peak, sampled every 30 s.
BASE=$(free -m | awk '/^Mem:/{print $3}')
echo "$BASE $BASE" > "$OUT/peakmem.txt"
( peak=$BASE
  while sleep 30; do
      used=$(free -m | awk '/^Mem:/{print $3}')
      if (( used > peak )); then peak=$used; echo "$BASE $peak" > "$OUT/peakmem.txt"; fi
  done ) &
MON=$!

echo "== $CFG: PATHS=$PATHS ENS_SIZE=$ENS_SIZE ITRS=$ITRS DISTRIBUTED=$DISTRIBUTED ${NPROCS:+NPROCS=$NPROCS PROC_THREADS=$PROC_THREADS}, start $(date -Is)"
STATUS=0
JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK srun --export=ALL --mpi=none "$JULIA" --project="$PROJ" optimize/g5run.jl || STATUS=$?
kill $MON
read -r b p < "$OUT/peakmem.txt"
echo "== $CFG end $(date -Is), status $STATUS; node memory in use: start $((b / 1024)) GiB, peak $((p / 1024)) GiB"
exit $STATUS
