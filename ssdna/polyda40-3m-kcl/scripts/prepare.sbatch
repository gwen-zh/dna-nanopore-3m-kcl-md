#!/bin/bash
#SBATCH --job-name=pda_prepare
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --output=../results/prepare.%j.out
#SBATCH --error=../results/prepare.%j.err

set -euo pipefail
SCRIPT_DIR="${SSDNA_DA_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_DIR="$SCRIPT_DIR/../results"
[[ -f "$SCRIPT_DIR/min.namd" && -d "$RUN_DIR" ]] || { echo "Submit from the scripts directory" >&2; exit 2; }
cd "$RUN_DIR"
module load namd/3.0/3.0.2
for stage in min heat eq1 eq2 eq3; do
    namd3 +p8 +devices 0 "$SCRIPT_DIR/${stage}.namd" > "${stage}.log" 2>&1
    grep -q "End of program" ${stage}.log
    ! grep -q "FATAL ERROR" ${stage}.log
done
