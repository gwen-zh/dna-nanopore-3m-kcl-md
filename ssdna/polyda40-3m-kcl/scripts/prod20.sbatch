#!/bin/bash
#SBATCH --job-name=pda_prod20
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00
#SBATCH --output=../results/prod20.%j.out
#SBATCH --error=../results/prod20.%j.err

set -euo pipefail
SCRIPT_DIR="${SSDNA_DA_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
[[ -f "$SCRIPT_DIR/prod20.namd" && -d "$SCRIPT_DIR/../results" ]] || { echo "Submit from the scripts directory" >&2; exit 2; }
cd "$SCRIPT_DIR/../results"
module load namd/3.0/3.0.2
namd3 +p8 +devices 0 "$SCRIPT_DIR/prod20.namd"
