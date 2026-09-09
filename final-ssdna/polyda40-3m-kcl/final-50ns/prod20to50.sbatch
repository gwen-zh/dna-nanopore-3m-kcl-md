#!/bin/bash
#SBATCH --job-name=pda_20to50
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=12G
#SBATCH --time=2-00:00:00
#SBATCH --nice=1000
#SBATCH --output=prod20to50.%j.out
#SBATCH --error=prod20to50.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load namd/3.0/3.0.2
namd3 +p8 +devices 0 prod20to50.namd
