#!/bin/bash
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --nice=0

set -euo pipefail
cd "$SYSTEM_DIR"
module load namd/3.0/3.0.2
namd3 +p8 +devices 0 ../../common/prod5.namd
