#!/bin/bash
#SBATCH --job-name=pda_formeq
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --output=prepare.%j.out
#SBATCH --error=prepare.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
module load namd/3.0/3.0.2
for stage in min heat eq1 eq2 eq3; do
    namd3 +p8 +devices 0 ${stage}.namd > ${stage}.log 2>&1
    grep -q "End of program" ${stage}.log
    ! grep -q "FATAL ERROR" ${stage}.log
done
