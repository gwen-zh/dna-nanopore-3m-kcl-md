#!/bin/bash
#SBATCH --partition=dept_gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --nice=1000

set -euo pipefail
cd "$SYSTEM_DIR"
module load namd/3.0/3.0.2
vmd -dispdev text -e ../../common/make_restraints.tcl \
    > make_restraints.log 2>&1
for stage in min heat eq1 eq2 eq3; do
    namd3 +p8 +devices 0 ../../common/${stage}.namd \
        > ${stage}.log 2>&1
    grep -q "End of program" ${stage}.log
    ! grep -q "FATAL ERROR" ${stage}.log
done
