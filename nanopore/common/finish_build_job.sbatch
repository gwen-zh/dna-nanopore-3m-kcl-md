#!/bin/bash
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -euo pipefail
cd "$SYSTEM_DIR"
export VMDFORCECPUCOUNT=1 VMDCUDA=0
vmd -dispdev text -e ../../common/finish_ionization.tcl
test -s system.psf
test -s system.pdb
python3 ../../common/validate_bonds.py
