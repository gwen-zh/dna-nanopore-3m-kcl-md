#!/bin/bash
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=12:00:00

set -euo pipefail
cd "$SYSTEM_DIR"
vmd -dispdev text -e ../../common/build_formal_pore.tcl \
    -args "$PORE_PSF" "$PORE_PDB" "$DONOR_PDB" "$DNA_BASE" \
    "$REMOVE_OLD_DNA" "$TARGET_Z" "$MIN_X" "$MIN_Y" "$MIN_Z" \
    "$MAX_X" "$MAX_Y" "$MAX_Z" "$SYSTEM_LABEL"
test -s system.psf
test -s system.pdb
grep -q "POT" system.pdb
grep -q "CLA" system.pdb
