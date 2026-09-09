#!/bin/bash
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00

set -euo pipefail
cd "$SYSTEM_DIR"
vmd -dispdev text -e ../../common/analyze_pore.tcl \
    -args system.psf system.pdb prod5.dcd "$ANALYSIS_PREFIX"
vmd -dispdev text -e ../../../final-ssdna/common/analyze_base_flipping.tcl \
    -args system.psf system.pdb prod5.dcd "$DNA_BASE" \
    "$ANALYSIS_PREFIX" 0.0 0.0 10 0.25

awk -F, 'NR>1{n++;rg+=$3;ee+=$4;lz+=$5;cz+=$6;ct+=$9;if(n==1){first=$5}last=$5}END{printf "frames %d\nmean_Rg_A %.4f\nmean_EE_A %.4f\nlead_z_initial_A %.4f\nlead_z_final_A %.4f\nmean_DNA_COM_z_A %.4f\nmean_DNA_protein_contacts_3p5A %.4f\n",n,rg/n,ee/n,first,last,cz/n,ct/n}' \
    ${ANALYSIS_PREFIX}_pore_metrics.csv > summary_pore.txt
awk -F, 'NR>1{n++;sa+=$2;st+=$3;un+=$4;sy+=$5}END{printf "mean_base_SASA_A2 %.4f\nmean_stacked_pairs %.4f\nmean_unstacked_fraction %.4f\nmean_syn_like_fraction %.4f\n",sa/n,st/n,un/n,sy/n}' \
    ${ANALYSIS_PREFIX}_base_flipping_time.csv >> summary_pore.txt
