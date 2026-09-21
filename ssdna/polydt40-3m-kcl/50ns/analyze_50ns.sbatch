#!/bin/bash
#SBATCH --job-name=pdt_analyze_50ns
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=analysis.%j.out
#SBATCH --error=analysis.%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

vmd -dispdev text -e ../../scripts/analyze_conformation.tcl \
    -args system.psf system.pdb prod20.dcd .segment20 10.0 0.0
vmd -dispdev text -e ../../scripts/analyze_conformation.tcl \
    -args system.psf system.pdb prod20to50.dcd .segment30 40.0 20.0
vmd -dispdev text -e ../../scripts/analyze_base_flipping.tcl \
    -args system.psf system.pdb prod20.dcd THY .segment20 0.0 10.0
vmd -dispdev text -e ../../scripts/analyze_base_flipping.tcl \
    -args system.psf system.pdb prod20to50.dcd THY .segment30 20.0 40.0

awk 'FNR == 1 && NR != 1 {next} {print}' \
    .segment20_metrics.csv .segment30_metrics.csv \
    > analysis_metrics_50ns.csv
awk 'FNR == 1 && NR != 1 {next} {print}' \
    .segment20_base_flipping_time.csv .segment30_base_flipping_time.csv \
    > analysis_base_flipping_time_50ns.csv
mv .segment30_base_flipping_residue.csv \
    analysis_base_flipping_residue_last10ns.csv
mv .segment30_representative_dna.pdb representative_last10ns_dna.pdb
mv .segment30_representative_full.pdb representative_last10ns_full.pdb
mv .segment30_selected.txt representative_last10ns_selection.txt

awk -F, 'NR>1{n++;rg+=$3;ee+=$4;if($2>=40){m++;rg2+=$3;ee2+=$4}}END{printf "all50ns_Rg_A %.4f\nall50ns_EE_A %.4f\nlast10ns_Rg_A %.4f\nlast10ns_EE_A %.4f\n",rg/n,ee/n,rg2/m,ee2/m}' \
    analysis_metrics_50ns.csv > summary.txt
awk -F, 'NR>1 && $1>=40{n++;sa+=$2;st+=$3;un+=$4;sy+=$5}END{printf "last10ns_base_SASA_A2 %.4f\nlast10ns_stacked_pairs %.4f\nlast10ns_unstacked_fraction %.4f\nlast10ns_syn_like_fraction %.4f\n",sa/n,st/n,un/n,sy/n}' \
    analysis_base_flipping_time_50ns.csv >> summary.txt

rm -f .segment20_* .segment30_metrics.csv .segment30_base_flipping_time.csv
