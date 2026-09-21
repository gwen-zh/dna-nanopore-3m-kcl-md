#!/bin/bash
#SBATCH --job-name=pda_analyze_50ns
#SBATCH --partition=dept_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=../results/analysis.%j.out
#SBATCH --error=../results/analysis.%j.err

set -euo pipefail
SCRIPT_DIR="${SSDNA_DA_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
RUN_DIR="$SCRIPT_DIR/../results"
ANALYSIS_DIR="$RUN_DIR/analysis"
[[ -f "$SCRIPT_DIR/analyze_conformation.tcl" && -d "$RUN_DIR" ]] || { echo "Submit from the scripts directory" >&2; exit 2; }
mkdir -p "$ANALYSIS_DIR"
cd "$RUN_DIR"

vmd -dispdev text -e "$SCRIPT_DIR/analyze_conformation.tcl" \
    -args "$SCRIPT_DIR/input/system.psf" "$SCRIPT_DIR/input/system.pdb" prod20.dcd "$ANALYSIS_DIR/.segment20" 10.0 0.0
vmd -dispdev text -e "$SCRIPT_DIR/analyze_conformation.tcl" \
    -args "$SCRIPT_DIR/input/system.psf" "$SCRIPT_DIR/input/system.pdb" prod20to50.dcd "$ANALYSIS_DIR/.segment30" 40.0 20.0
vmd -dispdev text -e "$SCRIPT_DIR/analyze_base_flipping.tcl" \
    -args "$SCRIPT_DIR/input/system.psf" "$SCRIPT_DIR/input/system.pdb" prod20.dcd ADE "$ANALYSIS_DIR/.segment20" 0.0 10.0
vmd -dispdev text -e "$SCRIPT_DIR/analyze_base_flipping.tcl" \
    -args "$SCRIPT_DIR/input/system.psf" "$SCRIPT_DIR/input/system.pdb" prod20to50.dcd ADE "$ANALYSIS_DIR/.segment30" 20.0 40.0

awk 'FNR == 1 && NR != 1 {next} {print}' \
    "$ANALYSIS_DIR/.segment20_metrics.csv" "$ANALYSIS_DIR/.segment30_metrics.csv" \
    > "$ANALYSIS_DIR/analysis_metrics_50ns.csv"
awk 'FNR == 1 && NR != 1 {next} {print}' \
    "$ANALYSIS_DIR/.segment20_base_flipping_time.csv" "$ANALYSIS_DIR/.segment30_base_flipping_time.csv" \
    > "$ANALYSIS_DIR/analysis_base_flipping_time_50ns.csv"
mv "$ANALYSIS_DIR/.segment30_base_flipping_residue.csv" \
    "$ANALYSIS_DIR/analysis_base_flipping_residue_last10ns.csv"
mv "$ANALYSIS_DIR/.segment30_representative_dna.pdb" "$ANALYSIS_DIR/representative_last10ns_dna.pdb"
mv "$ANALYSIS_DIR/.segment30_representative_full.pdb" "$ANALYSIS_DIR/representative_last10ns_full.pdb"
mv "$ANALYSIS_DIR/.segment30_selected.txt" "$ANALYSIS_DIR/representative_last10ns_selection.txt"

awk -F, 'NR>1{n++;rg+=$3;ee+=$4;if($2>=40){m++;rg2+=$3;ee2+=$4}}END{printf "all50ns_Rg_A %.4f\nall50ns_EE_A %.4f\nlast10ns_Rg_A %.4f\nlast10ns_EE_A %.4f\n",rg/n,ee/n,rg2/m,ee2/m}' \
    "$ANALYSIS_DIR/analysis_metrics_50ns.csv" > "$ANALYSIS_DIR/summary.txt"
awk -F, 'NR>1 && $1>=40{n++;sa+=$2;st+=$3;un+=$4;sy+=$5}END{printf "last10ns_base_SASA_A2 %.4f\nlast10ns_stacked_pairs %.4f\nlast10ns_unstacked_fraction %.4f\nlast10ns_syn_like_fraction %.4f\n",sa/n,st/n,un/n,sy/n}' \
    "$ANALYSIS_DIR/analysis_base_flipping_time_50ns.csv" >> "$ANALYSIS_DIR/summary.txt"

rm -f "$ANALYSIS_DIR"/.segment20_* "$ANALYSIS_DIR/.segment30_metrics.csv" "$ANALYSIS_DIR/.segment30_base_flipping_time.csv"
