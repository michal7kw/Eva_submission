#!/bin/bash

#===============================================================================
# SCRIPT: heatmap_tead_family.sh
# PURPOSE: Generate a TEAD-family (TEAD1-4) expression heatmap, Mock vs TES
#
# DESCRIPTION:
# Wraps heatmap_tead_family.R. Subsets the DESeq2 normalized counts to
# TEAD1/TEAD2/TEAD3/TEAD4 and draws a row z-scored heatmap (Mock = GFP control
# vs TES, 3 replicates each), reusing the project's publication heatmap style.
#
# EXPECTED INPUTS:
# - DESeq2 results:    results/05_deseq2/deseq2_results_TES_vs_GFP.txt
# - Normalized counts: results/05_deseq2/normalized_counts.txt
#
# EXPECTED OUTPUTS (PNG + PDF, labelled and nolabel variants):
# - results/05_deseq2/plots/heatmap_tead_family.{png,pdf}
# - results/05_deseq2/plots/heatmap_tead_family_nolabel.{png,pdf}
#
# USAGE:
# sbatch scripts/heatmap_tead_family.sh
#===============================================================================

#SBATCH --job-name=heatmap_tead
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/heatmap_tead_family.err"
#SBATCH --output="./logs/heatmap_tead_family.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
SCRIPT_DIR="${BASE_DIR}/scripts"
DESEQ_RESULTS="${BASE_DIR}/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
NORM_COUNTS="${BASE_DIR}/results/05_deseq2/normalized_counts.txt"
OUTPUT_DIR="${BASE_DIR}/results/05_deseq2/plots"

mkdir -p "${OUTPUT_DIR}"

echo "=== TEAD-family Expression Heatmap ==="
echo "Timestamp: $(date)"

for file in "${DESEQ_RESULTS}" "${NORM_COUNTS}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done

echo "Running R script..."
Rscript "${SCRIPT_DIR}/heatmap_tead_family.R"

echo ""
echo "=== Checking output files ==="
for name in heatmap_tead_family heatmap_tead_family_nolabel; do
    for ext in png pdf; do
        file="${OUTPUT_DIR}/${name}.${ext}"
        if [[ -f "$file" ]]; then
            size=$(ls -lh "$file" | awk '{print $5}')
            echo "SUCCESS: ${name}.${ext} ($size)"
        else
            echo "MISSING: ${name}.${ext}"
        fi
    done
done

echo ""
echo "=== Heatmap generation complete ==="
echo "Timestamp: $(date)"
