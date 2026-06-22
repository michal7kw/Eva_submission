#!/bin/bash

#===============================================================================
# SCRIPT: rerun_deseq2.sh
# PURPOSE: Run the UPDATED deseq2_analysis.R directly
#
# DESCRIPTION:
# This wrapper runs the updated deseq2_analysis.R (with improved filtering
# and apeglm LFC shrinkage) WITHOUT the 5_deseq2.sh heredoc that would
# overwrite the updated R script with the old version.
#
# KEY CHANGES IN UPDATED deseq2_analysis.R:
# 1. Improved filtering: rowSums(counts(dds) >= 10) >= 3
#    (old: rowSums(counts(dds)) >= 10)
# 2. Added apeglm LFC shrinkage for GSEA-compatible fold changes
# 3. Proper alpha=0.05 for optimized independent filtering
#
# EXPECTED INPUTS:
# - Count matrix: results/04_quantified/count_matrix.txt
# - Sample metadata: results/04_quantified/sample_metadata.txt
#
# EXPECTED OUTPUTS:
# - results/05_deseq2/deseq2_results_TES_vs_GFP.txt (with shrunken LFC)
# - results/05_deseq2/significant_genes_TES_vs_GFP.txt
# - results/05_deseq2/normalized_counts.txt
# - results/05_deseq2/plots/*.png
#
# USAGE:
# sbatch scripts/rerun_deseq2.sh
#===============================================================================

#SBATCH --job-name=rerun_deseq2
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/rerun_deseq2.err"
#SBATCH --output="./logs/rerun_deseq2.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
OUTPUT_DIR="${BASE_DIR}/results/05_deseq2"

# Create plots directory
mkdir -p ${OUTPUT_DIR}/plots

echo "=== Re-running DESeq2 with updated analysis script ==="
echo "Timestamp: $(date)"
echo "Working directory: ${OUTPUT_DIR}"
echo ""
echo "Changes in updated script:"
echo "  - Improved gene filtering: >= 10 counts in >= 3 samples"
echo "  - apeglm LFC shrinkage for better fold change estimates"
echo "  - alpha=0.05 for optimized independent filtering"
echo ""

# Run the UPDATED R script directly (not the heredoc version from 5_deseq2.sh)
Rscript ${OUTPUT_DIR}/deseq2_analysis.R

echo ""
echo "=== Checking output files ==="

if [[ -f "${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt" ]]; then
    RESULTS_LINES=$(wc -l < "${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt")
    echo "SUCCESS: DESeq2 results file created (${RESULTS_LINES} lines)"
else
    echo "ERROR: DESeq2 results file not created!"
    exit 1
fi

if [[ -f "${OUTPUT_DIR}/significant_genes_TES_vs_GFP.txt" ]]; then
    SIG_LINES=$(tail -n +2 "${OUTPUT_DIR}/significant_genes_TES_vs_GFP.txt" | wc -l)
    echo "SUCCESS: Significant genes file created (${SIG_LINES} genes)"
fi

if [[ -d "${OUTPUT_DIR}/plots" ]]; then
    PLOT_COUNT=$(find "${OUTPUT_DIR}/plots" -name "*.png" | wc -l)
    echo "SUCCESS: ${PLOT_COUNT} plots generated"
fi

echo ""
echo "=== DESeq2 re-analysis complete ==="
echo "Timestamp: $(date)"
