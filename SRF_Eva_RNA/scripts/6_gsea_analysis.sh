#!/bin/bash

#===============================================================================
# SCRIPT: gsea_analysis.sh
# PURPOSE: Gene Set Enrichment Analysis (GSEA) with publication-ready plots
#
# DESCRIPTION:
# Performs GSEA using ranked gene list from DESeq2 results.
# Creates classic enrichment plots showing:
#   - Running enrichment score
#   - Gene hit markers
#   - Ranked list gradient
#   - NES and FDR annotations
#
# EXPECTED INPUTS:
# - DESeq2 results: results/05_deseq2/deseq2_results_TES_vs_GFP.txt
#
# EXPECTED OUTPUTS:
# - results/06_gsea/gsea_results_*.csv
# - results/06_gsea/plots/gsea_enrichment_*.pdf/png
# - results/06_gsea/plots/gsea_dotplot_*.pdf/png
#
# USAGE:
# sbatch scripts/gsea_analysis.sh
#===============================================================================

#SBATCH --job-name=gsea_analysis
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/gsea_analysis.err"
#SBATCH --output="./logs/gsea_analysis.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate analysis3_env

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
INPUT_FILE="${BASE_DIR}/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
OUTPUT_DIR="${BASE_DIR}/results/06_gsea"

# Create output directories
mkdir -p ${OUTPUT_DIR}/plots

echo "=== Gene Set Enrichment Analysis (GSEA) ==="
echo "Timestamp: $(date)"
echo "Input: ${INPUT_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if input file exists
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "ERROR: DESeq2 results file not found: ${INPUT_FILE}"
    echo "Please run 5_deseq2.sh first."
    exit 1
fi

# Run the R script
echo "Running GSEA analysis in R..."
echo "This may take 15-30 minutes depending on the number of gene sets..."
Rscript ${BASE_DIR}/scripts/gsea_analysis.R

# Check outputs
echo ""
echo "=== Checking output files ==="

# Check results files
for collection in hallmark GO_BP reactome KEGG; do
    file="${OUTPUT_DIR}/gsea_results_${collection}.csv"
    if [[ -f "$file" ]]; then
        n_lines=$(wc -l < "$file")
        echo "SUCCESS: gsea_results_${collection}.csv ($n_lines pathways)"
    fi
done

# Check plots
n_plots=$(find ${OUTPUT_DIR}/plots -name "*.png" 2>/dev/null | wc -l)
echo "Total plots created: $n_plots"

# Check summary
if [[ -f "${OUTPUT_DIR}/gsea_significant_summary.csv" ]]; then
    n_sig=$(tail -n +2 "${OUTPUT_DIR}/gsea_significant_summary.csv" | wc -l)
    echo "Significant pathways in summary: $n_sig"
fi

echo ""
echo "=== GSEA analysis complete ==="
echo "Timestamp: $(date)"
