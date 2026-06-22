#!/bin/bash

#===============================================================================
# SCRIPT: heatmap_publication.sh
# PURPOSE: Generate publication-ready heatmaps from DESeq2 results
#
# DESCRIPTION:
# Creates heatmaps in three versions:
#   1. Top 50 DEGs (by adjusted p-value)
#   2. Genes from custom gene list (TES_degs.txt)
#   3. All significant DEGs
#
# Style: Dark blue to intensive red color scale, z-score normalized, no dendrograms
# Gene filtering: Genes without valid gene symbol mapping are excluded
#
# EXPECTED INPUTS:
# - DESeq2 results: results/05_deseq2/deseq2_results_TES_vs_GFP.txt
# - Normalized counts: results/05_deseq2/normalized_counts.txt
# - Gene list: SRF_Eva_integrated_analysis/data/TES_degs.txt
#
# EXPECTED OUTPUTS (PNG only, with labeled and nolabel versions):
# - results/05_deseq2/plots/heatmap_top50_degs.png (and _nolabel.png)
# - results/05_deseq2/plots/heatmap_genelist.png (and _nolabel.png)
# - results/05_deseq2/plots/heatmap_all_degs.png
# - results/05_deseq2/plots/heatmap_all_degs_labeled.png (and _nolabel.png)
#
# USAGE:
# sbatch scripts/heatmap_publication.sh
#===============================================================================

#SBATCH --job-name=heatmap_pub
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/heatmap_publication.err"
#SBATCH --output="./logs/heatmap_publication.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
DESEQ_RESULTS="${BASE_DIR}/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
NORM_COUNTS="${BASE_DIR}/results/05_deseq2/normalized_counts.txt"
GENE_LIST="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/data/TES_degs.txt"
OUTPUT_DIR="${BASE_DIR}/results/05_deseq2/plots"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

echo "=== Generating Publication-Ready Heatmaps ==="
echo "Timestamp: $(date)"
echo "DESeq2 results: ${DESEQ_RESULTS}"
echo "Normalized counts: ${NORM_COUNTS}"
echo "Gene list: ${GENE_LIST}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if input files exist
for file in "${DESEQ_RESULTS}" "${NORM_COUNTS}" "${GENE_LIST}"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: File not found: $file"
        exit 1
    fi
done

# Run the R script
echo "Running R script..."
Rscript ${OUTPUT_DIR}/heatmap_publication.R

# Check outputs
echo ""
echo "=== Checking output files ==="

for name in top50_degs top50_degs_nolabel genelist genelist_nolabel all_degs all_degs_labeled all_degs_labeled_nolabel; do
    file="${OUTPUT_DIR}/heatmap_${name}.png"
    if [[ -f "$file" ]]; then
        size=$(ls -lh "$file" | awk '{print $5}')
        echo "SUCCESS: heatmap_${name}.png ($size)"
    fi
done

echo ""
echo "=== Heatmap generation complete ==="
echo "Timestamp: $(date)"
