#!/bin/bash

#===============================================================================
# SCRIPT: generate_supplementary_tables.sh
# PURPOSE: Generate publication-ready supplementary tables (S2, S3, S5)
#
# DESCRIPTION:
# Compiles data from RNA-seq DEGs, GO/GSEA enrichment, and CUT&Run peaks
# into formatted CSV files for supplementary materials.
#
# EXPECTED INPUTS:
# - DESeq2 results: SRF_Eva_RNA/results/05_deseq2/significant_genes_TES_vs_GFP.txt
# - GSEA results: SRF_Eva_RNA/results/06_gsea/gsea_results_*.csv
# - GO results: SRF_Eva_RNA/results/06_go_enrichment/*/GO_*.csv
# - Peak files: SRF_Eva_Cut_and_Tag/results/07_analysis_narrow/*_peaks_annotated.csv
#
# EXPECTED OUTPUTS:
# - supplementary_tables/Table_S2_DEGs_RNA_seq.csv
# - supplementary_tables/Table_S3_GO_GSEA.csv
# - supplementary_tables/Table_S5a_TES_peaks.csv
# - supplementary_tables/Table_S5b_TEAD1_peaks.csv
#
# USAGE:
# sbatch scripts/generate_supplementary_tables.sh
#===============================================================================

#SBATCH --job-name=supp_tables
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/generate_supplementary_tables.err"
#SBATCH --output="./logs/generate_supplementary_tables.out"

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission"
SCRIPT_DIR="${BASE_DIR}/scripts"
OUTPUT_DIR="${BASE_DIR}/supplementary_tables"

echo "==============================================================================="
echo "Generating Supplementary Tables"
echo "==============================================================================="
echo "Timestamp: $(date)"
echo "Working directory: ${BASE_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Check required input files
echo "=== Checking input files ==="

# Table S2 source
DESEQ_FILE="${BASE_DIR}/SRF_Eva_RNA/results/05_deseq2/significant_genes_TES_vs_GFP.txt"
if [[ -f "${DESEQ_FILE}" ]]; then
    echo "OK: DESeq2 results found"
else
    echo "ERROR: DESeq2 results not found: ${DESEQ_FILE}"
    exit 1
fi

# Table S3 sources (GSEA)
GSEA_DIR="${BASE_DIR}/SRF_Eva_RNA/results/06_gsea"
if [[ -d "${GSEA_DIR}" ]]; then
    n_gsea=$(ls ${GSEA_DIR}/gsea_results_*.csv 2>/dev/null | wc -l)
    echo "OK: GSEA directory found (${n_gsea} result files)"
else
    echo "WARNING: GSEA directory not found: ${GSEA_DIR}"
fi

# Table S3 sources (GO)
GO_DIR="${BASE_DIR}/SRF_Eva_RNA/results/06_go_enrichment"
if [[ -d "${GO_DIR}" ]]; then
    echo "OK: GO enrichment directory found"
else
    echo "WARNING: GO enrichment directory not found: ${GO_DIR}"
fi

# Table S5 sources
PEAKS_DIR="${BASE_DIR}/SRF_Eva_Cut_and_Tag/results/07_analysis_narrow"
if [[ -f "${PEAKS_DIR}/TES_peaks_annotated.csv" ]]; then
    echo "OK: TES peaks annotated file found"
else
    echo "ERROR: TES peaks not found: ${PEAKS_DIR}/TES_peaks_annotated.csv"
    exit 1
fi

if [[ -f "${PEAKS_DIR}/TEAD1_peaks_annotated.csv" ]]; then
    echo "OK: TEAD1 peaks annotated file found"
else
    echo "ERROR: TEAD1 peaks not found: ${PEAKS_DIR}/TEAD1_peaks_annotated.csv"
    exit 1
fi

echo ""
echo "=== Running R script ==="
Rscript ${SCRIPT_DIR}/generate_supplementary_tables.R

# Check outputs
echo ""
echo "=== Checking output files ==="

for table in "Table_S2_DEGs_RNA_seq.csv" "Table_S3_GO_GSEA.csv" "Table_S5a_TES_peaks.csv" "Table_S5b_TEAD1_peaks.csv"; do
    file="${OUTPUT_DIR}/${table}"
    if [[ -f "$file" ]]; then
        n_rows=$(tail -n +2 "$file" | wc -l)
        size=$(du -h "$file" | cut -f1)
        echo "SUCCESS: ${table} (${n_rows} rows, ${size})"
    else
        echo "ERROR: ${table} not created"
    fi
done

echo ""
echo "==============================================================================="
echo "Supplementary table generation complete"
echo "Timestamp: $(date)"
echo "==============================================================================="
