#!/bin/bash

#===============================================================================
# SCRIPT: 9_diff_bind.sh
# PURPOSE: Differential binding analysis for Cut&Tag data using DiffBind
#
# DESCRIPTION:
# This script performs statistical differential binding analysis between
# experimental conditions using the DiffBind R package. It identifies
# regions with significantly different binding patterns between TES, TESmut,
# and TEAD1 conditions, providing statistical validation of binding differences.
#
# KEY OPERATIONS:
# 1. Sample correlation analysis and quality assessment
# 2. Peak consensus set generation across replicates
# 3. Read count normalization and statistical modeling
# 4. Differential binding site identification
# 5. Generation of MA plots, PCA plots, and heatmaps
# 6. Export of significantly different binding sites
#
# METHODOLOGY:
# - Uses DiffBind R package for comprehensive analysis
# - Employs DESeq2 or edgeR for statistical testing
# - Performs TMM normalization for count data
# - Applies FDR correction for multiple testing
# - Generates consensus peak sets from replicates
#
# IMPORTANT PARAMETERS:
# - Memory: 32GB (required for large matrix operations)
# - Time: 4 hours (sufficient for statistical analysis)
# - Threads: 16 (parallel processing for count operations)
# - FDR threshold: typically 0.05 (set in R script)
# - Fold change threshold: typically 2-fold (set in R script)
#
# INPUTS:
# - Sample sheet with BAM files and peak files (diffbind_samplesheet.csv)
# - Filtered BAM files from alignment pipeline
# - Peak files from MACS2 analysis
# - Configuration parameters in R script
#
# OUTPUTS:
# - Differential binding results (CSV files)
# - Quality control plots (PCA, correlation heatmaps)
# - MA plots showing fold changes and significance
# - Venn diagrams of overlapping differential sites
# - Normalized count matrices
#
# DEPENDENCIES:
# - R with DiffBind, DESeq2, and visualization packages
# - Conda environment: diffbind_analysis
# - Sample sheet configuration file
#
# USAGE:
# sbatch 9_diff_bind.sh
#
# NOTES:
# - Requires properly formatted sample sheet
# - Generates publication-ready statistical plots
# - Provides rigorous statistical validation of binding differences
# - Results inform biological interpretation of Cut&Tag data
#===============================================================================

#SBATCH --job-name=9_diffbind
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/9_diffbind.err"
#SBATCH --output="./logs/9_diffbind.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"

echo "=========================================="
echo "Starting DiffBind analysis"
echo "Date: $(date)"
echo "=========================================="

# IMPORTANT: Cut&Tag-specific DiffBind considerations
# The R script should use the following settings for Cut&Tag data:
#
# 1. Use consensus peaks from step 11 (not raw MACS2 peaks)
#    - Preferably use: *_consensus_peaks_scored.bed
#    - Alternative: *_idr_union.bed for IDR-based peaks
#
# 2. Normalization method:
#    - Recommended: "lib" (library size) or "RLE" (DESeq2 default)
#    - Avoid: "background" normalization (Cut&Tag has low background)
#
# 3. Fragment length:
#    - Cut&Tag fragments are ~150bp (nucleosome-free)
#    - Set fragmentSize=150 in DiffBind
#
# 4. Filtering:
#    - minCount threshold should be adjusted for Cut&Tag's lower read counts
#    - Consider minCount=10 instead of default 15
#
# 5. Summit usage:
#    - summits=TRUE to use peak summits (Cut&Tag has precise binding)
#
# Example DiffBind settings for Cut&Tag:
# dba <- dba.count(dba, fragmentSize=150, summits=TRUE, minCount=10)
# dba <- dba.normalize(dba, method=DBA_NORM_LIB)  # or DBA_NORM_RLE

echo "Running DiffBind analysis (ensure R script uses Cut&Tag-optimized settings)..."
Rscript ${BASE_DIR}/scripts/diffbind_analysis.R

if [ $? -eq 0 ]; then
    echo "Differential binding analysis complete!"
    echo "Results should be in: ${BASE_DIR}/results/diffbind_analysis/"
else
    echo "ERROR: DiffBind analysis failed. Check R script and logs."
    exit 1
fi

echo "=========================================="
echo "DiffBind analysis completed"
echo "Date: $(date)"
echo "=========================================="