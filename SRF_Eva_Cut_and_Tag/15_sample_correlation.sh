#!/bin/bash

#===============================================================================
# SCRIPT: 15_sample_correlation.sh
# PURPOSE: Generate sample correlation/distance heatmap for Cut&Tag data
#
# DESCRIPTION:
# This script creates publication-ready sample correlation heatmaps showing
# sample clustering patterns similar to RNA-seq distance plots. It uses
# DiffBind to count reads in peaks and generates:
# 1. Euclidean distance heatmap with hierarchical clustering
# 2. Pearson correlation heatmap
# Both include annotation bars for condition and replicate.
#
# METHODOLOGY:
# - Uses DiffBind to count reads in consensus peak regions
# - Log2 transforms counts with pseudocount
# - Calculates Euclidean distance and Pearson correlation matrices
# - Uses pheatmap for visualization with hierarchical clustering
#
# INPUTS:
# - Filtered BAM files (results/04_filtered/)
# - Narrow peak files (results/05_peaks_narrow/)
#
# OUTPUTS:
# - sample_distance_heatmap.pdf/png (Euclidean distance)
# - sample_correlation_heatmap.pdf/png (Pearson correlation)
# - diffbind_object.RData (cached DiffBind object for reuse)
#
# DEPENDENCIES:
# - R with DiffBind, pheatmap, RColorBrewer
# - Conda environment: diffbind_analysis
#
# USAGE:
# sbatch 15_sample_correlation.sh
#===============================================================================

#SBATCH --job-name=15_sample_corr
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/15_sample_correlation.err"
#SBATCH --output="./logs/15_sample_correlation.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
cd ${BASE_DIR}

echo "=========================================="
echo "Starting Sample Correlation Analysis"
echo "Date: $(date)"
echo "=========================================="

# Create logs directory if it doesn't exist
mkdir -p ./logs

# Run the R script
echo "Running sample correlation heatmap script..."
Rscript ${BASE_DIR}/scripts/sample_correlation_heatmap.R

if [ $? -eq 0 ]; then
    echo "Sample correlation analysis complete!"
    echo "Results saved to: ${BASE_DIR}/results/07_analysis_narrow/"
    ls -la ${BASE_DIR}/results/07_analysis_narrow/sample_*_heatmap.*
else
    echo "ERROR: Sample correlation analysis failed. Check R script and logs."
    exit 1
fi

echo "=========================================="
echo "Sample Correlation Analysis completed"
echo "Date: $(date)"
echo "=========================================="
