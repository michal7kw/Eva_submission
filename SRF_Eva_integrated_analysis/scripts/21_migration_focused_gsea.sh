#!/bin/bash
#SBATCH --job-name=a1_21_migration_gsea
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/21_migration_gsea.out
#SBATCH --error=logs/21_migration_gsea.err

# ==============================================================================
# Migration-Focused GSEA Analysis
# ==============================================================================
# This script runs GSEA focused on migration, motility, EMT, and cell adhesion
# pathways with publication-quality visualizations.
#
# Input: DESeq2 results from SRF_Eva_RNA
# Output: output/21_migration_focused_gsea/
#
# Usage: sbatch 21_migration_focused_gsea.sh
# ==============================================================================

echo "=========================================="
echo "Migration-Focused GSEA Analysis"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# Set working directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1

# Activate conda environment
echo "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Check R is available
echo "R version:"
R --version | head -1
echo ""

# Run the analysis
echo "Starting migration-focused GSEA analysis..."
echo ""

Rscript 21_migration_focused_gsea.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    echo "End time: $(date)"
    echo ""
    echo "Output files:"
    ls -la output/21_migration_focused_gsea/
    echo ""
    echo "Results:"
    ls -la output/21_migration_focused_gsea/results/
    echo ""
    echo "Plots:"
    ls -la output/21_migration_focused_gsea/plots/
else
    echo ""
    echo "=========================================="
    echo "ERROR: Analysis failed!"
    echo "=========================================="
    echo "Check logs/21_migration_gsea.err for details"
    exit 1
fi
