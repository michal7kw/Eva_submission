#!/bin/bash
#SBATCH --job-name=reactome_overview
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --output=logs/19_reactome_overview_plot_%j.out
#SBATCH --error=logs/19_reactome_overview_plot_%j.err

# ============================================================================
# 19_reactome_overview_plot.sh
# Creates GSEA overview bar chart for Reactome pathways
# ============================================================================

echo "=========================================="
echo "Reactome Pathway GSEA Overview Plot"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

# Setup environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Create logs directory if needed
mkdir -p logs

# Change to script directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1

# Run R script
echo "Running R script..."
Rscript 19_reactome_overview_plot.R

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS"
    echo "=========================================="
    echo "Output files:"
    ls -la output/12_msigdb_by_collection/06_Reactome_Pathways/plot_05_gsea_overview.*
else
    echo ""
    echo "=========================================="
    echo "FAILED"
    echo "=========================================="
    exit 1
fi

echo ""
echo "End: $(date)"
