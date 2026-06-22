#!/bin/bash
#SBATCH --job-name=go_gsea_viz
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/17_go_gsea_combined_visualization_%j.out
#SBATCH --error=logs/17_go_gsea_combined_visualization_%j.err

# GO and GSEA Combined Visualization (Simplified)

echo "GO and GSEA Combined Visualization"
echo "Start: $(date)"
echo ""

mkdir -p logs

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 17_go_gsea_combined_visualization.R

echo ""
echo "End: $(date)"
