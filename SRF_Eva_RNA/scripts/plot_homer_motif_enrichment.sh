#!/bin/bash

#===============================================================================
# SCRIPT: plot_homer_motif_enrichment.sh
# PURPOSE: Create motif enrichment scatter plots highlighting TEAD motifs
#
# DESCRIPTION:
# Creates scatter plots showing % enrichment vs log p-value for all motifs,
# with TEAD family motifs labeled and highlighted to show their dominance
# in regulating differentially expressed genes.
#
# OUTPUTS:
# - Individual plots for upregulated, downregulated, and all DEGs
# - Combined comparison plot
# - TEAD-focused summary plot and table
#
# USAGE:
# sbatch scripts/plot_homer_motif_enrichment.sh
#===============================================================================

#SBATCH --job-name=homer_plots
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/homer_enrichment_plots.err"
#SBATCH --output="logs/homer_enrichment_plots.out"

set -e

# Base directory
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
cd "${BASE_DIR}"
mkdir -p logs

echo "=============================================="
echo "HOMER Motif Enrichment Visualization"
echo "=============================================="
echo "Date: $(date)"
echo "Working directory: ${BASE_DIR}"
echo ""

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

# Check if required packages are available
echo "Checking R packages..."
Rscript -e "
required <- c('ggplot2', 'ggrepel', 'dplyr', 'tidyr', 'viridis', 'patchwork')
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
    cat('Missing packages:', paste(missing, collapse = ', '), '\n')
    cat('Installing missing packages...\n')
    install.packages(missing, repos = 'https://cran.r-project.org', quiet = TRUE)
}
cat('All required packages available\n')
"

# Check for MASS package (for density calculation)
Rscript -e "
if (!requireNamespace('MASS', quietly = TRUE)) {
    install.packages('MASS', repos = 'https://cran.r-project.org', quiet = TRUE)
}
"

echo ""
echo "Running R script..."
echo ""

# Run the R script
Rscript "${BASE_DIR}/scripts/plot_homer_motif_enrichment.R"

EXIT_CODE=$?

echo ""
echo "=============================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Visualization complete!"
    echo ""
    echo "Output files:"
    ls -la "${BASE_DIR}/results/07_homer_motifs/plots/"
else
    echo "Visualization failed with exit code ${EXIT_CODE}"
fi
echo "=============================================="
echo "Completed: $(date)"

exit $EXIT_CODE
