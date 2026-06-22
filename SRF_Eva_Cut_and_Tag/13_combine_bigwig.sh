#!/bin/bash

#===============================================================================
# SCRIPT: 13_combine_bigwig.sh
# PURPOSE: Combine replicate BigWig files into group-level averages.
#
# **IMPORTANT NOTE**: This script is now OPTIONAL/LEGACY
# Step 11 (11_combine_replicates_narrow.sh) already creates average BigWig
# tracks using bamCoverage on merged BAM files, which is the preferred method
# for Cut&Tag data. This script creates averages from individual BigWig files,
# which is an alternative approach but less optimal.
#
# **RECOMMENDATION**: Use the BigWig files from step 11 instead:
#   results/11_combined_replicates_narrow/bigwig/*_average.bw
#
# DESCRIPTION:
# This script takes multiple BigWig files from individual replicates of a
# ChIP-seq or Cut&Tag experiment and computes an average signal track for each
# experimental group. The output is a single, representative BigWig file
# per group, which is useful for visualization in a genome browser.
#
# KEY OPERATIONS:
# 1. Activate a conda environment containing the necessary tools (deepTools).
# 2. For each experimental group (TES, TEAD1), it identifies the
#    replicate BigWig files.
# 3. It uses the 'bigwigAverage' tool from the deepTools suite to calculate
#    the mean signal across all replicates for that group.
# 4. The resulting average signal is saved as a new, combined BigWig file.
#
# METHODOLOGY:
# The script employs 'bigwigAverage' which calculates the average score
# (e.g., read coverage) for genomic bins of a specified size. This approach
# smooths out replicate-specific noise and provides a robust representation
# of the signal for an entire experimental condition.
#
# IMPORTANT PARAMETERS:
# - --bigwigs: A space-separated list of input BigWig files to be averaged.
# - --outFileName: The path for the output combined BigWig file.
# - --binSize: The size of the genomic window (in base pairs) to average
#   the signal over. A smaller bin size provides higher resolution.
# - --numberOfProcessors: The number of CPU cores to use for the computation.
# - SLURM parameters (--mem, --time, etc.) define computational resources.
#
# INPUTS:
# - Individual replicate BigWig files located in ${BASE_DIR}/results/06_bigwig/.
#   - e.g., TES-1_CPM.bw, TES-2_CPM.bw, TES-3_CPM.bw
#
# OUTPUTS:
# - Combined, group-level BigWig files in ${BASE_DIR}/results/06_bigwig/.
#   - TES_comb.bw
#   - TEAD1_comb.bw
#   NOTE: TESmut removed - failed sample
#
# DEPENDENCIES:
# - deepTools (specifically the 'bigwigAverage' command)
# - A conda environment named 'bigwig' where deepTools is installed.
#===============================================================================
#SBATCH --job-name=13_combine_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/13_combine_bigwig.err"
#SBATCH --output="./logs/13_combine_bigwig.out"

# Set up conda environment with deepTools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BIGWIG_DIR="${BASE_DIR}/results/06_bigwig"
OUTPUT_DIR="${BASE_DIR}/results/06_bigwig"

echo "=========================================="
echo "WARNING: This script is OPTIONAL/LEGACY"
echo "=========================================="
echo "Step 11 already creates average BigWig files from merged BAMs."
echo "Those files are preferred for Cut&Tag data analysis."
echo ""
echo "Step 11 output: results/11_combined_replicates_narrow/bigwig/*_average.bw"
echo ""
echo "This script creates averages from individual BigWig files (alternative method)."
echo "Continuing anyway..."
echo ""

echo "Starting bigwig combination for TES and TEAD1 groups (TESmut removed - failed sample)"

# Create combined TES bigwig (average of TES-1, TES-2, TES-3)
echo "Combining TES replicates..."
bigwigAverage \
    --bigwigs ${BIGWIG_DIR}/TES-1_CPM.bw ${BIGWIG_DIR}/TES-2_CPM.bw ${BIGWIG_DIR}/TES-3_CPM.bw \
    --outFileName ${OUTPUT_DIR}/TES_comb.bw \
    --binSize 10 \
    --numberOfProcessors 8

# TESmut removed - failed sample

# Create combined TEAD1 bigwig (average of TEAD1-1, TEAD1-2, TEAD1-3)
echo "Combining TEAD1 replicates..."
bigwigAverage \
    --bigwigs ${BIGWIG_DIR}/TEAD1-1_CPM.bw ${BIGWIG_DIR}/TEAD1-2_CPM.bw ${BIGWIG_DIR}/TEAD1-3_CPM.bw \
    --outFileName ${OUTPUT_DIR}/TEAD1_comb.bw \
    --binSize 10 \
    --numberOfProcessors 8

echo "BigWig combination complete. Created:"
echo "  - ${OUTPUT_DIR}/TES_comb.bw"
echo "  - ${OUTPUT_DIR}/TEAD1_comb.bw"
