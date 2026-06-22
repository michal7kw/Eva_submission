#!/bin/bash

#===============================================================================
# SCRIPT: 6_bigwig.sh
# PURPOSE: Generate normalized BigWig tracks for genome browser visualization
#
# DESCRIPTION:
# This script converts filtered BAM files to normalized BigWig format for
# visualization in genome browsers (IGV, UCSC Genome Browser, etc.). BigWig
# files provide continuous coverage tracks that show signal intensity across
# the genome, enabling visual inspection of Cut&Tag binding patterns and
# comparison between samples.
#
# KEY OPERATIONS:
# 1. Array job processing for parallel BigWig generation
# 2. CPM (Counts Per Million) normalization for sample comparison
# 3. Read extension for better signal representation
# 4. High-resolution binning (10bp) for detailed visualization
# 5. Multi-threaded processing for efficiency
#
# METHODOLOGY:
# - Uses deepTools bamCoverage for BigWig generation
# - Applies CPM normalization to enable cross-sample comparison
# - Extends reads to estimated fragment length for better signal
# - Uses 10bp bins for high-resolution tracks
# - Processes samples in parallel using SLURM array jobs
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB per job (sufficient for BigWig generation)
# - Time: 2 hours (adequate for individual samples)
# - Threads: 8 per job (bamCoverage supports multi-threading)
# - Array: 0-10 (processes 11 samples in parallel)
# - Bin size: 10bp (high resolution for detailed visualization)
# - Normalization: CPM (Counts Per Million for cross-sample comparison)
#
# EXPECTED INPUTS:
# - Filtered BAM files from step 4 (04_filtered directory)
# - samples.txt configuration file listing all sample names
# - BAM files should be properly indexed
#
# EXPECTED OUTPUTS:
# - *_CPM.bw: Normalized BigWig files for each sample
# - Files suitable for loading into genome browsers
# - Tracks show normalized signal intensity across the genome
#
# DEPENDENCIES:
# - deepTools (bamCoverage tool)
# - Conda environment: bigwig
# - samples.txt configuration file
#
# USAGE:
# sbatch 6_bigwig.sh
#
# NOTES:
# - Uses array jobs for efficient parallel processing
# - CPM normalization enables direct comparison between samples
# - Read extension improves signal representation for Cut&Tag data
# - High-resolution bins provide detailed visualization
# - Output files can be directly loaded into genome browsers
#===============================================================================

#SBATCH --job-name=6_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/6_bigwig_%a.err"
#SBATCH --output="./logs/6_bigwig_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/06_bigwig"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Creating BigWig for ${SAMPLE}"

# Check fragment size first
echo "Analyzing fragment size for ${SAMPLE}..."
FRAGMENT_SIZE=$(samtools view -f 2 ${BAM_DIR}/${SAMPLE}_filtered.bam | head -10000 | awk '$9 > 0 {sum+=$9; count++} END {print int(sum/count)}')
echo "Estimated fragment size: ${FRAGMENT_SIZE}bp"

# Cut&Tag-optimized BigWig generation
# Key differences from ChIP-seq:
# - NO --centerReads (Cut&Tag fragments represent actual binding, not extended reads)
# - NO --smoothLength (Cut&Tag has precise signal, don't over-smooth)
# - Use --extendReads for better signal representation
# - Keep binSize at 10bp for high resolution

echo "Generating normalized BigWig files (CPM normalization)..."
bamCoverage \
    -b ${BAM_DIR}/${SAMPLE}_filtered.bam \
    -o ${OUTPUT_DIR}/${SAMPLE}_CPM.bw \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8 \
    --extendReads

# Also create RPGC normalized version (recommended for Cut&Tag cross-sample comparison)
echo "Generating RPGC normalized BigWig..."
bamCoverage \
    -b ${BAM_DIR}/${SAMPLE}_filtered.bam \
    -o ${OUTPUT_DIR}/${SAMPLE}_RPGC.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398 \
    --binSize 10 \
    --numberOfProcessors 8 \
    --extendReads

echo "BigWig generation complete for ${SAMPLE}"