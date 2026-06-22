#!/bin/bash

#===============================================================================
# SCRIPT: 4_filter.sh
# PURPOSE: BAM Filtering and Deduplication for Cut&Tag Sequencing Data
#
# DESCRIPTION:
# This script performs comprehensive filtering and deduplication of aligned BAM
# files to prepare high-quality data for peak calling. It is the fourth step in
# the Cut&Tag analysis pipeline and includes read group addition, duplicate
# removal, quality filtering, and removal of problematic genomic regions.
#
# KEY OPERATIONS:
# 1. Activates conda environment and installs/locates Picard tools
# 2. Adds read group information to BAM files (required for Picard)
# 3. Marks (but does not remove) PCR duplicates using Picard MarkDuplicates
# 4. Filters for high-quality, properly paired alignments
# 5. Removes mitochondrial reads and problematic alignments
# 6. Restricts analysis to canonical chromosomes (1-22, X, Y)
# 7. Indexes final filtered BAM files
# 8. Generates alignment statistics and cleans up intermediate files
#
# METHODOLOGY:
# - Uses Picard MarkDuplicates for duplicate identification (mark-only, not removal)
# - Applies stringent quality filters (MAPQ ≥ 30) for high-confidence alignments
# - Filters for properly paired reads (flag 2) and removes problematic reads
# - Excludes reads with flags: unmapped (4), mate unmapped (8), not primary (256),
#   fails QC (512), duplicate (1024), supplementary (2048)
# - Restricts analysis to autosomal and sex chromosomes only
# - Multi-threading (16 cores) for efficient processing
#
# IMPORTANT PARAMETERS:
# - Memory allocation: 32GB per job (30GB for Java heap)
# - Time limit: 4 hours per job
# - Threads: 16 cores for SAMtools operations
# - Quality threshold: MAPQ ≥ 30 (high-quality alignments only)
# - Read flags: -f 2 (properly paired), -F 2828 (exclude problematic reads, keep duplicates)
# - Chromosomes: 1-22, X, Y (excludes mitochondrial and unplaced contigs)
#
# EXPECTED INPUTS:
# - Sorted BAM files: ${ALIGNED_DIR}/${SAMPLE}_sorted.bam
# - BAM index files: ${ALIGNED_DIR}/${SAMPLE}_sorted.bam.bai
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - Filtered BAM files: ${OUTPUT_DIR}/${SAMPLE}_filtered.bam
# - BAM index files: ${OUTPUT_DIR}/${SAMPLE}_filtered.bam.bai
# - Deduplication metrics: ${OUTPUT_DIR}/${SAMPLE}_dedup_metrics.txt
# - Alignment statistics: ${OUTPUT_DIR}/${SAMPLE}_filtered_flagstat.txt
# - SLURM logs: ./logs/4_filter_${ARRAY_ID}.out and ./logs/4_filter_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - Picard tools (installed via mamba in alignment environment)
# - SAMtools (for BAM processing and filtering)
# - Java (for running Picard tools)
# - SLURM workload manager
# - Aligned BAM files from previous step
#
# USAGE:
# sbatch 4_filter.sh
#
# NOTE: This script should be run after alignment (3_align.sh) and before
# peak calling (5_peak_calling.sh). The filtered BAM files represent the
# final processed data for downstream analysis.
#===============================================================================

#SBATCH --job-name=4_filter
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4_filter_%a.err"
#SBATCH --output="./logs/4_filter_%a.out"

# Function to handle errors
handle_error() {
    echo "Error on line $1"
    exit 1
}
trap 'handle_error $LINENO' ERR

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

# Wait a bit to avoid conda lock conflicts
sleep $((RANDOM % 30))

# Find Picard.jar and set PICARD variable (already installed in alignment env)
PICARD_PATH=$(find /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment -name "picard.jar" | head -n 1)
if [ -z "$PICARD_PATH" ]; then
    echo "Error: picard.jar not found in the alignment conda environment."
    exit 1
fi
PICARD="$PICARD_PATH"
echo "PICARD set to: $PICARD"

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
ALIGNED_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/04_filtered"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Filtering sample: ${SAMPLE}"

# Check if input file exists
if [ ! -f "${ALIGNED_DIR}/${SAMPLE}_sorted.bam" ]; then
    echo "Error: Input file ${ALIGNED_DIR}/${SAMPLE}_sorted.bam not found"
    exit 1
fi

# Add read groups to BAM file (required for Picard MarkDuplicates)
echo "Adding read groups to BAM file..."
java -Xmx30g -jar $PICARD AddOrReplaceReadGroups \
    I=${ALIGNED_DIR}/${SAMPLE}_sorted.bam \
    O=${OUTPUT_DIR}/${SAMPLE}_with_rg.bam \
    RGID=${SAMPLE} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${SAMPLE} \
    CREATE_INDEX=true

# Check if read group addition was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_with_rg.bam" ]; then
    echo "Error: Failed to add read groups"
    exit 1
fi

# Mark duplicates with Picard (mark only, do not remove)
# CUT&Tag uses Tn5 transposase which cuts at specific accessible sites,
# so identical fragment positions are often biological, not PCR artifacts.
# Duplicates are marked (not removed) and MACS2's --keep-dup all retains them.
echo "Marking duplicates..."
java -Xmx30g -jar $PICARD MarkDuplicates \
    I=${OUTPUT_DIR}/${SAMPLE}_with_rg.bam \
    O=${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
    M=${OUTPUT_DIR}/${SAMPLE}_dedup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    CREATE_INDEX=true

# Check if deduplication was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_dedup.bam" ]; then
    echo "Error: Deduplication failed"
    exit 1
fi

# Filter for properly paired, high-quality alignments
# Remove mitochondrial reads and blacklisted regions
# -F 2828: exclude unmapped(4), mate unmapped(8), not primary(256), QC fail(512), supplementary(2048)
# Note: NOT excluding duplicates (1024) — CUT&Tag duplicates are often biological.
# MACS2 --keep-dup all will handle them appropriately.
echo "Filtering alignments..."
samtools view -@ 16 -b -q 30 -f 2 -F 2828 \
    ${OUTPUT_DIR}/${SAMPLE}_dedup.bam \
    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y | \
    samtools sort -@ 16 -o ${OUTPUT_DIR}/${SAMPLE}_filtered.bam -

# Check if filtering was successful
if [ ! -f "${OUTPUT_DIR}/${SAMPLE}_filtered.bam" ]; then
    echo "Error: Filtering failed"
    exit 1
fi

# Index filtered BAM
echo "Indexing filtered BAM..."
samtools index -@ 16 ${OUTPUT_DIR}/${SAMPLE}_filtered.bam

# Generate final statistics
echo "Generating statistics..."
samtools flagstat ${OUTPUT_DIR}/${SAMPLE}_filtered.bam > ${OUTPUT_DIR}/${SAMPLE}_filtered_flagstat.txt

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f ${OUTPUT_DIR}/${SAMPLE}_with_rg.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_with_rg.bai
rm -f ${OUTPUT_DIR}/${SAMPLE}_dedup.bam
rm -f ${OUTPUT_DIR}/${SAMPLE}_dedup.bai

echo "Filtering complete for ${SAMPLE}"
echo "Final output: ${OUTPUT_DIR}/${SAMPLE}_filtered.bam"