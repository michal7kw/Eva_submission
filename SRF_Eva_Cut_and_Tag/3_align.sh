#!/bin/bash

#===============================================================================
# SCRIPT: 3_align.sh
# PURPOSE: Read Alignment for Cut&Tag Sequencing Data
#
# DESCRIPTION:
# This script performs alignment of trimmed paired-end FASTQ files to the
# reference genome using Bowtie2. It is the third step in the Cut&Tag analysis
# pipeline and generates sorted, indexed BAM files with comprehensive alignment
# statistics. The script uses parameters optimized for Cut&Tag data characteristics.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with alignment tools (Bowtie2, SAMtools)
# 2. Validates input trimmed FASTQ files exist
# 3. Performs Bowtie2 alignment with Cut&Tag-optimized parameters
# 4. Converts SAM output to BAM format using SAMtools
# 5. Sorts BAM files by genomic coordinates
# 6. Indexes sorted BAM files for downstream analysis
# 7. Generates alignment statistics using flagstat
# 8. Validates output files and reports read counts
#
# METHODOLOGY:
# - Uses Bowtie2 with --very-sensitive-local for Cut&Tag (soft-clips adapter remnants)
# - Enables --dovetail for short fragments where mates extend past each other
# - Enforces proper paired-end alignment (--no-mixed --no-discordant)
# - Sets insert size range (10-700 bp) appropriate for Cut&Tag fragments
# - Includes read group information for sample tracking
# - Employs multi-threading (32 cores) for optimal performance
# - Pipes alignment directly to sorting to minimize I/O overhead
#
# IMPORTANT PARAMETERS:
# - Alignment sensitivity: --very-sensitive-local (local alignment, soft-clips adapters)
# - Insert size range: -I 10 -X 700 (minimum 10bp, maximum 700bp)
# - Quality encoding: --phred33 (standard Illumina encoding)
# - Threads: 32 cores for Bowtie2 and SAMtools operations
# - Memory allocation: 64GB per job
# - Time limit: 8 hours per job
# - Proper pairs only: --no-mixed --no-discordant
#
# EXPECTED INPUTS:
# - Trimmed FASTQ files: ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz
# - Trimmed FASTQ files: ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz
# - Reference genome index: ${GENOME_INDEX} (Bowtie2 index files)
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - Sorted BAM files: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam
# - BAM index files: ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai
# - Alignment logs: ${OUTPUT_DIR}/${SAMPLE}_bowtie2.log
# - Alignment statistics: ${OUTPUT_DIR}/${SAMPLE}_flagstat.txt
# - SLURM logs: ./logs/3_align_${ARRAY_ID}.out and ./logs/3_align_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - Bowtie2 (via conda environment 'alignment')
# - SAMtools (for BAM processing and indexing)
# - Reference genome index (GRCh38)
# - SLURM workload manager
# - Trimmed FASTQ files from previous step
#
# USAGE:
# sbatch 3_align.sh
#
# NOTE: This script should be run after adapter trimming (2_adapter_trim.sh)
# and before BAM filtering (4_filter.sh). The resulting BAM files will be
# used for all downstream peak calling and analysis steps.
#===============================================================================

#SBATCH --job-name=3_align
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-7
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3_align_%a.err"
#SBATCH --output="./logs/3_align_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
TRIMMED_DIR="${BASE_DIR}/results/02_trimmed"
OUTPUT_DIR="${BASE_DIR}/results/03_aligned"

# Reference genome - adjust path as needed
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/GRCh38"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting alignment for sample: ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input files:"
echo "  R1: ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
echo "  R2: ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"
echo "Output directory: ${OUTPUT_DIR}"
echo "Genome index: ${GENOME_INDEX}"
echo ""

# Check if input files exist
if [[ ! -f "${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz" ]]; then
    echo "ERROR: R1 file not found: ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
    exit 1
fi
if [[ ! -f "${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz" ]]; then
    echo "ERROR: R2 file not found: ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"
    exit 1
fi

echo "=== Step 1: Running Bowtie2 alignment ==="
echo "Timestamp: $(date)"
# Bowtie2 alignment with parameters optimized for Cut&Tag (Henikoff protocol)
# --very-sensitive-local: local alignment soft-clips residual Tn5 adapter sequences
# --dovetail: allows mates to extend past each other (common with short CUT&Tag fragments)
# --no-mixed --no-discordant: require proper pairs only
bowtie2 \
    -p 32 \
    --very-sensitive-local \
    --dovetail \
    --no-mixed \
    --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    --rg-id "${SAMPLE}" \
    --rg "SM:${SAMPLE}" \
    --rg "LB:${SAMPLE}" \
    --rg "PL:ILLUMINA" \
    --rg "PU:${SAMPLE}" \
    -x ${GENOME_INDEX} \
    -1 ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz \
    -2 ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz \
    2> ${OUTPUT_DIR}/${SAMPLE}_bowtie2.log | \
    samtools view -@ 32 -bS - | \
    samtools sort -@ 32 -o ${OUTPUT_DIR}/${SAMPLE}_sorted.bam -

echo "=== Step 2: Indexing BAM file ==="
echo "Timestamp: $(date)"
# Index the BAM file
samtools index -@ 32 ${OUTPUT_DIR}/${SAMPLE}_sorted.bam

echo "=== Step 3: Generating alignment statistics ==="
echo "Timestamp: $(date)"
# Generate alignment statistics
samtools flagstat ${OUTPUT_DIR}/${SAMPLE}_sorted.bam > ${OUTPUT_DIR}/${SAMPLE}_flagstat.txt

echo "=== Step 4: Checking output file ==="
if [[ -f "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" ]]; then
    BAM_SIZE=$(ls -lh "${OUTPUT_DIR}/${SAMPLE}_sorted.bam" | awk '{print $5}')
    echo "SUCCESS: BAM file created - ${OUTPUT_DIR}/${SAMPLE}_sorted.bam (${BAM_SIZE})"
    
    # Quick read count check
    READ_COUNT=$(samtools view -c "${OUTPUT_DIR}/${SAMPLE}_sorted.bam")
    echo "Total reads in BAM: ${READ_COUNT}"
else
    echo "ERROR: BAM file not created!"
    exit 1
fi

echo "=== Alignment complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"