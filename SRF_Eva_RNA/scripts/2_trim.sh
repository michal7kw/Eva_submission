#!/bin/bash

#===============================================================================
# SCRIPT: 2_trim.sh
# PURPOSE: Adapter Trimming for RNA-seq Sequencing Data
#
# DESCRIPTION:
# This script performs adapter trimming and quality filtering of paired-end
# FASTQ files using Trim Galore. It is the second step in the RNA-seq analysis
# pipeline and generates high-quality, adapter-free reads for alignment.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with Trim Galore and Cutadapt
# 2. Validates input FASTQ files exist
# 3. Runs Trim Galore with RNA-seq optimized parameters
# 4. Removes adapters and low-quality bases
# 5. Generates trimming reports and statistics
#
# METHODOLOGY:
# - Uses Trim Galore (wrapper for Cutadapt) for RNA-seq data
# - Automatically detects and removes Illumina adapters
# - Quality trimming with Phred score threshold of 20
# - Minimum read length filter of 20 bp
# - Paired-end mode with length validation
#
# EXPECTED INPUTS:
# - Raw FASTQ files: ${RAW_DIR}/${SAMPLE}/${SAMPLE}_R1_001.fastq.gz
# - Raw FASTQ files: ${RAW_DIR}/${SAMPLE}/${SAMPLE}_R2_001.fastq.gz
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - Trimmed FASTQ files: ${OUTPUT_DIR}/${SAMPLE}_R1_001_val_1.fq.gz
# - Trimmed FASTQ files: ${OUTPUT_DIR}/${SAMPLE}_R2_001_val_2.fq.gz
# - Trimming reports: ${OUTPUT_DIR}/${SAMPLE}_R1_001.fastq.gz_trimming_report.txt
# - SLURM logs: ./logs/2_trim_${ARRAY_ID}.out and ./logs/2_trim_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - Trim Galore (via conda environment)
# - Cutadapt (dependency of Trim Galore)
# - FastQC (for post-trimming QC)
# - SLURM workload manager
#
# USAGE:
# sbatch 2_trim.sh
#===============================================================================

#SBATCH --job-name=2_trim
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/2_trim_%a.err"
#SBATCH --output="./logs/2_trim_%a.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/trim

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
RAW_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/download/TES_RNA/TES"
OUTPUT_DIR="${BASE_DIR}/results/02_trimmed"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting adapter trimming for sample: ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input directory: ${RAW_DIR}/${SAMPLE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if input files exist
R1_FILE="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R1_001.fastq.gz"
R2_FILE="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R2_001.fastq.gz"

if [[ ! -f "${R1_FILE}" ]]; then
    echo "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi
if [[ ! -f "${R2_FILE}" ]]; then
    echo "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

echo "=== Running Trim Galore ==="
echo "Timestamp: $(date)"

# Run Trim Galore with RNA-seq optimized parameters
trim_galore \
    --paired \
    --quality 20 \
    --phred33 \
    --length 20 \
    --stringency 3 \
    --fastqc \
    --cores 4 \
    --output_dir ${OUTPUT_DIR} \
    ${R1_FILE} ${R2_FILE}

echo "=== Checking output files ==="
TRIMMED_R1="${OUTPUT_DIR}/${SAMPLE}_val_1.fq.gz"
TRIMMED_R2="${OUTPUT_DIR}/${SAMPLE}_val_2.fq.gz"

if [[ -f "${TRIMMED_R1}" && -f "${TRIMMED_R2}" ]]; then
    # Get file sizes
    R1_SIZE=$(ls -lh "${TRIMMED_R1}" | awk '{print $5}')
    R2_SIZE=$(ls -lh "${TRIMMED_R2}" | awk '{print $5}')

    echo "SUCCESS: Trimmed files created for ${SAMPLE}"
    echo "  R1: ${TRIMMED_R1} (${R1_SIZE})"
    echo "  R2: ${TRIMMED_R2} (${R2_SIZE})"

    # Count reads in trimmed files
    R1_READS=$(zcat "${TRIMMED_R1}" | wc -l | awk '{print $1/4}')
    R2_READS=$(zcat "${TRIMMED_R2}" | wc -l | awk '{print $1/4}')

    echo "  R1 reads after trimming: ${R1_READS}"
    echo "  R2 reads after trimming: ${R2_READS}"

    if [[ "${R1_READS}" != "${R2_READS}" ]]; then
        echo "WARNING: Unequal number of reads in R1 and R2 after trimming!"
    fi

else
    echo "ERROR: Trimmed files not created!"
    exit 1
fi

# Check for trimming report
REPORT_R1="${OUTPUT_DIR}/${SAMPLE}_R1_001.fastq.gz_trimming_report.txt"
REPORT_R2="${OUTPUT_DIR}/${SAMPLE}_R2_001.fastq.gz_trimming_report.txt"

if [[ -f "${REPORT_R1}" && -f "${REPORT_R2}" ]]; then
    echo "SUCCESS: Trimming reports generated"
    echo "  R1 report: ${REPORT_R1}"
    echo "  R2 report: ${REPORT_R2}"
else
    echo "WARNING: Trimming reports not found"
fi

echo "=== Trimming complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"