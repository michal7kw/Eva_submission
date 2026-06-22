#!/bin/bash

#===============================================================================
# SCRIPT: 1_fastqc.sh
# PURPOSE: Quality Control for RNA-seq Sequencing Data
#
# DESCRIPTION:
# This script performs quality control analysis of raw paired-end FASTQ files
# using FastQC. It is the first step in the RNA-seq analysis pipeline and
# generates comprehensive quality reports for all sequencing samples.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with FastQC
# 2. Validates input FASTQ files exist
# 3. Runs FastQC on both R1 and R2 reads
# 4. Generates HTML and ZIP quality reports
# 5. Creates a MultiQC summary report (optional)
#
# METHODOLOGY:
# - Processes both forward (R1) and reverse (R2) reads
# - Uses multi-threading for faster processing
# - Generates per-base quality scores, GC content, adapter content analysis
# - Provides overrepresented sequences detection
#
# EXPECTED INPUTS:
# - Raw FASTQ files: ${RAW_DIR}/${SAMPLE}/${SAMPLE}_R1_001.fastq.gz
# - Raw FASTQ files: ${RAW_DIR}/${SAMPLE}/${SAMPLE}_R2_001.fastq.gz
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - FastQC reports: ${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html/zip
# - FastQC reports: ${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html/zip
# - SLURM logs: ./logs/1_fastqc_${ARRAY_ID}.out and ./logs/1_fastqc_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - FastQC (via conda environment)
# - SLURM workload manager
#
# USAGE:
# sbatch 1_fastqc.sh
#===============================================================================

#SBATCH --job-name=1_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/1_fastqc_%a.err"
#SBATCH --output="./logs/1_fastqc_%a.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/quality

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
RAW_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/download/TES_RNA/TES"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting FastQC for sample: ${SAMPLE} ==="
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

echo "=== Running FastQC ==="
echo "Timestamp: $(date)"

# Run FastQC on both R1 and R2 files
fastqc \
    --threads 8 \
    --outdir ${OUTPUT_DIR} \
    --format fastq \
    ${R1_FILE} ${R2_FILE}

echo "=== Checking output files ==="
if [[ -f "${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html" && -f "${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html" ]]; then
    echo "SUCCESS: FastQC reports generated for ${SAMPLE}"
    echo "  R1 report: ${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html"
    echo "  R2 report: ${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html"
else
    echo "ERROR: FastQC reports not generated!"
    exit 1
fi

echo "=== FastQC complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"