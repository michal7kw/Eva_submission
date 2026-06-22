#!/bin/bash

#===============================================================================
# SCRIPT: 1_qc.sh
# PURPOSE: Quality Control Analysis for Cut&Tag Sequencing Data
#
# DESCRIPTION:
# This script performs initial quality control assessment of raw FASTQ files
# using FastQC. It is the first step in the Cut&Tag analysis pipeline and
# generates comprehensive quality reports for both forward (R1) and reverse (R2)
# reads from paired-end sequencing data.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with FastQC and quality control tools
# 2. Processes samples in parallel using SLURM job arrays
# 3. Runs FastQC analysis on both R1 and R2 FASTQ files for each sample
# 4. Generates HTML and ZIP reports containing quality metrics
#
# METHODOLOGY:
# - Uses SLURM job arrays to process multiple samples in parallel
# - Each array task processes one sample from the samples.txt configuration file
# - FastQC analyzes sequence quality, GC content, adapter content, and other metrics
# - Multi-threading (4 cores) is used to accelerate processing
#
# IMPORTANT PARAMETERS:
# - SLURM_ARRAY_TASK_ID: Array index used to select sample from samples list
# - Memory allocation: 8GB per job
# - Time limit: 2 hours per job
# - Threads: 4 cores for FastQC processing
# - Array range: 0-10 (processes 11 samples)
#
# EXPECTED INPUTS:
# - Raw FASTQ files: ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz
# - Raw FASTQ files: ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - FastQC HTML reports: ${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.html
# - FastQC HTML reports: ${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.html
# - FastQC ZIP archives: ${OUTPUT_DIR}/${SAMPLE}_R1_001_fastqc.zip
# - FastQC ZIP archives: ${OUTPUT_DIR}/${SAMPLE}_R2_001_fastqc.zip
# - Log files: ./logs/1_fastqc_${ARRAY_ID}.out and ./logs/1_fastqc_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - FastQC (via conda environment 'quality')
# - SLURM workload manager
# - Input FASTQ files in gzipped format
#
# USAGE:
# sbatch 1_qc.sh
#
# NOTE: This script should be run before any preprocessing steps to assess
# the quality of raw sequencing data and identify potential issues.
#===============================================================================

#SBATCH --job-name=1_fastqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/1_fastqc_%a.err"
#SBATCH --output="./logs/1_fastqc_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate quality

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
FASTQ_DIR="${BASE_DIR}/90-1222471453/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/01_fastqc"

# Array of samples
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: ${SAMPLE}"

# Run FastQC on both R1 and R2
fastqc -t 4 -o ${OUTPUT_DIR} \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

echo "FastQC complete for ${SAMPLE}"