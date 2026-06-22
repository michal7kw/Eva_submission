#!/bin/bash

#===============================================================================
# SCRIPT: 2_adapter_trim.sh
# PURPOSE: Adapter Trimming and Quality Filtering for Cut&Tag Sequencing Data
#
# DESCRIPTION:
# This script performs adapter trimming and quality filtering of raw FASTQ files
# using Trim Galore. It is the second step in the Cut&Tag analysis pipeline and
# removes adapter sequences, low-quality bases, and short reads to improve
# downstream alignment and analysis quality.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with Trim Galore and trimming tools
# 2. Processes paired-end FASTQ files for each sample in parallel
# 3. Removes Illumina adapter sequences automatically
# 4. Trims low-quality bases (Phred score < 20)
# 5. Filters out reads shorter than 20 bp after trimming
# 6. Generates FastQC reports for trimmed reads
# 7. Outputs gzip-compressed trimmed FASTQ files
#
# METHODOLOGY:
# - Uses Trim Galore wrapper around Cutadapt for robust adapter removal
# - Employs stringent quality trimming (Q20) suitable for Cut&Tag data
# - Maintains paired-end read relationships during trimming
# - Automatic adapter detection for standard Illumina sequencing
# - Multi-threading (8 cores) for accelerated processing
#
# IMPORTANT PARAMETERS:
# - Quality threshold: 20 (removes bases with Phred score < 20)
# - Stringency: 5 (minimum overlap for adapter detection)
# - Minimum length: 20 bp (removes shorter reads after trimming)
# - Cores: 8 (parallel processing threads)
# - Memory allocation: 16GB per job
# - Time limit: 4 hours per job
#
# EXPECTED INPUTS:
# - Raw FASTQ files: ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz
# - Raw FASTQ files: ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - Trimmed FASTQ files: ${OUTPUT_DIR}/${SAMPLE}_R1_001_val_1.fq.gz
# - Trimmed FASTQ files: ${OUTPUT_DIR}/${SAMPLE}_R2_001_val_2.fq.gz
# - Trimming reports: ${OUTPUT_DIR}/${SAMPLE}_R1_001.fastq.gz_trimming_report.txt
# - Trimming reports: ${OUTPUT_DIR}/${SAMPLE}_R2_001.fastq.gz_trimming_report.txt
# - FastQC reports: ${OUTPUT_DIR}/${SAMPLE}_R1_001_val_1_fastqc.html
# - FastQC reports: ${OUTPUT_DIR}/${SAMPLE}_R2_001_val_2_fastqc.html
# - Log files: ./logs/2_trim_${ARRAY_ID}.out and ./logs/2_trim_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - Trim Galore (via conda environment 'trim')
# - Cutadapt (dependency of Trim Galore)
# - FastQC (for post-trimming quality assessment)
# - SLURM workload manager
#
# USAGE:
# sbatch 2_adapter_trim.sh
#
# NOTE: This script should be run after quality control assessment (1_qc.sh)
# and before alignment (3_align.sh). The trimmed reads will be used for all
# downstream analyses.
#===============================================================================

#SBATCH --job-name=2_trim
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/2_trim_%a.err"
#SBATCH --output="./logs/2_trim_%a.out"

conda activate  /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/trim

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
FASTQ_DIR="${BASE_DIR}/90-1222471453/00_fastq"
OUTPUT_DIR="${BASE_DIR}/results/02_trimmed"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Trimming sample: ${SAMPLE}"

# Trim Galore for paired-end Cut&Tag data
# Using standard Illumina adapters and stringent quality trimming
trim_galore \
    --paired \
    --cores 8 \
    --quality 20 \
    --stringency 5 \
    --length 20 \
    --fastqc \
    --gzip \
    --output_dir ${OUTPUT_DIR} \
    ${FASTQ_DIR}/${SAMPLE}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${SAMPLE}_R2_001.fastq.gz

echo "Trimming complete for ${SAMPLE}"