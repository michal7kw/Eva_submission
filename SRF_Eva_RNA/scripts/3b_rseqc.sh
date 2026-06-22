#!/bin/bash

#===============================================================================
# SCRIPT: 3b_rseqc.sh
# PURPOSE: RNA-seq QC metrics using RSeQC and related tools
#
# DESCRIPTION:
# Generates comprehensive RNA-seq quality metrics that reviewers increasingly
# expect: gene body coverage, read distribution (exonic/intronic), 5'/3' bias,
# inner distance, strandedness validation, and junction saturation.
#
# USAGE:
# sbatch scripts/3b_rseqc.sh
#
# INPUTS:
# - Sorted BAM files from STAR alignment (results/03_aligned/)
# - BED12 annotation file (gene model)
#
# OUTPUTS:
# - results/03b_rseqc/ (per-sample QC reports)
# - results/03b_rseqc/summary/ (aggregated summaries)
#===============================================================================

#SBATCH --job-name=3b_rseqc
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3b_rseqc_%a.err"
#SBATCH --output="./logs/3b_rseqc_%a.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate quality

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
ALIGNED_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/03b_rseqc"
GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# BED12 annotation (convert from GTF if needed)
BED_ANNOTATION="${BASE_DIR}/results/03b_rseqc/gencode.v44.bed12"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/summary
mkdir -p ${BASE_DIR}/logs

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Find BAM file (STAR output naming convention)
BAM_FILE="${ALIGNED_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
if [[ ! -f "${BAM_FILE}" ]]; then
    # Try alternative naming
    BAM_FILE="${ALIGNED_DIR}/${SAMPLE}_sorted.bam"
fi

if [[ ! -f "${BAM_FILE}" ]]; then
    echo "ERROR: BAM file not found for ${SAMPLE}"
    echo "  Tried: ${ALIGNED_DIR}/${SAMPLE}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    echo "  Tried: ${ALIGNED_DIR}/${SAMPLE}_sorted.bam"
    exit 1
fi

echo "=== RSeQC Quality Control: ${SAMPLE} ==="
echo "BAM file: ${BAM_FILE}"
echo "Timestamp: $(date)"

#===============================================================================
# Step 0: Convert GTF to BED12 if needed (only first array job does this)
#===============================================================================
if [[ ! -f "${BED_ANNOTATION}" ]]; then
    echo "Converting GTF to BED12 format..."
    # Check if gtfToGenePred and genePredToBed are available
    if command -v gtfToGenePred &> /dev/null && command -v genePredToBed &> /dev/null; then
        gtfToGenePred ${GTF} /dev/stdout | genePredToBed /dev/stdin ${BED_ANNOTATION}
    else
        # Alternative: use awk-based conversion or install UCSC tools
        echo "UCSC tools (gtfToGenePred) not found."
        echo "Install via: conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed"
        echo "Or download from: https://hgdownload.soe.ucsc.edu/admin/exe/"
        echo ""
        echo "Attempting to proceed with limited QC (strandedness only)..."
    fi
fi

# Index BAM if needed
if [[ ! -f "${BAM_FILE}.bai" ]]; then
    echo "Indexing BAM file..."
    samtools index -@ 8 ${BAM_FILE}
fi

#===============================================================================
# Step 1: Infer strandedness (infer_experiment.py)
#===============================================================================
echo ""
echo "Step 1: Inferring strandedness..."
if [[ -f "${BED_ANNOTATION}" ]]; then
    infer_experiment.py \
        -i ${BAM_FILE} \
        -r ${BED_ANNOTATION} \
        -s 1000000 \
        > ${OUTPUT_DIR}/${SAMPLE}_strandedness.txt 2>&1
    cat ${OUTPUT_DIR}/${SAMPLE}_strandedness.txt
else
    echo "  Skipped: BED12 annotation not available"
fi

#===============================================================================
# Step 2: Read distribution (read_distribution.py)
#===============================================================================
echo ""
echo "Step 2: Analyzing read distribution..."
if [[ -f "${BED_ANNOTATION}" ]]; then
    read_distribution.py \
        -i ${BAM_FILE} \
        -r ${BED_ANNOTATION} \
        > ${OUTPUT_DIR}/${SAMPLE}_read_distribution.txt 2>&1
    cat ${OUTPUT_DIR}/${SAMPLE}_read_distribution.txt
else
    echo "  Skipped: BED12 annotation not available"
fi

#===============================================================================
# Step 3: Gene body coverage (geneBody_coverage.py)
#===============================================================================
echo ""
echo "Step 3: Gene body coverage analysis..."
if [[ -f "${BED_ANNOTATION}" ]]; then
    geneBody_coverage.py \
        -i ${BAM_FILE} \
        -r ${BED_ANNOTATION} \
        -o ${OUTPUT_DIR}/${SAMPLE}_genebody
    echo "  Gene body coverage plot saved"
else
    echo "  Skipped: BED12 annotation not available"
fi

#===============================================================================
# Step 4: Inner distance (inner_distance.py)
#===============================================================================
echo ""
echo "Step 4: Inner distance analysis..."
if [[ -f "${BED_ANNOTATION}" ]]; then
    inner_distance.py \
        -i ${BAM_FILE} \
        -r ${BED_ANNOTATION} \
        -o ${OUTPUT_DIR}/${SAMPLE}_inner_distance
    echo "  Inner distance plot saved"
else
    echo "  Skipped: BED12 annotation not available"
fi

#===============================================================================
# Step 5: Junction saturation (junction_saturation.py)
#===============================================================================
echo ""
echo "Step 5: Junction saturation analysis..."
if [[ -f "${BED_ANNOTATION}" ]]; then
    junction_saturation.py \
        -i ${BAM_FILE} \
        -r ${BED_ANNOTATION} \
        -o ${OUTPUT_DIR}/${SAMPLE}_junction_saturation
    echo "  Junction saturation plot saved"
else
    echo "  Skipped: BED12 annotation not available"
fi

#===============================================================================
# Step 6: BAM stat (bam_stat.py)
#===============================================================================
echo ""
echo "Step 6: BAM statistics..."
bam_stat.py \
    -i ${BAM_FILE} \
    > ${OUTPUT_DIR}/${SAMPLE}_bam_stat.txt 2>&1
cat ${OUTPUT_DIR}/${SAMPLE}_bam_stat.txt

echo ""
echo "=== RSeQC QC complete for ${SAMPLE} ==="
echo "Results: ${OUTPUT_DIR}/${SAMPLE}_*"
echo "Timestamp: $(date)"
