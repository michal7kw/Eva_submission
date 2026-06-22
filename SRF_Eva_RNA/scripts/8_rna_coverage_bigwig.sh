#!/bin/bash

#===============================================================================
# SCRIPT: 8_rna_coverage_bigwig.sh
# PURPOSE: Generate normalized RNA-seq coverage BigWig tracks for genome-browser
#          / pyGenomeTracks visualization (e.g. the TEAD1-locus snapshot).
#
# DESCRIPTION:
# The RNA-seq pipeline (steps 1-5) produces STAR BAMs but no coverage tracks.
# This step converts each aligned BAM to a CPM-normalized BigWig so the 6 samples
# (Mock = GFP, and TES; 3 replicates each) can be compared on a shared y-axis.
#
# NOTE (vs Cut&Tag 6_bigwig.sh): RNA reads are spliced, so we do NOT use
# --extendReads (extending across introns would misrepresent coverage). CPM
# normalization is used so tracks are directly comparable between samples.
#
# EXPECTED INPUTS:
# - Sorted, indexed BAMs: results/03_aligned/${SAMPLE}_sorted.bam (+ .bai)
# - Sample list:          config/samples.txt (6 samples -> array 0-5)
#
# EXPECTED OUTPUTS:
# - results/08_bigwig/${SAMPLE}.bw  (CPM-normalized coverage)
#
# DEPENDENCIES:
# - deepTools (bamCoverage), conda environment: bigwig
#
# USAGE:
# sbatch 8_rna_coverage_bigwig.sh
#===============================================================================

#SBATCH --job-name=8_rna_bigwig
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8_rna_bigwig_%a.err"
#SBATCH --output="./logs/8_rna_bigwig_%a.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate bigwig

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
BAM_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/08_bigwig"

mkdir -p "${OUTPUT_DIR}"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

BAM="${BAM_DIR}/${SAMPLE}_sorted.bam"

echo "=== RNA-seq BigWig for ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input BAM: ${BAM}"

if [[ ! -f "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}"
    exit 1
fi

# Ensure the BAM is indexed (3_align.sh copies the .bai, but be safe)
if [[ ! -f "${BAM}.bai" ]]; then
    echo "Index not found, creating with samtools index..."
    samtools index -@ 8 "${BAM}"
fi

echo "Generating CPM-normalized BigWig..."
bamCoverage \
    -b "${BAM}" \
    -o "${OUTPUT_DIR}/${SAMPLE}.bw" \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8

if [[ -f "${OUTPUT_DIR}/${SAMPLE}.bw" ]]; then
    echo "SUCCESS: ${OUTPUT_DIR}/${SAMPLE}.bw"
else
    echo "ERROR: BigWig not created for ${SAMPLE}"
    exit 1
fi

echo "=== Done for ${SAMPLE} ==="
echo "Timestamp: $(date)"
