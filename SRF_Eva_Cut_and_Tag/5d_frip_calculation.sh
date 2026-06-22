#!/bin/bash

#===============================================================================
# SCRIPT: 5d_frip_calculation.sh
# PURPOSE: Calculate Fraction of Reads in Peaks (FRiP) per sample
#
# DESCRIPTION:
# FRiP is a standard quality metric for ChIP-seq/CUT&Tag experiments measuring
# signal-to-noise ratio. It is defined as:
#   FRiP = reads_overlapping_peaks / total_mapped_reads
#
# ENCODE quality thresholds:
#   Transcription factors: FRiP > 5% is acceptable, >20% is excellent
#   CUT&Tag typically achieves higher FRiP than ChIP-seq due to lower background.
#
# USAGE:
# sbatch 5d_frip_calculation.sh
#
# INPUTS:
# - Filtered BAM files from step 4
# - NarrowPeak files from step 5 (individual replicates)
#
# OUTPUTS:
# - results/05_peaks_narrow/frip_scores.txt (per-sample FRiP values)
#===============================================================================

#SBATCH --job-name=5d_frip
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5d_frip.err"
#SBATCH --output="./logs/5d_frip.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BAM_DIR="${BASE_DIR}/results/04_filtered"
PEAK_DIR="${BASE_DIR}/results/05_peaks_narrow"
OUTPUT_FILE="${PEAK_DIR}/frip_scores.txt"

echo "=========================================="
echo "FRiP (Fraction of Reads in Peaks) Calculation"
echo "Date: $(date)"
echo "=========================================="

printf "%-15s %15s %15s %10s %10s\n" "Sample" "Total_Reads" "Reads_in_Peaks" "FRiP" "Quality" > ${OUTPUT_FILE}
printf "%-15s %15s %15s %10s %10s\n" "------" "-----------" "--------------" "----" "-------" >> ${OUTPUT_FILE}

for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_filtered.bam"
    PEAK_FILE="${PEAK_DIR}/${SAMPLE}_peaks.narrowPeak"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "WARNING: BAM not found for ${SAMPLE}, skipping"
        continue
    fi
    if [[ ! -f "${PEAK_FILE}" ]]; then
        echo "WARNING: Peak file not found for ${SAMPLE}, skipping"
        continue
    fi

    echo "Processing ${SAMPLE}..."

    # Total mapped reads (properly paired)
    TOTAL_READS=$(samtools view -c -F 4 ${BAM_FILE})

    # Reads overlapping peaks
    READS_IN_PEAKS=$(bedtools intersect -u -a ${BAM_FILE} -b ${PEAK_FILE} -ubam | samtools view -c -)

    # Calculate FRiP
    if [[ ${TOTAL_READS} -gt 0 ]]; then
        FRIP=$(awk "BEGIN {printf \"%.4f\", ${READS_IN_PEAKS}/${TOTAL_READS}}")
        FRIP_PCT=$(awk "BEGIN {printf \"%.1f%%\", ${READS_IN_PEAKS}/${TOTAL_READS}*100}")
    else
        FRIP="0.0000"
        FRIP_PCT="0.0%"
    fi

    # Quality assessment
    FRIP_NUM=$(awk "BEGIN {print ${READS_IN_PEAKS}/${TOTAL_READS}+0}")
    if (( $(awk "BEGIN {print (${FRIP_NUM} >= 0.20)}") )); then
        QUALITY="Excellent"
    elif (( $(awk "BEGIN {print (${FRIP_NUM} >= 0.05)}") )); then
        QUALITY="Good"
    elif (( $(awk "BEGIN {print (${FRIP_NUM} >= 0.01)}") )); then
        QUALITY="Acceptable"
    else
        QUALITY="Low"
    fi

    printf "%-15s %15d %15d %10s %10s\n" "${SAMPLE}" "${TOTAL_READS}" "${READS_IN_PEAKS}" "${FRIP_PCT}" "${QUALITY}" >> ${OUTPUT_FILE}
    echo "  ${SAMPLE}: FRiP = ${FRIP_PCT} (${QUALITY})"
done

echo ""
echo "=== FRiP Summary ==="
cat ${OUTPUT_FILE}

echo ""
echo "=========================================="
echo "FRiP calculation complete"
echo "Results: ${OUTPUT_FILE}"
echo "Date: $(date)"
echo "=========================================="
