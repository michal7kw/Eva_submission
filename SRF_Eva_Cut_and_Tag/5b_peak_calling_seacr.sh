#!/bin/bash

#===============================================================================
# SCRIPT: 5b_peak_calling_seacr.sh
# PURPOSE: Complementary peak calling using SEACR (CUT&Tag-optimized)
#
# DESCRIPTION:
# SEACR (Sparse Enrichment Analysis for CUT&RUN) is specifically designed for
# the sparse background of CUT&Tag data. This provides complementary peak calls
# to MACS2 (5_peak_calling_narrow.sh). Peaks found by both callers represent
# the highest-confidence binding sites.
#
# SEACR uses the IgG control signal distribution directly to set an empirical
# threshold, rather than modeling background from the treatment sample.
#
# USAGE:
# sbatch 5b_peak_calling_seacr.sh
#
# INPUTS:
# - Filtered BAM files from step 4
# - IgG control BAM files
#
# OUTPUTS:
# - results/05_peaks_seacr/*.stringent.bed (stringent threshold)
# - results/05_peaks_seacr/*.relaxed.bed (relaxed threshold)
# - results/05_peaks_seacr/seacr_vs_macs2_comparison.txt
#===============================================================================

#SBATCH --job-name=5b_seacr
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5b_seacr.err"
#SBATCH --output="./logs/5b_seacr.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/05_peaks_seacr"
MACS2_DIR="${BASE_DIR}/results/05_peaks_narrow"

mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "SEACR Peak Calling for CUT&Tag Data"
echo "Date: $(date)"
echo "=========================================="

#===============================================================================
# Check if SEACR is available
#===============================================================================
SEACR_PATH=$(which SEACR_1.3.sh 2>/dev/null || which SEACR.sh 2>/dev/null || echo "")

if [[ -z "${SEACR_PATH}" ]]; then
    echo "SEACR not found in PATH."
    echo "Install: pip install seacr or download from https://github.com/FredHutch/SEACR"
    echo ""
    echo "To install SEACR:"
    echo "  git clone https://github.com/FredHutch/SEACR.git"
    echo "  chmod +x SEACR/SEACR_1.3.sh"
    echo "  export PATH=\$PATH:/path/to/SEACR"
    echo ""
    echo "Skipping SEACR peak calling. MACS2 results remain the primary peak calls."
    echo "Justify in Methods: MACS2 with BAMPE mode and IgG control is a well-validated"
    echo "approach for CUT&Tag data, used by numerous published studies."
    exit 0
fi

echo "Using SEACR: ${SEACR_PATH}"

#===============================================================================
# Convert BAM to bedgraph (SEACR input format)
#===============================================================================
echo ""
echo "Step 1: Converting BAM files to bedgraph format..."

GENOME_SIZE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes"

convert_bam_to_bedgraph() {
    local bam=$1
    local name=$2
    local bdg="${OUTPUT_DIR}/${name}.bedgraph"

    if [[ -f "${bdg}" ]]; then
        echo "  ${name}: bedgraph already exists, skipping"
        return
    fi

    echo "  Converting ${name}..."
    # Use bedtools genomecov to create bedgraph
    bedtools genomecov -ibam ${bam} -bg -pc | \
        sort -k1,1 -k2,2n > ${bdg}

    echo "  ${name}: $(wc -l < ${bdg}) regions"
}

# Convert treatment samples
for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    convert_bam_to_bedgraph "${BAM_DIR}/${SAMPLE}_filtered.bam" "${SAMPLE}"
done

# Convert controls
convert_bam_to_bedgraph "${BAM_DIR}/IggMs_filtered.bam" "IggMs"
convert_bam_to_bedgraph "${BAM_DIR}/IggRb_filtered.bam" "IggRb"

#===============================================================================
# Run SEACR peak calling
#===============================================================================
echo ""
echo "Step 2: Running SEACR peak calling..."

run_seacr() {
    local treatment=$1
    local control=$2
    local output_name=$3
    local threshold=$4  # "stringent" or "relaxed"

    echo "  ${output_name} (${threshold})..."

    ${SEACR_PATH} \
        ${OUTPUT_DIR}/${treatment}.bedgraph \
        ${OUTPUT_DIR}/${control}.bedgraph \
        non \
        ${threshold} \
        ${OUTPUT_DIR}/${output_name}_${threshold}

    if [[ -f "${OUTPUT_DIR}/${output_name}_${threshold}.stringent.bed" ]] || \
       [[ -f "${OUTPUT_DIR}/${output_name}_${threshold}.bed" ]]; then
        local peak_file=$(ls ${OUTPUT_DIR}/${output_name}_${threshold}*.bed 2>/dev/null | head -1)
        local count=$(wc -l < "${peak_file}" 2>/dev/null || echo 0)
        echo "    Found ${count} peaks"
    else
        echo "    WARNING: No peak file generated"
    fi
}

# Individual replicates
for SAMPLE in TES-1 TES-2 TES-3; do
    run_seacr "${SAMPLE}" "IggMs" "${SAMPLE}" "stringent"
    run_seacr "${SAMPLE}" "IggMs" "${SAMPLE}" "relaxed"
done

for SAMPLE in TEAD1-1 TEAD1-2 TEAD1-3; do
    run_seacr "${SAMPLE}" "IggRb" "${SAMPLE}" "stringent"
    run_seacr "${SAMPLE}" "IggRb" "${SAMPLE}" "relaxed"
done

#===============================================================================
# Step 3: Compare SEACR vs MACS2 peak calls
#===============================================================================
echo ""
echo "Step 3: Comparing SEACR vs MACS2 peak calls..."

COMPARISON_FILE="${OUTPUT_DIR}/seacr_vs_macs2_comparison.txt"

echo "=== SEACR vs MACS2 Peak Call Comparison ===" > ${COMPARISON_FILE}
echo "Date: $(date)" >> ${COMPARISON_FILE}
echo "" >> ${COMPARISON_FILE}

printf "%-15s %10s %10s %10s %10s\n" "Sample" "MACS2" "SEACR_str" "Overlap" "Pct_overlap" >> ${COMPARISON_FILE}
printf "%-15s %10s %10s %10s %10s\n" "------" "-----" "---------" "-------" "-----------" >> ${COMPARISON_FILE}

for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    MACS2_FILE="${MACS2_DIR}/${SAMPLE}_peaks.narrowPeak"
    SEACR_FILE=$(ls ${OUTPUT_DIR}/${SAMPLE}_stringent*.bed 2>/dev/null | head -1)

    if [[ -f "${MACS2_FILE}" && -f "${SEACR_FILE}" ]]; then
        MACS2_COUNT=$(wc -l < ${MACS2_FILE})
        SEACR_COUNT=$(wc -l < ${SEACR_FILE})

        # Count MACS2 peaks overlapping SEACR peaks
        OVERLAP_COUNT=$(bedtools intersect -u -a ${MACS2_FILE} -b ${SEACR_FILE} | wc -l)

        if [[ ${MACS2_COUNT} -gt 0 ]]; then
            OVERLAP_PCT=$((OVERLAP_COUNT * 100 / MACS2_COUNT))
        else
            OVERLAP_PCT=0
        fi

        printf "%-15s %10d %10d %10d %9d%%\n" "${SAMPLE}" "${MACS2_COUNT}" "${SEACR_COUNT}" "${OVERLAP_COUNT}" "${OVERLAP_PCT}" >> ${COMPARISON_FILE}
    else
        printf "%-15s %10s %10s %10s %10s\n" "${SAMPLE}" "N/A" "N/A" "N/A" "N/A" >> ${COMPARISON_FILE}
    fi
done

echo ""
cat ${COMPARISON_FILE}

echo ""
echo "=========================================="
echo "SEACR peak calling complete"
echo "Results: ${OUTPUT_DIR}/"
echo "Comparison: ${COMPARISON_FILE}"
echo "Date: $(date)"
echo "=========================================="
