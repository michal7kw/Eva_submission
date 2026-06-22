#!/bin/bash

#===============================================================================
# SCRIPT: 4b_fragment_stratification.sh
# PURPOSE: Stratify CUT&Tag reads by fragment size for TF-specific analysis
#
# DESCRIPTION:
# CUT&Tag produces fragments of distinct sizes reflecting different chromatin
# states. Sub-nucleosomal fragments (<120bp) represent direct TF binding sites
# and are most informative for transcription factor peak calling. This script
# separates fragments by size class and generates fragment size distribution
# plots, enabling TF-enriched peak calling as a complementary analysis.
#
# Fragment size classes (Henikoff lab CUT&Tag protocol):
#   <120bp:    Sub-nucleosomal (TF footprints) — most informative for TF binding
#   150-300bp: Mono-nucleosomal (flanking nucleosomes)
#   >300bp:    Multi-nucleosomal (broader chromatin)
#
# USAGE:
# sbatch 4b_fragment_stratification.sh
#
# INPUTS:
# - Filtered BAM files from step 4
#
# OUTPUTS:
# - results/04b_fragment_stratification/<sample>_subnucleosomal.bam (<120bp)
# - results/04b_fragment_stratification/<sample>_mononucleosomal.bam (150-300bp)
# - results/04b_fragment_stratification/fragment_size_summary.txt
#===============================================================================

#SBATCH --job-name=4b_frag_strat
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4b_frag_strat.err"
#SBATCH --output="./logs/4b_frag_strat.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate cutntag

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/04b_fragment_stratification"

mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "Fragment Size Stratification for CUT&Tag"
echo "Date: $(date)"
echo "=========================================="

# Size thresholds
SUBNUC_MAX=120
MONONUC_MIN=150
MONONUC_MAX=300

SUMMARY_FILE="${OUTPUT_DIR}/fragment_size_summary.txt"
printf "%-15s %12s %12s %12s %12s %10s\n" \
    "Sample" "Total" "SubNuc(<120)" "MonoNuc" "MultiNuc" "SubNuc_Pct" > ${SUMMARY_FILE}
printf "%-15s %12s %12s %12s %12s %10s\n" \
    "------" "-----" "------------" "--------" "--------" "----------" >> ${SUMMARY_FILE}

for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_filtered.bam"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "WARNING: BAM not found for ${SAMPLE}, skipping"
        continue
    fi

    echo "Processing ${SAMPLE}..."

    # Index if needed
    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        samtools index -@ 8 ${BAM_FILE}
    fi

    # Extract sub-nucleosomal fragments (<120bp)
    # -f 2: properly paired; fragment size filtered via awk on TLEN field
    echo "  Extracting sub-nucleosomal fragments (<${SUBNUC_MAX}bp)..."
    samtools view -h -f 2 ${BAM_FILE} | \
        awk -v max=${SUBNUC_MAX} 'BEGIN {OFS="\t"} /^@/ {print; next} {if ($9 > 0 && $9 < max) print; else if ($9 < 0 && $9 > -max) print}' | \
        samtools view -bS -@ 8 - | \
        samtools sort -@ 8 -o ${OUTPUT_DIR}/${SAMPLE}_subnucleosomal.bam -
    samtools index -@ 8 ${OUTPUT_DIR}/${SAMPLE}_subnucleosomal.bam

    # Extract mono-nucleosomal fragments (150-300bp)
    echo "  Extracting mono-nucleosomal fragments (${MONONUC_MIN}-${MONONUC_MAX}bp)..."
    samtools view -h -f 2 ${BAM_FILE} | \
        awk -v min=${MONONUC_MIN} -v max=${MONONUC_MAX} 'BEGIN {OFS="\t"} /^@/ {print; next} {tlen=($9>0?$9:-$9); if (tlen >= min && tlen <= max) print}' | \
        samtools view -bS -@ 8 - | \
        samtools sort -@ 8 -o ${OUTPUT_DIR}/${SAMPLE}_mononucleosomal.bam -
    samtools index -@ 8 ${OUTPUT_DIR}/${SAMPLE}_mononucleosomal.bam

    # Count reads in each fraction
    TOTAL=$(samtools view -c -f 2 ${BAM_FILE})
    SUBNUC=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE}_subnucleosomal.bam)
    MONONUC=$(samtools view -c ${OUTPUT_DIR}/${SAMPLE}_mononucleosomal.bam)
    MULTINUC=$((TOTAL - SUBNUC - MONONUC))

    if [[ ${TOTAL} -gt 0 ]]; then
        SUBNUC_PCT=$(awk "BEGIN {printf \"%.1f%%\", ${SUBNUC}/${TOTAL}*100}")
    else
        SUBNUC_PCT="0.0%"
    fi

    printf "%-15s %12d %12d %12d %12d %10s\n" \
        "${SAMPLE}" "${TOTAL}" "${SUBNUC}" "${MONONUC}" "${MULTINUC}" "${SUBNUC_PCT}" >> ${SUMMARY_FILE}

    echo "  Total: ${TOTAL}, SubNuc: ${SUBNUC} (${SUBNUC_PCT}), MonoNuc: ${MONONUC}, MultiNuc: ${MULTINUC}"
done

echo ""
echo "=== Fragment Size Summary ==="
cat ${SUMMARY_FILE}

# Generate fragment size distribution using samtools
echo ""
echo "Generating fragment size distributions..."
for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    BAM_FILE="${BAM_DIR}/${SAMPLE}_filtered.bam"
    if [[ -f "${BAM_FILE}" ]]; then
        echo "  ${SAMPLE}: extracting insert sizes..."
        samtools view -f 2 ${BAM_FILE} | \
            awk '$9 > 0 {print $9}' | \
            sort -n | uniq -c | \
            awk '{print $2"\t"$1}' > ${OUTPUT_DIR}/${SAMPLE}_fragment_sizes.txt
    fi
done

echo ""
echo "=========================================="
echo "Fragment stratification complete"
echo "Sub-nucleosomal BAMs: ${OUTPUT_DIR}/*_subnucleosomal.bam"
echo "  (Use these for TF-enriched peak calling with MACS2)"
echo "Summary: ${SUMMARY_FILE}"
echo ""
echo "To call peaks on sub-nucleosomal fragments:"
echo "  macs2 callpeak -t ${OUTPUT_DIR}/<sample>_subnucleosomal.bam \\"
echo "    -c <control>.bam -f BAMPE -g hs -q 0.05 --keep-dup all"
echo "Date: $(date)"
echo "=========================================="
