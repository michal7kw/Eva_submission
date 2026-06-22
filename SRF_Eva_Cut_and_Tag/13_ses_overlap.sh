#!/bin/bash

#===============================================================================
# SCRIPT: 13_ses_overlap.sh
# PURPOSE: To quantify the overlap between experimental ChIP-seq peaks and a
#          pre-defined set of Super-Enhancer Signature (SES) regions.
#
# DESCRIPTION:
# This script performs an intersection analysis to determine the extent to which
# ChIP-seq peaks (both narrow and broad) from different experimental conditions
# co-localize with known SES consensus sites. It handles potential mismatches in
# chromosome naming conventions (e.g., '1' vs 'chr1') and produces detailed
# statistical reports and BED files of the overlapping regions.
#
# KEY OPERATIONS:
# 1. Iterates through all narrow and broad peak files for each sample.
# 2. Standardizes chromosome names in peak files to the 'chr' format (e.g., 'chr1')
#    to ensure compatibility with the SES consensus peak file.
# 3. Uses 'bedtools intersect' to find genomic regions common to both the
#    sample's peaks and the SES peaks.
# 4. Calculates key statistics for each sample:
#    - The number and percentage of the sample's peaks that overlap with SES regions.
#    - The number and percentage of SES regions that are covered by the sample's peaks.
# 5. Generates individual statistics files and a comprehensive summary report.
# 6. Creates combined BED files containing all overlapping regions for downstream
#    analysis or visualization.
#
# METHODOLOGY:
# The core of the analysis relies on 'bedtools intersect' for genomic coordinate
# comparisons. Chromosome name standardization is performed on-the-fly using 'awk'.
# Statistical calculations are done using shell commands like 'wc', 'cut', 'sort',
# and 'bc' for floating-point arithmetic.
#
# IMPORTANT PARAMETERS:
# - SES_CONSENSUS: Path to the BED file containing SES consensus peaks.
# - NARROW_PEAKS_DIR: Directory containing narrow peak files (*.narrowPeak).
# - BROAD_PEAKS_DIR: Directory containing broad peak files (*.broadPeak).
# - OUTPUT_DIR: Directory where all results will be saved.
# - SLURM parameters (--mem, --time, etc.) define computational resources.
#
# INPUTS:
# - SES consensus peaks: ${BASE_DIR}/DATA/SES_consensus.broadPeak
# - Individual sample narrow peaks: ${NARROW_PEAKS_DIR}/*_peaks.narrowPeak
# - Combined sample broad peaks: ${BROAD_PEAKS_DIR}/*_peaks.broadPeak
#
# OUTPUTS:
# - Overlap BED files: ${OUTPUT_DIR}/*_ses_overlap.bed
# - Statistics summary files: ${OUTPUT_DIR}/*_ses_stats.txt
# - A global summary report: ${OUTPUT_DIR}/ses_overlap_summary.txt
# - Merged BED files of unique overlapping regions:
#   - ${OUTPUT_DIR}/unique_narrow_overlap_regions.bed
#   - ${OUTPUT_DIR}/unique_broad_overlap_regions.bed
#
# DEPENDENCIES:
# - bedtools
# - bc (for calculations)
# - A conda environment named 'overlap' where dependencies are installed.
#===============================================================================

#SBATCH --job-name=13_ses_overlap
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/13_ses_overlap.err"
#SBATCH --output="./logs/13_ses_overlap.out"

# Initialize conda for the current shell
echo "Initializing conda..."
eval "$(conda shell.bash hook)"

# Activate conda environment
echo "Activating conda environment..."
conda activate overlap

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
cd $BASE_DIR

# Define paths
SES_CONSENSUS="${BASE_DIR}/DATA/SES_consensus.broadPeak"
NARROW_PEAKS_DIR="${BASE_DIR}/results/05_peaks_narrow"
BROAD_PEAKS_DIR="${BASE_DIR}/results/05_peaks_broad"
CONSENSUS_NARROW_DIR="${BASE_DIR}/results/11_combined_replicates_narrow"
OUTPUT_DIR="${BASE_DIR}/results/13_ses_overlap"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "SES Consensus Overlap Analysis (Improved)"
echo "=========================================="
echo "This script will use high-quality consensus peaks when available"
echo "Priority: Consensus > IDR > MACS2 combined peaks"
echo ""
echo "SES consensus peaks: $(wc -l < ${SES_CONSENSUS})"
echo ""

# Function to convert chromosome names and perform overlap analysis
perform_overlap() {
    local peak_file=$1
    local peak_type=$2
    local sample_name=$(basename ${peak_file} | sed 's/_peaks\.(narrowPeak|broadPeak)//')

    echo "Processing ${sample_name} (${peak_type})..."

    # Create temporary file with fixed chromosome names
    local temp_peak_file="${OUTPUT_DIR}/${sample_name}_${peak_type}_chr.bed"

    # Convert chromosome names: 1 -> chr1, 2 -> chr2, etc.
    # Skip lines that already have 'chr' prefix or are not standard chromosomes
    awk 'BEGIN{OFS="\t"} {
        if ($1 ~ /^[0-9]+$/ || $1 == "X" || $1 == "Y" || $1 == "MT" || $1 == "M") {
            $1 = "chr" $1
        }
        print
    }' ${peak_file} > ${temp_peak_file}

    # Create output files
    local overlap_file="${OUTPUT_DIR}/${sample_name}_${peak_type}_ses_overlap.bed"
    local stats_file="${OUTPUT_DIR}/${sample_name}_${peak_type}_ses_stats.txt"

    # Find overlaps using bedtools intersect
    bedtools intersect -a ${temp_peak_file} -b ${SES_CONSENSUS} -wa -wb > ${overlap_file}

    # Calculate statistics
    local total_peaks=$(wc -l < ${peak_file})
    local overlapping_peaks=$(cut -f1-3 ${overlap_file} | sort -u | wc -l)

    # Dynamically determine SES peak columns based on input file column count
    # SES columns start after the peak file columns
    local peak_cols=$(head -1 ${temp_peak_file} | awk '{print NF}')
    local ses_start_col=$((peak_cols + 1))
    local ses_end_col=$((peak_cols + 3))
    local ses_peaks_overlapped=$(cut -f${ses_start_col}-${ses_end_col} ${overlap_file} | sort -u | wc -l)
    local total_ses_peaks=$(wc -l < ${SES_CONSENSUS})

    # Calculate percentages
    local overlap_percentage=0
    local ses_coverage_percentage=0

    if [ ${total_peaks} -gt 0 ]; then
        overlap_percentage=$(echo "scale=2; ${overlapping_peaks} * 100 / ${total_peaks}" | bc -l)
    fi

    if [ ${total_ses_peaks} -gt 0 ]; then
        ses_coverage_percentage=$(echo "scale=2; ${ses_peaks_overlapped} * 100 / ${total_ses_peaks}" | bc -l)
    fi

    # Write statistics
    cat > ${stats_file} << EOF
Sample: ${sample_name}
Peak Type: ${peak_type}
Total ${sample_name} peaks: ${total_peaks}
${sample_name} peaks overlapping SES: ${overlapping_peaks}
Percentage of ${sample_name} peaks overlapping SES: ${overlap_percentage}%
SES peaks overlapped by ${sample_name}: ${ses_peaks_overlapped}
Total SES consensus peaks: ${total_ses_peaks}
Percentage of SES peaks covered by ${sample_name}: ${ses_coverage_percentage}%
EOF

    echo "  ${sample_name} results:"
    echo "    ${overlapping_peaks}/${total_peaks} (${overlap_percentage}%) peaks overlap with SES"
    echo "    Covers ${ses_peaks_overlapped}/${total_ses_peaks} (${ses_coverage_percentage}%) SES peaks"

    # Clean up temporary file
    rm ${temp_peak_file}
}

# Function to select best narrow peak file
select_narrow_peak() {
    local condition=$1
    local consensus_scored="${CONSENSUS_NARROW_DIR}/peaks/${condition}_consensus_peaks_scored.bed"
    local consensus_full="${CONSENSUS_NARROW_DIR}/peaks/${condition}_consensus_peaks.bed"
    local idr_union="${CONSENSUS_NARROW_DIR}/idr/${condition}_idr_union.bed"
    local macs2_file="${NARROW_PEAKS_DIR}/${condition}_peaks.narrowPeak"

    if [ -f "$consensus_scored" ]; then
        echo "$consensus_scored"
    elif [ -f "$consensus_full" ]; then
        echo "$consensus_full"
    elif [ -f "$idr_union" ]; then
        echo "$idr_union"
    elif [ -f "$macs2_file" ]; then
        echo "$macs2_file"
    else
        echo ""
    fi
}

# Process narrow peaks (use consensus when available)
echo ""
echo "=== Analyzing Narrow Peaks (Consensus/IDR Preferred) ==="

# Check for consensus peaks for combined conditions
for condition in TES TEAD1; do  # TESmut removed - failed sample
    selected_file=$(select_narrow_peak "$condition")
    if [ -n "$selected_file" ]; then
        echo "Selected for $condition: $(basename $selected_file)"
        perform_overlap "$selected_file" "narrow"
    fi
done

# Also process individual replicate peaks if they exist (for comparison)
echo ""
echo "=== Analyzing Individual Replicate Narrow Peaks (if available) ==="
for peak_file in ${NARROW_PEAKS_DIR}/*-[123]_peaks.narrowPeak; do
    if [ -f "$peak_file" ]; then
        echo "Processing individual replicate: $(basename $peak_file)"
        perform_overlap "$peak_file" "narrow_replicate"
    fi
done

# Process broad peaks (combined samples)
echo ""
echo "=== Analyzing Broad Peaks ==="
for peak_file in ${BROAD_PEAKS_DIR}/*_peaks.broadPeak; do
    if [ -f "$peak_file" ]; then
        perform_overlap "$peak_file" "broad"
    fi
done

# Create summary report
echo ""
echo "=== Creating Summary Report ==="
summary_file="${OUTPUT_DIR}/ses_overlap_summary.txt"

cat > ${summary_file} << EOF
SES Consensus Overlap Analysis Summary (Chromosome Names Fixed)
===============================================================
Date: $(date)
SES Consensus File: ${SES_CONSENSUS}
Total SES Consensus Peaks: $(wc -l < ${SES_CONSENSUS})

Note: Chromosome names were converted from numeric (1,2,3...) to chr format (chr1,chr2,chr3...)
to match SES consensus peak naming convention.

EOF

echo "NARROW PEAKS ANALYSIS:" >> ${summary_file}
echo "======================" >> ${summary_file}
for stats_file in ${OUTPUT_DIR}/*_narrow_ses_stats.txt; do
    if [ -f "$stats_file" ]; then
        echo "" >> ${summary_file}
        cat "$stats_file" >> ${summary_file}
    fi
done

echo "" >> ${summary_file}
echo "BROAD PEAKS ANALYSIS:" >> ${summary_file}
echo "=====================" >> ${summary_file}
for stats_file in ${OUTPUT_DIR}/*_broad_ses_stats.txt; do
    if [ -f "$stats_file" ]; then
        echo "" >> ${summary_file}
        cat "$stats_file" >> ${summary_file}
    fi
done

# Create combined overlap BED files for visualization
echo ""
echo "=== Creating Combined Overlap Files ==="

# Combine all narrow peak overlaps
if ls ${OUTPUT_DIR}/*_narrow_ses_overlap.bed 1> /dev/null 2>&1; then
    cat ${OUTPUT_DIR}/*_narrow_ses_overlap.bed > ${OUTPUT_DIR}/all_narrow_ses_overlaps.bed
    # Create unique overlapping regions
    cut -f1-3 ${OUTPUT_DIR}/all_narrow_ses_overlaps.bed | sort -k1,1 -k2,2n | bedtools merge > ${OUTPUT_DIR}/unique_narrow_overlap_regions.bed
fi

# Combine all broad peak overlaps
if ls ${OUTPUT_DIR}/*_broad_ses_overlap.bed 1> /dev/null 2>&1; then
    cat ${OUTPUT_DIR}/*_broad_ses_overlap.bed > ${OUTPUT_DIR}/all_broad_ses_overlaps.bed
    # Create unique overlapping regions
    cut -f1-3 ${OUTPUT_DIR}/all_broad_ses_overlaps.bed | sort -k1,1 -k2,2n | bedtools merge > ${OUTPUT_DIR}/unique_broad_overlap_regions.bed
fi

echo ""
echo "SES overlap analysis complete (with chromosome name correction)!"
echo "Results saved in: ${OUTPUT_DIR}"
echo "Summary report: ${summary_file}"
echo ""
echo "Key output files:"
echo "  - Individual overlap files: *_ses_overlap.bed"
echo "  - Statistics files: *_ses_stats.txt"
echo "  - Combined overlaps: all_*_ses_overlaps.bed"
echo "  - Unique overlap regions: unique_*_overlap_regions.bed"
echo "  - Summary report: ses_overlap_summary.txt"
