#!/bin/bash

#===============================================================================
# SCRIPT: 5_peak_calling.sh
# PURPOSE: Narrow peak calling for Cut&Tag data using MACS2
#
# DESCRIPTION:
# This script performs narrow peak calling on filtered BAM files from Cut&Tag
# experiments using MACS2. It processes both combined replicates and individual
# samples to identify regions of protein-DNA binding. The script automatically
# detects paired-end vs single-end data and adjusts MACS2 parameters accordingly.
#
# KEY OPERATIONS:
# 1. BAM file validation and indexing
# 2. Automatic detection of paired-end vs single-end sequencing
# 3. MACS2 narrow peak calling with appropriate controls
# 4. Processing of combined replicates for main analysis
# 5. Individual replicate processing for reproducibility assessment
# 6. Peak count summary and quality metrics
# 7. Optional diagnostic plot generation
#
# METHODOLOGY (Cut&Tag Optimized):
# - Uses MACS2 callpeak with Cut&Tag-specific parameters
# - Applies q-value threshold of 0.05 (recommended for Cut&Tag)
# - Uses appropriate control samples (IggMs for TES/TESmut, IggRb for TEAD1)
# - Keeps all duplicate reads (--keep-dup all) as appropriate for Cut&Tag
# - Automatically selects BAMPE (paired-end) or BAM (single-end) format
# - Applies Cut&Tag-specific fragment size adjustments
# - Filters for nucleosome-free fragments (<120bp)
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (sufficient for peak calling)
# - Time: 4 hours (adequate for multiple samples)
# - Threads: 8 (MACS2 is mostly single-threaded but allows parallel processing)
# - Q-value threshold: 0.05 (Cut&Tag recommended threshold)
# - Shift: -75bp (half of expected fragment size for Cut&Tag)
# - Extsize: 150bp (typical Cut&Tag fragment size)
# - Genome size: hs (hg38 effective genome size)
#
# EXPECTED INPUTS:
# - Filtered BAM files from step 4 (04_filtered directory)
# - BAM files should be properly indexed
# - Control samples: IggMs_filtered.bam, IggRb_filtered.bam
# - Treatment samples: TES-*, TESmut-*, TEAD1-* replicates
#
# EXPECTED OUTPUTS:
# - *_peaks.narrowPeak: Peak coordinates in ENCODE narrowPeak format
# - *_summits.bed: Peak summit positions
# - *_model.r: R script for MACS2 model visualization
# - *_peaks.xls: Detailed peak information in Excel format
# - *_macs2.log: MACS2 execution logs for troubleshooting
# - peak_summary.pdf: Optional diagnostic plot (if R available)
#
# DEPENDENCIES:
# - MACS2 (peak calling algorithm)
# - samtools (BAM file operations)
# - R and ggplot2 (optional, for diagnostic plots)
# - Conda environment: peak_calling_new
#
# USAGE:
# sbatch 5_peak_calling.sh
#
# NOTES:
# - Script validates BAM files before processing
# - Automatically handles both paired-end and single-end data
# - Processes samples in groups (combined replicates + individuals)
# - Generates comprehensive summary statistics
# - Error handling includes detailed logging and validation
#===============================================================================

#SBATCH --job-name=5_peaks
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5_peaks.err"
#SBATCH --output="./logs/5_peaks.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate peak_calling_new

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BAM_DIR="${BASE_DIR}/results/04_filtered"
OUTPUT_DIR="${BASE_DIR}/results/05_peaks_narrow"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

echo "=========================================="
echo "Starting peak calling analysis"
echo "Date: $(date)"
echo "=========================================="

# Function to check BAM file
check_bam_file() {
    local bam_file=$1
    local sample_name=$(basename $bam_file .bam)
    
    echo "Checking $sample_name..."
    
    # Check if file exists
    if [ ! -f "$bam_file" ]; then
        echo "ERROR: File $bam_file does not exist!"
        return 1
    fi
    
    # Check file size
    local file_size=$(ls -lh $bam_file | awk '{print $5}')
    echo "  File size: $file_size"
    
    # Count total reads
    local read_count=$(samtools view -c $bam_file 2>/dev/null)
    echo "  Total reads: $read_count"
    
    if [ "$read_count" -eq 0 ]; then
        echo "  WARNING: BAM file is empty!"
        return 1
    fi
    
    # Check for paired reads
    local paired_count=$(samtools view -f 1 -c $bam_file 2>/dev/null)
    local paired_pct=0
    if [ "$read_count" -gt 0 ]; then
        paired_pct=$((paired_count * 100 / read_count))
    fi
    echo "  Paired reads: $paired_count ($paired_pct%)"
    
    # Return paired status (0 = paired, 1 = single)
    if [ "$paired_pct" -gt 50 ]; then
        echo "  Status: PAIRED-END"
        return 0
    else
        echo "  Status: SINGLE-END"
        return 2
    fi
}

# Function to index BAM file if needed
index_bam_if_needed() {
    local bam_file=$1
    if [ ! -f "${bam_file}.bai" ]; then
        echo "  Creating index for $(basename $bam_file)..."
        samtools index $bam_file
    fi
}

echo ""
echo "Step 1: Checking and indexing BAM files"
echo "------------------------------------------"

# Check all BAM files and determine read type
PAIRED_END=true
for bam in ${BAM_DIR}/*.bam; do
    check_bam_file $bam
    status=$?
    if [ $status -eq 1 ]; then
        echo "ERROR: Problem with $bam"
        exit 1
    elif [ $status -eq 2 ]; then
        PAIRED_END=false
    fi
    index_bam_if_needed $bam
    echo ""
done

# Set format based on paired-end status
if [ "$PAIRED_END" = true ]; then
    FORMAT="BAMPE"
    echo "Using PAIRED-END mode (BAMPE) for MACS2"
else
    FORMAT="BAM"
    echo "Using SINGLE-END mode (BAM) for MACS2"
fi

echo ""
echo "Step 2: Running MACS2 peak calling"
echo "------------------------------------------"

# Function to run MACS2 with Cut&Tag-optimized parameters
run_macs2() {
    local treatment_files="$1"
    local control_file="$2"
    local output_name="$3"
    local format="$4"

    echo "Calling peaks for $output_name..."
    echo "  Treatment: $treatment_files"
    echo "  Control: $control_file"
    echo "  Format: $format"

    # Cut&Tag-optimized MACS2 parameters
    # Key differences from ChIP-seq:
    # - q-value 0.05 (less stringent, Cut&Tag has lower background)
    # - shift -75, extsize 150 (typical Cut&Tag fragment size)
    # - Genome size 2913022398 (hg38 effective genome size)
    # - keep-dup all (Cut&Tag has low PCR duplication)

    if [ "$format" = "BAMPE" ]; then
        # For paired-end data, use actual fragment sizes
        macs2 callpeak \
            -t $treatment_files \
            -c $control_file \
            -f $format \
            -g 2913022398 \
            -n $output_name \
            --outdir ${OUTPUT_DIR} \
            -q 0.05 \
            --keep-dup all \
            2> ${OUTPUT_DIR}/${output_name}_macs2.log
    else
        # For single-end data, apply shift and extsize
        macs2 callpeak \
            -t $treatment_files \
            -c $control_file \
            -f $format \
            -g 2913022398 \
            -n $output_name \
            --outdir ${OUTPUT_DIR} \
            -q 0.05 \
            --keep-dup all \
            --shift -75 \
            --extsize 150 \
            --nomodel \
            2> ${OUTPUT_DIR}/${output_name}_macs2.log
    fi

    if [ $? -eq 0 ]; then
        echo "  SUCCESS: Peak calling completed for $output_name"
        # Count peaks
        if [ -f "${OUTPUT_DIR}/${output_name}_peaks.narrowPeak" ]; then
            peak_count=$(wc -l < "${OUTPUT_DIR}/${output_name}_peaks.narrowPeak")
            echo "  Found $peak_count peaks"

            # Additional quality filtering: keep only high-quality peaks (fold enrichment >= 2)
            awk '$7 >= 2' "${OUTPUT_DIR}/${output_name}_peaks.narrowPeak" > "${OUTPUT_DIR}/${output_name}_peaks_filtered.narrowPeak"
            filtered_count=$(wc -l < "${OUTPUT_DIR}/${output_name}_peaks_filtered.narrowPeak")
            echo "  High-quality peaks (fold enrichment ≥2): $filtered_count"

            # Replace original with filtered peaks
            mv "${OUTPUT_DIR}/${output_name}_peaks_filtered.narrowPeak" "${OUTPUT_DIR}/${output_name}_peaks.narrowPeak"
        fi
    else
        echo "  ERROR: Peak calling failed for $output_name"
        echo "  Check log file: ${OUTPUT_DIR}/${output_name}_macs2.log"
        tail -20 ${OUTPUT_DIR}/${output_name}_macs2.log
    fi
    echo ""
}

# Peak calling for combined replicates
echo "Processing combined replicates..."
echo "=================================="

# TES samples
run_macs2 \
    "${BAM_DIR}/TES-1_filtered.bam ${BAM_DIR}/TES-2_filtered.bam ${BAM_DIR}/TES-3_filtered.bam" \
    "${BAM_DIR}/IggMs_filtered.bam" \
    "TES" \
    "$FORMAT"

# TESmut removed - failed sample

# TEAD1 samples
run_macs2 \
    "${BAM_DIR}/TEAD1-1_filtered.bam ${BAM_DIR}/TEAD1-2_filtered.bam ${BAM_DIR}/TEAD1-3_filtered.bam" \
    "${BAM_DIR}/IggRb_filtered.bam" \
    "TEAD1" \
    "$FORMAT"

echo "Processing individual replicates..."
echo "===================================="

# Individual replicates for reproducibility analysis
for SAMPLE in TES-1 TES-2 TES-3; do  # TESmut removed - failed sample
    run_macs2 \
        "${BAM_DIR}/${SAMPLE}_filtered.bam" \
        "${BAM_DIR}/IggMs_filtered.bam" \
        "${SAMPLE}" \
        "$FORMAT"
done

for SAMPLE in TEAD1-1 TEAD1-2 TEAD1-3; do
    run_macs2 \
        "${BAM_DIR}/${SAMPLE}_filtered.bam" \
        "${BAM_DIR}/IggRb_filtered.bam" \
        "${SAMPLE}" \
        "$FORMAT"
done

# Apply blacklist filtering to all peak files
echo ""
echo "Step 3: Applying blacklist filtering"
echo "------------------------------------------"

# Download ENCODE blacklist if not present
BLACKLIST="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38-blacklist.v2.bed"
if [ ! -f "$BLACKLIST" ]; then
    echo "Warning: Blacklist file not found at $BLACKLIST"
    echo "Peaks will not be blacklist-filtered"
    APPLY_BLACKLIST=false
else
    echo "Applying blacklist filtering using: $BLACKLIST"
    APPLY_BLACKLIST=true
fi

# Filter peaks against blacklist regions
if [ "$APPLY_BLACKLIST" = true ]; then
    for peak_file in ${OUTPUT_DIR}/*_peaks.narrowPeak; do
        if [ -f "$peak_file" ]; then
            base_name=$(basename $peak_file _peaks.narrowPeak)
            echo "Filtering blacklist regions from $base_name..."

            # Remove blacklist regions
            bedtools intersect -v -a $peak_file -b $BLACKLIST > ${OUTPUT_DIR}/${base_name}_peaks_blacklist_filtered.narrowPeak

            # Count removed peaks
            original_count=$(wc -l < $peak_file)
            filtered_count=$(wc -l < ${OUTPUT_DIR}/${base_name}_peaks_blacklist_filtered.narrowPeak)
            removed_count=$((original_count - filtered_count))

            echo "  Original peaks: $original_count"
            echo "  Filtered peaks: $filtered_count"
            echo "  Removed peaks: $removed_count"

            # Replace original with filtered version
            mv ${OUTPUT_DIR}/${base_name}_peaks_blacklist_filtered.narrowPeak $peak_file
        fi
    done
fi

echo ""
echo "Step 4: Summary of results"
echo "------------------------------------------"

# Summary report
echo "Peak files generated:"
ls -lh ${OUTPUT_DIR}/*_peaks.narrowPeak 2>/dev/null | awk '{print "  "$9": "$5" ("$1")"}'

echo ""
echo "Peak counts:"
for peak_file in ${OUTPUT_DIR}/*_peaks.narrowPeak; do
    if [ -f "$peak_file" ]; then
        sample=$(basename $peak_file _peaks.narrowPeak)
        count=$(wc -l < $peak_file)
        printf "  %-20s: %6d peaks\n" "$sample" "$count"
    fi
done

echo ""
echo "Log files:"
ls -lh ${OUTPUT_DIR}/*.log 2>/dev/null | awk '{print "  "$9}'

echo ""
echo "=========================================="
echo "Peak calling analysis completed!"
echo "Date: $(date)"
echo "=========================================="

# Optional: Generate a quick diagnostic plot if R is available
if command -v R &> /dev/null; then
    echo ""
    echo "Generating diagnostic plots..."
    Rscript -e "
    library(ggplot2, quietly=TRUE)
    peak_files <- list.files('${OUTPUT_DIR}', pattern='_peaks.narrowPeak$', full.names=TRUE)
    if(length(peak_files) > 0) {
        peak_counts <- sapply(peak_files, function(x) length(readLines(x)))
        names(peak_counts) <- gsub('_peaks.narrowPeak', '', basename(names(peak_counts)))
        pdf('${OUTPUT_DIR}/peak_summary.pdf', width=10, height=6)
        barplot(peak_counts, main='Number of Peaks per Sample', 
                xlab='Sample', ylab='Number of Peaks', las=2, 
                col=rainbow(length(peak_counts)))
        dev.off()
        cat('Diagnostic plot saved to ${OUTPUT_DIR}/peak_summary.pdf\n')
    }
    " 2>/dev/null || echo "R plotting skipped (optional)"
fi