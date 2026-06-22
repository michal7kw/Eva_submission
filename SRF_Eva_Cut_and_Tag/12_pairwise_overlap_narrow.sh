#!/bin/bash

#===============================================================================
# SCRIPT: 12_pairwise_overlap_narrow.sh
# PURPOSE: Comprehensive pairwise binding site overlap analysis for NARROW PEAKS
#
# DESCRIPTION:
# This script performs detailed pairwise overlap analysis between different
# experimental conditions (TES, TEAD1) using narrow peaks to identify
# shared and condition-specific binding sites. It generates Venn diagrams,
# overlap statistics, and comprehensive visualizations.
# NOTE: TESmut samples removed - failed sample
#
# KEY OPERATIONS:
# 1. Pairwise overlap analysis using bedtools intersect on narrow peaks
# 2. Identification of unique and shared binding sites
# 3. Venn diagram generation for visual representation
# 4. Statistical analysis of overlap percentages
# 5. Comprehensive summary report generation
#
# METHODOLOGY:
# - Uses narrow peaks from MACS2 peak calling
# - Uses bedtools for genomic interval overlap calculations
# - Employs R/VennDiagram for publication-quality visualizations
# - Calculates multiple overlap metrics (intersection, unique regions)
# - Generates both numerical and visual summaries
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (for large peak file operations)
# - Time: 2 hours (adequate for comprehensive analysis)
# - Threads: 8 (parallel processing for efficiency)
# - Overlap threshold: Minimum 1bp overlap
# - Peak format: narrowPeak format
#
# INPUTS:
# - Narrow peak files from step 5 (05_peaks_narrow)
# - Three experimental conditions: TES, TESmut, TEAD1
#
# OUTPUTS:
# - Pairwise overlap BED files for each comparison
# - Condition-specific unique peak sets
# - Venn diagrams for visual overlap representation
# - Comprehensive overlap statistics
# - Detailed summary report with biological interpretation
#
# DEPENDENCIES:
# - bedtools for genomic interval operations
# - R with VennDiagram, ggplot2, dplyr packages
# - Conda environment: overlap
# - Narrow peak calling results from step 5
#
# USAGE:
# sbatch 12_pairwise_overlap_narrow.sh
#
# NOTES:
# - Uses narrow peaks for precise binding site overlap analysis
# - Generates publication-ready visualizations
# - Essential for understanding binding specificity
# - Includes biological interpretation and next steps
#===============================================================================

#SBATCH --job-name=pairwise_overlap_narrow
#SBATCH --account=kubacki.michal
#SBATCH --output=./logs/12_pairwise_overlap_narrow.out
#SBATCH --error=./logs/12_pairwise_overlap_narrow.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

echo "=== Pairwise Binding Site Overlap Analysis (Narrow Peaks) ==="
echo "Start time: $(date)"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate overlap

# Set up directories
OUTDIR="results/12_pairwise_overlap_narrow"
mkdir -p "$OUTDIR"/{overlaps,unique,venn,matrices,plots,stats}

# Define peak directories
PEAK_DIR="results/05_peaks_narrow"
CONSENSUS_DIR="results/11_combined_replicates_narrow"
ANALYSIS_TYPE="narrow"

# Function to select best available peak file
select_best_peak_file() {
    local condition=$1
    local consensus_scored="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks_scored.bed"
    local consensus_full="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks.bed"
    local idr_union="${CONSENSUS_DIR}/idr/${condition}_idr_union.bed"
    local macs2_combined="${PEAK_DIR}/${condition}_peaks.narrowPeak"

    if [ -f "$consensus_scored" ]; then
        echo "$consensus_scored"
        echo "  Using quality-scored consensus peaks for $condition" >&2
    elif [ -f "$consensus_full" ]; then
        echo "$consensus_full"
        echo "  Using full consensus peaks for $condition" >&2
    elif [ -f "$idr_union" ]; then
        echo "$idr_union"
        echo "  Using IDR union peaks for $condition" >&2
    elif [ -f "$macs2_combined" ]; then
        echo "$macs2_combined"
        echo "  Using MACS2 combined peaks for $condition (consensus not available)" >&2
    else
        echo ""
        echo "  ERROR: No peak file found for $condition" >&2
    fi
}

# Define conditions and select peak files (TES and TEAD1 only - TESmut removed)
echo "=== Selecting optimal peak sets for pairwise comparison ==="
CONDITIONS=("TES" "TEAD1")
declare -A PEAK_FILES

for condition in "${CONDITIONS[@]}"; do
    selected_file=$(select_best_peak_file "$condition")
    if [ -z "$selected_file" ]; then
        echo "Error: No peak file found for $condition"
        exit 1
    fi
    PEAK_FILES[$condition]="$selected_file"
done

# Check if peak files exist
echo "=== Checking narrow peak files ==="
for condition in "${CONDITIONS[@]}"; do
    if [ -f "${PEAK_FILES[$condition]}" ]; then
        peak_count=$(wc -l < "${PEAK_FILES[$condition]}")
        echo "$condition: $peak_count peaks (${PEAK_FILES[$condition]})"
    else
        echo "Error: Peak file not found for $condition: ${PEAK_FILES[$condition]}"
        exit 1
    fi
done

echo ""
echo "=== Performing pairwise overlap analysis (narrow peaks) ==="

# Function to perform pairwise comparison with multiple overlap methods
# For narrow peaks, we use several approaches:
# 1. Direct overlap (any overlap)
# 2. Summit-based proximity (peaks within distance threshold)
# 3. Reciprocal overlap (minimum fraction requirement)
perform_pairwise_comparison() {
    local cond1=$1
    local cond2=$2
    local file1="${PEAK_FILES[$cond1]}"
    local file2="${PEAK_FILES[$cond2]}"

    # Send informational output to STDERR so it displays but doesn't interfere with parsing
    echo "Analyzing $cond1 vs $cond2 using narrow peaks..." >&2

    # Create comparison directory
    mkdir -p "$OUTDIR/overlaps/${cond1}_vs_${cond2}"

    # Convert to BED3 format for consistent comparison, removing potential headers
    grep -vE '^#|^track' "$file1" | cut -f1-3 | sort -k1,1 -k2,2n > "$OUTDIR/${cond1}_temp.bed"
    grep -vE '^#|^track' "$file2" | cut -f1-3 | sort -k1,1 -k2,2n > "$OUTDIR/${cond2}_temp.bed"

    # METHOD 1: Direct overlap (minimum 1bp overlap)
    echo "  Finding overlapping regions (direct overlap)..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" -wa -u > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed"
    bedtools intersect -a "$OUTDIR/${cond2}_temp.bed" -b "$OUTDIR/${cond1}_temp.bed" -wa -u > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_overlapping_${cond1}.bed"

    # Find actual overlapping regions (intersection geometry)
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_intersection.bed"

    # METHOD 2: Summit-based proximity for narrow peaks
    # Create summit (center) positions and check if within 200bp of each other
    echo "  Finding summit-based overlaps (within 200bp)..." >&2
    awk -v OFS='\t' '{center=int(($2+$3)/2); print $1, center, center+1}' "$OUTDIR/${cond1}_temp.bed" > "$OUTDIR/${cond1}_summits.bed"
    awk -v OFS='\t' '{center=int(($2+$3)/2); print $1, center, center+1}' "$OUTDIR/${cond2}_temp.bed" > "$OUTDIR/${cond2}_summits.bed"

    # Expand summits by 200bp and find overlaps (peaks within 200bp of each other)
    bedtools slop -i "$OUTDIR/${cond1}_summits.bed" -g /beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes -b 200 2>/dev/null | \
    bedtools intersect -a - -b "$OUTDIR/${cond2}_summits.bed" -wa -u > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_summit_overlap_${cond2}.bed" 2>/dev/null || true

    bedtools slop -i "$OUTDIR/${cond2}_summits.bed" -g /beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/hg38.chrom.sizes -b 200 2>/dev/null | \
    bedtools intersect -a - -b "$OUTDIR/${cond1}_summits.bed" -wa -u > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_summit_overlap_${cond1}.bed" 2>/dev/null || true

    # METHOD 3: Reciprocal overlap (at least 50% of smaller peak must overlap)
    echo "  Finding reciprocal overlaps (>=50% overlap)..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" -wa -wb -f 0.5 -r > "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_reciprocal_50pct.bed" 2>/dev/null || true

    # Find unique peaks (no overlap)
    echo "  Finding unique regions..." >&2
    bedtools intersect -a "$OUTDIR/${cond1}_temp.bed" -b "$OUTDIR/${cond2}_temp.bed" -v > "$OUTDIR/unique/${cond1}_unique_vs_${cond2}.bed"
    bedtools intersect -a "$OUTDIR/${cond2}_temp.bed" -b "$OUTDIR/${cond1}_temp.bed" -v > "$OUTDIR/unique/${cond2}_unique_vs_${cond1}.bed"

    # Count results
    overlap_cond1=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_overlapping_${cond2}.bed")
    overlap_cond2=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_overlapping_${cond1}.bed")
    intersection_count=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_intersection.bed")
    unique_cond1=$(wc -l < "$OUTDIR/unique/${cond1}_unique_vs_${cond2}.bed")
    unique_cond2=$(wc -l < "$OUTDIR/unique/${cond2}_unique_vs_${cond1}.bed")
    total_cond1=$(wc -l < "$OUTDIR/${cond1}_temp.bed")
    total_cond2=$(wc -l < "$OUTDIR/${cond2}_temp.bed")

    # Summit-based overlap counts
    summit_overlap_cond1=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_summit_overlap_${cond2}.bed" 2>/dev/null || echo 0)
    summit_overlap_cond2=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond2}_summit_overlap_${cond1}.bed" 2>/dev/null || echo 0)

    # Reciprocal overlap count
    reciprocal_count=$(wc -l < "$OUTDIR/overlaps/${cond1}_vs_${cond2}/${cond1}_${cond2}_reciprocal_50pct.bed" 2>/dev/null || echo 0)

    # Calculate statistics
    overlap_percent_cond1=$(awk "BEGIN {printf \"%.1f\", ($overlap_cond1/$total_cond1)*100}")
    overlap_percent_cond2=$(awk "BEGIN {printf \"%.1f\", ($overlap_cond2/$total_cond2)*100}")

    # Calculate Jaccard index (intersection / union)
    union_count=$((total_cond1 + total_cond2 - overlap_cond1))
    jaccard=$(awk "BEGIN {printf \"%.3f\", ($overlap_cond1/$union_count)}")

    echo "  Results:" >&2
    echo "    $cond1 peaks overlapping with $cond2: $overlap_cond1 / $total_cond1 (${overlap_percent_cond1}%)" >&2
    echo "    $cond2 peaks overlapping with $cond1: $overlap_cond2 / $total_cond2 (${overlap_percent_cond2}%)" >&2
    echo "    Summit-based overlap (within 200bp): $summit_overlap_cond1 ($cond1), $summit_overlap_cond2 ($cond2)" >&2
    echo "    Reciprocal overlap (>=50%): $reciprocal_count pairs" >&2
    echo "    Jaccard index: $jaccard" >&2
    echo "    $cond1 unique peaks: $unique_cond1" >&2
    echo "    $cond2 unique peaks: $unique_cond2" >&2

    # Clean up temp files
    rm -f "$OUTDIR/${cond1}_temp.bed" "$OUTDIR/${cond2}_temp.bed"
    rm -f "$OUTDIR/${cond1}_summits.bed" "$OUTDIR/${cond2}_summits.bed"

    # Return statistics for summary (ONLY this line goes to STDOUT for parsing)
    echo "${cond1}_vs_${cond2}:${overlap_cond1}:${overlap_cond2}:${unique_cond1}:${unique_cond2}:${total_cond1}:${total_cond2}:${intersection_count}:${summit_overlap_cond1}:${summit_overlap_cond2}:${reciprocal_count}:${jaccard}"
}

# Perform all pairwise comparisons
echo ""
declare -A TOTAL_PEAKS
declare -A INTERSECTION_COUNTS
declare -A UNIQUE_PEAKS_1
declare -A UNIQUE_PEAKS_2

for condition in "${CONDITIONS[@]}"; do
    TOTAL_PEAKS[$condition]=$(wc -l < "${PEAK_FILES[$condition]}")
done

comparisons=("TES:TEAD1")  # TESmut removed - failed sample

# Additional arrays for new metrics
declare -A OVERLAP_COUNTS_1
declare -A OVERLAP_COUNTS_2
declare -A SUMMIT_OVERLAPS_1
declare -A SUMMIT_OVERLAPS_2
declare -A RECIPROCAL_COUNTS
declare -A JACCARD_INDICES

for comparison in "${comparisons[@]}"; do
    IFS=':' read -r cond1 cond2 <<< "$comparison"
    echo "Processing comparison: ${cond1} vs ${cond2}"
    stats_output=$(perform_pairwise_comparison "$cond1" "$cond2")

    # Parse the stats output (now only the stats line goes to STDOUT)
    IFS=':' read -r comp_name overlap1 overlap2 unique1 unique2 total1 total2 intersection_count summit1 summit2 reciprocal jaccard <<< "$stats_output"

    OVERLAP_COUNTS_1["${cond1}_${cond2}"]="$overlap1"
    OVERLAP_COUNTS_2["${cond1}_${cond2}"]="$overlap2"
    INTERSECTION_COUNTS["${cond1}_${cond2}"]="$intersection_count"
    UNIQUE_PEAKS_1["${cond1}_${cond2}"]="$unique1"
    UNIQUE_PEAKS_2["${cond1}_${cond2}"]="$unique2"
    SUMMIT_OVERLAPS_1["${cond1}_${cond2}"]="$summit1"
    SUMMIT_OVERLAPS_2["${cond1}_${cond2}"]="$summit2"
    RECIPROCAL_COUNTS["${cond1}_${cond2}"]="$reciprocal"
    JACCARD_INDICES["${cond1}_${cond2}"]="$jaccard"

    echo ""
done

echo "=== Creating Venn diagram data ==="

# Create R script for Venn diagrams and statistics
cat > "$OUTDIR/create_venn_diagrams.R" << EOF
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Peak counts from shell script (TESmut removed - failed sample)
tes_total=${TOTAL_PEAKS[TES]:-0}
tead1_total=${TOTAL_PEAKS[TEAD1]:-0}

# Direct overlap counts (peaks with any overlap)
tes_overlap_tead1=${OVERLAP_COUNTS_1[TES_TEAD1]:-0}
tead1_overlap_tes=${OVERLAP_COUNTS_2[TES_TEAD1]:-0}

# For Venn diagrams, use overlap counts (peaks that have any overlap)
tes_tead1_overlap=${OVERLAP_COUNTS_1[TES_TEAD1]:-0}

tes_unique_tead1=${UNIQUE_PEAKS_1[TES_TEAD1]:-0}
tead1_unique_tes=${UNIQUE_PEAKS_2[TES_TEAD1]:-0}

# Summit-based overlaps
summit_tes_tead1=${SUMMIT_OVERLAPS_1[TES_TEAD1]:-0}

# Reciprocal overlaps (>=50%)
reciprocal_tes_tead1=${RECIPROCAL_COUNTS[TES_TEAD1]:-0}

# Jaccard indices
jaccard_tes_tead1=${JACCARD_INDICES[TES_TEAD1]:-0}

cat(paste("Narrow Peak counts (TESmut removed - failed sample):\\n"))
cat(paste("TES:", tes_total, "\\n"))
cat(paste("TEAD1:", tead1_total, "\\n\\n"))

cat(paste("Direct Overlaps:\\n"))
cat(paste("TES overlapping TEAD1:", tes_overlap_tead1, "\\n\\n"))

# Create Venn diagrams for pairwise comparisons (TESmut removed - failed sample)
output_dir <- "results/12_pairwise_overlap_narrow/venn"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# TES vs TEAD1 (only comparison now - TESmut removed)
cat("Creating TES vs TEAD1 Venn diagram...\\n")
png(file.path(output_dir, "TES_vs_TEAD1_narrow_venn.png"),
    width = 3000, height = 3000, res = 300)
venn.plot <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tead1_total,
    cross.area = tes_tead1_overlap,
    category = c("TES", "TEAD1"),
    col = "transparent",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TEAD1 Narrow Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot)
dev.off()

# Create comprehensive summary statistics table (TESmut removed)
stats_table <- data.frame(
    Comparison = c("TES vs TEAD1"),
    Direct_Overlap = c(tes_tead1_overlap),
    Summit_Overlap = c(summit_tes_tead1),
    Reciprocal_50pct = c(reciprocal_tes_tead1),
    TES_Only = c(tes_unique_tead1),
    TEAD1_Only = c(tead1_unique_tes),
    TES_Total = c(tes_total),
    TEAD1_Total = c(tead1_total),
    Overlap_Pct_TES = round(c(tes_tead1_overlap/tes_total*100), 1),
    Overlap_Pct_TEAD1 = round(c(tes_tead1_overlap/tead1_total*100), 1),
    Jaccard_Index = c(jaccard_tes_tead1)
)

# Write statistics to file
write.table(stats_table, file.path(output_dir, "narrow_peak_overlap_statistics.txt"),
            sep="\\t", row.names=FALSE, quote=FALSE)

# Create bar plot comparing overlap methods
overlap_comparison <- data.frame(
    Comparison = rep(c("TES vs TEAD1"), 3),
    Method = c("Direct Overlap", "Summit-based (200bp)", "Reciprocal (50%)"),
    Count = c(tes_tead1_overlap, summit_tes_tead1, reciprocal_tes_tead1)
)

p_methods <- ggplot(overlap_comparison, aes(x = Comparison, y = Count, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparison of Overlap Methods for Narrow Peaks",
         y = "Number of Overlapping Peaks", x = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")

ggsave(file.path(output_dir, "overlap_methods_comparison.pdf"), p_methods, width = 10, height = 6)
ggsave(file.path(output_dir, "overlap_methods_comparison.png"), p_methods, width = 10, height = 6, dpi = 300)

cat("Narrow peak Venn diagram analysis completed (TESmut excluded - failed sample).\\n")
cat("\\nSummary Statistics:\\n")
print(stats_table)

cat("\\nInterpretation Guide:\\n")
cat("- Direct Overlap: Any genomic overlap (most permissive)\\n")
cat("- Summit-based: Peak centers within 200bp (appropriate for TF binding)\\n")
cat("- Reciprocal 50%: At least 50% mutual overlap (most stringent)\\n")
cat("- Jaccard Index: Overlap similarity score (0-1, higher = more similar)\\n")
cat("- NOTE: TESmut samples removed from analysis (failed sample)\\n")
EOF

# Run R script for Venn diagrams
echo "=== Creating Venn diagrams with R ==="
Rscript "$OUTDIR/create_venn_diagrams.R"

echo ""
echo "=== Creating comprehensive summary ==="

# Create comprehensive summary file
SUMMARY_FILE="$OUTDIR/PAIRWISE_OVERLAP_SUMMARY.md"

cat > "$SUMMARY_FILE" << EOF
# Pairwise Binding Site Overlap Analysis Summary (Narrow Peaks - Improved)

**Generated on:** $(date)
**Peak Type Used:** High-quality narrow peaks (consensus/IDR when available)
**Peak Files Used:**
EOF

# Add actual peak files used
for condition in "${CONDITIONS[@]}"; do
    echo "- **$condition**: ${PEAK_FILES[$condition]}" >> "$SUMMARY_FILE"
done

cat >> "$SUMMARY_FILE" << EOF

## Overview

This analysis compares binding sites between TES and TEAD1 conditions using pairwise overlap analysis.
NOTE: TESmut samples have been excluded from this analysis (failed sample).
The analysis prioritizes high-quality consensus peaks (from step 11) when available, falling back to MACS2
combined peaks if consensus analysis hasn't been run. This approach provides more reproducible and
statistically robust overlap estimates.

## Peak Counts by Condition

EOF

for condition in "${CONDITIONS[@]}"; do
    if [ -f "${PEAK_FILES[$condition]}" ]; then
        peak_count=$(wc -l < "${PEAK_FILES[$condition]}")
        echo "- **$condition**: $peak_count peaks" >> "$SUMMARY_FILE"
    fi
done

cat >> "$SUMMARY_FILE" << EOF

## Pairwise Overlap Results

NOTE: TESmut samples have been excluded from this analysis (failed sample).

### TES vs TEAD1
EOF

if [ -f "$OUTDIR/overlaps/TES_vs_TEAD1/TES_overlapping_TEAD1.bed" ]; then
    tes_overlap_tead1=$(wc -l < "$OUTDIR/overlaps/TES_vs_TEAD1/TES_overlapping_TEAD1.bed")
    tead1_overlap_tes=$(wc -l < "$OUTDIR/overlaps/TES_vs_TEAD1/TEAD1_overlapping_TES.bed")
    tes_unique_tead1=$(wc -l < "$OUTDIR/unique/TES_unique_vs_TEAD1.bed")
    tead1_unique_tes=$(wc -l < "$OUTDIR/unique/TEAD1_unique_vs_TES.bed")

    cat >> "$SUMMARY_FILE" << EOF
- TES peaks overlapping with TEAD1: **$tes_overlap_tead1**
- TEAD1 peaks overlapping with TES: **$tead1_overlap_tes**
- TES-specific peaks (not in TEAD1): **$tes_unique_tead1**
- TEAD1-specific peaks (not in TES): **$tead1_unique_tes**

EOF
fi

cat >> "$SUMMARY_FILE" << EOF

## Files Generated

### Overlap Files
- \`overlaps/\`: Peak regions that overlap between conditions
- \`unique/\`: Condition-specific peaks with no overlap

### Visualization Files
- \`venn/\`: Venn diagrams showing overlap relationships

### Statistics Files
- \`venn/narrow_peak_overlap_statistics.txt\`: Detailed numerical overlap statistics
- \`PAIRWISE_OVERLAP_SUMMARY.md\`: This comprehensive summary

## Interpretation Notes

1. **TES vs TEAD1**: Reveals similarity between TES and endogenous TEAD1 binding
2. **TESmut samples excluded**: Failed sample, not included in this analysis

## Next Steps

1. Review Venn diagrams for visual overlap assessment
2. Compare with broad peak overlap analysis for sensitivity differences
3. Perform functional annotation of condition-specific peak sets
4. Consider motif analysis in overlap vs unique regions

EOF

echo "Files created:"
echo "- Overlap regions: $OUTDIR/overlaps/"
echo "- Unique regions: $OUTDIR/unique/"
echo "- Venn diagrams: $OUTDIR/venn/"
echo "- Summary report: $SUMMARY_FILE"
echo ""
echo "=== Narrow peak pairwise overlap analysis completed ==="
echo "Results saved to: $OUTDIR"
echo "End time: $(date)"