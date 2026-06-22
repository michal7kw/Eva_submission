#!/usr/bin/env Rscript
#
# ENCODE ENHANCER METHYLATION ANALYSIS - DATA PREPARATION
# =============================================================================
#
# Purpose: Prepare ENCODE enhancer regions stratified by TES/TEAD1 binding status.
#
# Input:
#   - ENCODE dELS (distal enhancer-like signatures)
#   - TES peaks (narrowPeak)
#   - TEAD1 peaks (narrowPeak)
#
# Output:
#   - TES_bound_enhancers.bed (enhancers overlapping TES peaks)
#   - TEAD1_bound_enhancers.bed (enhancers overlapping TEAD1 peaks)
#   - TES_only_enhancers.bed (TES but not TEAD1)
#   - TEAD1_only_enhancers.bed (TEAD1 but not TES)
#   - Both_bound_enhancers.bed (both TES and TEAD1)
#   - Unbound_enhancers.bed (no TES or TEAD1 overlap)
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("ENCODE ENHANCER PREPARATION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# ENCODE enhancers (from download step)
ENCODE_ENHANCERS <- "output/32_encode_enhancer/ENCODE_distal_enhancers.bed"

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# Output directory
OUTPUT_DIR <- "output/32_encode_enhancer"

# =============================================================================
# PHASE 1: LOAD ENCODE ENHANCERS
# =============================================================================

cat("=== PHASE 1: Loading ENCODE Enhancers ===\n")

# ENCODE cCRE BED format (from UCSC encodeCcreCombined):
# V1: chr, V2: start, V3: end, V4: accession, V5: score, V6: strand,
# V7: thickStart, V8: thickEnd, V9: itemRgb, V10: type,CTCF, V11: simplified type
# V12: zscore, V13: category, V14: shortID, V15: description

# Read with basic columns first
encode_raw <- read.table(ENCODE_ENHANCERS, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")

cat(sprintf("  Raw ENCODE enhancers loaded: %d\n", nrow(encode_raw)))

# Create GRanges - using first 3 columns (chr, start, end)
encode_gr <- GRanges(
    seqnames = encode_raw$V1,
    ranges = IRanges(start = encode_raw$V2 + 1, end = encode_raw$V3)  # Convert to 1-based
)

# Add metadata
encode_gr$accession <- encode_raw$V4
encode_gr$cre_type <- encode_raw$V11

cat(sprintf("  ENCODE dELS enhancers: %d\n", length(encode_gr)))

# =============================================================================
# PHASE 2: LOAD TES AND TEAD1 PEAKS
# =============================================================================

cat("\n=== PHASE 2: Loading TES/TEAD1 Peaks ===\n")

# Load peaks function
load_peaks <- function(peak_file) {
    peaks <- read.table(peak_file,
        header = FALSE, stringsAsFactors = FALSE,
        col.names = c(
            "chr", "start", "end", "name", "score",
            "strand", "signalValue", "pValue", "qValue", "peak"
        )
    )
    # Fix chr names
    if (!grepl("^chr", peaks$chr[1])) peaks$chr <- paste0("chr", peaks$chr)

    GRanges(
        seqnames = peaks$chr,
        ranges = IRanges(peaks$start + 1, peaks$end),
        signal = peaks$signalValue
    )
}

tes_gr <- load_peaks(TES_PEAKS)
tead1_gr <- load_peaks(TEAD1_PEAKS)

cat(sprintf("  TES peaks: %d\n", length(tes_gr)))
cat(sprintf("  TEAD1 peaks: %d\n", length(tead1_gr)))

# =============================================================================
# PHASE 3: OVERLAP ENHANCERS WITH PEAKS
# =============================================================================

cat("\n=== PHASE 3: Overlapping Enhancers with Peaks ===\n")

# Find overlaps
tes_overlaps <- findOverlaps(encode_gr, tes_gr)
tead1_overlaps <- findOverlaps(encode_gr, tead1_gr)

# Get indices of enhancers with overlaps
tes_bound_idx <- unique(queryHits(tes_overlaps))
tead1_bound_idx <- unique(queryHits(tead1_overlaps))

cat(sprintf("  Enhancers overlapping TES peaks: %d (%.1f%%)\n",
    length(tes_bound_idx), 100*length(tes_bound_idx)/length(encode_gr)))
cat(sprintf("  Enhancers overlapping TEAD1 peaks: %d (%.1f%%)\n",
    length(tead1_bound_idx), 100*length(tead1_bound_idx)/length(encode_gr)))

# Create categories
both_idx <- intersect(tes_bound_idx, tead1_bound_idx)
tes_only_idx <- setdiff(tes_bound_idx, tead1_bound_idx)
tead1_only_idx <- setdiff(tead1_bound_idx, tes_bound_idx)
unbound_idx <- setdiff(1:length(encode_gr), union(tes_bound_idx, tead1_bound_idx))

cat(sprintf("\n  Category breakdown:\n"))
cat(sprintf("    TES + TEAD1 (both): %d (%.1f%%)\n",
    length(both_idx), 100*length(both_idx)/length(encode_gr)))
cat(sprintf("    TES only: %d (%.1f%%)\n",
    length(tes_only_idx), 100*length(tes_only_idx)/length(encode_gr)))
cat(sprintf("    TEAD1 only: %d (%.1f%%)\n",
    length(tead1_only_idx), 100*length(tead1_only_idx)/length(encode_gr)))
cat(sprintf("    Unbound: %d (%.1f%%)\n",
    length(unbound_idx), 100*length(unbound_idx)/length(encode_gr)))

# =============================================================================
# PHASE 4: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 4: Exporting BED Files ===\n")

# Function to write BED (without scientific notation)
write_bed <- function(gr, filename) {
    # Disable scientific notation
    options(scipen = 999)

    bed_df <- data.frame(
        chr = as.character(seqnames(gr)),
        start = as.integer(start(gr) - 1),  # Convert to 0-based
        end = as.integer(end(gr)),
        name = if (!is.null(gr$accession)) gr$accession else paste0("enh_", 1:length(gr)),
        score = 0,
        strand = "."
    )
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d regions)\n", basename(filename), nrow(bed_df)))
}

# Export all categories
write_bed(encode_gr[tes_bound_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers.bed"))
write_bed(encode_gr[tead1_bound_idx], file.path(OUTPUT_DIR, "TEAD1_bound_enhancers.bed"))
write_bed(encode_gr[tes_only_idx], file.path(OUTPUT_DIR, "TES_only_enhancers.bed"))
write_bed(encode_gr[tead1_only_idx], file.path(OUTPUT_DIR, "TEAD1_only_enhancers.bed"))
write_bed(encode_gr[both_idx], file.path(OUTPUT_DIR, "Both_bound_enhancers.bed"))
write_bed(encode_gr[unbound_idx], file.path(OUTPUT_DIR, "Unbound_enhancers.bed"))

# =============================================================================
# PHASE 5: SUMMARY STATISTICS
# =============================================================================

cat("\n=== PHASE 5: Summary Statistics ===\n")

summary_df <- data.frame(
    Category = c("TES_bound", "TEAD1_bound", "TES_only", "TEAD1_only", "Both_TES_TEAD1", "Unbound", "Total"),
    N_enhancers = c(
        length(tes_bound_idx),
        length(tead1_bound_idx),
        length(tes_only_idx),
        length(tead1_only_idx),
        length(both_idx),
        length(unbound_idx),
        length(encode_gr)
    ),
    Percentage = c(
        100*length(tes_bound_idx)/length(encode_gr),
        100*length(tead1_bound_idx)/length(encode_gr),
        100*length(tes_only_idx)/length(encode_gr),
        100*length(tead1_only_idx)/length(encode_gr),
        100*length(both_idx)/length(encode_gr),
        100*length(unbound_idx)/length(encode_gr),
        100
    )
)

print(summary_df)

write.csv(summary_df, file.path(OUTPUT_DIR, "enhancer_binding_summary.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("PREPARATION COMPLETE\n")
cat("==========================================================\n")
