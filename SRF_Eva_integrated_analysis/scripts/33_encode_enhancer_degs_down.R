#!/usr/bin/env Rscript
#
# ENCODE ENHANCERS OF DEGs DOWN - METHYLATION ANALYSIS
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE enhancers specifically associated
#          with downregulated DEGs, stratified by TES/TEAD1 binding status.
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ChIPseeker)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(dplyr)
})

# Disable scientific notation for BED output
options(scipen = 999)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("ENCODE ENHANCERS OF DEGs DOWN - PREPARATION\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

# ENCODE enhancers (from script 32)
ENCODE_ENHANCERS <- "output/32_encode_enhancer/ENCODE_distal_enhancers.bed"

# Peak files
TES_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/05_peaks_narrow/TES_peaks.narrowPeak"
TEAD1_PEAKS <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/05_peaks_narrow/TEAD1_peaks.narrowPeak"

# DESeq2 results
DESEQ2_FILE <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"

# Output directory
OUTPUT_DIR <- "output/33_encode_enhancer_degs_down"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD ENCODE ENHANCERS
# =============================================================================

cat("=== PHASE 1: Loading ENCODE Enhancers ===\n")

# Read ENCODE enhancers (already filtered for dELS)
encode_raw <- read.table(ENCODE_ENHANCERS, header = FALSE, stringsAsFactors = FALSE, fill = TRUE, sep = "\t")

cat(sprintf("  Raw ENCODE enhancers loaded: %d\n", nrow(encode_raw)))

# Create GRanges
encode_gr <- GRanges(
    seqnames = encode_raw$V1,
    ranges = IRanges(start = encode_raw$V2 + 1, end = encode_raw$V3)
)
encode_gr$accession <- encode_raw$V4

cat(sprintf("  ENCODE dELS enhancers: %d\n", length(encode_gr)))

# =============================================================================
# PHASE 2: ANNOTATE ENHANCERS TO NEAREST GENES
# =============================================================================

cat("\n=== PHASE 2: Annotating Enhancers to Nearest Genes ===\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate to nearest genes
cat("  Running ChIPseeker annotation (this may take a few minutes)...\n")
peak_anno <- annotatePeak(encode_gr, TxDb = txdb, annoDb = "org.Hs.eg.db", verbose = FALSE)
peak_anno_df <- as.data.frame(peak_anno)

cat(sprintf("  Annotated enhancers: %d\n", nrow(peak_anno_df)))

# Check annotation categories
cat("  Annotation categories:\n")
print(table(gsub(" \\(.*", "", peak_anno_df$annotation)))

# =============================================================================
# PHASE 3: LOAD DESeq2 RESULTS AND IDENTIFY DEGs DOWN
# =============================================================================

cat("\n=== PHASE 3: Loading DESeq2 Results ===\n")

deseq2 <- read.delim(DESEQ2_FILE, stringsAsFactors = FALSE)
cat(sprintf("  Total genes in DESeq2: %d\n", nrow(deseq2)))

# DEGs DOWN: padj < 0.05 AND log2FoldChange < 0
degs_down_genes <- deseq2 %>%
    filter(padj < 0.05, log2FoldChange < 0) %>%
    pull(gene_symbol)

cat(sprintf("  DEGs DOWN genes (padj < 0.05, log2FC < 0): %d\n", length(degs_down_genes)))

# Control genes: expressed but NOT differentially expressed
control_genes <- deseq2 %>%
    filter(!is.na(padj), padj >= 0.05) %>%
    pull(gene_symbol)

cat(sprintf("  Control genes (non-DE): %d\n", length(control_genes)))

# =============================================================================
# PHASE 4: FILTER ENHANCERS FOR DEGs DOWN
# =============================================================================

cat("\n=== PHASE 4: Filtering Enhancers for DEGs DOWN ===\n")

# Enhancers of DEGs DOWN
enhancers_degs_down <- peak_anno_df %>%
    filter(SYMBOL %in% degs_down_genes)

cat(sprintf("  Enhancers linked to DEGs DOWN: %d\n", nrow(enhancers_degs_down)))

# Enhancers of Control genes
enhancers_control_pool <- peak_anno_df %>%
    filter(SYMBOL %in% control_genes)

cat(sprintf("  Enhancers linked to Control genes (pool): %d\n", nrow(enhancers_control_pool)))

# =============================================================================
# PHASE 5: LOAD TES/TEAD1 PEAKS AND OVERLAP WITH ENHANCERS
# =============================================================================

cat("\n=== PHASE 5: Overlapping with TES/TEAD1 Peaks ===\n")

# Load peaks function
load_peaks <- function(peak_file) {
    peaks <- read.table(peak_file,
        header = FALSE, stringsAsFactors = FALSE,
        col.names = c(
            "chr", "start", "end", "name", "score",
            "strand", "signalValue", "pValue", "qValue", "peak"
        )
    )
    if (!grepl("^chr", peaks$chr[1])) peaks$chr <- paste0("chr", peaks$chr)
    GRanges(seqnames = peaks$chr, ranges = IRanges(peaks$start + 1, peaks$end))
}

tes_gr <- load_peaks(TES_PEAKS)
tead1_gr <- load_peaks(TEAD1_PEAKS)

cat(sprintf("  TES peaks: %d\n", length(tes_gr)))
cat(sprintf("  TEAD1 peaks: %d\n", length(tead1_gr)))

# Convert enhancers to GRanges
enhancers_degs_down_gr <- GRanges(
    seqnames = enhancers_degs_down$seqnames,
    ranges = IRanges(start = enhancers_degs_down$start, end = enhancers_degs_down$end)
)
enhancers_degs_down_gr$SYMBOL <- enhancers_degs_down$SYMBOL
enhancers_degs_down_gr$accession <- enhancers_degs_down$accession

# Find overlaps
tes_overlaps <- findOverlaps(enhancers_degs_down_gr, tes_gr)
tead1_overlaps <- findOverlaps(enhancers_degs_down_gr, tead1_gr)

tes_bound_idx <- unique(queryHits(tes_overlaps))
tead1_bound_idx <- unique(queryHits(tead1_overlaps))

# Create categories
both_idx <- intersect(tes_bound_idx, tead1_bound_idx)
tes_only_idx <- setdiff(tes_bound_idx, tead1_bound_idx)
tead1_only_idx <- setdiff(tead1_bound_idx, tes_bound_idx)
unbound_idx <- setdiff(1:length(enhancers_degs_down_gr), union(tes_bound_idx, tead1_bound_idx))

cat(sprintf("\n  DEGs DOWN Enhancers by binding status:\n"))
cat(sprintf("    TES-bound: %d (%.1f%%)\n", length(tes_bound_idx), 100*length(tes_bound_idx)/length(enhancers_degs_down_gr)))
cat(sprintf("    TEAD1-bound: %d (%.1f%%)\n", length(tead1_bound_idx), 100*length(tead1_bound_idx)/length(enhancers_degs_down_gr)))
cat(sprintf("    TES-only: %d (%.1f%%)\n", length(tes_only_idx), 100*length(tes_only_idx)/length(enhancers_degs_down_gr)))
cat(sprintf("    TEAD1-only: %d (%.1f%%)\n", length(tead1_only_idx), 100*length(tead1_only_idx)/length(enhancers_degs_down_gr)))
cat(sprintf("    Both: %d (%.1f%%)\n", length(both_idx), 100*length(both_idx)/length(enhancers_degs_down_gr)))
cat(sprintf("    Unbound: %d (%.1f%%)\n", length(unbound_idx), 100*length(unbound_idx)/length(enhancers_degs_down_gr)))

# =============================================================================
# PHASE 6: CREATE CONTROL SET (MATCHED SIZE)
# =============================================================================

cat("\n=== PHASE 6: Creating Control Set ===\n")

# Match size to TES-bound DEGs DOWN enhancers
n_target <- length(tes_bound_idx)
set.seed(42)

if (nrow(enhancers_control_pool) > n_target) {
    control_sample_idx <- sample(1:nrow(enhancers_control_pool), n_target)
    enhancers_control <- enhancers_control_pool[control_sample_idx, ]
    cat(sprintf("  Subsampled Control enhancers: %d\n", nrow(enhancers_control)))
} else {
    enhancers_control <- enhancers_control_pool
    cat(sprintf("  Using all Control enhancers: %d\n", nrow(enhancers_control)))
}

# =============================================================================
# PHASE 7: EXPORT BED FILES
# =============================================================================

cat("\n=== PHASE 7: Exporting BED Files ===\n")

# Function to write BED
write_bed <- function(gr, filename) {
    bed_df <- data.frame(
        chr = as.character(seqnames(gr)),
        start = as.integer(start(gr) - 1),
        end = as.integer(end(gr)),
        name = if (!is.null(gr$SYMBOL)) gr$SYMBOL else if (!is.null(gr$accession)) gr$accession else paste0("enh_", 1:length(gr)),
        score = 0,
        strand = "."
    )
    bed_df <- bed_df[order(bed_df$chr, bed_df$start), ]
    write.table(bed_df, filename, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    cat(sprintf("  Created: %s (%d regions)\n", basename(filename), nrow(bed_df)))
}

# Write DEGs DOWN enhancer subsets
write_bed(enhancers_degs_down_gr[tes_bound_idx], file.path(OUTPUT_DIR, "TES_bound_enhancers_DEGs_DOWN.bed"))
write_bed(enhancers_degs_down_gr[tead1_only_idx], file.path(OUTPUT_DIR, "TEAD1_only_enhancers_DEGs_DOWN.bed"))
write_bed(enhancers_degs_down_gr[unbound_idx], file.path(OUTPUT_DIR, "Unbound_enhancers_DEGs_DOWN.bed"))

# Write Control enhancers
control_gr <- GRanges(
    seqnames = enhancers_control$seqnames,
    ranges = IRanges(start = enhancers_control$start, end = enhancers_control$end)
)
control_gr$SYMBOL <- enhancers_control$SYMBOL
write_bed(control_gr, file.path(OUTPUT_DIR, "Control_enhancers.bed"))

# =============================================================================
# PHASE 8: SUMMARY STATISTICS
# =============================================================================

cat("\n=== PHASE 8: Summary Statistics ===\n")

summary_df <- data.frame(
    Category = c("TES_bound_DEGs_DOWN", "TEAD1_only_DEGs_DOWN", "Unbound_DEGs_DOWN", "Control"),
    N_enhancers = c(length(tes_bound_idx), length(tead1_only_idx), length(unbound_idx), nrow(enhancers_control)),
    Description = c(
        "ENCODE enhancers of DEGs DOWN, overlapping TES peaks",
        "ENCODE enhancers of DEGs DOWN, overlapping TEAD1 only",
        "ENCODE enhancers of DEGs DOWN, no TF overlap",
        "ENCODE enhancers of non-DE genes (matched control)"
    )
)

print(summary_df)

write.csv(summary_df, file.path(OUTPUT_DIR, "enhancer_summary.csv"), row.names = FALSE)

cat("\n")
cat("==========================================================\n")
cat("PREPARATION COMPLETE\n")
cat("==========================================================\n")
