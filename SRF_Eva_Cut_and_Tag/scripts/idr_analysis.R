#!/usr/bin/env Rscript

# ===============================================================================
# SCRIPT: idr_analysis.R
# PURPOSE: Irreproducible Discovery Rate (IDR) analysis for Cut&Tag replicates
#
# DESCRIPTION:
# This script performs IDR analysis to assess reproducibility between
# biological replicates. IDR is crucial for ensuring that peaks are
# consistently detected across replicates before downstream analysis.
# ===============================================================================

# Load required packages with graceful handling
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(ggplot2)
    library(RColorBrewer)
    library(dplyr)
})

# Try to load corrplot, use base R heatmap if not available
use_corrplot <- requireNamespace("corrplot", quietly = TRUE)
if (use_corrplot) {
    library(corrplot)
} else {
    cat("Note: corrplot package not available, using base R heatmap instead\n")
}

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Create output directory
dir.create("results/idr_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to read narrowPeak files
read_narrowpeak <- function(file_path) {
    if (!file.exists(file_path)) {
        cat("Warning: File not found:", file_path, "\n")
        return(NULL)
    }

    peaks <- read.table(file_path, sep = "\t", stringsAsFactors = FALSE)
    colnames(peaks) <- c(
        "chr", "start", "end", "name", "score", "strand",
        "signalValue", "pValue", "qValue", "peak"
    )

    # Create GRanges object
    gr <- GRanges(
        seqnames = peaks$chr,
        ranges = IRanges(start = peaks$start + 1, end = peaks$end),
        score = peaks$score,
        signalValue = peaks$signalValue,
        pValue = peaks$pValue,
        qValue = peaks$qValue
    )

    return(gr)
}

# Function to calculate overlap between two peak sets
calculate_overlap <- function(peaks1, peaks2, min_overlap = 0.5) {
    if (is.null(peaks1) || is.null(peaks2)) {
        return(0)
    }

    overlaps <- findOverlaps(peaks1, peaks2)

    # Calculate reciprocal overlap
    hits1 <- queryHits(overlaps)
    hits2 <- subjectHits(overlaps)

    # Calculate overlap percentages
    overlap_width <- width(pintersect(peaks1[hits1], peaks2[hits2]))
    width1 <- width(peaks1[hits1])
    width2 <- width(peaks2[hits2])

    # Reciprocal overlap: min(overlap/width1, overlap/width2) >= threshold
    reciprocal_overlap <- pmin(overlap_width / width1, overlap_width / width2) >= min_overlap

    # Count high-quality overlaps
    n_overlaps <- sum(reciprocal_overlap)

    return(list(
        n_overlaps = n_overlaps,
        total_peaks1 = length(peaks1),
        total_peaks2 = length(peaks2),
        overlap_rate1 = n_overlaps / length(peaks1),
        overlap_rate2 = n_overlaps / length(peaks2)
    ))
}

# Define sample groups and replicates (TES and TEAD1 only - TESmut removed)
groups <- list(
    TES = c("TES-1", "TES-2", "TES-3"),
    TEAD1 = c("TEAD1-1", "TEAD1-2", "TEAD1-3")
)

peaks_dir <- "results/05_peaks_narrow"

# Read all peak files
all_peaks <- list()
for (group_name in names(groups)) {
    for (sample in groups[[group_name]]) {
        peak_file <- file.path(peaks_dir, paste0(sample, "_peaks.narrowPeak"))
        all_peaks[[sample]] <- read_narrowpeak(peak_file)
    }
}

# Remove NULL entries
all_peaks <- all_peaks[!sapply(all_peaks, is.null)]

# Calculate pairwise overlaps within groups
overlap_results <- list()
correlation_matrix <- matrix(NA, nrow = length(all_peaks), ncol = length(all_peaks))
rownames(correlation_matrix) <- names(all_peaks)
colnames(correlation_matrix) <- names(all_peaks)

for (group_name in names(groups)) {
    group_samples <- groups[[group_name]]
    group_samples <- group_samples[group_samples %in% names(all_peaks)]

    if (length(group_samples) < 2) next

    cat("Analyzing group:", group_name, "\n")

    # Pairwise comparisons within group
    for (i in 1:(length(group_samples) - 1)) {
        for (j in (i + 1):length(group_samples)) {
            sample1 <- group_samples[i]
            sample2 <- group_samples[j]

            overlap_res <- calculate_overlap(all_peaks[[sample1]], all_peaks[[sample2]])

            # Store results
            result_key <- paste(sample1, sample2, sep = "_vs_")
            overlap_results[[result_key]] <- overlap_res

            # Calculate correlation based on overlap rates
            correlation <- sqrt(overlap_res$overlap_rate1 * overlap_res$overlap_rate2)
            correlation_matrix[sample1, sample2] <- correlation
            correlation_matrix[sample2, sample1] <- correlation

            cat(sprintf(
                "%s vs %s: %d overlaps, %.2f%% and %.2f%% overlap rates, correlation: %.3f\n",
                sample1, sample2, overlap_res$n_overlaps,
                overlap_res$overlap_rate1 * 100, overlap_res$overlap_rate2 * 100,
                correlation
            ))
        }
    }
    cat("\n")
}

# Set diagonal to 1
diag(correlation_matrix) <- 1

# Create summary statistics
overlap_summary <- data.frame(
    comparison = names(overlap_results),
    group = sapply(names(overlap_results), function(x) {
        sample1 <- strsplit(x, "_vs_")[[1]][1]
        if (sample1 %in% groups$TES) {
            return("TES")
        }
        if (sample1 %in% groups$TEAD1) {
            return("TEAD1")
        }
        return("Unknown")
    }),
    n_overlaps = sapply(overlap_results, function(x) x$n_overlaps),
    total_peaks1 = sapply(overlap_results, function(x) x$total_peaks1),
    total_peaks2 = sapply(overlap_results, function(x) x$total_peaks2),
    overlap_rate1 = sapply(overlap_results, function(x) x$overlap_rate1),
    overlap_rate2 = sapply(overlap_results, function(x) x$overlap_rate2),
    avg_overlap_rate = sapply(overlap_results, function(x) (x$overlap_rate1 + x$overlap_rate2) / 2),
    correlation = sapply(overlap_results, function(x) sqrt(x$overlap_rate1 * x$overlap_rate2)),
    stringsAsFactors = FALSE
)

# Save results
write.csv(overlap_summary, "results/idr_analysis/replicate_overlap_summary.csv", row.names = FALSE)

# Create visualizations
pdf("results/idr_analysis/correlation_matrix.pdf", width = 10, height = 8)
if (use_corrplot) {
    corrplot(correlation_matrix,
        method = "color", type = "upper",
        order = "hclust", tl.cex = 0.8, tl.col = "black",
        col = brewer.pal(n = 10, name = "RdYlBu")
    )
    title("Replicate Correlation Matrix\n(Based on Peak Overlap)")
} else {
    # Use base R heatmap as fallback
    heatmap(correlation_matrix,
        col = colorRampPalette(brewer.pal(9, "RdYlBu"))(100),
        scale = "none",
        main = "Replicate Correlation Matrix\n(Based on Peak Overlap)")
}
dev.off()

# Overlap rate boxplot
p1 <- ggplot(overlap_summary, aes(x = group, y = avg_overlap_rate, fill = group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.8, size = 2) +
    labs(
        title = "Average Overlap Rate Between Replicates",
        x = "Experimental Group", y = "Average Overlap Rate"
    ) +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme(legend.position = "none")

ggsave("results/idr_analysis/overlap_rates_by_group.pdf", p1, width = 8, height = 6)
ggsave("results/idr_analysis/overlap_rates_by_group.png", p1, width = 8, height = 6, dpi = 300)

# Correlation distribution
p2 <- ggplot(overlap_summary, aes(x = correlation, fill = group)) +
    geom_histogram(bins = 20, alpha = 0.7, position = "dodge") +
    labs(
        title = "Distribution of Replicate Correlations",
        x = "Correlation (sqrt of product of overlap rates)", y = "Count"
    ) +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set1")

ggsave("results/idr_analysis/correlation_distribution.pdf", p2, width = 10, height = 6)
ggsave("results/idr_analysis/correlation_distribution.png", p2, width = 10, height = 6, dpi = 300)

# Quality assessment
cat("\n=== REPLICATE REPRODUCIBILITY ASSESSMENT ===\n")
cat("Quality thresholds:\n")
cat("- Excellent reproducibility: >70% overlap\n")
cat("- Good reproducibility: 50-70% overlap\n")
cat("- Moderate reproducibility: 30-50% overlap\n")
cat("- Poor reproducibility: <30% overlap\n\n")

# Group-wise summary
group_summary <- overlap_summary %>%
    group_by(group) %>%
    summarise(
        n_comparisons = n(),
        mean_overlap = mean(avg_overlap_rate, na.rm = TRUE),
        min_overlap = min(avg_overlap_rate, na.rm = TRUE),
        max_overlap = max(avg_overlap_rate, na.rm = TRUE),
        .groups = "drop"
    )

print(group_summary)

# Flag quality issues
cat("\nQuality assessment by group:\n")
for (i in 1:nrow(group_summary)) {
    group <- group_summary$group[i]
    mean_overlap <- group_summary$mean_overlap[i]
    min_overlap <- group_summary$min_overlap[i]

    if (mean_overlap > 0.7) {
        cat("EXCELLENT:", group, "- Mean overlap:", sprintf("%.1f%%", mean_overlap * 100), "\n")
    } else if (mean_overlap > 0.5) {
        cat("GOOD:", group, "- Mean overlap:", sprintf("%.1f%%", mean_overlap * 100), "\n")
    } else if (mean_overlap > 0.3) {
        cat("MODERATE:", group, "- Mean overlap:", sprintf("%.1f%%", mean_overlap * 100), "\n")
        if (min_overlap < 0.3) {
            cat("  WARNING: Some replicate pairs have <30% overlap\n")
        }
    } else {
        cat("POOR:", group, "- Mean overlap:", sprintf("%.1f%%", mean_overlap * 100), "\n")
        cat("  CRITICAL: Consider excluding low-quality replicates\n")
    }
}

cat("\nIDR analysis complete! Check results/idr_analysis/ for detailed results.\n")
cat("Consider replicates with <30% overlap for exclusion from downstream analysis.\n")
