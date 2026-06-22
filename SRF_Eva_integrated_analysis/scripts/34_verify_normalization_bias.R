#!/usr/bin/env Rscript
# meDIP NOTE: reads DNA-methylation (meDIP) data from the ARCHIVED tree SRF_Eva_top/Archive/meDIP, which is NOT self-contained in Eva_submission. If the Archive is moved or removed, repoint or copy the meDIP inputs before running.
#
# VERIFY NORMALIZATION BIAS IN meDIP ANALYSIS
# =============================================================================
#
# Purpose: Statistically verify whether observed hypomethylation at TES-bound
#          enhancers could be an artifact of INPUT normalization
#
# Three verification tasks:
# 1. Quantify INPUT signal differences at TES-bound vs Unbound enhancers
# 2. Analyze INPUT-subtracted methylation data
# 3. Compare MEDIPS IP-only vs INPUT-normalized DMR results
#
# =============================================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(dplyr)
    library(ggplot2)
    library(VennDiagram)
})

options(scipen = 999)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==========================================================\n")
cat("NORMALIZATION BIAS VERIFICATION - STATISTICAL ANALYSIS\n")
cat("==========================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

OUTPUT_DIR <- "output/34_verify_normalization_bias"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Matrix files created by the shell script
INPUT_MATRIX <- file.path(OUTPUT_DIR, "input_signal_matrix.gz")
INPUT_SUB_MATRIX <- file.path(OUTPUT_DIR, "input_subtracted_enhancer_matrix.gz")

# MEDIPS results
MEDIPS_IP_ONLY <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/07_differential_MEDIPS"
MEDIPS_INPUT_NORM <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/07_differential_MEDIPS_INPUT_normalized"

# =============================================================================
# HELPER FUNCTION: Parse deepTools matrix
# =============================================================================

parse_deeptools_matrix <- function(matrix_file) {
    # Read header line
    con <- gzfile(matrix_file, "r")
    header_line <- readLines(con, n = 1)
    close(con)

    # Parse JSON header
    header_json <- gsub("^@", "", header_line)
    header <- jsonlite::fromJSON(header_json)

    # Read data (skip header)
    mat_data <- read.table(gzfile(matrix_file), skip = 1, header = FALSE)

    # Extract sample and group info
    sample_labels <- header$sample_labels
    group_labels <- header$group_labels
    group_boundaries <- header$group_boundaries

    # Number of bins
    n_bins <- (header$upstream + header$downstream) / header$`bin size`

    # Matrix starts at column 7 (1-6 are chr, start, end, name, score, strand)
    mat_values <- as.matrix(mat_data[, 7:ncol(mat_data)])

    return(list(
        matrix = mat_values,
        samples = sample_labels,
        groups = group_labels,
        boundaries = group_boundaries,
        n_bins = n_bins,
        n_samples = length(sample_labels)
    ))
}

# =============================================================================
# TASK 1: QUANTIFY INPUT SIGNAL AT ENHANCERS
# =============================================================================

cat("=== TASK 1: INPUT Signal Analysis ===\n\n")

if (file.exists(INPUT_MATRIX)) {
    mat_info <- parse_deeptools_matrix(INPUT_MATRIX)

    cat(sprintf("Matrix info:\n"))
    cat(sprintf("  Samples: %s\n", paste(mat_info$samples, collapse = ", ")))
    cat(sprintf("  Groups: %s\n", paste(mat_info$groups, collapse = ", ")))
    cat(sprintf("  Group boundaries: %s\n", paste(mat_info$boundaries, collapse = ", ")))
    cat("\n")

    # Extract center region (±500bp around center)
    center_idx <- round(mat_info$n_bins / 2)
    center_range <- (center_idx - 10):(center_idx + 10)  # ±500bp with 50bp bins

    # Group indices
    n_regions <- nrow(mat_info$matrix)
    g1_end <- mat_info$boundaries[2]

    group1_idx <- 1:g1_end  # TES-bound enhancers
    group2_idx <- (g1_end + 1):n_regions  # Unbound enhancers

    cat(sprintf("  Group 1 (TES-bound): %d enhancers\n", length(group1_idx)))
    cat(sprintf("  Group 2 (Unbound): %d enhancers\n", length(group2_idx)))
    cat("\n")

    # Calculate mean INPUT signal at center for each group
    input_col_idx <- 1:mat_info$n_bins  # First sample's bins

    group1_center <- mat_info$matrix[group1_idx, input_col_idx[center_range]]
    group2_center <- mat_info$matrix[group2_idx, input_col_idx[center_range]]

    # Mean per region
    group1_mean <- rowMeans(group1_center, na.rm = TRUE)
    group2_mean <- rowMeans(group2_center, na.rm = TRUE)

    cat("========================================\n")
    cat("INPUT SIGNAL AT ENHANCERS (Center ±500bp)\n")
    cat("========================================\n\n")

    cat("TES-bound enhancers:\n")
    cat(sprintf("  Mean INPUT signal: %.4f\n", mean(group1_mean, na.rm = TRUE)))
    cat(sprintf("  Median INPUT signal: %.4f\n", median(group1_mean, na.rm = TRUE)))
    cat(sprintf("  SD: %.4f\n", sd(group1_mean, na.rm = TRUE)))
    cat("\n")

    cat("Unbound enhancers:\n")
    cat(sprintf("  Mean INPUT signal: %.4f\n", mean(group2_mean, na.rm = TRUE)))
    cat(sprintf("  Median INPUT signal: %.4f\n", median(group2_mean, na.rm = TRUE)))
    cat(sprintf("  SD: %.4f\n", sd(group2_mean, na.rm = TRUE)))
    cat("\n")

    # Statistical test
    wilcox_result <- wilcox.test(group1_mean, group2_mean)
    t_result <- t.test(group1_mean, group2_mean)

    cat("========================================\n")
    cat("STATISTICAL COMPARISON\n")
    cat("========================================\n\n")

    cat(sprintf("Wilcoxon test p-value: %.4e\n", wilcox_result$p.value))
    cat(sprintf("T-test p-value: %.4e\n", t_result$p.value))

    diff_mean <- mean(group1_mean, na.rm = TRUE) - mean(group2_mean, na.rm = TRUE)

    if (abs(diff_mean) < 0.1 && wilcox_result$p.value > 0.05) {
        cat("\n*** CONCLUSION: INPUT signal is SIMILAR at TES-bound and Unbound enhancers\n")
        cat("*** This suggests the hypomethylation pattern is NOT an INPUT bias artifact\n")
        input_bias <- "NO_BIAS"
    } else if (diff_mean < 0) {
        cat(sprintf("\n*** CONCLUSION: INPUT signal is LOWER at TES-bound enhancers (diff = %.4f)\n", diff_mean))
        cat("*** This could contribute to apparent hypomethylation (POTENTIAL BIAS)\n")
        input_bias <- "POTENTIAL_BIAS_LOWER"
    } else {
        cat(sprintf("\n*** CONCLUSION: INPUT signal is HIGHER at TES-bound enhancers (diff = %.4f)\n", diff_mean))
        cat("*** This would actually MASK hypomethylation - the real effect may be stronger!\n")
        input_bias <- "POTENTIAL_BIAS_HIGHER"
    }
    cat("\n")

    # Save INPUT signal statistics
    input_stats <- data.frame(
        Group = c("TES_bound", "Unbound"),
        Mean_INPUT = c(mean(group1_mean, na.rm = TRUE), mean(group2_mean, na.rm = TRUE)),
        Median_INPUT = c(median(group1_mean, na.rm = TRUE), median(group2_mean, na.rm = TRUE)),
        SD_INPUT = c(sd(group1_mean, na.rm = TRUE), sd(group2_mean, na.rm = TRUE)),
        N = c(length(group1_idx), length(group2_idx))
    )
    write.csv(input_stats, file.path(OUTPUT_DIR, "input_signal_statistics.csv"), row.names = FALSE)

} else {
    cat("INPUT matrix file not found - skipping Task 1\n")
    input_bias <- "NOT_TESTED"
}

# =============================================================================
# TASK 2: ANALYZE INPUT-SUBTRACTED DATA
# =============================================================================

cat("\n=== TASK 2: INPUT-Subtracted Analysis ===\n\n")

if (file.exists(INPUT_SUB_MATRIX)) {
    mat_info2 <- parse_deeptools_matrix(INPUT_SUB_MATRIX)

    cat(sprintf("Matrix info:\n"))
    cat(sprintf("  Samples: %s\n", paste(mat_info2$samples, collapse = ", ")))
    cat(sprintf("  Groups: %s\n", paste(mat_info2$groups, collapse = ", ")))
    cat("\n")

    # Extract center region
    center_idx <- round(mat_info2$n_bins / 2)
    center_range <- (center_idx - 10):(center_idx + 10)

    # Group indices
    n_regions <- nrow(mat_info2$matrix)
    g1_end <- mat_info2$boundaries[2]

    group1_idx <- 1:g1_end
    group2_idx <- (g1_end + 1):n_regions

    # Sample 1 = TES-INPUT, Sample 2 = GFP-INPUT
    tes_cols <- 1:mat_info2$n_bins
    gfp_cols <- (mat_info2$n_bins + 1):(2 * mat_info2$n_bins)

    # TES-bound enhancers
    g1_tes <- rowMeans(mat_info2$matrix[group1_idx, tes_cols[center_range]], na.rm = TRUE)
    g1_gfp <- rowMeans(mat_info2$matrix[group1_idx, gfp_cols[center_range]], na.rm = TRUE)
    g1_diff <- g1_tes - g1_gfp

    # Unbound enhancers
    g2_tes <- rowMeans(mat_info2$matrix[group2_idx, tes_cols[center_range]], na.rm = TRUE)
    g2_gfp <- rowMeans(mat_info2$matrix[group2_idx, gfp_cols[center_range]], na.rm = TRUE)
    g2_diff <- g2_tes - g2_gfp

    cat("========================================\n")
    cat("INPUT-SUBTRACTED METHYLATION DIFFERENCE\n")
    cat("(TES-INPUT) - (GFP-INPUT) at Center ±500bp\n")
    cat("========================================\n\n")

    cat("TES-bound enhancers:\n")
    cat(sprintf("  Mean diff:   %.4f\n", mean(g1_diff, na.rm = TRUE)))
    cat(sprintf("  Median diff: %.4f\n", median(g1_diff, na.rm = TRUE)))
    cat(sprintf("  %% hypermethylated: %.1f%%\n", 100 * mean(g1_diff > 0, na.rm = TRUE)))
    cat(sprintf("  %% hypomethylated:  %.1f%%\n", 100 * mean(g1_diff < 0, na.rm = TRUE)))
    cat("\n")

    cat("Unbound enhancers:\n")
    cat(sprintf("  Mean diff:   %.4f\n", mean(g2_diff, na.rm = TRUE)))
    cat(sprintf("  Median diff: %.4f\n", median(g2_diff, na.rm = TRUE)))
    cat(sprintf("  %% hypermethylated: %.1f%%\n", 100 * mean(g2_diff > 0, na.rm = TRUE)))
    cat(sprintf("  %% hypomethylated:  %.1f%%\n", 100 * mean(g2_diff < 0, na.rm = TRUE)))
    cat("\n")

    # Statistical test
    wilcox_sub <- wilcox.test(g1_diff, g2_diff)

    cat("========================================\n")
    cat("STATISTICAL TEST (INPUT-SUBTRACTED)\n")
    cat("========================================\n\n")

    cat(sprintf("TES-bound vs Unbound Wilcoxon p-value: %.4e\n", wilcox_sub$p.value))

    if (mean(g1_diff, na.rm = TRUE) < mean(g2_diff, na.rm = TRUE)) {
        cat("Direction: TES-bound shows LESS methylation (same as non-subtracted)\n")
        input_sub_result <- "HYPOMETHYLATION_CONFIRMED"
    } else {
        cat("Direction: TES-bound shows MORE methylation (DIFFERENT from non-subtracted)\n")
        input_sub_result <- "PATTERN_CHANGED"
    }
    cat("\n")

    # Save statistics
    input_sub_stats <- data.frame(
        Group = c("TES_bound", "Unbound"),
        Mean_diff = c(mean(g1_diff, na.rm = TRUE), mean(g2_diff, na.rm = TRUE)),
        Median_diff = c(median(g1_diff, na.rm = TRUE), median(g2_diff, na.rm = TRUE)),
        Pct_hyper = c(100 * mean(g1_diff > 0, na.rm = TRUE), 100 * mean(g2_diff > 0, na.rm = TRUE)),
        Pct_hypo = c(100 * mean(g1_diff < 0, na.rm = TRUE), 100 * mean(g2_diff < 0, na.rm = TRUE)),
        N = c(length(g1_diff), length(g2_diff))
    )
    write.csv(input_sub_stats, file.path(OUTPUT_DIR, "input_subtracted_statistics.csv"), row.names = FALSE)

} else {
    cat("INPUT-subtracted matrix file not found - skipping Task 2\n")
    input_sub_result <- "NOT_TESTED"
}

# =============================================================================
# TASK 3: COMPARE MEDIPS IP-ONLY vs INPUT-NORMALIZED
# =============================================================================

cat("\n=== TASK 3: MEDIPS Comparison ===\n\n")

# Look for DMR files
ip_only_files <- list.files(MEDIPS_IP_ONLY, pattern = "DMRs.*\\.csv$|DMRs.*\\.txt$", full.names = TRUE)
input_norm_files <- list.files(MEDIPS_INPUT_NORM, pattern = "DMRs.*\\.csv$|DMRs.*\\.txt$", full.names = TRUE)

cat(sprintf("IP-only DMR files found: %d\n", length(ip_only_files)))
cat(sprintf("INPUT-normalized DMR files found: %d\n", length(input_norm_files)))
cat("\n")

if (length(ip_only_files) > 0 && length(input_norm_files) > 0) {

    # Read the first DMR file from each
    ip_only_dmrs <- tryCatch({
        f <- ip_only_files[grep("TES_vs_GFP", ip_only_files)][1]
        if (grepl("\\.csv$", f)) {
            read.csv(f, stringsAsFactors = FALSE)
        } else {
            read.delim(f, stringsAsFactors = FALSE)
        }
    }, error = function(e) NULL)

    input_norm_dmrs <- tryCatch({
        f <- input_norm_files[grep("TES_vs_GFP", input_norm_files)][1]
        if (grepl("\\.csv$", f)) {
            read.csv(f, stringsAsFactors = FALSE)
        } else {
            read.delim(f, stringsAsFactors = FALSE)
        }
    }, error = function(e) NULL)

    if (!is.null(ip_only_dmrs) && !is.null(input_norm_dmrs)) {
        cat("========================================\n")
        cat("MEDIPS DMR COMPARISON\n")
        cat("========================================\n\n")

        cat(sprintf("IP-only DMRs: %d\n", nrow(ip_only_dmrs)))
        cat(sprintf("INPUT-normalized DMRs: %d\n", nrow(input_norm_dmrs)))
        cat("\n")

        # Create genomic ranges for overlap
        if ("chr" %in% colnames(ip_only_dmrs) && "start" %in% colnames(ip_only_dmrs)) {
            ip_gr <- GRanges(
                seqnames = ip_only_dmrs$chr,
                ranges = IRanges(start = ip_only_dmrs$start, end = ip_only_dmrs$end)
            )
        } else if ("seqnames" %in% colnames(ip_only_dmrs)) {
            ip_gr <- GRanges(
                seqnames = ip_only_dmrs$seqnames,
                ranges = IRanges(start = ip_only_dmrs$start, end = ip_only_dmrs$end)
            )
        } else {
            ip_gr <- NULL
        }

        if ("chr" %in% colnames(input_norm_dmrs) && "start" %in% colnames(input_norm_dmrs)) {
            input_gr <- GRanges(
                seqnames = input_norm_dmrs$chr,
                ranges = IRanges(start = input_norm_dmrs$start, end = input_norm_dmrs$end)
            )
        } else if ("seqnames" %in% colnames(input_norm_dmrs)) {
            input_gr <- GRanges(
                seqnames = input_norm_dmrs$seqnames,
                ranges = IRanges(start = input_norm_dmrs$start, end = input_norm_dmrs$end)
            )
        } else {
            input_gr <- NULL
        }

        if (!is.null(ip_gr) && !is.null(input_gr)) {
            # Find overlaps
            overlaps <- findOverlaps(ip_gr, input_gr)
            n_overlap <- length(unique(queryHits(overlaps)))

            cat(sprintf("DMRs overlapping: %d (%.1f%% of IP-only, %.1f%% of INPUT-norm)\n",
                        n_overlap,
                        100 * n_overlap / length(ip_gr),
                        100 * n_overlap / length(input_gr)))
            cat("\n")

            # Check direction consistency
            if ("logFC" %in% colnames(ip_only_dmrs) && "logFC" %in% colnames(input_norm_dmrs)) {
                ip_hyper <- sum(ip_only_dmrs$logFC > 0, na.rm = TRUE)
                ip_hypo <- sum(ip_only_dmrs$logFC < 0, na.rm = TRUE)
                input_hyper <- sum(input_norm_dmrs$logFC > 0, na.rm = TRUE)
                input_hypo <- sum(input_norm_dmrs$logFC < 0, na.rm = TRUE)

                cat("Direction breakdown:\n")
                cat(sprintf("  IP-only: %d hyper, %d hypo (%.1f%% hyper)\n",
                            ip_hyper, ip_hypo, 100 * ip_hyper / (ip_hyper + ip_hypo)))
                cat(sprintf("  INPUT-norm: %d hyper, %d hypo (%.1f%% hyper)\n",
                            input_hyper, input_hypo, 100 * input_hyper / (input_hyper + input_hypo)))
                cat("\n")
            }

            # Create Venn diagram
            pdf(file.path(OUTPUT_DIR, "DMR_overlap_venn.pdf"), width = 6, height = 6)
            grid.newpage()
            venn <- draw.pairwise.venn(
                area1 = length(ip_gr),
                area2 = length(input_gr),
                cross.area = n_overlap,
                category = c("IP-only", "INPUT-norm"),
                fill = c("#E31A1C", "#1F78B4"),
                alpha = 0.5,
                cat.pos = c(-20, 20)
            )
            dev.off()
            cat("  Created: DMR_overlap_venn.pdf\n")

            medips_comparison <- "COMPLETED"
        } else {
            cat("Could not create genomic ranges from DMR files\n")
            medips_comparison <- "FAILED"
        }
    } else {
        cat("Could not read DMR files\n")
        medips_comparison <- "FAILED"
    }
} else {
    cat("DMR files not found in MEDIPS directories\n")
    cat(sprintf("  IP-only path: %s\n", MEDIPS_IP_ONLY))
    cat(sprintf("  INPUT-norm path: %s\n", MEDIPS_INPUT_NORM))
    medips_comparison <- "FILES_NOT_FOUND"
}

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("==========================================================\n")
cat("VERIFICATION SUMMARY\n")
cat("==========================================================\n\n")

cat("TASK 1 - INPUT Signal Check:\n")
if (exists("input_bias")) {
    cat(sprintf("  Result: %s\n", input_bias))
} else {
    cat("  Result: NOT_TESTED\n")
}
cat("\n")

cat("TASK 2 - INPUT-Subtracted Analysis:\n")
if (exists("input_sub_result")) {
    cat(sprintf("  Result: %s\n", input_sub_result))
} else {
    cat("  Result: NOT_TESTED\n")
}
cat("\n")

cat("TASK 3 - MEDIPS Comparison:\n")
if (exists("medips_comparison")) {
    cat(sprintf("  Result: %s\n", medips_comparison))
} else {
    cat("  Result: NOT_TESTED\n")
}
cat("\n")

cat("========================================\n")
cat("OVERALL CONCLUSION\n")
cat("========================================\n\n")

# Determine overall conclusion
if (exists("input_bias") && input_bias == "NO_BIAS" &&
    exists("input_sub_result") && input_sub_result == "HYPOMETHYLATION_CONFIRMED") {
    cat("The hypomethylation at TES-bound enhancers is REAL and NOT an artifact of\n")
    cat("INPUT normalization. The pattern persists after INPUT subtraction and INPUT\n")
    cat("signal is similar at TES-bound and unbound enhancers.\n")
} else if (exists("input_bias") && input_bias != "NO_BIAS") {
    cat("CAUTION: INPUT signal differs between TES-bound and unbound enhancers.\n")
    cat("This could contribute to the observed methylation differences.\n")
    cat("Consider using INPUT-normalized values for final conclusions.\n")
} else if (exists("input_sub_result") && input_sub_result == "PATTERN_CHANGED") {
    cat("WARNING: The methylation pattern CHANGES after INPUT subtraction.\n")
    cat("The original hypomethylation may be at least partially an artifact.\n")
} else {
    cat("Verification incomplete. Check individual task results above.\n")
}

cat("\n")
cat("==========================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==========================================================\n")
cat("Finished:", as.character(Sys.time()), "\n")
