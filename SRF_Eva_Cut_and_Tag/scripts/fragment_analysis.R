#!/usr/bin/env Rscript

# ===============================================================================
# SCRIPT: fragment_analysis.R
# PURPOSE: Fragment size analysis for Cut&Tag quality assessment
#
# DESCRIPTION:
# This script performs fragment size distribution analysis which is crucial
# for Cut&Tag data quality assessment. It analyzes insert size distributions
# to ensure proper tagmentation and library preparation.
# ===============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Create output directory
dir.create("results/fragment_analysis", recursive = TRUE, showWarnings = FALSE)

# Function to extract fragment sizes from BAM file
extract_fragment_sizes <- function(bam_file, sample_name) {
    cat("Analyzing fragment sizes for:", sample_name, "\n")

    # Use samtools to extract fragment sizes
    cmd <- paste("samtools view -f 2", bam_file, "| awk '$9 > 0 {print $9}' | head -100000")
    fragment_sizes <- as.numeric(system(cmd, intern = TRUE))

    # Filter reasonable fragment sizes (10-1000bp for Cut&Tag)
    fragment_sizes <- fragment_sizes[fragment_sizes >= 10 & fragment_sizes <= 1000]

    return(data.frame(
        sample = sample_name,
        fragment_size = fragment_sizes,
        stringsAsFactors = FALSE
    ))
}

# Sample list (TESmut removed - failed sample)
samples <- c(
    "TES-1", "TES-2", "TES-3",
    "TEAD1-1", "TEAD1-2", "TEAD1-3", "IggMs", "IggRb"
)

# Extract fragment sizes for all samples
fragment_data_list <- list()
bam_dir <- "results/04_filtered"

for (sample in samples) {
    bam_file <- file.path(bam_dir, paste0(sample, "_filtered.bam"))
    if (file.exists(bam_file)) {
        fragment_data_list[[sample]] <- extract_fragment_sizes(bam_file, sample)
    } else {
        cat("Warning: BAM file not found for", sample, "\n")
    }
}

# Combine all fragment data
fragment_data <- do.call(rbind, fragment_data_list)

# Add experimental groups
fragment_data$group <- case_when(
    grepl("^TES-[0-9]", fragment_data$sample) ~ "TES",
    grepl("^TESmut-[0-9]", fragment_data$sample) ~ "TESmut",
    grepl("^TEAD1-[0-9]", fragment_data$sample) ~ "TEAD1",
    grepl("^Igg", fragment_data$sample) ~ "Control",
    TRUE ~ "Other"
)

# Calculate summary statistics
fragment_stats <- fragment_data %>%
    group_by(sample, group) %>%
    summarise(
        median_size = median(fragment_size, na.rm = TRUE),
        mean_size = mean(fragment_size, na.rm = TRUE),
        q25 = quantile(fragment_size, 0.25, na.rm = TRUE),
        q75 = quantile(fragment_size, 0.75, na.rm = TRUE),
        nucleosome_peak = fragment_size[which.max(table(cut(fragment_size, breaks = seq(140, 180, by = 2))))],
        .groups = "drop"
    )

# Save statistics
write.csv(fragment_stats, "results/fragment_analysis/fragment_statistics.csv", row.names = FALSE)

# Create comprehensive plots
p1 <- ggplot(fragment_data, aes(x = fragment_size, fill = group)) +
    geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
    facet_wrap(~sample, scales = "free_y") +
    labs(
        title = "Fragment Size Distribution by Sample",
        x = "Fragment Size (bp)", y = "Count"
    ) +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set1")

# Overlay plot comparing groups
p2 <- ggplot(fragment_data, aes(x = fragment_size, color = group)) +
    geom_density(size = 1.2, alpha = 0.8) +
    labs(
        title = "Fragment Size Distribution by Group",
        x = "Fragment Size (bp)", y = "Density"
    ) +
    theme_minimal() +
    scale_color_brewer(type = "qual", palette = "Set1") +
    xlim(0, 500)

# Nucleosome pattern analysis (focused on 100-300bp range)
p3 <- ggplot(
    filter(fragment_data, fragment_size >= 100 & fragment_size <= 300),
    aes(x = fragment_size, color = group)
) +
    geom_density(size = 1.2, alpha = 0.8) +
    labs(
        title = "Nucleosome Pattern Analysis (100-300bp)",
        x = "Fragment Size (bp)", y = "Density"
    ) +
    theme_minimal() +
    scale_color_brewer(type = "qual", palette = "Set1") +
    geom_vline(xintercept = c(147, 294), linetype = "dashed", alpha = 0.5)

# Box plot showing median fragment sizes
p4 <- ggplot(fragment_data, aes(x = sample, y = fragment_size, fill = group)) +
    geom_boxplot(outlier.size = 0.1) +
    labs(
        title = "Fragment Size Distribution by Sample",
        x = "Sample", y = "Fragment Size (bp)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ylim(0, 500)

# Save plots
ggsave("results/fragment_analysis/fragment_size_by_sample.pdf", p1, width = 16, height = 12)
ggsave("results/fragment_analysis/fragment_size_by_group.pdf", p2, width = 10, height = 6)
ggsave("results/fragment_analysis/nucleosome_pattern.pdf", p3, width = 10, height = 6)
ggsave("results/fragment_analysis/fragment_size_boxplot.pdf", p4, width = 12, height = 8)

# PNG versions
ggsave("results/fragment_analysis/fragment_size_by_sample.png", p1, width = 16, height = 12, dpi = 300)
ggsave("results/fragment_analysis/fragment_size_by_group.png", p2, width = 10, height = 6, dpi = 300)
ggsave("results/fragment_analysis/nucleosome_pattern.png", p3, width = 10, height = 6, dpi = 300)
ggsave("results/fragment_analysis/fragment_size_boxplot.png", p4, width = 12, height = 8, dpi = 300)

# Quality assessment report
cat("\n=== CUT&TAG FRAGMENT SIZE ANALYSIS REPORT ===\n")
cat("Expected patterns for good Cut&Tag data:\n")
cat("1. Peak around 147bp (mononucleosome)\n")
cat("2. Smaller peak around 294bp (dinucleosome)\n")
cat("3. Minimal fragments < 100bp (over-tagmentation)\n")
cat("4. Minimal fragments > 500bp (under-tagmentation)\n\n")

cat("Sample-wise median fragment sizes:\n")
print(fragment_stats[, c("sample", "group", "median_size", "mean_size")])

# Flag potential issues
cat("\nQuality assessment:\n")
for (i in 1:nrow(fragment_stats)) {
    sample_name <- fragment_stats$sample[i]
    median_size <- fragment_stats$median_size[i]

    if (median_size < 120) {
        cat("WARNING:", sample_name, "- Small median fragment size (", median_size, "bp) suggests over-tagmentation\n")
    } else if (median_size > 200) {
        cat("WARNING:", sample_name, "- Large median fragment size (", median_size, "bp) suggests under-tagmentation\n")
    } else {
        cat("GOOD:", sample_name, "- Fragment size distribution looks normal (", median_size, "bp)\n")
    }
}

cat("\nFragment analysis complete! Check results/fragment_analysis/ for plots and statistics.\n")
