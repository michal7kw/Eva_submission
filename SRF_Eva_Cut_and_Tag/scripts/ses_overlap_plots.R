#!/usr/bin/env Rscript

# SES Overlap Analysis - Graphical Summary
# Generates comprehensive plots for TES/TEAD1 vs SES consensus overlap analysis

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(gridExtra)
library(readr)
library(tibble)

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Create output directory
dir.create("results/13_ses_overlap/plots", showWarnings = FALSE, recursive = TRUE)

# Function to parse stats files and extract overlap data
parse_stats_files <- function(stats_dir = "results/13_ses_overlap") {
  stats_files <- list.files(stats_dir, pattern = "*_ses_stats.txt", full.names = TRUE)

  overlap_data <- data.frame()

  for (file in stats_files) {
    lines <- readLines(file)

    # Extract information
    sample <- gsub("Sample: ", "", lines[grep("^Sample:", lines)])
    peak_type <- gsub("Peak Type: ", "", lines[grep("^Peak Type:", lines)])
    total_peaks <- as.numeric(gsub(".*: ", "", lines[grep("^Total.*peaks:", lines)]))
    overlapping_peaks <- as.numeric(gsub(".*: ", "", lines[grep(".*overlapping SES:", lines)]))
    overlap_percentage <- as.numeric(gsub(".*: ([0-9.]+)%", "\\1", lines[grep("Percentage.*overlapping SES:", lines)]))
    ses_covered <- as.numeric(gsub(".*: ", "", lines[grep("SES peaks overlapped by", lines)]))
    ses_coverage_percentage <- as.numeric(gsub(".*: ([0-9.]+)%", "\\1", lines[grep("Percentage of SES peaks covered", lines)]))

    # Clean sample name
    sample_clean <- gsub("_peaks\\.(narrowPeak|broadPeak)", "", sample)

    # Determine group and replicate (TESmut removed - failed sample)
    if (grepl("^TES-[123]$", sample_clean)) {
      group <- "TES"
      replicate <- gsub("TES-", "", sample_clean)
    } else if (grepl("^TEAD1-[123]$", sample_clean)) {
      group <- "TEAD1"
      replicate <- gsub("TEAD1-", "", sample_clean)
    } else if (sample_clean == "TES") {
      group <- "TES"
      replicate <- "combined"
    } else if (sample_clean == "TEAD1") {
      group <- "TEAD1"
      replicate <- "combined"
    } else if (grepl("TESmut", sample_clean)) {
      # Skip TESmut samples - failed sample
      next
    } else {
      group <- sample_clean
      replicate <- "unknown"
    }

    overlap_data <- rbind(overlap_data, data.frame(
      sample = sample_clean,
      group = group,
      replicate = replicate,
      peak_type = peak_type,
      total_peaks = total_peaks,
      overlapping_peaks = overlapping_peaks,
      overlap_percentage = overlap_percentage,
      ses_covered = ses_covered,
      ses_coverage_percentage = ses_coverage_percentage
    ))
  }

  return(overlap_data)
}

# Parse the overlap data
cat("Parsing overlap statistics...\n")
overlap_data <- parse_stats_files()

# Print summary
cat("Data summary:\n")
print(overlap_data)

# Define color palette (TESmut removed - failed sample)
group_colors <- c("TES" = "#2E86AB", "TEAD1" = "#F18F01")

# 1. Bar plot: Overlap percentages by sample and peak type
cat("Creating overlap percentage bar plots...\n")

# Filter individual replicates for clarity
individual_data <- overlap_data[overlap_data$replicate != "combined", ]
combined_data <- overlap_data[overlap_data$replicate == "combined", ]

p1 <- ggplot(individual_data, aes(x = sample, y = overlap_percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~peak_type, scales = "free_x") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Percentage of Peaks Overlapping with SES Consensus",
    subtitle = "Individual Replicates",
    x = "Sample", y = "Overlap Percentage (%)",
    fill = "Group"
  ) +
  ylim(0, 50)

# Combined samples plot
p2 <- ggplot(combined_data, aes(x = sample, y = overlap_percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~peak_type) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Percentage of Peaks Overlapping with SES Consensus",
    subtitle = "Combined Samples",
    x = "Sample", y = "Overlap Percentage (%)",
    fill = "Group"
  ) +
  ylim(0, 50)

# Save plots
pdf("results/13_ses_overlap/plots/overlap_percentage_barplots.pdf", width = 14, height = 10)
grid.arrange(p1, p2, nrow = 2)
dev.off()

# 2. SES Coverage plot
cat("Creating SES coverage plots...\n")

p3 <- ggplot(individual_data, aes(x = sample, y = ses_coverage_percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~peak_type, scales = "free_x") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Percentage of SES Consensus Peaks Covered",
    subtitle = "Individual Replicates",
    x = "Sample", y = "SES Coverage (%)",
    fill = "Group"
  ) +
  ylim(0, 25)

p4 <- ggplot(combined_data, aes(x = sample, y = ses_coverage_percentage, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~peak_type) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Percentage of SES Consensus Peaks Covered",
    subtitle = "Combined Samples",
    x = "Sample", y = "SES Coverage (%)",
    fill = "Group"
  ) +
  ylim(0, 25)

pdf("results/13_ses_overlap/plots/ses_coverage_barplots.pdf", width = 14, height = 10)
grid.arrange(p3, p4, nrow = 2)
dev.off()

# 3. Scatter plot: Overlap vs SES Coverage
cat("Creating scatter plot...\n")

p5 <- ggplot(overlap_data, aes(
  x = overlap_percentage, y = ses_coverage_percentage,
  color = group, shape = peak_type, size = total_peaks
)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_size_continuous(range = c(2, 8), name = "Total Peaks") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Peak Overlap vs SES Coverage",
    subtitle = "Bubble size represents total number of peaks",
    x = "Percentage of Peaks Overlapping SES (%)",
    y = "Percentage of SES Peaks Covered (%)",
    color = "Group", shape = "Peak Type"
  ) +
  xlim(0, 50) +
  ylim(0, 25)

ggsave("results/13_ses_overlap/plots/overlap_vs_coverage_scatter.pdf", p5, width = 10, height = 8)

# 4. Heatmap of overlap percentages
cat("Creating heatmap...\n")

# Prepare matrix for heatmap
heatmap_data <- overlap_data %>%
  select(sample, peak_type, overlap_percentage) %>%
  pivot_wider(
    names_from = peak_type, values_from = overlap_percentage,
    values_fn = mean
  ) %>%
  column_to_rownames("sample")

# Handle missing values
heatmap_data[is.na(heatmap_data)] <- 0

pdf("results/13_ses_overlap/plots/overlap_heatmap.pdf", width = 8, height = 10)
pheatmap(as.matrix(heatmap_data),
  color = colorRampPalette(c("white", "lightblue", "darkblue"))(50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  main = "Peak Overlap with SES Consensus (%)",
  cellwidth = 60,
  cellheight = 20,
  fontsize = 10,
  fontsize_number = 8
)
dev.off()

# 5. Summary statistics table
cat("Creating summary statistics...\n")

summary_stats <- overlap_data %>%
  group_by(group, peak_type) %>%
  summarise(
    n_samples = n(),
    mean_overlap = round(mean(overlap_percentage), 2),
    sd_overlap = round(sd(overlap_percentage), 2),
    mean_ses_coverage = round(mean(ses_coverage_percentage), 2),
    sd_ses_coverage = round(sd(ses_coverage_percentage), 2),
    .groups = "drop"
  )

write.csv(summary_stats, "results/13_ses_overlap/plots/summary_statistics.csv", row.names = FALSE)

# 6. Box plots comparing groups
cat("Creating box plots...\n")

p6 <- ggplot(individual_data, aes(x = group, y = overlap_percentage, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~peak_type) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Distribution of Overlap Percentages by Group",
    subtitle = "Individual replicates with jittered points",
    x = "Group", y = "Overlap Percentage (%)",
    fill = "Group"
  ) +
  ylim(0, 50)

p7 <- ggplot(individual_data, aes(x = group, y = ses_coverage_percentage, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~peak_type) +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "Distribution of SES Coverage by Group",
    subtitle = "Individual replicates with jittered points",
    x = "Group", y = "SES Coverage (%)",
    fill = "Group"
  ) +
  ylim(0, 25)

pdf("results/13_ses_overlap/plots/group_comparison_boxplots.pdf", width = 12, height = 8)
grid.arrange(p6, p7, nrow = 2)
dev.off()

# 7. Absolute numbers plot
cat("Creating absolute numbers plot...\n")

p8 <- ggplot(overlap_data, aes(x = sample, y = overlapping_peaks, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~peak_type, scales = "free") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Absolute Number of Overlapping Peaks",
    x = "Sample", y = "Number of Overlapping Peaks",
    fill = "Group"
  )

ggsave("results/13_ses_overlap/plots/absolute_overlap_numbers.pdf", p8, width = 14, height = 8)

# Create a comprehensive summary report
cat("Creating summary report...\n")

summary_text <- paste(
  "SES Consensus Overlap Analysis - Graphical Summary",
  "==================================================",
  "",
  paste("Analysis Date:", Sys.Date()),
  paste("Total SES Consensus Peaks:", 44726),
  "NOTE: TESmut samples excluded from analysis (failed sample)",
  "",
  "Key Findings:",
  "",
  "1. OVERLAP PERCENTAGES:",
  sprintf(
    "   - TES samples: %.1f%% ± %.1f%% (narrow)",
    mean(individual_data[individual_data$group == "TES" & individual_data$peak_type == "narrow", "overlap_percentage"], na.rm = TRUE),
    sd(individual_data[individual_data$group == "TES" & individual_data$peak_type == "narrow", "overlap_percentage"], na.rm = TRUE)
  ),
  sprintf(
    "   - TEAD1 samples: %.1f%% ± %.1f%% (narrow)",
    mean(individual_data[individual_data$group == "TEAD1" & individual_data$peak_type == "narrow", "overlap_percentage"], na.rm = TRUE),
    sd(individual_data[individual_data$group == "TEAD1" & individual_data$peak_type == "narrow", "overlap_percentage"], na.rm = TRUE)
  ),
  "",
  "2. SES COVERAGE:",
  sprintf(
    "   - Highest coverage by TEAD1 peaks: %.1f%%",
    max(overlap_data[overlap_data$group == "TEAD1", "ses_coverage_percentage"], na.rm = TRUE)
  ),
  "",
  "3. PEAK TYPE COMPARISON:",
  "   - Broad peaks generally show higher SES coverage than narrow peaks",
  "",
  "Generated plots:",
  "   - overlap_percentage_barplots.pdf: Bar plots of overlap percentages",
  "   - ses_coverage_barplots.pdf: Bar plots of SES coverage",
  "   - overlap_vs_coverage_scatter.pdf: Scatter plot of overlap vs coverage",
  "   - overlap_heatmap.pdf: Heatmap of overlap percentages",
  "   - group_comparison_boxplots.pdf: Box plots comparing groups",
  "   - absolute_overlap_numbers.pdf: Absolute numbers of overlapping peaks",
  "",
  sep = "\n"
)

writeLines(summary_text, "results/13_ses_overlap/plots/analysis_summary.txt")

cat("\nSES overlap graphical analysis complete!\n")
cat("Results saved in: results/13_ses_overlap/plots/\n")
cat("\nGenerated files:\n")
cat("  - overlap_percentage_barplots.pdf\n")
cat("  - ses_coverage_barplots.pdf\n")
cat("  - overlap_vs_coverage_scatter.pdf\n")
cat("  - overlap_heatmap.pdf\n")
cat("  - group_comparison_boxplots.pdf\n")
cat("  - absolute_overlap_numbers.pdf\n")
cat("  - summary_statistics.csv\n")
cat("  - analysis_summary.txt\n")
