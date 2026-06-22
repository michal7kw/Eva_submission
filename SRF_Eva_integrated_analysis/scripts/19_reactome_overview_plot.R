#!/usr/bin/env Rscript
# =============================================================================
# Reactome Pathway GSEA Plots
# =============================================================================
# Creates:
# 1. GSEA overview bar chart (plot_05_gsea_overview)
# 2. Enhanced volcano plot matching 16_enhanced_plots style (plot_06_volcano_enhanced)
#
# Usage: Rscript 19_reactome_overview_plot.R
# =============================================================================

# -----------------------------------------------------------------------------
# Load Libraries
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

cat("===========================================\n")
cat("REACTOME PATHWAY GSEA PLOTS\n")
cat("===========================================\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
INPUT_FILE <- "output/12_msigdb_by_collection/06_Reactome_Pathways/all_pathways.csv"
OUTPUT_DIR <- "output/12_msigdb_by_collection/06_Reactome_Pathways"

# Check input exists
if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE, "\n",
       "Run 12_msigdb_gsea_by_collection.R first to generate Reactome GSEA results.")
}

# -----------------------------------------------------------------------------
# Publication Theme (matching 16_enhanced_gsea_visualizations.R)
# -----------------------------------------------------------------------------
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(0.9), hjust = 0.5),
      panel.border = element_rect(colour = "black", linewidth = 1),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(size = rel(0.9)),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(colour = "black", fill = "#f7f7f7"),
      strip.text = element_text(face = "bold")
    )
}

# -----------------------------------------------------------------------------
# Read Data
# -----------------------------------------------------------------------------
cat("Reading Reactome GSEA results from:", INPUT_FILE, "\n")
reactome_results <- read.csv(INPUT_FILE)

cat("  Total pathways loaded:", nrow(reactome_results), "\n")

# -----------------------------------------------------------------------------
# Calculate Summary Statistics
# -----------------------------------------------------------------------------
cat("\nCalculating summary statistics...\n")

total_tested <- nrow(reactome_results)
significant <- sum(reactome_results$padj < 0.05, na.rm = TRUE)
downregulated <- sum(reactome_results$padj < 0.05 & reactome_results$NES < 0, na.rm = TRUE)
upregulated <- sum(reactome_results$padj < 0.05 & reactome_results$NES > 0, na.rm = TRUE)

summary_data <- data.frame(
  Category = c("Total Tested", "Significant (FDR<0.05)",
               "Downregulated (NES<0)", "Upregulated (NES>0)"),
  Count = c(total_tested, significant, downregulated, upregulated)
)

# Print summary
cat("\nSummary Statistics:\n")
cat("  Total Tested:", total_tested, "\n")
cat("  Significant (FDR<0.05):", significant, "\n")
cat("  Downregulated (NES<0):", downregulated, "\n")
cat("  Upregulated (NES>0):", upregulated, "\n")

# =============================================================================
# PLOT 1: GSEA OVERVIEW BAR CHART
# =============================================================================
cat("\n1. Creating GSEA overview bar chart...\n")

p_overview <- ggplot(summary_data, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 5, fontface = "bold") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Reactome Pathway GSEA Overview",
    subtitle = "MSigDB C2:CP:REACTOME Collection | TES vs GFP | FDR < 0.05",
    x = NULL,
    y = "Number of Pathways"
  ) +
  theme_publication() +
  theme(legend.position = "none") +
  ylim(0, max(summary_data$Count) * 1.15)

# Save overview plot
ggsave(file.path(OUTPUT_DIR, "plot_05_gsea_overview.pdf"), p_overview,
       width = 10, height = 5, device = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "plot_05_gsea_overview.png"), p_overview,
       width = 10, height = 5, dpi = 300)
cat("  Saved: plot_05_gsea_overview.pdf/png\n")

# =============================================================================
# PLOT 2: ENHANCED VOLCANO PLOT (matching 16_enhanced_plots style)
# =============================================================================
cat("\n2. Creating enhanced volcano plot...\n")

# Prepare data for volcano plot
volcano_data <- reactome_results %>%
  mutate(
    neg_log10_padj = -log10(padj),
    # Handle infinite values from very small p-values
    neg_log10_padj = ifelse(is.infinite(neg_log10_padj), max(neg_log10_padj[!is.infinite(neg_log10_padj)]) * 1.1, neg_log10_padj)
  )

# Select top 10 pathways to label (by adjusted p-value) - reduced from 20 to avoid overlap
top_pathways <- volcano_data %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(10)

# Create label column - shorten pathway names for readability
volcano_data$label <- ifelse(
  volcano_data$pathway %in% top_pathways$pathway,
  gsub("^REACTOME_", "", volcano_data$pathway),  # Remove REACTOME_ prefix
  ""
)

# Further shorten long labels
volcano_data$label <- ifelse(
  nchar(volcano_data$label) > 40,
  paste0(substr(volcano_data$label, 1, 37), "..."),
  volcano_data$label
)

# Create the enhanced volcano plot
p_volcano <- ggplot(volcano_data, aes(x = NES, y = neg_log10_padj)) +
  # All points with NES-based color gradient and size by gene set size
  geom_point(aes(color = NES, size = size), alpha = 0.6) +
  # Highlight labeled points with ring
  geom_point(
    data = subset(volcano_data, label != ""),
    aes(color = NES, size = size),
    alpha = 0.9, shape = 21, stroke = 1.5, fill = NA
  ) +
  # Add labels with ggrepel - adjusted parameters for better separation
  geom_text_repel(
    aes(label = label),
    size = 2.8,
    max.overlaps = 15,
    box.padding = 0.8,
    point.padding = 0.5,
    force = 3,
    force_pull = 0.5,
    segment.color = "grey50",
    segment.size = 0.3,
    min.segment.length = 0.2,
    max.time = 2,
    max.iter = 20000
  ) +
  # Threshold lines
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  # Color scale: blue (negative) to red (positive)
  scale_color_gradient2(
    low = "#2166AC", mid = "grey90", high = "#B2182B",
    midpoint = 0, name = "NES"
  ) +
  # Size scale
  scale_size_continuous(range = c(1, 8), name = "Gene Set Size") +
  # Labels
  labs(
    title = "Reactome Pathway Enrichment Landscape",
    subtitle = "MSigDB C2:CP:REACTOME | TES vs GFP | Top 10 labeled | Thresholds: |NES| = 1.5, FDR = 0.01",
    x = "Normalized Enrichment Score (NES)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_publication()

# Save enhanced volcano plot
ggsave(file.path(OUTPUT_DIR, "plot_06_volcano_enhanced.pdf"), p_volcano,
       width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(OUTPUT_DIR, "plot_06_volcano_enhanced.png"), p_volcano,
       width = 14, height = 10, dpi = 300)
cat("  Saved: plot_06_volcano_enhanced.pdf/png\n")

# -----------------------------------------------------------------------------
# Completion
# -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("REACTOME PLOTS COMPLETE\n")
cat("===========================================\n")
cat("\nGenerated files:\n")
cat("  - plot_05_gsea_overview.pdf/png - GSEA summary bar chart\n")
cat("  - plot_06_volcano_enhanced.pdf/png - Enhanced volcano plot\n")
cat("\nStatistics:\n")
cat(sprintf("  Total Tested: %d\n", total_tested))
cat(sprintf("  Significant: %d\n", significant))
cat(sprintf("  Downregulated: %d\n", downregulated))
cat(sprintf("  Upregulated: %d\n", upregulated))
