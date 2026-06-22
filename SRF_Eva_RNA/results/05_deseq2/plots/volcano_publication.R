#!/usr/bin/env Rscript

#===============================================================================
# Publication-Ready Volcano Plot Script
# Creates two versions: labeled and clean (no labels)
#===============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(scales)
})

cat("=== Publication-Ready Volcano Plot Generation ===\n")
cat("Timestamp:", as.character(Sys.time()), "\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/05_deseq2")

#===============================================================================
# Configuration - Customize these parameters
#===============================================================================

# Thresholds
PADJ_CUTOFF <- 0.05        # Adjusted p-value cutoff for significance
FC_CUTOFF <- 1.0           # Log2 fold change cutoff (absolute value)

# Labeling (reduced for cleaner plot)
N_TOP_GENES <- 10          # Number of top genes to label (by padj)
N_TOP_UP <- 5              # Additional top upregulated genes to label
N_TOP_DOWN <- 5            # Additional top downregulated genes to label

# Plot aesthetics (PUBLICATION-READY - larger fonts)
POINT_SIZE <- 2.5          # Size of points
POINT_ALPHA <- 0.7         # Transparency of points
LABEL_SIZE <- 5.0          # Size of gene labels
AXIS_TITLE_SIZE <- 20      # Axis title font size
AXIS_TEXT_SIZE <- 16       # Axis text font size
TITLE_SIZE <- 22           # Main title font size
SUBTITLE_SIZE <- 16        # Subtitle font size
LEGEND_TITLE_SIZE <- 16    # Legend title font size
LEGEND_TEXT_SIZE <- 14     # Legend text font size

# Colors (colorblind-friendly)
COLOR_UP <- "#D73027"      # Upregulated (red)
COLOR_DOWN <- "#4575B4"    # Downregulated (blue)
COLOR_NS <- "#808080"      # Not significant (grey)

# Plot dimensions (larger for publication)
PLOT_WIDTH <- 10           # inches
PLOT_HEIGHT <- 9           # inches
DPI <- 300                 # Resolution for PNG

#===============================================================================
# Load and prepare data
#===============================================================================

cat("Loading DESeq2 results...\n")
res <- read.table("deseq2_results_TES_vs_GFP.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat(sprintf("Loaded %d genes\n", nrow(res)))

# Remove rows with NA padj values
res <- res[!is.na(res$padj), ]
cat(sprintf("After removing NA padj: %d genes\n", nrow(res)))

#===============================================================================
# Classify genes and prepare for plotting
#===============================================================================

# Classify genes based on significance and fold change
res$significance <- "Not Significant"
res$significance[res$padj < PADJ_CUTOFF & res$log2FoldChange > FC_CUTOFF] <- "Upregulated"
res$significance[res$padj < PADJ_CUTOFF & res$log2FoldChange < -FC_CUTOFF] <- "Downregulated"

# Convert to factor with specific order for legend
res$significance <- factor(res$significance,
                           levels = c("Upregulated", "Downregulated", "Not Significant"))

# Calculate -log10(padj) for y-axis
res$neg_log10_padj <- -log10(res$padj)

# Handle infinite values (padj = 0)
max_y <- max(res$neg_log10_padj[is.finite(res$neg_log10_padj)], na.rm = TRUE)
res$neg_log10_padj[!is.finite(res$neg_log10_padj)] <- max_y * 1.05

# Count genes in each category
n_up <- sum(res$significance == "Upregulated", na.rm = TRUE)
n_down <- sum(res$significance == "Downregulated", na.rm = TRUE)
n_ns <- sum(res$significance == "Not Significant", na.rm = TRUE)

cat(sprintf("\nGene classification:\n"))
cat(sprintf("  Upregulated: %d\n", n_up))
cat(sprintf("  Downregulated: %d\n", n_down))
cat(sprintf("  Not Significant: %d\n", n_ns))

#===============================================================================
# Select genes to label
#===============================================================================

# Top genes by adjusted p-value (most significant overall)
top_by_padj <- res %>%
    filter(significance != "Not Significant") %>%
    arrange(padj) %>%
    head(N_TOP_GENES)

# Top upregulated genes
top_up <- res %>%
    filter(significance == "Upregulated") %>%
    arrange(padj) %>%
    head(N_TOP_UP)

# Top downregulated genes
top_down <- res %>%
    filter(significance == "Downregulated") %>%
    arrange(padj) %>%
    head(N_TOP_DOWN)

# Combine and remove duplicates
genes_to_label <- unique(c(top_by_padj$gene_symbol,
                           top_up$gene_symbol,
                           top_down$gene_symbol))

# Create label column
res$label <- ifelse(res$gene_symbol %in% genes_to_label, res$gene_symbol, "")

cat(sprintf("\nGenes selected for labeling: %d\n", length(genes_to_label)))

#===============================================================================
# Define common theme for publication
#===============================================================================

theme_publication <- function(base_size = 14) {
    theme_classic(base_size = base_size) +
    theme(
        # Title
        plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5, color = "grey40"),

        # Axis
        axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
        axis.text = element_text(size = AXIS_TEXT_SIZE, color = "black"),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length = unit(0.25, "cm"),

        # Legend (larger and more prominent)
        legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.8, "cm"),
        legend.spacing.y = unit(0.2, "cm"),

        # Panel
        panel.grid = element_blank(),
        panel.border = element_blank(),

        # Margins
        plot.margin = margin(20, 20, 20, 20)
    )
}

#===============================================================================
# Create base volcano plot (shared elements)
#===============================================================================

create_volcano_base <- function(data) {
    # Calculate axis limits
    x_max <- max(abs(data$log2FoldChange), na.rm = TRUE) * 1.1
    y_max <- max(data$neg_log10_padj, na.rm = TRUE) * 1.05

    p <- ggplot(data, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
        # Points
        geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +

        # Threshold lines (thicker for visibility)
        geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
                   linetype = "dashed", color = "grey40", linewidth = 0.8) +
        geom_hline(yintercept = -log10(PADJ_CUTOFF),
                   linetype = "dashed", color = "grey40", linewidth = 0.8) +

        # Colors
        scale_color_manual(
            values = c("Upregulated" = COLOR_UP,
                      "Downregulated" = COLOR_DOWN,
                      "Not Significant" = COLOR_NS),
            labels = c(paste0("Up (", n_up, ")"),
                      paste0("Down (", n_down, ")"),
                      paste0("NS (", n_ns, ")"))
        ) +

        # Axis labels
        labs(
            x = expression(bold(log[2]~"Fold Change (TES / GFP)")),
            y = expression(bold(-log[10]~"(adjusted p-value)")),
            color = "Expression"
        ) +

        # Axis limits (symmetric x-axis)
        coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +

        # Theme
        theme_publication()

    return(p)
}

#===============================================================================
# Version 1: Volcano plot WITH labels
#===============================================================================

cat("\nCreating labeled volcano plot...\n")

p_labeled <- create_volcano_base(res) +
    # Gene labels with repel to avoid overlap (larger for publication)
    geom_text_repel(
        aes(label = label),
        size = LABEL_SIZE,
        fontface = "italic",
        max.overlaps = 35,
        box.padding = 0.6,
        point.padding = 0.4,
        segment.color = "grey40",
        segment.size = 0.4,
        segment.alpha = 0.7,
        min.segment.length = 0.3,
        force = 3,
        force_pull = 0.5,
        show.legend = FALSE
    ) +
    labs(
        title = "Differential Gene Expression: TES vs GFP",
        subtitle = sprintf("padj < %.2f, |log2FC| > %.1f", PADJ_CUTOFF, FC_CUTOFF)
    )

# Save labeled version
ggsave("plots/volcano_publication_labeled.pdf", p_labeled,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, device = cairo_pdf)
ggsave("plots/volcano_publication_labeled.png", p_labeled,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI)

cat("  Saved: volcano_publication_labeled.pdf\n")
cat("  Saved: volcano_publication_labeled.png\n")

#===============================================================================
# Version 2: Volcano plot WITHOUT labels (clean)
#===============================================================================

cat("\nCreating clean volcano plot (no labels)...\n")

p_clean <- create_volcano_base(res) +
    labs(
        title = "Differential Gene Expression: TES vs GFP",
        subtitle = sprintf("padj < %.2f, |log2FC| > %.1f", PADJ_CUTOFF, FC_CUTOFF)
    )

# Save clean version
ggsave("plots/volcano_publication_clean.pdf", p_clean,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, device = cairo_pdf)
ggsave("plots/volcano_publication_clean.png", p_clean,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI)

cat("  Saved: volcano_publication_clean.pdf\n")
cat("  Saved: volcano_publication_clean.png\n")

#===============================================================================
# Version 3: Minimal version for supplementary (small, clean)
#===============================================================================

cat("\nCreating minimal volcano plot...\n")

# Simplified for smaller figures (but still readable)
p_minimal <- create_volcano_base(res) +
    theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = LEGEND_TEXT_SIZE - 2),
        legend.title = element_text(size = LEGEND_TITLE_SIZE - 2),
        legend.box.margin = margin(10, 0, 0, 0),
        legend.margin = margin(0, 0, 0, 0),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.margin = margin(20, 20, 30, 20)  # Extra bottom margin for legend
    ) +
    guides(color = guide_legend(nrow = 1, title.position = "left"))

ggsave("plots/volcano_publication_minimal.pdf", p_minimal,
       width = 10, height = 8, device = cairo_pdf)
ggsave("plots/volcano_publication_minimal.png", p_minimal,
       width = 10, height = 8, dpi = DPI)

cat("  Saved: volcano_publication_minimal.pdf\n")
cat("  Saved: volcano_publication_minimal.png\n")

#===============================================================================
# Summary
#===============================================================================

cat("\n=== Volcano Plot Generation Complete ===\n")
cat(sprintf("Total genes plotted: %d\n", nrow(res)))
cat(sprintf("Upregulated genes: %d\n", n_up))
cat(sprintf("Downregulated genes: %d\n", n_down))
cat(sprintf("\nOutput files:\n"))
cat("  - volcano_publication_labeled.pdf/png (with gene labels)\n")
cat("  - volcano_publication_clean.pdf/png (no labels)\n")
cat("  - volcano_publication_minimal.pdf/png (compact version)\n")

