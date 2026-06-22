#!/usr/bin/env Rscript

# GO Enrichment Analysis for TES vs GFP DEGs
# Separates upregulated and downregulated genes for pathway analysis

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(readr)
library(dplyr)
library(stringr)

#===============================================================================
# PLOT CONFIGURATION - Adjust these parameters to change text sizes
#===============================================================================

# Text size scaling factor (increase this to make all text bigger)
# Default is 1.0, use 1.5 for 50% larger text, 2.0 for double size, etc.
TEXT_SCALE <- 1.5

# Output settings
DPI <- 300

# Base font sizes (will be multiplied by TEXT_SCALE)
BASE_SIZE <- 14 * TEXT_SCALE        # Base font size for theme
TITLE_SIZE <- 18 * TEXT_SCALE       # Plot title size
SUBTITLE_SIZE <- 14 * TEXT_SCALE    # Plot subtitle size
AXIS_TITLE_SIZE <- 16 * TEXT_SCALE  # Axis title size
AXIS_TEXT_SIZE <- 14 * TEXT_SCALE   # Axis text size
AXIS_TEXT_Y_SIZE <- 13 * TEXT_SCALE # Y-axis text (GO terms)
LEGEND_TITLE_SIZE <- 14 * TEXT_SCALE # Legend title size
LEGEND_TEXT_SIZE <- 12 * TEXT_SCALE  # Legend text size
LABEL_SIZE <- 5 * TEXT_SCALE         # Bar labels size
DOTPLOT_FONT_SIZE <- 13 * TEXT_SCALE # Font size for dotplot/barplot

#===============================================================================
# SET UP DIRECTORIES
#===============================================================================

# Set up directories
input_file <- "results/05_deseq2/significant_genes_TES_vs_GFP.txt"
all_genes_file <- "results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
output_dir <- "results/06_go_enrichment"
upregulated_dir <- file.path(output_dir, "upregulated")
downregulated_dir <- file.path(output_dir, "downregulated")
plots_dir <- file.path(output_dir, "plots")

# Create output directories if they don't exist
dir.create(upregulated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(downregulated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

cat("Reading DEGs data...\n")
degs <- read_tsv(input_file, show_col_types = FALSE)

# Load ALL tested genes for proper background universe
cat("Reading all tested genes for background universe...\n")
all_genes <- read_tsv(all_genes_file, show_col_types = FALSE)
cat("Total genes in DESeq2 results:", nrow(all_genes), "\n")

# Convert ENSG IDs to gene symbols and Entrez IDs
cat("Converting gene IDs...\n")
# Strip version numbers from Ensembl IDs (e.g., ENSG00000060718.22 -> ENSG00000060718)
degs <- degs %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\.\\d+$"))

degs <- degs %>%
  mutate(
    gene_symbol = mapIds(org.Hs.eg.db,
                         keys = gene_id_clean,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first"),
    entrez_id = mapIds(org.Hs.eg.db,
                       keys = gene_id_clean,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  )
# Remove genes with missing mappings
degs_filtered <- degs %>%
  filter(!is.na(gene_symbol), !is.na(entrez_id))

cat("Total significant genes after ID mapping:", nrow(degs_filtered), "\n")

# Build background universe from ALL tested genes (not just significant ones)
all_genes <- all_genes %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\.\\d+$"))
all_genes <- all_genes %>%
  mutate(
    entrez_id = mapIds(org.Hs.eg.db,
                       keys = gene_id_clean,
                       column = "ENTREZID",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  )
background_genes <- unique(na.omit(all_genes$entrez_id))
cat("Background universe size (all tested genes with Entrez IDs):", length(background_genes), "\n")

# Separate upregulated and downregulated genes
# Upregulated: log2FoldChange > 0 (higher in TES)
# Downregulated: log2FoldChange < 0 (lower in TES)
upregulated_genes <- degs_filtered %>%
  filter(log2FoldChange > 0) %>%
  pull(entrez_id)

downregulated_genes <- degs_filtered %>%
  filter(log2FoldChange < 0) %>%
  pull(entrez_id)

cat("Upregulated genes in TES:", length(upregulated_genes), "\n")
cat("Downregulated genes in TES:", length(downregulated_genes), "\n")

# Function to perform GO enrichment analysis
perform_go_enrichment <- function(gene_list, background, direction, ontology = "BP") {
  cat(paste("Performing GO", ontology, "enrichment for", direction, "genes...\n"))
  
  ego <- enrichGO(
    gene = gene_list,
    universe = background,
    OrgDb = org.Hs.eg.db,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  return(ego)
}

# Function to save results and plots
save_results <- function(ego, direction, ontology) {
  # Save results table
  results_file <- file.path(output_dir, tolower(direction), 
                           paste0("GO_", ontology, "_", direction, "_results.csv"))
  write.csv(as.data.frame(ego), results_file, row.names = FALSE)
  
  # Create summary
  summary_file <- file.path(output_dir, tolower(direction),
                           paste0("GO_", ontology, "_", direction, "_summary.txt"))
  sink(summary_file)
  cat("GO", ontology, "Enrichment Summary -", direction, "Genes\n")
  cat("=" , rep("=", 50), "\n", sep = "")
  
  # Access gene universe size from the enrichResult object properly
  ego_df <- as.data.frame(ego)
  cat("Significant terms found:", nrow(ego_df), "\n")
  cat("P-value cutoff: 0.05\n")
  cat("Q-value cutoff: 0.2\n\n")
  
  if (nrow(ego_df) > 0) {
    cat("Top 10 enriched terms:\n")
    print(head(ego_df[order(ego_df$pvalue), ], 10))
  }
  sink()
  
  cat("Results saved to:", results_file, "\n")
  cat("Summary saved to:", summary_file, "\n")
}

# Publication-ready theme for all plots (uses TEXT_SCALE parameters)
pub_theme <- theme_bw(base_size = BASE_SIZE) +
  theme(
    # Title styling
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5, margin = margin(b = 15)),
    plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5, margin = margin(b = 10)),
    # Axis styling
    axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT_SIZE, color = "black"),
    axis.text.y = element_text(size = AXIS_TEXT_Y_SIZE),
    axis.line = element_line(color = "black", linewidth = 0.8),
    # Legend styling
    legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
    legend.text = element_text(size = LEGEND_TEXT_SIZE),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.8 * TEXT_SCALE, "cm"),
    # Panel styling
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # Margin
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt")
  )

# Function to create plots
create_plots <- function(ego, direction, ontology) {
  ego_df <- as.data.frame(ego)
  if (nrow(ego_df) == 0) {
    cat("No significant terms found for plotting\n")
    return()
  }

  # Map ontology codes to full names
  ont_names <- c("BP" = "Biological Process", "MF" = "Molecular Function", "CC" = "Cellular Component")
  ont_full <- ont_names[ontology]

  # Direction labels for titles
  dir_label <- ifelse(direction == "Upregulated", "Upregulated in TES", "Downregulated in TES")

  # Dot plot - Publication ready
  p1 <- dotplot(ego, showCategory = 20,
                title = paste0("GO ", ont_full, "\n", dir_label),
                font.size = DOTPLOT_FONT_SIZE) +
    pub_theme +
    theme(
      axis.text.y = element_text(size = AXIS_TEXT_Y_SIZE, lineheight = 0.9),
      plot.title = element_text(size = TITLE_SIZE * 0.9, face = "bold", hjust = 0.5, lineheight = 1.2)
    ) +
    scale_color_gradient(low = "#E64B35", high = "#4DBBD5",
                         name = "Adjusted\np-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3 * TEXT_SCALE, 10 * TEXT_SCALE))

  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_dotplot.pdf")),
         plot = p1, width = 12, height = 10, device = cairo_pdf)

  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_dotplot.png")),
         plot = p1, width = 12, height = 10, dpi = DPI)

  # Enrichment map - Publication ready
  if (nrow(ego_df) >= 5) {
    # Pairwise term similarity for emapplot
    ego_sim <- pairwise_termsim(ego)

    p2 <- emapplot(ego_sim, showCategory = 30,
                   cex_label_category = 1.2 * TEXT_SCALE,
                   cex_category = 1.5 * TEXT_SCALE) +
      ggtitle(paste0("GO ", ont_full, " Network\n", dir_label)) +
      pub_theme +
      theme(
        plot.title = element_text(size = TITLE_SIZE * 0.9, face = "bold", hjust = 0.5, lineheight = 1.2),
        legend.position = "right"
      )

    ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_emap.pdf")),
           plot = p2, width = 12, height = 10, device = cairo_pdf)

    ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_emap.png")),
           plot = p2, width = 12, height = 10, dpi = DPI)
  }

  # Bar plot - Publication ready
  p3 <- barplot(ego, showCategory = 20,
                title = paste0("GO ", ont_full, "\n", dir_label),
                font.size = DOTPLOT_FONT_SIZE) +
    pub_theme +
    theme(
      axis.text.y = element_text(size = AXIS_TEXT_Y_SIZE, lineheight = 0.9),
      plot.title = element_text(size = TITLE_SIZE * 0.9, face = "bold", hjust = 0.5, lineheight = 1.2)
    ) +
    scale_fill_gradient(low = "#4DBBD5", high = "#E64B35",
                        name = "Adjusted\np-value")

  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_barplot.pdf")),
         plot = p3, width = 12, height = 10, device = cairo_pdf)

  ggsave(filename = file.path(plots_dir, paste0("GO_", ontology, "_", direction, "_barplot.png")),
         plot = p3, width = 12, height = 10, dpi = DPI)
  cat("Plots saved for", direction, ontology, "\n")
}

# Perform GO enrichment analysis for different ontologies
ontologies <- c("BP", "MF", "CC")  # Biological Process, Molecular Function, Cellular Component

for (direction in c("Upregulated", "Downregulated")) {
  if (direction == "Upregulated") {
    gene_list <- upregulated_genes
  } else {
    gene_list <- downregulated_genes
  }
  
  cat("\n=== Analyzing", direction, "Genes ===\n")
  
  for (ont in ontologies) {
    ego <- perform_go_enrichment(gene_list, background_genes, direction, ont)
    save_results(ego, direction, ont)
    create_plots(ego, direction, ont)
  }
}

# Create comparative summary
cat("\nCreating comparative summary...\n")

summary_data <- data.frame(
  Direction = character(),
  Ontology = character(),
  Significant_Terms = integer(),
  Genes_in_Terms = integer(),
  stringsAsFactors = FALSE
)

for (direction in c("Upregulated", "Downregulated")) {
  for (ont in ontologies) {
    results_file <- file.path(output_dir, tolower(direction),
                             paste0("GO_", ont, "_", direction, "_results.csv"))
    if (file.exists(results_file)) {
      results <- read.csv(results_file)
      if (nrow(results) > 0) {
        summary_data <- rbind(summary_data, data.frame(
          Direction = direction,
          Ontology = ont,
          Significant_Terms = nrow(results),
          Genes_in_Terms = sum(str_count(results$geneID, "/") + 1),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Save comparative summary
summary_file <- file.path(output_dir, "GO_enrichment_comparative_summary.csv")
write.csv(summary_data, summary_file, row.names = FALSE)

# Create comparative plot - Publication ready
if (nrow(summary_data) > 0) {
  # Map ontology codes to full names for x-axis
  summary_data <- summary_data %>%
    mutate(Ontology_Full = case_when(
      Ontology == "BP" ~ "Biological\nProcess",
      Ontology == "MF" ~ "Molecular\nFunction",
      Ontology == "CC" ~ "Cellular\nComponent",
      TRUE ~ Ontology
    ))

  p_comparison <- ggplot(summary_data, aes(x = Ontology_Full, y = Significant_Terms, fill = Direction)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.5) +
    geom_text(aes(label = Significant_Terms),
              position = position_dodge(width = 0.8),
              vjust = -0.5, size = LABEL_SIZE, fontface = "bold") +
    scale_fill_manual(values = c("Upregulated" = "#E64B35", "Downregulated" = "#4DBBD5"),
                      labels = c("Upregulated" = "Upregulated in TES",
                                 "Downregulated" = "Downregulated in TES")) +
    labs(title = "GO Enrichment Comparison",
         subtitle = "Number of Significant Terms by Direction",
         x = "GO Ontology",
         y = "Number of Significant Terms",
         fill = "Gene Set") +
    pub_theme +
    theme(
      plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = SUBTITLE_SIZE, hjust = 0.5, color = "grey40"),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE, lineheight = 0.9),
      legend.position = "top",
      legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
      legend.text = element_text(size = LEGEND_TEXT_SIZE)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

  ggsave(filename = file.path(plots_dir, "GO_enrichment_comparison.pdf"),
         plot = p_comparison, width = 10, height = 7, device = cairo_pdf)

  ggsave(filename = file.path(plots_dir, "GO_enrichment_comparison.png"),
         plot = p_comparison, width = 10, height = 7, dpi = DPI)
}

# Create gene lists for reference
upregulated_df <- degs_filtered %>%
  filter(log2FoldChange > 0) %>%
  select(gene_id, gene_symbol, entrez_id, log2FoldChange, padj) %>%
  arrange(padj)

downregulated_df <- degs_filtered %>%
  filter(log2FoldChange < 0) %>%
  select(gene_id, gene_symbol, entrez_id, log2FoldChange, padj) %>%
  arrange(padj)

write.csv(upregulated_df, file.path(upregulated_dir, "upregulated_genes_list.csv"), row.names = FALSE)
write.csv(downregulated_df, file.path(downregulated_dir, "downregulated_genes_list.csv"), row.names = FALSE)

cat("\n=== GO Enrichment Analysis Complete ===\n")
cat("Results saved in:", output_dir, "\n")
cat("Upregulated results:", upregulated_dir, "\n")
cat("Downregulated results:", downregulated_dir, "\n")
cat("Plots saved in:", plots_dir, "\n")
cat("Comparative summary:", summary_file, "\n")

# Print final statistics
cat("\nFinal Statistics:\n")
cat("- Total significant genes analyzed:", nrow(degs_filtered), "\n")
cat("- Upregulated genes:", length(upregulated_genes), "\n")
cat("- Downregulated genes:", length(downregulated_genes), "\n")

if (file.exists(summary_file)) {
  summary_stats <- read.csv(summary_file)
  cat("\nSignificant GO terms found:\n")
  print(summary_stats)
}

cat("\nAnalysis completed successfully!\n")