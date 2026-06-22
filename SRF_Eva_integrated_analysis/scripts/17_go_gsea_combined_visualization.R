#!/usr/bin/env Rscript
#
# GO AND GSEA COMBINED VISUALIZATION (SIMPLIFIED)
# Publication-quality comparison plots for GO ORA and GSEA results
#

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== GO AND GSEA COMBINED VISUALIZATION ===\n\n")

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

GSEA_INPUT <- "output/01_true_gsea_analysis"
GO_INPUT <- "output/02_directional_go_enrichment"
OUTPUT_BASE <- file.path("output", "17_go_gsea_combined")
dir.create(OUTPUT_BASE, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PUBLICATION-QUALITY THEME
# =============================================================================

theme_publication <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(size = rel(0.9)),
      legend.title = element_text(face = "bold", size = rel(1.0)),
      legend.text = element_text(size = rel(0.9)),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold", size = rel(1.0))
    )
}

# Color palette
colors_direction <- c("Upregulated" = "#D6604D", "Downregulated" = "#4393C3")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

# Load GSEA
gsea_sig <- tryCatch(
  read.csv(file.path(GSEA_INPUT, "fgsea_significant_pathways.csv")),
  error = function(e) {
    cat("  ERROR: Could not load GSEA results\n")
    return(NULL)
  }
)

# Load GO
go_up <- tryCatch(
  read.csv(file.path(GO_INPUT, "upregulated_GO_BP.csv")),
  error = function(e) NULL
)
go_down <- tryCatch(
  read.csv(file.path(GO_INPUT, "downregulated_GO_BP.csv")),
  error = function(e) NULL
)

# Check data
if (is.null(gsea_sig) || nrow(gsea_sig) == 0) {
  cat("ERROR: No GSEA results found. Exiting.\n")
  quit(status = 1)
}

if ((is.null(go_up) || nrow(go_up) == 0) && (is.null(go_down) || nrow(go_down) == 0)) {
  cat("ERROR: No GO results found. Exiting.\n")
  quit(status = 1)
}

# Process GO
if (!is.null(go_up) && nrow(go_up) > 0) go_up$direction <- "Upregulated"
if (!is.null(go_down) && nrow(go_down) > 0) go_down$direction <- "Downregulated"

go_combined <- bind_rows(go_up, go_down) %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    term_short = str_trunc(Description, 50),
    log10_padj = -log10(p.adjust)
  )

# Process GSEA
gsea_sig <- gsea_sig %>%
  mutate(
    direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
    term_short = str_trunc(pathway, 50),
    log10_padj = -log10(padj)
  )

cat(sprintf("  GO terms: %d\n", nrow(go_combined)))
cat(sprintf("  GSEA pathways: %d\n", nrow(gsea_sig)))

# Prepare plot data
go_plot <- go_combined %>%
  select(
    term = Description, term_short,
    effect_size = FoldEnrichment,
    padj = p.adjust, log10_padj,
    direction, count = Count
  ) %>%
  mutate(
    method = "GO ORA",
    effect_size_directed = ifelse(direction == "Downregulated", -effect_size, effect_size)
  )

gsea_plot <- gsea_sig %>%
  select(
    term = pathway, term_short,
    effect_size = NES,
    padj, log10_padj,
    direction, count = size
  ) %>%
  mutate(
    method = "GSEA",
    effect_size_directed = effect_size
  )

# Labels for top terms
top_go <- go_plot %>% arrange(padj) %>% head(10) %>% pull(term)
top_gsea <- gsea_plot %>% arrange(padj) %>% head(10) %>% pull(term)
go_plot$label <- ifelse(go_plot$term %in% top_go, go_plot$term_short, "")
gsea_plot$label <- ifelse(gsea_plot$term %in% top_gsea, gsea_plot$term_short, "")

# =============================================================================
# PLOT 1: FACETED COMPARISON (Combined + Separate)
# =============================================================================

cat("\n1. Creating volcano plots...\n")

# GO panel
p_go <- ggplot(go_plot, aes(x = effect_size_directed, y = log10_padj)) +
  geom_point(aes(color = direction, size = count), alpha = 0.7, shape = 17) +
  geom_text_repel(
    aes(label = label), size = 2.5, max.overlaps = 15,
    box.padding = 0.3, segment.color = "grey50", segment.size = 0.3
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40", alpha = 0.5) +
  scale_color_manual(values = colors_direction, name = "Direction") +
  scale_size_continuous(range = c(2, 8), name = "Gene Count") +
  labs(
    title = "GO Over-Representation Analysis",
    subtitle = "X-axis: Fold Enrichment (directional)",
    x = "Fold Enrichment (- = downreg, + = upreg)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_publication() +
  theme(legend.position = "right")

# GSEA panel
p_gsea <- ggplot(gsea_plot, aes(x = effect_size_directed, y = log10_padj)) +
  geom_point(aes(color = direction, size = count), alpha = 0.7, shape = 16) +
  geom_text_repel(
    aes(label = label), size = 2.5, max.overlaps = 15,
    box.padding = 0.3, segment.color = "grey50", segment.size = 0.3
  ) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey40", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40", alpha = 0.5) +
  scale_color_manual(values = colors_direction, name = "Direction") +
  scale_size_continuous(range = c(2, 8), name = "Gene Set Size") +
  labs(
    title = "Gene Set Enrichment Analysis",
    subtitle = "X-axis: Normalized Enrichment Score",
    x = "Normalized Enrichment Score (NES)",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_publication() +
  theme(legend.position = "right")

# Combined side-by-side (wider)
p_combined <- p_go + p_gsea +
  plot_annotation(
    title = "Comparison of GO ORA and GSEA Enrichment Results",
    subtitle = "Note: X-axes represent different statistics and should not be directly compared",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30")
    )
  )

# Save combined (wider)
ggsave(file.path(OUTPUT_BASE, "01_faceted_comparison.png"), p_combined,
       width = 20, height = 9, dpi = 300)
ggsave(file.path(OUTPUT_BASE, "01_faceted_comparison.pdf"), p_combined,
       width = 20, height = 9, device = cairo_pdf)

# Save separate GO plot
ggsave(file.path(OUTPUT_BASE, "01a_GO_volcano.png"), p_go,
       width = 12, height = 9, dpi = 300)
ggsave(file.path(OUTPUT_BASE, "01a_GO_volcano.pdf"), p_go,
       width = 12, height = 9, device = cairo_pdf)

# Save separate GSEA plot
ggsave(file.path(OUTPUT_BASE, "01b_GSEA_volcano.png"), p_gsea,
       width = 12, height = 9, dpi = 300)
ggsave(file.path(OUTPUT_BASE, "01b_GSEA_volcano.pdf"), p_gsea,
       width = 12, height = 9, device = cairo_pdf)

# =============================================================================
# PLOT 2: SUMMARY (Panels A and C only)
# =============================================================================

cat("2. Creating summary figure...\n")

# Combine data for summary
combined_data <- bind_rows(
  go_plot %>% mutate(method = "GO ORA"),
  gsea_plot %>% mutate(method = "GSEA")
)

# Panel A: Significant terms by method
method_summary <- data.frame(
  Method = c("GO ORA", "GO ORA", "GSEA", "GSEA"),
  Direction = c("Upregulated", "Downregulated", "Upregulated", "Downregulated"),
  Count = c(
    sum(go_plot$direction == "Upregulated"),
    sum(go_plot$direction == "Downregulated"),
    sum(gsea_plot$direction == "Upregulated"),
    sum(gsea_plot$direction == "Downregulated")
  )
)

pA <- ggplot(method_summary, aes(x = Method, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.9),
    vjust = -0.5, size = 4.5, fontface = "bold"
  ) +
  scale_fill_manual(values = colors_direction) +
  labs(
    title = "A) Significant Terms by Method",
    x = NULL,
    y = "Number of Terms"
  ) +
  theme_publication() +
  theme(legend.position = "bottom") +
  ylim(0, max(method_summary$Count) * 1.15)

# Panel C: Direction breakdown (percentage)
direction_summary <- combined_data %>%
  group_by(method, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(method) %>%
  mutate(
    pct = n / sum(n) * 100,
    label = sprintf("%.0f%%", pct)
  )

pC <- ggplot(direction_summary, aes(x = method, y = pct, fill = direction)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = colors_direction, name = "Direction") +
  labs(
    title = "B) Direction Breakdown",
    x = NULL,
    y = "Percentage"
  ) +
  theme_publication() +
  theme(legend.position = "bottom")

# Combine A and C side by side
p_summary <- pA + pC +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "GO ORA and GSEA Enrichment Summary",
    subtitle = "TES vs GFP Differential Expression",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(file.path(OUTPUT_BASE, "02_summary.png"), p_summary,
       width = 14, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_BASE, "02_summary.pdf"), p_summary,
       width = 14, height = 7, device = cairo_pdf)

# =============================================================================
# COMPLETION
# =============================================================================

cat("\n========================================\n")
cat("VISUALIZATION COMPLETE\n")
cat("========================================\n")
cat(sprintf("Output: %s\n", OUTPUT_BASE))
cat("\nGenerated files:\n")
cat("  01_faceted_comparison.png/pdf - Combined GO + GSEA volcano plots\n")
cat("  01a_GO_volcano.png/pdf        - GO ORA volcano plot (separate)\n")
cat("  01b_GSEA_volcano.png/pdf      - GSEA volcano plot (separate)\n")
cat("  02_summary.png/pdf            - Summary with counts and direction breakdown\n")
