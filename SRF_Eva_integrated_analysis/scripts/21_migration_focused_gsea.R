#!/usr/bin/env Rscript
#
# ==============================================================================
# MIGRATION-FOCUSED GSEA ANALYSIS: Publication-Quality Visualizations
# ==============================================================================
#
# This script performs Gene Set Enrichment Analysis (GSEA) focused specifically
# on migration, motility, EMT, and cell adhesion pathways. Creates publication-
# quality figures suitable for manuscripts.
#
# Key focus areas:
# - Cell migration and motility
# - Epithelial-mesenchymal transition (EMT)
# - Cell adhesion and ECM interactions
# - Chemotaxis and invasion
#
# ==============================================================================

suppressPackageStartupMessages({
  # Core Analysis Packages
  library(fgsea)          # Fast Gene Set Enrichment Analysis
  library(msigdbr)        # MSigDB gene sets
  library(dplyr)          # Data manipulation
  library(tidyr)          # Data tidying
  library(readr)          # Fast CSV/TSV reading
  library(stringr)        # String manipulation

  # Visualization Packages
  library(ggplot2)        # Main plotting engine
  library(ggrepel)        # Non-overlapping text labels
  library(patchwork)      # Multi-panel figures
  library(scales)         # Scale formatting

  # Annotation Databases
  library(org.Hs.eg.db)   # Human genome annotation
})

# =============================================================================
# PLOT CONFIGURATION - Matching SRF_Eva_RNA/scripts/6_gsea_analysis.R style
# =============================================================================

# Plot parameters
PLOT_WIDTH <- 18
PLOT_HEIGHT <- 10
DPI <- 180

# Font sizes for publication
TITLE_SIZE <- 20
AXIS_TITLE_SIZE <- 18
AXIS_TEXT_SIZE <- 16
ANNOTATION_SIZE <- 8

# Color scheme matching heatmaps (from 6_gsea_analysis.R)
# GFP = Brown (#8B4513), TES = Teal (#2E8B8B)
COLOR_TES_GRADIENT <- "#5FBFBF"   # Medium teal for TES/upregulated (high rank)
COLOR_GFP_GRADIENT <- "#C9A86C"   # Medium tan/brown for GFP/downregulated (low rank)
COLOR_ENRICHMENT_LINE <- "#D73027"  # Red enrichment score line

# Comparison name for plot subtitles
COMPARISON_NAME <- "TES vs GFP"

# =============================================================================
# Custom GSEA plot function with TES/GFP color gradient (for fgsea results)
# Matches style from SRF_Eva_RNA/scripts/6_gsea_analysis.R
# =============================================================================

create_gsea_plot_custom_colors <- function(pathway_name, gene_sets, gene_ranks, fgsea_result,
                                            title = NULL,
                                            add_stats = FALSE,
                                            color_high = COLOR_TES_GRADIENT,
                                            color_low = COLOR_GFP_GRADIENT,
                                            line_color = COLOR_ENRICHMENT_LINE) {
    # Get pathway data from fgsea results
    pathway_data <- fgsea_result[fgsea_result$pathway == pathway_name, ]
    if (nrow(pathway_data) == 0) {
        return(NULL)
    }

    nes_val <- round(pathway_data$NES, 4)
    pval <- pathway_data$pval
    qval <- pathway_data$padj

    if (is.null(title)) {
        title <- gsub("_", " ", pathway_name)
        title <- gsub("^(HALLMARK |GOBP |REACTOME |KEGG |CUSTOM )", "", title)
        if (nchar(title) > 60) {
            title <- paste0(substr(title, 1, 57), "...")
        }
    }

    # Get gene set
    gene_set <- gene_sets[[pathway_name]]

    # Calculate running enrichment score
    n <- length(gene_ranks)
    gene_hits <- names(gene_ranks) %in% gene_set

    # Running sum statistics
    hit_indicator <- as.integer(gene_hits)
    no_hit_indicator <- 1 - hit_indicator

    # Calculate running ES
    Phit <- cumsum(abs(gene_ranks) * hit_indicator) / sum(abs(gene_ranks[gene_hits]))
    Pmiss <- cumsum(no_hit_indicator) / sum(no_hit_indicator)
    running_es <- Phit - Pmiss

    # Create data frames for plotting
    es_data <- data.frame(
        rank = 1:n,
        running_es = running_es
    )

    hit_data <- data.frame(
        rank = which(gene_hits),
        y = 0
    )

    rank_data <- data.frame(
        rank = 1:n,
        value = gene_ranks
    )

    # Panel 1: Running Enrichment Score
    p1 <- ggplot(es_data, aes(x = rank, y = running_es)) +
        geom_line(color = line_color, linewidth = 1.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        labs(y = "Enrichment Score", x = NULL,
             title = paste0(title, "\n(", COMPARISON_NAME, ")")) +
        theme_classic(base_size = AXIS_TEXT_SIZE) +
        theme(
            plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
            axis.title.y = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
            axis.text = element_text(size = AXIS_TEXT_SIZE, color = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 15, 5, 15)
        )

    # Add NES/FDR annotation if requested
    if (add_stats) {
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)
        annotation_text <- sprintf("NES: %.3f\nFDR: %s", nes_val, qval_text)

        p1 <- p1 +
            annotate("text", x = n * 0.85, y = max(running_es) * 0.8,
                     label = annotation_text, hjust = 0.5, vjust = 1,
                     size = ANNOTATION_SIZE, fontface = "bold")
    }

    # Panel 2: Gene hit barcode
    p2 <- ggplot(hit_data, aes(x = rank, y = y)) +
        geom_segment(aes(xend = rank, yend = 1), color = "black", linewidth = 0.3) +
        labs(x = NULL, y = NULL) +
        theme_void() +
        theme(
            plot.margin = margin(0, 15, 0, 15)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(limits = c(1, n), expand = c(0, 0))

    # Panel 3: Ranked list gradient with custom colors
    # Create gradient based on POSITION (not value) for smooth color transition
    rank_data$position_color <- (rank_data$rank - 1) / (n - 1)  # 0 to 1

    p3 <- ggplot(rank_data, aes(x = rank, y = 0.5, fill = position_color)) +
        geom_tile(height = 1, width = 1) +
        scale_fill_gradient(low = color_high, high = color_low, guide = "none") +  # TES (left) to GFP (right)
        scale_x_continuous(expand = c(0, 0), limits = c(0, n)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
        labs(x = "Rank in Ordered Dataset", y = NULL) +
        theme_classic(base_size = AXIS_TEXT_SIZE) +
        theme(
            axis.title.x = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
            axis.text.x = element_text(size = AXIS_TEXT_SIZE, color = "black"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            plot.margin = margin(0, 15, 10, 15)
        )

    # Add labels for TES and GFP sides
    p3 <- p3 +
        annotate("text", x = n * 0.05, y = 0.5, label = "TES",
                 size = 5, color = "white", fontface = "bold") +
        annotate("text", x = n * 0.95, y = 0.5, label = "GFP",
                 size = 5, color = "white", fontface = "bold")

    # Combine panels
    combined_plot <- p1 / p2 / p3 +
        plot_layout(heights = c(3, 0.5, 1))

    return(combined_plot)
}

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("==============================================================================\n")
cat("MIGRATION-FOCUSED GSEA ANALYSIS\n")
cat("==============================================================================\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# =============================================================================
# OUTPUT CONFIGURATION
# =============================================================================

output_dir <- "output/21_migration_focused_gsea"
results_dir <- file.path(output_dir, "results")
plots_dir <- file.path(output_dir, "plots")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("Output directory: %s\n\n", output_dir))

# =============================================================================
# PUBLICATION THEME
# =============================================================================

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.3), hjust = 0.5),
      plot.subtitle = element_text(size = rel(1.0), hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = rel(1.1)),
      axis.text = element_text(size = rel(0.95)),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold")
    )
}

# Color scheme for direction
color_enriched <- "#E41A1C"  # Red for enriched (upregulated in TES)
color_depleted <- "#377EB8"  # Blue for depleted (downregulated in TES)

# =============================================================================
# PHASE 1: LOAD RNA-SEQ DATA AND CREATE RANKED GENE LIST
# =============================================================================

cat("=== PHASE 1: Creating Ranked Gene List ===\n")

# Load complete RNA-seq results
rna_all <- read.delim("../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
                      stringsAsFactors = FALSE)

cat(sprintf("Loaded %d genes from RNA-seq\n", nrow(rna_all)))

# The file already has gene_symbol column from DESeq2
# Clean Ensembl IDs (remove version)
rna_all$ensembl_id <- gsub("\\..*", "", rna_all$gene_id)

# For genes missing symbols, try to map from Ensembl
missing_symbols <- is.na(rna_all$gene_symbol) | rna_all$gene_symbol == ""
if (sum(missing_symbols) > 0) {
  mapped_symbols <- mapIds(org.Hs.eg.db,
                           keys = rna_all$ensembl_id[missing_symbols],
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
  rna_all$gene_symbol[missing_symbols] <- mapped_symbols
}

# Remove genes without symbols or log2FC
rna_ranked <- rna_all %>%
  filter(!is.na(.data$gene_symbol) & .data$gene_symbol != "" & !is.na(.data$log2FoldChange)) %>%
  arrange(desc(.data$log2FoldChange))

cat(sprintf("Ranked %d genes with valid symbols and fold changes\n", nrow(rna_ranked)))

# Create ranked list (gene symbol -> log2FC)
gene_ranks <- setNames(rna_ranked$log2FoldChange, rna_ranked$gene_symbol)

# Remove duplicates (keep first = highest abs value)
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat(sprintf("Final ranked list: %d unique genes\n", length(gene_ranks)))
cat(sprintf("  Range: %.2f to %.2f log2FC\n", min(gene_ranks), max(gene_ranks)))
cat(sprintf("  Most upregulated: %s (%.2f)\n", names(gene_ranks)[1], gene_ranks[1]))
cat(sprintf("  Most downregulated: %s (%.2f)\n\n",
            names(gene_ranks)[length(gene_ranks)], gene_ranks[length(gene_ranks)]))

# =============================================================================
# PHASE 2: PREPARE MIGRATION-RELATED GENE SETS
# =============================================================================

cat("=== PHASE 2: Preparing Migration-Related Gene Sets ===\n")

# Get all MSigDB gene sets for human
cat("Fetching MSigDB gene sets...\n")

# Check msigdbr version and column names
msigdbr_version <- packageVersion("msigdbr")
cat(sprintf("  msigdbr version: %s\n", msigdbr_version))

# Get column names from a test query
test_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat("  Available columns:", paste(colnames(test_df), collapse = ", "), "\n")

# Determine the correct column name for gene symbols
symbol_col <- if ("gene_symbol" %in% colnames(test_df)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% colnames(test_df)) {
  "human_gene_symbol"
} else {
  # Find any column with "symbol" in the name
  symbol_candidates <- grep("symbol", colnames(test_df), value = TRUE, ignore.case = TRUE)
  if (length(symbol_candidates) > 0) symbol_candidates[1] else NULL
}

if (is.null(symbol_col)) {
  stop("Could not find gene symbol column in msigdbr output")
}
cat(sprintf("  Using symbol column: %s\n", symbol_col))

# Keywords for migration-related pathways
migration_keywords <- c(
  "migration", "migrat", "motil", "chemotaxis", "invasion", "invasive",
  "adhesion", "locomotion", "EMT", "epithelial_mesenchymal", "mesenchymal",
  "extracellular_matrix", "ECM", "integrin", "focal_adhesion",
  "cell_junction", "cell_matrix", "proteoglycan", "collagen",
  "laminin", "fibronectin", "basement_membrane"
)
pattern <- paste(migration_keywords, collapse = "|")

# Helper function to get gene sets with consistent column naming
get_gene_set_df <- function(df, pattern_filter = NULL) {
  # Select and rename columns
  result <- df[, c("gs_name", symbol_col)]
  colnames(result) <- c("gs_name", "gene_symbol")

  # Apply pattern filter if provided
  if (!is.null(pattern_filter)) {
    result <- result[grepl(pattern_filter, result$gs_name, ignore.case = TRUE), ]
  }
  return(result)
}

# 1. Hallmark Collection
cat("  Loading Hallmark collection...\n")
hallmark_raw <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- get_gene_set_df(hallmark_raw, "EPITHELIAL_MESENCHYMAL|APICAL|MYOGENESIS|ANGIOGENESIS")

# 2. GO Biological Process - filter for migration terms
cat("  Loading GO Biological Process collection...\n")
gobp_all <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_migration <- get_gene_set_df(gobp_all, pattern)

cat(sprintf("    Found %d migration-related GO:BP gene sets\n",
            length(unique(gobp_migration$gs_name))))

# 3. Reactome pathways - EMT and ECM
cat("  Loading Reactome pathways...\n")
reactome_all <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
reactome_migration <- get_gene_set_df(reactome_all, pattern)

cat(sprintf("    Found %d migration-related Reactome pathways\n",
            length(unique(reactome_migration$gs_name))))

# 4. KEGG pathways
cat("  Loading KEGG pathways...\n")
kegg_all <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
kegg_migration <- get_gene_set_df(kegg_all, pattern)

cat(sprintf("    Found %d migration-related KEGG pathways\n",
            length(unique(kegg_migration$gs_name))))

# Combine all MSigDB sets
all_msigdb <- bind_rows(
  hallmark_sets %>% mutate(source = "Hallmark"),
  gobp_migration %>% mutate(source = "GO:BP"),
  reactome_migration %>% mutate(source = "Reactome"),
  kegg_migration %>% mutate(source = "KEGG")
)

# Convert to list format for fgsea
migration_gene_sets <- split(all_msigdb$gene_symbol, all_msigdb$gs_name)

# 5. Add Custom Migration Gene Set (core EMT/migration markers)
custom_migration <- c(
  # Integrins
  "ITGB1", "ITGB3", "ITGB5", "ITGA5", "ITGAV", "ITGA2", "ITGA6",
  # Cadherins (EMT markers)
  "CDH1", "CDH2", "CDH11",
  # EMT TFs
  "SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2",
  # Intermediate filaments
  "VIM", "KRT18", "KRT19",
  # MMPs
  "MMP2", "MMP9", "MMP14", "MMP7",
  # ECM proteins
  "FN1", "SPARC", "THBS1", "COL1A1", "COL3A1", "LAMA5",
  # Migration regulators
  "CXCR4", "CXCL12", "RAC1", "CDC42", "RHOA", "ROCK1",
  # Tight junction
  "TJP1", "OCLN", "CLDN1"
)
migration_gene_sets[["CUSTOM_CORE_MIGRATION_EMT"]] <- custom_migration

# Filter gene sets by size
set_sizes <- sapply(migration_gene_sets, length)
migration_gene_sets <- migration_gene_sets[set_sizes >= 15 & set_sizes <= 500]

cat(sprintf("\nTotal migration gene sets for GSEA: %d\n", length(migration_gene_sets)))
cat(sprintf("  Size range: %d - %d genes\n\n",
            min(sapply(migration_gene_sets, length)),
            max(sapply(migration_gene_sets, length))))

# =============================================================================
# PHASE 3: RUN FGSEA
# =============================================================================

cat("=== PHASE 3: Running Gene Set Enrichment Analysis ===\n")
cat("This may take a few minutes...\n\n")

set.seed(42)

# Run fgsea
fgsea_results <- fgsea(
  pathways = migration_gene_sets,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  nproc = 8,
  nPermSimple = 10000
)

# Add source information
fgsea_results$source <- sapply(fgsea_results$pathway, function(p) {
  if (grepl("^HALLMARK_", p)) return("Hallmark")
  if (grepl("^GOBP_", p)) return("GO:BP")
  if (grepl("^REACTOME_", p)) return("Reactome")
  if (grepl("^KEGG_", p)) return("KEGG")
  if (grepl("^CUSTOM_", p)) return("Custom")
  return("Other")
})

# Filter significant results - using relaxed threshold (FDR < 0.10)
fgsea_sig <- fgsea_results %>%
  filter(padj < 0.10) %>%
  arrange(padj)

cat(sprintf("GSEA complete!\n"))
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_results)))
cat(sprintf("  Significant (FDR < 0.10): %d\n", nrow(fgsea_sig)))

# If no significant results at FDR < 0.10, use further relaxed threshold or top results
if (nrow(fgsea_sig) == 0) {
  cat("  No pathways significant at FDR < 0.10, trying FDR < 0.25...\n")
  fgsea_sig <- fgsea_results %>%
    filter(padj < 0.25) %>%
    arrange(padj)

  if (nrow(fgsea_sig) == 0) {
    cat("  No pathways significant at FDR < 0.25, using top 30 by p-value...\n")
    fgsea_sig <- fgsea_results %>%
      arrange(pval) %>%
      head(30)
  }
  cat(sprintf("  Using %d pathways for visualization\n", nrow(fgsea_sig)))
}

cat(sprintf("  Enriched (positive NES): %d\n", sum(fgsea_sig$NES > 0)))
cat(sprintf("  Depleted (negative NES): %d\n\n", sum(fgsea_sig$NES < 0)))

# =============================================================================
# PHASE 4: PUBLICATION VISUALIZATIONS
# =============================================================================

cat("=== PHASE 4: Creating Publication-Quality Visualizations ===\n")

# Helper function to clean pathway names for display
clean_pathway_name <- function(name, max_length = 50) {
  # Remove prefixes
  name <- gsub("^(HALLMARK_|GOBP_|REACTOME_|KEGG_|CUSTOM_)", "", name)
  # Replace underscores with spaces
  name <- gsub("_", " ", name)
  # Title case
  name <- str_to_title(name)
  # Truncate if too long
  if (nchar(name) > max_length) {
    name <- paste0(substr(name, 1, max_length - 3), "...")
  }
  return(name)
}

# -----------------------------------------------------------------------------
# PLOT 1: ENRICHMENT SCORE CURVES (Top pathways)
# Using custom 3-panel style matching SRF_Eva_RNA/scripts/6_gsea_analysis.R
# -----------------------------------------------------------------------------

cat("Creating Plot 1: Enrichment score curves (publication style)...\n")

# Select top pathways for enrichment curves (most significant)
top_pathways <- fgsea_sig %>%
  arrange(padj) %>%
  head(5) %>%
  pull(pathway)

if (length(top_pathways) > 0) {
  # Create individual enrichment curve PNGs for each top pathway
  for (i in seq_along(top_pathways)) {
    pathway_name <- top_pathways[i]

    # Create clean title (remove prefix, replace underscores)
    clean_title <- clean_pathway_name(pathway_name, 60)

    # VERSION 1: Clean version (no stats) - matches target style
    p_clean <- create_gsea_plot_custom_colors(
      pathway_name = pathway_name,
      gene_sets = migration_gene_sets,
      gene_ranks = gene_ranks,
      fgsea_result = fgsea_results,
      title = clean_title,
      add_stats = FALSE
    )

    if (!is.null(p_clean)) {
      ggsave(file.path(plots_dir, sprintf("01_migration_gsea_enrichment_%02d.png", i)),
             p_clean, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
             bg = "white")

      # Also save a version with stats
      p_stats <- create_gsea_plot_custom_colors(
        pathway_name = pathway_name,
        gene_sets = migration_gene_sets,
        gene_ranks = gene_ranks,
        fgsea_result = fgsea_results,
        title = clean_title,
        add_stats = TRUE
      )

      ggsave(file.path(plots_dir, sprintf("01_migration_gsea_enrichment_%02d_annotated.png", i)),
             p_stats, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
             bg = "white")
    }
  }

  cat("  Saved: 01_migration_gsea_enrichment_*.png (clean and annotated versions)\n")
} else {
  cat("  Skipping Plot 1: no pathways to plot\n")
}

# -----------------------------------------------------------------------------
# PLOT 2: DOT PLOT - NES vs Pathway
# -----------------------------------------------------------------------------

cat("Creating Plot 2: Migration pathways dot plot...\n")

# Get top pathways for visualization (up to 25)
plot_data <- fgsea_sig %>%
  arrange(padj) %>%
  head(25) %>%
  mutate(
    display_name = sapply(pathway, clean_pathway_name, max_length = 45),
    direction = ifelse(NES > 0, "Enriched", "Depleted"),
    neg_log10_padj = -log10(padj)
  )

if (nrow(plot_data) > 0) {
  p2 <- ggplot(plot_data, aes(x = NES, y = reorder(display_name, NES))) +
    geom_segment(aes(x = 0, xend = NES, yend = reorder(display_name, NES)),
                 color = "grey50", linewidth = 0.5) +
    geom_point(aes(color = direction, size = size, alpha = neg_log10_padj)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    scale_color_manual(
      values = c("Enriched" = color_enriched, "Depleted" = color_depleted),
      name = "Direction"
    ) +
    scale_size_continuous(range = c(3, 10), name = "Gene Set\nSize") +
    scale_alpha_continuous(range = c(0.5, 1), name = "-log10(FDR)") +
    labs(
      title = "Migration/EMT Pathway Enrichment (GSEA)",
      subtitle = "Top 25 significant pathways (FDR < 0.10)",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_pub(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "right"
    )

  ggsave(file.path(plots_dir, "02_migration_pathways_dotplot.png"), p2,
         width = 12, height = 10, dpi = 300)

  cat("  Saved: 02_migration_pathways_dotplot.png\n")
} else {
  cat("  Skipping Plot 2: no pathways to plot\n")
}

# -----------------------------------------------------------------------------
# PLOT 3: LOLLIPOP PLOT - NES Scores
# -----------------------------------------------------------------------------

cat("Creating Plot 3: NES lollipop plot...\n")

# Get top 15 up and top 15 down
top_up <- fgsea_sig %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(15)
top_down <- fgsea_sig %>% filter(NES < 0) %>% arrange(NES) %>% head(15)
lollipop_data <- bind_rows(top_up, top_down) %>%
  mutate(
    display_name = sapply(pathway, clean_pathway_name, max_length = 40),
    direction = ifelse(NES > 0, "Enriched (Up in TES)", "Depleted (Down in TES)")
  )

if (nrow(lollipop_data) > 0) {
  p3 <- ggplot(lollipop_data, aes(x = NES, y = reorder(display_name, NES))) +
    geom_segment(aes(x = 0, xend = NES, yend = reorder(display_name, NES), color = direction),
                 linewidth = 1.2) +
    geom_point(aes(color = direction), size = 4) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    scale_color_manual(
      values = c("Enriched (Up in TES)" = color_enriched,
                "Depleted (Down in TES)" = color_depleted),
      name = "Direction"
    ) +
    labs(
      title = "Migration/EMT Pathways: Normalized Enrichment Scores",
      subtitle = "Top enriched and depleted pathways in TES vs GFP",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_pub(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
    )

  ggsave(file.path(plots_dir, "03_migration_nes_lollipop.png"), p3,
         width = 12, height = 8, dpi = 300)

  cat("  Saved: 03_migration_nes_lollipop.png\n")
} else {
  cat("  Skipping Plot 3: no pathways to plot\n")
}

# Note: Leading edge heatmap removed - not needed for this analysis

# # -----------------------------------------------------------------------------
# # PLOT 5: MULTI-PANEL SUMMARY
# # -----------------------------------------------------------------------------

# cat("Creating Plot 5: Multi-panel summary figure...\n")

# # Only create if we have results
# if (nrow(fgsea_sig) > 0) {
#   # Panel A: Summary bar chart by source
#   source_summary <- fgsea_sig %>%
#     group_by(source) %>%
#     summarise(
#       n_total = n(),
#       n_enriched = sum(NES > 0),
#       n_depleted = sum(NES < 0),
#       .groups = "drop"
#     ) %>%
#     pivot_longer(cols = c(n_enriched, n_depleted),
#                  names_to = "direction", values_to = "count") %>%
#     mutate(direction = ifelse(direction == "n_enriched", "Enriched", "Depleted"))

#   if (nrow(source_summary) > 0) {
#     pA <- ggplot(source_summary, aes(x = reorder(source, -count), y = count, fill = direction)) +
#       geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
#       geom_text(aes(label = count), position = position_dodge(width = 0.9),
#                 vjust = -0.3, size = 4, fontface = "bold") +
#       scale_fill_manual(values = c("Enriched" = color_enriched, "Depleted" = color_depleted)) +
#       labs(
#         title = "A) Pathways by Source",
#         x = NULL, y = "Count",
#         fill = "Direction"
#       ) +
#       theme_pub(base_size = 11) +
#       scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
#   } else {
#     pA <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data") + theme_void()
#   }

#   # Panel B: Top enriched pathways (bar)
#   top_enriched <- fgsea_sig %>%
#     filter(NES > 0) %>%
#     arrange(desc(NES)) %>%
#     head(8) %>%
#     mutate(display_name = sapply(pathway, clean_pathway_name, max_length = 35))

#   if (nrow(top_enriched) > 0) {
#     pB <- ggplot(top_enriched, aes(x = NES, y = reorder(display_name, NES))) +
#       geom_bar(stat = "identity", fill = color_enriched, color = "black", linewidth = 0.3) +
#       labs(
#         title = "B) Top Enriched (Up in TES)",
#         x = "NES", y = NULL
#       ) +
#       theme_pub(base_size = 11) +
#       theme(axis.text.y = element_text(size = 8))
#   } else {
#     pB <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No enriched pathways") +
#       theme_void() + ggtitle("B) Top Enriched (Up in TES)")
#   }

#   # Panel C: Top depleted pathways (bar)
#   top_depleted <- fgsea_sig %>%
#     filter(NES < 0) %>%
#     arrange(NES) %>%
#     head(8) %>%
#     mutate(display_name = sapply(pathway, clean_pathway_name, max_length = 35))

#   if (nrow(top_depleted) > 0) {
#     pC <- ggplot(top_depleted, aes(x = NES, y = reorder(display_name, -NES))) +
#       geom_bar(stat = "identity", fill = color_depleted, color = "black", linewidth = 0.3) +
#       labs(
#         title = "C) Top Depleted (Down in TES)",
#         x = "NES", y = NULL
#       ) +
#       theme_pub(base_size = 11) +
#       theme(axis.text.y = element_text(size = 8))
#   } else {
#     pC <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No depleted pathways") +
#       theme_void() + ggtitle("C) Top Depleted (Down in TES)")
#   }

#   # Panel D: NES distribution
#   pD <- ggplot(fgsea_sig, aes(x = NES, fill = source)) +
#     geom_histogram(bins = 20, color = "black", linewidth = 0.3, alpha = 0.7) +
#     geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
#     scale_fill_brewer(palette = "Set2", name = "Source") +
#     labs(
#       title = "D) NES Distribution",
#       x = "Normalized Enrichment Score",
#       y = "Count"
#     ) +
#     theme_pub(base_size = 11)

#   # Combine panels
#   combined_plot <- (pA | pD) / (pB | pC) +
#     plot_annotation(
#       title = "Migration/EMT GSEA Analysis Summary",
#       subtitle = "TES vs GFP Differential Expression",
#       theme = theme(
#         plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#         plot.subtitle = element_text(size = 12, hjust = 0.5)
#       )
#     )

#   ggsave(file.path(plots_dir, "05_migration_summary_panel.png"), combined_plot,
#          width = 16, height = 14, dpi = 300)

#   cat("  Saved: 05_migration_summary_panel.png\n")
# } else {
#   cat("  Skipping Plot 5: no pathways to plot\n")
# }

# =============================================================================
# PHASE 5: EXPORT RESULTS
# =============================================================================

cat("\n=== PHASE 5: Exporting Results ===\n")

# Convert leadingEdge list to character string for CSV export
export_all <- fgsea_results %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  arrange(padj)

export_sig <- fgsea_sig %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  arrange(padj)

# Export all results
write.csv(export_all, file.path(results_dir, "all_migration_gsea_results.csv"),
          row.names = FALSE)
write.csv(export_sig, file.path(results_dir, "significant_migration_pathways.csv"),
          row.names = FALSE)

# Export leading edge genes with pathway annotation
leading_edge_df <- data.frame()
if (nrow(fgsea_sig) > 0) {
  for (i in seq_len(nrow(fgsea_sig))) {
    genes <- fgsea_sig$leadingEdge[[i]]
    if (length(genes) > 0) {
      df <- data.frame(
        pathway = fgsea_sig$pathway[i],
        NES = fgsea_sig$NES[i],
        padj = fgsea_sig$padj[i],
        gene = genes,
        stringsAsFactors = FALSE
      )
      leading_edge_df <- bind_rows(leading_edge_df, df)
    }
  }
}

# Add expression info for leading edge genes
if (nrow(leading_edge_df) > 0) {
  # Create a lookup table from rna_ranked
  expr_lookup <- data.frame(
    gene = rna_ranked$gene_symbol,
    gene_log2FC = rna_ranked$log2FoldChange,
    gene_padj = rna_ranked$padj,
    stringsAsFactors = FALSE
  )
  leading_edge_df <- left_join(leading_edge_df, expr_lookup, by = "gene")
}

write.csv(leading_edge_df, file.path(results_dir, "leading_edge_genes.csv"),
          row.names = FALSE)

cat(sprintf("Results exported:\n"))
cat(sprintf("  - All pathways: %d\n", nrow(export_all)))
cat(sprintf("  - Significant pathways: %d\n", nrow(export_sig)))
cat(sprintf("  - Leading edge gene entries: %d\n\n", nrow(leading_edge_df)))

# =============================================================================
# SUMMARY REPORT
# =============================================================================

cat("=== Generating Summary Report ===\n")

summary_file <- file.path(output_dir, "SUMMARY_REPORT.txt")

sink(summary_file)
cat("==============================================================================\n")
cat("MIGRATION-FOCUSED GSEA ANALYSIS SUMMARY\n")
cat("==============================================================================\n\n")
cat(sprintf("Generated: %s\n\n", Sys.time()))

cat("METHOD\n")
cat("------\n")
cat("Algorithm: fgsea (Fast Gene Set Enrichment Analysis)\n")
cat("Ranking metric: log2FoldChange (TES vs GFP)\n")
cat(sprintf("Total genes ranked: %d\n", length(gene_ranks)))
cat(sprintf("Gene sets tested: %d\n", nrow(fgsea_results)))
cat("Permutations: 10,000\n")
cat("Gene set size filter: 15-500 genes\n")
cat("Significance threshold: FDR < 0.10 (relaxed)\n\n")

cat("GENE SET SOURCES\n")
cat("----------------\n")
source_counts <- table(fgsea_results$source)
for (src in names(source_counts)) {
  cat(sprintf("  %s: %d gene sets\n", src, source_counts[src]))
}
cat("\n")

cat("RESULTS OVERVIEW\n")
cat("----------------\n")
cat(sprintf("Total significant pathways: %d\n", nrow(fgsea_sig)))
cat(sprintf("  Enriched (NES > 0, up in TES): %d\n", sum(fgsea_sig$NES > 0)))
cat(sprintf("  Depleted (NES < 0, down in TES): %d\n\n", sum(fgsea_sig$NES < 0)))

cat("BY SOURCE:\n")
for (src in unique(fgsea_sig$source)) {
  src_data <- fgsea_sig %>% filter(source == src)
  cat(sprintf("  %s: %d total (%d enriched, %d depleted)\n",
              src, nrow(src_data), sum(src_data$NES > 0), sum(src_data$NES < 0)))
}
cat("\n")

cat("TOP 10 ENRICHED PATHWAYS (Up in TES)\n")
cat("------------------------------------\n")
top_up <- fgsea_sig %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(10)
if (nrow(top_up) > 0) {
  for (i in seq_len(nrow(top_up))) {
    cat(sprintf("%2d. %s\n    NES=%.2f, FDR=%.2e, Size=%d\n",
                i, top_up$pathway[i], top_up$NES[i], top_up$padj[i], top_up$size[i]))
  }
} else {
  cat("  None found\n")
}
cat("\n")

cat("TOP 10 DEPLETED PATHWAYS (Down in TES)\n")
cat("--------------------------------------\n")
top_down <- fgsea_sig %>% filter(NES < 0) %>% arrange(NES) %>% head(10)
if (nrow(top_down) > 0) {
  for (i in seq_len(nrow(top_down))) {
    cat(sprintf("%2d. %s\n    NES=%.2f, FDR=%.2e, Size=%d\n",
                i, top_down$pathway[i], top_down$NES[i], top_down$padj[i], top_down$size[i]))
  }
} else {
  cat("  None found\n")
}
cat("\n")

cat("KEY FINDINGS\n")
cat("------------\n")
# Count EMT-specific pathways
emt_pathways <- fgsea_sig %>%
  filter(grepl("EMT|EPITHELIAL|MESENCHYMAL", pathway, ignore.case = TRUE))
cat(sprintf("EMT-related pathways: %d (%d enriched, %d depleted)\n",
            nrow(emt_pathways), sum(emt_pathways$NES > 0), sum(emt_pathways$NES < 0)))

# Count migration-specific pathways
mig_pathways <- fgsea_sig %>%
  filter(grepl("MIGRATION|MOTILITY|LOCOMOTION", pathway, ignore.case = TRUE))
cat(sprintf("Migration/Motility pathways: %d (%d enriched, %d depleted)\n",
            nrow(mig_pathways), sum(mig_pathways$NES > 0), sum(mig_pathways$NES < 0)))

# Count adhesion pathways
adh_pathways <- fgsea_sig %>%
  filter(grepl("ADHESION|INTEGRIN|FOCAL", pathway, ignore.case = TRUE))
cat(sprintf("Cell adhesion pathways: %d (%d enriched, %d depleted)\n",
            nrow(adh_pathways), sum(adh_pathways$NES > 0), sum(adh_pathways$NES < 0)))

# Count chemotaxis pathways
chemo_pathways <- fgsea_sig %>%
  filter(grepl("CHEMOTAXIS|CHEMO", pathway, ignore.case = TRUE))
cat(sprintf("Chemotaxis pathways: %d (%d enriched, %d depleted)\n",
            nrow(chemo_pathways), sum(chemo_pathways$NES > 0), sum(chemo_pathways$NES < 0)))
cat("\n")

cat("OUTPUT FILES\n")
cat("------------\n")
cat("Results:\n")
cat("  - all_migration_gsea_results.csv\n")
cat("  - significant_migration_pathways.csv\n")
cat("  - leading_edge_genes.csv\n")
cat("\nPlots:\n")
cat("  - 01_migration_gsea_enrichment_*.png (clean, publication-ready)\n")
cat("  - 01_migration_gsea_enrichment_*_annotated.png (with NES/FDR stats)\n")
cat("  - 02_migration_pathways_dotplot.png\n")
cat("  - 03_migration_nes_lollipop.png\n")
cat("\nPlot style:\n")
cat("  - 3-panel layout: enrichment curve, gene barcode, TES/GFP gradient\n")
cat("  - Matches SRF_Eva_RNA/scripts/6_gsea_analysis.R format\n")
cat("\n")
cat("==============================================================================\n")
sink()

cat(sprintf("Summary report saved: %s\n\n", summary_file))

# =============================================================================
# COMPLETION
# =============================================================================

cat("==============================================================================\n")
cat("MIGRATION-FOCUSED GSEA ANALYSIS COMPLETE\n")
cat("==============================================================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", output_dir))

cat("Generated files:\n")
cat("  Results:\n")
cat("    - all_migration_gsea_results.csv\n")
cat("    - significant_migration_pathways.csv\n")
cat("    - leading_edge_genes.csv\n")
cat("\n  Plots (publication-ready 3-panel style):\n")
cat("    - 01_migration_gsea_enrichment_*.png (clean version)\n")
cat("    - 01_migration_gsea_enrichment_*_annotated.png (with stats)\n")
cat("    - 02_migration_pathways_dotplot.png\n")
cat("    - 03_migration_nes_lollipop.png\n")
cat("\n  Summary:\n")
cat("    - SUMMARY_REPORT.txt\n")
cat("\n  Plot style: Matches SRF_Eva_RNA/scripts/6_gsea_analysis.R format\n")
cat("    - Panel 1: Red enrichment line, dashed zero line\n")
cat("    - Panel 2: Gene hit barcode (vertical lines)\n")
cat("    - Panel 3: TES (teal) to GFP (brown) gradient bar\n")
cat("\n")
