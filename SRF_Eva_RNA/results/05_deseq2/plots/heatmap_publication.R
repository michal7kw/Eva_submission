#!/usr/bin/env Rscript

# ===============================================================================
# Publication-Ready Heatmap Script
# Creates three versions: Top DEGs, Gene List, All DEGs
# Style: Blue-white-red, z-score scaled, NO dendrograms
# ===============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(org.Hs.eg.db)
    library(grDevices)
    library(grid)
})

cat("=== Publication-Ready Heatmap Generation ===\n")
cat("Timestamp:", as.character(Sys.time()), "\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/05_deseq2")

# ===============================================================================
# Configuration
# ===============================================================================

# Significance thresholds
PADJ_CUTOFF <- 0.05
FC_CUTOFF <- 1.0 # |log2FC| > 1

# Number of top genes for top DEGs heatmap
N_TOP_GENES <- 50

# Plot dimensions (INCREASED for larger fonts)
PLOT_WIDTH <- 16 # Increased to give more space for legend
PLOT_HEIGHT <- 16
DPI <- 300

# Font sizes (publication-ready - INCREASED for readability)
FONTSIZE_MAIN <- 22 # Main title
FONTSIZE_ROW <- 14 # Row labels (gene names)
FONTSIZE_COL <- 18 # Column labels
FONTSIZE_LEGEND <- 18 # Legend text (BIGGER)
FONTSIZE_ANNOTATION <- 16 # Annotation labels

# Color palette: Dark blue to Intensive red (continuous gradient)
color_palette <- colorRampPalette(c(
    "#08306B", "#08519C", "#2171B5", "#4292C6",
    "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7",
    "white",
    "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
    "#EF3B2C", "#CB181D", "#A50F15", "#67000D"
))(100)

# Annotation colors (matching reference image)
annotation_colors <- list(
    condition = c("GFP" = "#921100", "TES" = "#069093"), # Dark red and Teal
    samples = c("1" = "#4D4D4D", "2" = "#FFB6C1", "3" = "#CD5C5C") # Grey, Pink, Red
)

# ===============================================================================
# Load data
# ===============================================================================

cat("Loading data...\n")

# Load DESeq2 results
res <- read.table("deseq2_results_TES_vs_GFP.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
cat(sprintf("  DESeq2 results: %d genes\n", nrow(res)))

# Load normalized counts
counts <- read.table("normalized_counts.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
    row.names = 1
)
cat(sprintf("  Normalized counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# Load custom gene list
gene_list_file <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_integrated_analysis/data/TES_degs.txt"
custom_genes <- unique(trimws(readLines(gene_list_file)))
custom_genes <- custom_genes[custom_genes != ""]
cat(sprintf("  Custom gene list: %d unique genes\n", length(custom_genes)))

# ===============================================================================
# Prepare sample metadata
# ===============================================================================

# Create sample metadata from column names
sample_names <- colnames(counts)
sample_data <- data.frame(
    sample = sample_names,
    condition = ifelse(grepl("GFP", sample_names), "GFP", "TES"),
    samples = gsub(".*([0-9])$", "\\1", sample_names),
    stringsAsFactors = FALSE
)
rownames(sample_data) <- sample_names

# Order samples: GFP first, then TES
sample_order <- c(
    sample_names[grepl("GFP", sample_names)],
    sample_names[grepl("TES", sample_names)]
)

cat("\nSample order:\n")
print(sample_data[sample_order, ])

# ===============================================================================
# Create gene symbol to Ensembl ID mapping
# ===============================================================================

cat("\nMapping gene symbols to Ensembl IDs...\n")

# Clean Ensembl IDs (remove version)
res$ensembl_clean <- gsub("\\..*", "", res$gene_id)

# Create mapping from symbol to Ensembl ID
symbol_to_ensembl <- setNames(res$gene_id, res$gene_symbol)

# Also create mapping from counts rownames
counts_ensembl_clean <- gsub("\\..*", "", rownames(counts))
ensembl_to_symbol <- setNames(
    res$gene_symbol[match(counts_ensembl_clean, res$ensembl_clean)],
    rownames(counts)
)

# ===============================================================================
# Helper function to create heatmap
# ===============================================================================

create_publication_heatmap <- function(gene_ids, title, filename_base,
                                       show_rownames = TRUE, max_genes = 100,
                                       create_nolabel_version = TRUE,
                                       update_title_with_count = TRUE) {
    cat(sprintf("\nCreating heatmap: %s\n", title))
    cat(sprintf("  Input genes: %d\n", length(gene_ids)))

    # Filter to genes present in counts
    available_genes <- gene_ids[gene_ids %in% rownames(counts)]
    cat(sprintf("  Genes found in counts: %d\n", length(available_genes)))

    if (length(available_genes) == 0) {
        cat("  WARNING: No genes found, skipping this heatmap\n")
        return(NULL)
    }

    # Get log2FC for each gene and sort by it (upregulated at top, downregulated at bottom)
    gene_lfc <- res$log2FoldChange[match(gsub("\\..*", "", available_genes), res$ensembl_clean)]
    available_genes <- available_genes[order(gene_lfc, decreasing = TRUE, na.last = TRUE)]
    gene_lfc <- gene_lfc[order(gene_lfc, decreasing = TRUE, na.last = TRUE)]

    cat(sprintf("  Genes sorted by log2FC (up to down)\n"))
    cat(sprintf("    Upregulated (FC > 0): %d\n", sum(gene_lfc > 0, na.rm = TRUE)))
    cat(sprintf("    Downregulated (FC < 0): %d\n", sum(gene_lfc < 0, na.rm = TRUE)))

    # Limit number of genes if too many (for readability)
    # Keep balanced: top upregulated and top downregulated
    if (length(available_genes) > max_genes && show_rownames) {
        cat(sprintf("  Limiting to %d genes for readability (balanced up/down)\n", max_genes))
        n_up <- sum(gene_lfc > 0, na.rm = TRUE)
        n_down <- sum(gene_lfc < 0, na.rm = TRUE)

        # Take proportional amounts from up and down
        take_up <- min(n_up, ceiling(max_genes * n_up / length(available_genes)))
        take_down <- min(n_down, max_genes - take_up)

        up_genes <- available_genes[gene_lfc > 0][1:take_up]
        down_genes <- available_genes[gene_lfc < 0][1:take_down]
        available_genes <- c(up_genes, down_genes)
        available_genes <- available_genes[!is.na(available_genes)]
    }

    # Extract count data
    mat <- as.matrix(counts[available_genes, sample_order])

    # Z-score normalize by row
    mat_scaled <- t(scale(t(mat)))

    # Cap extreme values for better visualization
    mat_scaled[mat_scaled > 2] <- 2
    mat_scaled[mat_scaled < -2] <- -2

    # Handle any NA/Inf values
    mat_scaled[is.na(mat_scaled)] <- 0
    mat_scaled[is.infinite(mat_scaled)] <- 0

    # Get gene symbols for row labels
    row_symbols <- ensembl_to_symbol[available_genes]

    # Filter out genes without valid gene symbol mapping (keep only genes with proper names)
    valid_symbol_mask <- !is.na(row_symbols) & row_symbols != "" & !grepl("^ENSG", row_symbols)

    if (sum(valid_symbol_mask) == 0) {
        cat("  WARNING: No genes with valid gene symbols found, skipping this heatmap\n")
        return(NULL)
    }

    if (sum(!valid_symbol_mask) > 0) {
        cat(sprintf("  Filtering out %d genes without valid gene symbol mapping\n", sum(!valid_symbol_mask)))
    }

    # Keep only genes with valid symbols
    mat_scaled <- mat_scaled[valid_symbol_mask, , drop = FALSE]
    row_symbols <- row_symbols[valid_symbol_mask]
    available_genes <- available_genes[valid_symbol_mask]

    rownames(mat_scaled) <- row_symbols
    cat(sprintf("  Final genes for heatmap: %d\n", length(row_symbols)))

    # Update title with actual gene count if requested
    if (update_title_with_count) {
        # Replace any number in parentheses with actual count, or append count
        if (grepl("\\(n=\\d+\\)", title)) {
            title <- gsub("\\(n=\\d+\\)", sprintf("(n=%d)", length(row_symbols)), title)
        } else if (grepl("Top \\d+", title)) {
            title <- gsub("Top \\d+", sprintf("Top %d", length(row_symbols)), title)
        }
    }

    # Prepare annotation
    annotation_col <- sample_data[sample_order, c("condition", "samples"), drop = FALSE]

    # Determine row font size based on number of genes (INCREASED for publication)
    fontsize_row <- if (length(available_genes) <= 30) {
        FONTSIZE_ROW
    } else if (length(available_genes) <= 50) {
        FONTSIZE_ROW - 2
    } else if (length(available_genes) <= 100) {
        FONTSIZE_ROW - 3
    } else {
        FONTSIZE_ROW - 4
    }

    # Adjust cell height based on number of genes
    cellheight <- if (length(available_genes) <= 30) {
        22
    } else if (length(available_genes) <= 50) {
        18
    } else if (length(available_genes) <= 100) {
        14
    } else if (length(available_genes) <= 200) {
        8
    } else if (length(available_genes) <= 500) {
        4
    } else {
        2
    } # Very small for large heatmaps (>500 genes)

    # Create heatmap (NO clustering/dendrograms)
    p <- pheatmap(
        mat_scaled,
        color = color_palette,
        breaks = seq(-2, 2, length.out = 101),
        cluster_rows = FALSE, # NO row dendrogram
        cluster_cols = FALSE, # NO column dendrogram
        show_rownames = show_rownames,
        show_colnames = FALSE,
        annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        annotation_names_col = FALSE, # Hide annotation names to avoid overlap
        annotation_legend = TRUE,
        legend = TRUE,
        fontsize = FONTSIZE_LEGEND,
        fontsize_row = fontsize_row,
        fontsize_col = FONTSIZE_COL,
        cellwidth = 50, # Cell width (reduced to give space for legend)
        cellheight = cellheight, # Adaptive cell height
        border_color = NA,
        main = title,
        angle_col = 0, # Horizontal column labels
        gaps_col = 3, # Gap between GFP and TES samples
        silent = TRUE
    )

    # Calculate dynamic plot height based on number of genes
    n_genes <- length(available_genes)
    # Cap height at 30 inches to keep plots reasonable
    dynamic_height <- min(30, max(PLOT_HEIGHT, (n_genes * cellheight / 72) + 5)) # Convert points to inches + margins for title/legend

    # Helper function to save PNG
    save_pheatmap_png <- function(x, filename, width, height, res) {
        png(filename, width = width, height = height, units = "in", res = res)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
    }

    # Save PNG (with labels if show_rownames=TRUE)
    png_file <- paste0("plots/", filename_base, ".png")
    save_pheatmap_png(p, png_file, PLOT_WIDTH, dynamic_height, DPI)
    cat(sprintf("  Saved: %s\n", png_file))

    # Create version WITHOUT labels if original has labels
    if (show_rownames && create_nolabel_version) {
        cat("  Creating no-label version...\n")

        # Create heatmap without row names
        p_nolabel <- pheatmap(
            mat_scaled,
            color = color_palette,
            breaks = seq(-2, 2, length.out = 101),
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            show_rownames = FALSE, # NO row labels
            show_colnames = FALSE,
            annotation_col = annotation_col,
            annotation_colors = annotation_colors,
            annotation_names_col = FALSE, # Hide annotation names to avoid overlap
            annotation_legend = TRUE,
            legend = TRUE,
            fontsize = FONTSIZE_LEGEND,
            fontsize_col = FONTSIZE_COL,
            cellwidth = 50,
            cellheight = cellheight,
            border_color = NA,
            main = title,
            angle_col = 0,
            gaps_col = 3,
            silent = TRUE
        )

        # Save no-label version
        png_nolabel <- paste0("plots/", filename_base, "_nolabel.png")
        save_pheatmap_png(p_nolabel, png_nolabel, PLOT_WIDTH, dynamic_height, DPI)
        cat(sprintf("  Saved: %s\n", png_nolabel))
    }

    return(length(available_genes))
}

# ===============================================================================
# Version 1: Top 50 DEGs by adjusted p-value
# ===============================================================================

cat("\n", strrep("=", 60), "\n")
cat("VERSION 1: Top 50 DEGs\n")
cat(strrep("=", 60), "\n")

# Get significant DEGs sorted by padj - take extra to account for filtering
# We want exactly N_TOP_GENES after filtering out genes without symbols
sig_genes_all <- res %>%
    filter(!is.na(padj), padj < PADJ_CUTOFF) %>%
    arrange(padj)

# Filter to genes with valid symbols first, then take top N
sig_genes_with_symbols <- sig_genes_all %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSG", gene_symbol)) %>%
    head(N_TOP_GENES)

cat(sprintf(
    "  Selected %d genes with valid symbols from %d significant genes\n",
    nrow(sig_genes_with_symbols), nrow(sig_genes_all)
))

top_gene_ids <- sig_genes_with_symbols$gene_id

# Match to counts rownames
top_counts_ids <- rownames(counts)[gsub("\\..*", "", rownames(counts)) %in%
    gsub("\\..*", "", top_gene_ids)]

n1 <- create_publication_heatmap(
    gene_ids = top_counts_ids,
    title = sprintf("Top %d DEGs (TES vs GFP)", N_TOP_GENES),
    filename_base = "heatmap_top50_degs",
    show_rownames = TRUE,
    update_title_with_count = FALSE
)

# ===============================================================================
# Version 2: Custom gene list (TES_degs.txt)
# ===============================================================================

cat("\n", strrep("=", 60), "\n")
cat("VERSION 2: Custom Gene List\n")
cat(strrep("=", 60), "\n")

# Find Ensembl IDs for custom genes
custom_ensembl_ids <- c()
for (gene in custom_genes) {
    # Try exact match in gene_symbol
    matched <- res$gene_id[res$gene_symbol == gene]
    if (length(matched) > 0) {
        custom_ensembl_ids <- c(custom_ensembl_ids, matched[1])
    }
}

cat(sprintf(
    "Matched %d of %d custom genes to DESeq2 results\n",
    length(custom_ensembl_ids), length(custom_genes)
))

# Match to counts rownames
custom_counts_ids <- rownames(counts)[gsub("\\..*", "", rownames(counts)) %in%
    gsub("\\..*", "", custom_ensembl_ids)]

n2 <- create_publication_heatmap(
    gene_ids = custom_counts_ids,
    title = "TES Target Genes Expression",
    filename_base = "heatmap_genelist",
    show_rownames = TRUE,
    max_genes = 150
)

# ===============================================================================
# Version 3: All significant DEGs
# ===============================================================================

# cat("\n", strrep("=", 60), "\n")
# cat("VERSION 3: All Significant DEGs\n")
# cat(strrep("=", 60), "\n")

# Get all significant DEGs with FC cutoff
all_sig <- res %>%
    filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > FC_CUTOFF) %>%
    arrange(padj)

# cat(sprintf("Total significant DEGs (padj < %.2f, |log2FC| > %.1f): %d\n",
#             PADJ_CUTOFF, FC_CUTOFF, nrow(all_sig)))

# all_sig_ids <- res$gene_id[res$gene_id %in% all_sig$gene_id]

# # Match to counts rownames
# all_counts_ids <- rownames(counts)[gsub("\\..*", "", rownames(counts)) %in%
#                                     gsub("\\..*", "", all_sig_ids)]

# # For all DEGs, don't show row names if too many
# show_names <- length(all_counts_ids) <= 100

# n3 <- create_publication_heatmap(
#     gene_ids = all_counts_ids,
#     title = sprintf("All Significant DEGs (n=%d)", length(all_counts_ids)),
#     filename_base = "heatmap_all_degs",
#     show_rownames = show_names,
#     max_genes = 500
# )

# ===============================================================================
# Version 4: Top 100 DEGs with row names (pre-filtered for valid symbols)
# ===============================================================================

cat("\n", strrep("=", 60), "\n")
cat("VERSION 4: Top 100 DEGs with labels\n")
cat(strrep("=", 60), "\n")

# Get top 100 significant DEGs that have valid gene symbols
top100_with_symbols <- all_sig %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSG", gene_symbol)) %>%
    head(100)

cat(sprintf("  Selected %d genes with valid symbols\n", nrow(top100_with_symbols)))

top100_ids <- top100_with_symbols$gene_id

# Match to counts rownames
top100_counts_ids <- rownames(counts)[gsub("\\..*", "", rownames(counts)) %in%
    gsub("\\..*", "", top100_ids)]

n4 <- create_publication_heatmap(
    gene_ids = top100_counts_ids,
    title = "Top 100 Significant DEGs",
    filename_base = "heatmap_all_degs_labeled",
    show_rownames = TRUE,
    max_genes = 100,
    update_title_with_count = FALSE
)

# ===============================================================================
# Summary
# ===============================================================================

cat("\n", strrep("=", 60), "\n")
cat("=== Heatmap Generation Complete ===\n")
cat(strrep("=", 60), "\n")
cat(sprintf("\nOutput files (PNG only):\n"))
cat(sprintf("  1. heatmap_top50_degs.png - Top %d DEGs (with labels)\n", N_TOP_GENES))
cat(sprintf("     heatmap_top50_degs_nolabel.png - Top %d DEGs (no labels)\n", N_TOP_GENES))
cat(sprintf(
    "  2. heatmap_genelist.png - Custom gene list (%d genes, with labels)\n",
    ifelse(is.null(n2), 0, n2)
))
cat(sprintf("     heatmap_genelist_nolabel.png - Custom gene list (no labels)\n"))
cat(sprintf("  3. heatmap_all_degs.png - All significant DEGs\n"))
cat(sprintf("  4. heatmap_all_degs_labeled.png - Top 100 DEGs (with labels)\n"))
cat(sprintf("     heatmap_all_degs_labeled_nolabel.png - Top 100 DEGs (no labels)\n"))
cat("\nStyle: Dark blue (downregulated) - White - Intensive red (upregulated)\n")
cat("Values: Z-score normalized (row-scaled)\n")
cat("Dendrograms: None (as requested)\n")
cat("Gene filtering: Genes without valid gene symbol mapping are excluded\n")
