#!/usr/bin/env Rscript

# ===============================================================================
# TEAD-family Expression Heatmap (TEAD1, TEAD2, TEAD3, TEAD4)
# Compares Mock (= GFP control) vs TES RNA-seq samples (3 replicates each).
#
# Adapted from results/05_deseq2/plots/heatmap_publication.R:
#   - same blue->white->red palette, row z-scoring (capped +/-2), no dendrograms
#   - only the gene-selection step differs (fixed TEAD1-4 set instead of DEGs)
#   - GFP columns are relabelled "Mock" to match the manuscript terminology
#
# Style: Dark blue (low) - White - Intensive red (high), row z-score, no clustering
# ===============================================================================

suppressPackageStartupMessages({
    library(pheatmap)
    library(RColorBrewer)
    library(dplyr)
    library(org.Hs.eg.db)
    library(grDevices)
    library(grid)
})

cat("=== TEAD-family Expression Heatmap ===\n")
cat("Timestamp:", as.character(Sys.time()), "\n\n")

# ===============================================================================
# Configuration
# ===============================================================================

# Work from the DESeq2 results directory (inputs + plots/ live here)
DESEQ_DIR <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/05_deseq2"
setwd(DESEQ_DIR)

TEAD_GENES <- c("TEAD1", "TEAD2", "TEAD3", "TEAD4") # fixed row order (top -> bottom)

# Wide canvas so the long row labels (gene + log2FC/padj), all 6 columns, and the
# condition legend all fit without clipping the title or edges.
PLOT_WIDTH <- 16
PLOT_HEIGHT <- 6.5
DPI <- 300

FONTSIZE_LEGEND <- 16
FONTSIZE_ROW <- 16
FONTSIZE_COL <- 16

# Color palette: identical to heatmap_publication.R (dark blue -> white -> red)
color_palette <- colorRampPalette(c(
    "#08306B", "#08519C", "#2171B5", "#4292C6",
    "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7",
    "white",
    "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
    "#EF3B2C", "#CB181D", "#A50F15", "#67000D"
))(100)

# Condition colors reuse the project convention (GFP/Mock red, TES teal)
annotation_colors <- list(
    condition = c("Mock" = "#921100", "TES" = "#069093")
)

# ===============================================================================
# Load data
# ===============================================================================

cat("Loading data...\n")

res <- read.table("deseq2_results_TES_vs_GFP.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
cat(sprintf("  DESeq2 results: %d genes\n", nrow(res)))

counts <- read.table("normalized_counts.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
    row.names = 1
)
cat(sprintf("  Normalized counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

# ===============================================================================
# Map TEAD symbols -> Ensembl IDs present in the counts matrix
# ===============================================================================

# Prefer the gene_symbol column already in the DESeq2 table; fall back to org.Hs.eg.db
res$ensembl_clean <- gsub("\\..*", "", res$gene_id)

symbol_for_tead <- res %>%
    filter(gene_symbol %in% TEAD_GENES) %>%
    distinct(gene_symbol, .keep_all = TRUE)

# Fallback for any TEAD gene missing a symbol in the DESeq2 table
missing <- setdiff(TEAD_GENES, symbol_for_tead$gene_symbol)
if (length(missing) > 0) {
    cat(sprintf("  Symbol not in DESeq2 table for: %s -- trying org.Hs.eg.db\n",
        paste(missing, collapse = ", ")))
    ens <- mapIds(org.Hs.eg.db, keys = missing, keytype = "SYMBOL",
        column = "ENSEMBL", multiVals = "first")
    ens <- ens[!is.na(ens)]
    if (length(ens) > 0) {
        extra <- data.frame(
            gene_symbol = names(ens),
            ensembl_clean = unname(ens),
            stringsAsFactors = FALSE
        )
        symbol_for_tead <- bind_rows(
            symbol_for_tead[, c("gene_symbol", "ensembl_clean")], extra
        )
    }
}

# Match TEAD Ensembl IDs to the (versioned) rownames of the counts matrix
counts_ensembl_clean <- gsub("\\..*", "", rownames(counts))
tead_rows <- list()
for (g in TEAD_GENES) {
    ec <- symbol_for_tead$ensembl_clean[symbol_for_tead$gene_symbol == g]
    if (length(ec) == 0) {
        cat(sprintf("  WARNING: no Ensembl ID resolved for %s -- skipping\n", g))
        next
    }
    hit <- rownames(counts)[counts_ensembl_clean %in% ec]
    if (length(hit) == 0) {
        cat(sprintf("  WARNING: %s (%s) not present in counts matrix -- skipping\n", g, ec[1]))
        next
    }
    tead_rows[[g]] <- hit[1]
}

if (length(tead_rows) == 0) {
    stop("None of TEAD1-4 were found in the normalized counts matrix.")
}
cat(sprintf("  TEAD genes found: %s\n", paste(names(tead_rows), collapse = ", ")))

# ===============================================================================
# Build sample metadata + ordering (Mock first, then TES)
# ===============================================================================

sample_names <- colnames(counts)
is_ctrl <- grepl("GFP", sample_names)
sample_order <- c(sample_names[is_ctrl], sample_names[!is_ctrl])

sample_data <- data.frame(
    condition = ifelse(grepl("GFP", sample_order), "Mock", "TES"),
    row.names = sample_order,
    stringsAsFactors = FALSE
)

# Pretty column labels: "Mock 1..3", "TES 1..3"
rep_num <- gsub(".*([0-9])$", "\\1", sample_order)
col_labels <- paste0(sample_data$condition, " ", rep_num)

# ===============================================================================
# Assemble matrix (genes x samples), z-score by row
# ===============================================================================

gene_ids <- unlist(tead_rows[TEAD_GENES[TEAD_GENES %in% names(tead_rows)]])
mat <- as.matrix(counts[gene_ids, sample_order])
rownames(mat) <- names(gene_ids) # TEAD symbols

# Row labels annotated with log2FC / padj so the up/down story is explicit
lfc <- res$log2FoldChange[match(gsub("\\..*", "", gene_ids), res$ensembl_clean)]
padj <- res$padj[match(gsub("\\..*", "", gene_ids), res$ensembl_clean)]
row_labels <- sprintf("%s  (log2FC %+.2f, padj %.1e)", rownames(mat), lfc, padj)

cat("\nTEAD expression summary (TES vs Mock):\n")
print(data.frame(gene = rownames(mat), log2FC = round(lfc, 2), padj = signif(padj, 3)))

# Row z-score (same approach as the template), capped at +/-2
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2
mat_scaled[is.na(mat_scaled)] <- 0
mat_scaled[is.infinite(mat_scaled)] <- 0

# ===============================================================================
# Draw + save (labelled and no-label variants, PNG + PDF)
# ===============================================================================

dir.create("plots", showWarnings = FALSE)

draw_heatmap <- function(labels_row, show_rownames) {
    pheatmap(
        mat_scaled,
        color = color_palette,
        breaks = seq(-2, 2, length.out = 101),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = show_rownames,
        labels_row = labels_row,
        show_colnames = TRUE,
        labels_col = col_labels,
        annotation_col = sample_data,
        annotation_colors = annotation_colors,
        annotation_names_col = FALSE,
        annotation_legend = TRUE,
        legend = TRUE,
        fontsize = FONTSIZE_LEGEND,
        fontsize_row = FONTSIZE_ROW,
        fontsize_col = FONTSIZE_COL,
        cellwidth = 42,
        cellheight = 42,
        border_color = "grey80",
        main = "TEAD family expression (TES vs Mock)",
        angle_col = 45,
        gaps_col = sum(is_ctrl),
        silent = TRUE
    )
}

save_png <- function(p, filename) {
    png(filename, width = PLOT_WIDTH, height = PLOT_HEIGHT, units = "in", res = DPI)
    grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
    cat(sprintf("  Saved: %s\n", filename))
}
save_pdf <- function(p, filename) {
    pdf(filename, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off()
    cat(sprintf("  Saved: %s\n", filename))
}

cat("\nWriting plots...\n")

p_lab <- draw_heatmap(row_labels, TRUE)
save_png(p_lab, "plots/heatmap_tead_family.png")
save_pdf(p_lab, "plots/heatmap_tead_family.pdf")

# No-label variant (for figure assembly), matches existing convention
p_nolab <- draw_heatmap(rep("", nrow(mat_scaled)), FALSE)
save_png(p_nolab, "plots/heatmap_tead_family_nolabel.png")
save_pdf(p_nolab, "plots/heatmap_tead_family_nolabel.pdf")

cat("\n=== Done ===\n")
cat("Style: row z-score, dark blue (low) -> white -> red (high), no clustering\n")
cat("Columns: Mock (GFP) 1-3 | TES 1-3\n")
