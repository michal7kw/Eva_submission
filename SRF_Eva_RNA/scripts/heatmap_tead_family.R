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

# Row labels are now just gene symbols (log2FC/padj moved to a side annotation),
# so the canvas no longer needs to be as wide as when labels were long strings.
# Height must fit the three stacked legends (condition + Significance + log2FC)
# plus the z-score colorbar without clipping the bottom swatch.
PLOT_WIDTH <- 11
PLOT_HEIGHT <- 6.8
DPI <- 300

FONTSIZE_LEGEND <- 14
FONTSIZE_ROW <- 17   # gene symbols, slightly larger now that they stand alone
FONTSIZE_COL <- 12   # small enough that horizontal "Mock 1".."TES 3" don't collide

# Color palette: identical to heatmap_publication.R (dark blue -> white -> red)
color_palette <- colorRampPalette(c(
    "#08306B", "#08519C", "#2171B5", "#4292C6",
    "#6BAED6", "#9ECAE1", "#C6DBEF", "#DEEBF7",
    "white",
    "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
    "#EF3B2C", "#CB181D", "#A50F15", "#67000D"
))(100)

# Annotation colors. All three annotation legends use *discrete* swatches so they
# read as one tidy stacked block, visually separate from the main red<->blue
# z-score colorbar:
#   - condition  : project convention (Mock red, TES teal)
#   - log2FC     : binned, purple (up in TES) <-> green (down) -- PRGn endpoints,
#                  deliberately NOT red/blue so it isn't confused with the z-score
#   - Significance: grey -> black ramp by padj
# Levels are ordered up->down / strong->weak so each legend reads top-to-bottom.
LOG2FC_LEVELS <- c(">= +1", "0 to +1", "-1 to 0", "<= -1")
SIG_LEVELS <- c("ns", "*", "**", "***")
annotation_colors <- list(
    condition = c("Mock" = "#921100", "TES" = "#069093"),
    log2FC = c(
        ">= +1"   = "#762A83",  # dark purple  (strong up)
        "0 to +1" = "#C2A5CF",  # light purple (up)
        "-1 to 0" = "#A6DBA0",  # light green  (mild down)
        "<= -1"   = "#1B7837"   # dark green   (strong down)
    ),
    Significance = c(
        "ns"  = "#E8E8E8",
        "*"   = "#BDBDBD",
        "**"  = "#737373",
        "***" = "#252525"
    )
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

# Pull log2FC / padj for each TEAD gene; these now ride in a side annotation bar
# (annotation_row) instead of being crammed into the row label text.
lfc <- res$log2FoldChange[match(gsub("\\..*", "", gene_ids), res$ensembl_clean)]
padj <- res$padj[match(gsub("\\..*", "", gene_ids), res$ensembl_clean)]

# padj -> significance stars (ns / * / ** / ***); NA padj treated as ns
sig_stars <- cut(padj,
    breaks = c(-Inf, 1e-3, 1e-2, 5e-2, Inf),
    labels = c("***", "**", "*", "ns")
)
sig_stars <- factor(as.character(sig_stars), levels = SIG_LEVELS)
sig_stars[is.na(sig_stars)] <- "ns"

# log2FC -> discrete bins so its legend is a tidy swatch block (not a gradient)
log2fc_bin <- cut(lfc,
    breaks = c(-Inf, -1, 0, 1, Inf),
    labels = c("<= -1", "-1 to 0", "0 to +1", ">= +1")
)
log2fc_bin <- factor(as.character(log2fc_bin), levels = LOG2FC_LEVELS)

# annotation_row is matched to the heatmap by rownames (= gene symbols).
# Order: log2FC (effect) then Significance (confidence).
annotation_row <- data.frame(
    log2FC = log2fc_bin,
    Significance = sig_stars,
    row.names = rownames(mat),
    stringsAsFactors = FALSE
)

cat("\nTEAD expression summary (TES vs Mock):\n")
print(data.frame(gene = rownames(mat), log2FC = round(lfc, 2),
    padj = signif(padj, 3), sig = as.character(sig_stars)))

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

# show_rownames + annotation_row are toggled together: the labelled variant shows
# gene symbols and the log2FC/Significance bars; the no-label variant is a clean
# cells-only heatmap (just the column condition bar) for figure assembly.
draw_heatmap <- function(show_rownames, annotation_row_arg) {
    pheatmap(
        mat_scaled,
        color = color_palette,
        breaks = seq(-2, 2, length.out = 101),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = show_rownames,   # rownames(mat_scaled) = gene symbols
        show_colnames = TRUE,
        labels_col = col_labels,
        annotation_col = sample_data,
        annotation_row = annotation_row_arg,
        annotation_colors = annotation_colors,
        annotation_names_col = FALSE,
        annotation_names_row = FALSE,   # names shown via legends; avoids colliding with col labels
        annotation_legend = TRUE,
        legend = TRUE,
        legend_breaks = c(-2, -1, 0, 1, 2),
        legend_labels = c("-2 (low)", "-1", "0", "+1", "+2 (high)"),
        fontsize = FONTSIZE_LEGEND,
        fontsize_row = FONTSIZE_ROW,
        fontsize_col = FONTSIZE_COL,
        cellwidth = 52,
        cellheight = 44,
        border_color = "white",          # thin clean gridlines
        main = "TEAD family expression  (row z-score, TES vs Mock)",
        angle_col = 0,                   # short labels read better horizontal
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

p_lab <- draw_heatmap(TRUE, annotation_row)
save_png(p_lab, "plots/heatmap_tead_family.png")
save_pdf(p_lab, "plots/heatmap_tead_family.pdf")

# No-label variant (for figure assembly): clean cells, no row text/annotation
p_nolab <- draw_heatmap(FALSE, NA)
save_png(p_nolab, "plots/heatmap_tead_family_nolabel.png")
save_pdf(p_nolab, "plots/heatmap_tead_family_nolabel.pdf")

cat("\n=== Done ===\n")
cat("Style: row z-score, dark blue (low) -> white -> red (high), no clustering\n")
cat("Columns: Mock (GFP) 1-3 | TES 1-3\n")
