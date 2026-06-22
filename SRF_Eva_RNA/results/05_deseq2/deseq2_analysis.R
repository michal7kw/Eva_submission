#!/usr/bin/env Rscript

# Load required libraries
cat("Loading required R packages...\n")

# Load required libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tibble)
library(ggrepel)
library(org.Hs.eg.db)

# Try to load EnhancedVolcano (optional)
if (!require(EnhancedVolcano, quietly = TRUE)) {
    cat("Warning: EnhancedVolcano not available, skipping volcano plot\n")
}

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/05_deseq2")

cat("=== DESeq2 Differential Expression Analysis ===\n")
cat("Timestamp:", as.character(Sys.time()), "\n\n")

# ============================================================================
# 1. Load data
# ============================================================================
cat("1. Loading data...\n")

# Load count matrix
count_data <- read.table("../04_quantified/count_matrix.txt",
                        header = TRUE, row.names = 1, sep = "\t")

# Load sample metadata
sample_data <- read.table("../04_quantified/sample_metadata.txt",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure sample order matches
sample_data <- sample_data[match(colnames(count_data), sample_data$sample), ]
rownames(sample_data) <- sample_data$sample

cat(sprintf("Loaded count data: %d genes x %d samples\n", nrow(count_data), ncol(count_data)))
cat("Sample metadata:\n")
print(sample_data)

# Convert condition to factor with appropriate levels
sample_data$condition <- factor(sample_data$condition, levels = c("GFP", "TES"))

# ============================================================================
# 2. Create DESeq2 object and run analysis
# ============================================================================
cat("\n2. Creating DESeq2 object and running analysis...\n")

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                             colData = sample_data,
                             design = ~ condition)

# Filter low count genes - require at least 10 counts in at least 3 samples
# (minimum group size), rather than total counts >= 10 across all samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

cat(sprintf("After filtering: %d genes retained\n", nrow(dds)))

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results for TES vs GFP comparison
# Set alpha = 0.05 to match our significance cutoff (optimizes independent filtering)
res <- results(dds, contrast = c("condition", "TES", "GFP"), alpha = 0.05)
res <- res[order(res$padj), ]

# Apply apeglm LFC shrinkage for more accurate fold change estimates
# Shrunken LFC reduces noise from low-count genes and is recommended for:
# - GSEA ranking, volcano plots, heatmaps, integration analyses
# Note: significance testing (padj) uses the unshrunken results above
cat("Applying apeglm LFC shrinkage...\n")
cat("Available coefficients:", paste(resultsNames(dds), collapse = ", "), "\n")
shrink_coef <- "condition_TES_vs_GFP"
if (!shrink_coef %in% resultsNames(dds)) {
    stop("Coefficient '", shrink_coef, "' not found. Available: ",
         paste(resultsNames(dds), collapse = ", "))
}
res_shrunk <- lfcShrink(dds, coef = shrink_coef, type = "apeglm")
res_shrunk <- res_shrunk[order(res_shrunk$padj), ]

cat("DESeq2 analysis completed.\n")

# ============================================================================
# 3. Generate summary statistics
# ============================================================================
cat("\n3. Generating summary statistics...\n")

# Summary of results
summary(res)

# Count significant genes
sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
total_sig <- sum(res$padj < 0.05, na.rm = TRUE)

cat(sprintf("\nDifferential Expression Summary (TES vs GFP):\n"))
cat(sprintf("Total genes tested: %d\n", sum(!is.na(res$padj))))
cat(sprintf("Significantly upregulated in TES: %d\n", sig_up))
cat(sprintf("Significantly downregulated in TES: %d\n", sig_down))
cat(sprintf("Total significantly changed: %d\n", total_sig))

# Save summary to file
summary_text <- capture.output({
    cat("Differential Expression Analysis Summary\n")
    cat("========================================\n")
    cat("Analysis date:", as.character(Sys.time()), "\n")
    cat("Comparison: TES vs GFP\n")
    cat("Significance threshold: padj < 0.05\n\n")
    cat(sprintf("Total genes in count matrix: %d\n", nrow(count_data)))
    cat(sprintf("Genes after filtering (counts >= 10): %d\n", nrow(dds)))
    cat(sprintf("Genes with valid p-values: %d\n", sum(!is.na(res$padj))))
    cat(sprintf("Significantly upregulated in TES: %d\n", sig_up))
    cat(sprintf("Significantly downregulated in TES: %d\n", sig_down))
    cat(sprintf("Total significantly changed genes: %d\n", total_sig))
})

writeLines(summary_text, "differential_expression_summary.txt")

# ============================================================================
# 4. Save results
# ============================================================================
cat("\n4. Saving results...\n")

# Convert results to data frame and add gene names
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Add shrunken LFC from apeglm shrinkage
res_shrunk_df <- as.data.frame(res_shrunk)
res_df$log2FoldChange_shrunk <- res_shrunk_df[rownames(res_df), "log2FoldChange"]
res_df$lfcSE_shrunk <- res_shrunk_df[rownames(res_df), "lfcSE"]

# Map Ensembl IDs to gene symbols
cat("Mapping Ensembl IDs to gene symbols...\n")
# Remove version numbers from Ensembl IDs for mapping
ensembl_ids_clean <- gsub("\\..*", "", res_df$gene_id)

# Get gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                      keys = ensembl_ids_clean,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Add gene symbols to results, use gene_id if symbol not found
res_df$gene_symbol <- ifelse(!is.na(gene_symbols), gene_symbols, res_df$gene_id)

# Reorder columns: include both unshrunken (for significance) and shrunken (for effect size) LFC
res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean",
                      "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
                      "log2FoldChange_shrunk", "lfcSE_shrunk")]

# Save results
write.table(res_df, "deseq2_results_TES_vs_GFP.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

# Save significant genes only
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
write.table(sig_genes, "significant_genes_TES_vs_GFP.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(normalized_counts), normalized_counts),
           "normalized_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================================
# 5. Generate visualizations
# ============================================================================
cat("\n5. Generating visualizations...\n")

# Dispersion plot - standard DESeq2 diagnostic
png("plots/dispersion_plot.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()

# Variance stabilizing transformation for visualization
vsd <- vst(dds, blind = FALSE)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = name), size = 3) +
    labs(title = "PCA Plot",
         x = paste0("PC1: ", percent_var[1], "% variance"),
         y = paste0("PC2: ", percent_var[2], "% variance")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("plots/pca_plot.png", p_pca, width = 8, height = 6, dpi = 300)

# Sample distance heatmap
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png("plots/sample_distance_heatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colors,
         main = "Sample Distance Heatmap")
dev.off()

# MA plot (shrunken LFC)
png("plots/ma_plot.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res_shrunk, main = "MA Plot (TES vs GFP, apeglm shrinkage)", ylim = c(-5, 5))
dev.off()

# Volcano plot using EnhancedVolcano (uses shrunken LFC for cleaner visualization)
if ("EnhancedVolcano" %in% rownames(installed.packages())) {
    # Get top 20 gene symbols for labeling
    top_genes_idx <- head(order(res_df$padj), 20)
    top_gene_symbols <- res_df$gene_symbol[top_genes_idx]

    p_volcano <- EnhancedVolcano(res_df,
                                lab = res_df$gene_symbol,
                                x = 'log2FoldChange_shrunk',
                                y = 'padj',
                                title = 'Volcano Plot: TES vs GFP',
                                subtitle = 'Differential Expression (apeglm shrunken LFC)',
                                pCutoff = 0.05,
                                FCcutoff = 1.0,
                                pointSize = 2.0,
                                labSize = 3.0,
                                selectLab = top_gene_symbols,
                                drawConnectors = TRUE,
                                widthConnectors = 0.5)

    ggsave("plots/volcano_plot.png", p_volcano, width = 10, height = 8, dpi = 300)
}

# Heatmap of top differentially expressed genes
top_genes <- head(rownames(res[order(res$padj), ]), 50)
top_counts <- assay(vsd)[top_genes, ]

# Map row names to gene symbols for the heatmap
top_genes_clean <- gsub("\\..*", "", top_genes)
top_gene_symbols <- mapIds(org.Hs.eg.db,
                           keys = top_genes_clean,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

# Use gene symbols as row names, fall back to Ensembl ID if not found
rownames(top_counts) <- ifelse(!is.na(top_gene_symbols), top_gene_symbols, top_genes)

png("plots/top_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(top_counts,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         annotation_col = sample_data["condition"],
         scale = "row",
         main = "Top 50 Differentially Expressed Genes")
dev.off()

# Count plot of significant genes
count_data_plot <- data.frame(
    Direction = c("Upregulated in TES", "Downregulated in TES"),
    Count = c(sig_up, sig_down)
)

p_counts <- ggplot(count_data_plot, aes(x = Direction, y = Count, fill = Direction)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5) +
    labs(title = "Differentially Expressed Genes",
         x = "", y = "Number of Genes") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("red", "blue"))

ggsave("plots/de_gene_counts.png", p_counts, width = 8, height = 6, dpi = 300)

cat("\nAnalysis completed successfully!\n")
cat("Results saved to:\n")
cat("  - deseq2_results_TES_vs_GFP.txt\n")
cat("  - significant_genes_TES_vs_GFP.txt\n")
cat("  - normalized_counts.txt\n")
cat("  - differential_expression_summary.txt\n")
cat("  - plots/ directory with visualizations\n")

