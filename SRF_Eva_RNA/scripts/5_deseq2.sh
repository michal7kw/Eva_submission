#!/bin/bash

#===============================================================================
# SCRIPT: 5_deseq2.sh
# PURPOSE: Differential Expression Analysis using DESeq2
#
# DESCRIPTION:
# This script performs differential expression analysis using DESeq2 in R.
# It is the fifth step in the RNA-seq analysis pipeline and identifies
# significantly differentially expressed genes between conditions (TES vs GFP).
#
# KEY OPERATIONS:
# 1. Activates R environment with DESeq2 and related packages
# 2. Loads count matrix and sample metadata
# 3. Performs DESeq2 differential expression analysis
# 4. Generates quality control plots (PCA, heatmaps)
# 5. Creates volcano plots and MA plots
# 6. Exports results tables and visualizations
# 7. Performs gene set enrichment analysis (optional)
#
# METHODOLOGY:
# - DESeq2 for robust differential expression analysis
# - Variance stabilizing transformation for visualization
# - Multiple testing correction (Benjamini-Hochberg)
# - Principal component analysis for sample clustering
# - Comprehensive visualization and reporting
#
# EXPECTED INPUTS:
# - Count matrix: ${QUANTIFIED_DIR}/count_matrix.txt
# - Sample metadata: ${QUANTIFIED_DIR}/sample_metadata.txt
#
# EXPECTED OUTPUTS:
# - DE results: ${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt
# - Normalized counts: ${OUTPUT_DIR}/normalized_counts.txt
# - QC plots: ${OUTPUT_DIR}/plots/*.png
# - Summary report: ${OUTPUT_DIR}/differential_expression_summary.txt
# - SLURM logs: ./logs/5_deseq2.out and ./logs/5_deseq2.err
#
# DEPENDENCIES:
# - R with DESeq2, ggplot2, pheatmap, EnhancedVolcano
# - SLURM workload manager
#
# USAGE:
# sbatch 5_deseq2.sh
#===============================================================================

#SBATCH --job-name=5_deseq2
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5_deseq2.err"
#SBATCH --output="./logs/5_deseq2.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
QUANTIFIED_DIR="${BASE_DIR}/results/04_quantified"
OUTPUT_DIR="${BASE_DIR}/results/05_deseq2"

# Create plots directory
mkdir -p ${OUTPUT_DIR}/plots

echo "=== Starting differential expression analysis ==="
echo "Timestamp: $(date)"
echo "Input directory: ${QUANTIFIED_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create comprehensive R script for DESeq2 analysis
cat > ${OUTPUT_DIR}/deseq2_analysis.R << 'EOF'
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

# Filter low count genes (optional but recommended)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

cat(sprintf("After filtering: %d genes retained\n", nrow(dds)))

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results for TES vs GFP comparison
res <- results(dds, contrast = c("condition", "TES", "GFP"))
res <- res[order(res$padj), ]

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

# Reorder columns to include gene_symbol
res_df <- res_df[, c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

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

# MA plot
png("plots/ma_plot.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res, main = "MA Plot (TES vs GFP)", ylim = c(-5, 5))
dev.off()

# Volcano plot using EnhancedVolcano
if ("EnhancedVolcano" %in% rownames(installed.packages())) {
    # Get top 20 gene symbols for labeling
    top_genes_idx <- head(order(res_df$padj), 20)
    top_gene_symbols <- res_df$gene_symbol[top_genes_idx]

    p_volcano <- EnhancedVolcano(res_df,
                                lab = res_df$gene_symbol,
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'Volcano Plot: TES vs GFP',
                                subtitle = 'Differential Expression Analysis',
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

EOF

echo "=== Running DESeq2 analysis in R ==="
echo "This may take several minutes..."

# Run the R script
Rscript ${OUTPUT_DIR}/deseq2_analysis.R

echo "=== Checking output files ==="

# Check if main results file was created
if [[ -f "${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt" ]]; then
    echo "SUCCESS: DESeq2 results file created"
    RESULTS_LINES=$(wc -l < "${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt")
    echo "  File: ${OUTPUT_DIR}/deseq2_results_TES_vs_GFP.txt"
    echo "  Lines: ${RESULTS_LINES}"
else
    echo "ERROR: DESeq2 results file not created!"
    exit 1
fi

# Check significant genes file
if [[ -f "${OUTPUT_DIR}/significant_genes_TES_vs_GFP.txt" ]]; then
    SIG_LINES=$(tail -n +2 "${OUTPUT_DIR}/significant_genes_TES_vs_GFP.txt" | wc -l)
    echo "SUCCESS: Significant genes file created (${SIG_LINES} genes)"
else
    echo "WARNING: Significant genes file not created!"
fi

# Check plots directory
if [[ -d "${OUTPUT_DIR}/plots" ]]; then
    PLOT_COUNT=$(find "${OUTPUT_DIR}/plots" -name "*.png" | wc -l)
    echo "SUCCESS: Plots directory created with ${PLOT_COUNT} plots"
    find "${OUTPUT_DIR}/plots" -name "*.png" | sed 's/^/  /'
else
    echo "WARNING: Plots directory not created!"
fi

# Display summary if available
if [[ -f "${OUTPUT_DIR}/differential_expression_summary.txt" ]]; then
    echo ""
    echo "=== Analysis Summary ==="
    cat "${OUTPUT_DIR}/differential_expression_summary.txt"
fi

echo ""
echo "=== Differential expression analysis complete ==="
echo "Timestamp: $(date)"
echo ""
echo "Analysis complete! Check the results in: ${OUTPUT_DIR}"