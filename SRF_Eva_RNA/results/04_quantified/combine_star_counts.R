#!/usr/bin/env Rscript

# Load required libraries
if (!require("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "https://cran.r-project.org/")
    library(tidyverse)
}

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/04_quantified")

# Read sample names
samples <- readLines("../../config/samples.txt")
aligned_dir <- "../03_aligned"

# Initialize count matrix
count_matrix <- NULL
gene_info <- NULL

# Read the strandedness column determined by the shell script
# 2 = unstranded, 3 = sense (forward-stranded), 4 = antisense (reverse-stranded)
strand_col <- as.integer(readLines("strand_column.txt"))
col_names <- c("gene_id", "unstranded", "sense", "antisense")
strand_label <- col_names[strand_col]
cat(sprintf("Using STAR count column %d (%s) based on strandedness check\n", strand_col, strand_label))

cat("Processing STAR count files...\n")

for (i in seq_along(samples)) {
    sample <- samples[i]
    count_file <- file.path(aligned_dir, paste0(sample, "_counts.tab"))

    cat(sprintf("Processing %s (%d/%d)...\n", sample, i, length(samples)))

    # Read STAR count file
    # Columns: gene_id, unstranded, sense(++,--), antisense(+-,-+)
    counts <- read.table(count_file, header = FALSE, stringsAsFactors = FALSE,
                        col.names = col_names)

    # Remove summary lines at the end
    counts <- counts[!grepl("^N_", counts$gene_id), ]

    # For first sample, initialize matrix with gene info
    if (is.null(count_matrix)) {
        gene_info <- counts[, "gene_id", drop = FALSE]
        count_matrix <- data.frame(gene_id = counts$gene_id)
    }

    # Add sample counts using the strand-appropriate column
    count_matrix[[sample]] <- counts[[strand_label]]
}

# Write count matrix
cat("Writing count matrix...\n")
write.table(count_matrix, "count_matrix.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("Count matrix saved to count_matrix.txt\n")
cat(sprintf("Dimensions: %d genes x %d samples\n", nrow(count_matrix), ncol(count_matrix)-1))

# Basic statistics
total_counts <- colSums(count_matrix[, -1])
cat("\nTotal counts per sample:\n")
for (i in 1:length(total_counts)) {
    cat(sprintf("  %s: %s\n", names(total_counts)[i], format(total_counts[i], big.mark = ",")))
}

# Save summary statistics
summary_stats <- data.frame(
    sample = names(total_counts),
    total_counts = total_counts,
    genes_detected = colSums(count_matrix[, -1] > 0)
)

write.table(summary_stats, "quantification_summary.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nQuantification complete!\n")
