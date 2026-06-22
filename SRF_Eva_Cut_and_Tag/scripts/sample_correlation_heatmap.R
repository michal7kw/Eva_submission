#!/usr/bin/env Rscript
#===============================================================================
# SCRIPT: sample_correlation_heatmap.R
# PURPOSE: Generate sample correlation/distance heatmap for Cut&Tag data
#
# DESCRIPTION:
# Creates a publication-ready sample correlation heatmap similar to RNA-seq
# Euclidean distance plots, showing sample clustering and similarity patterns.
# Uses DiffBind count matrix to compute sample distances.
#
# OUTPUT:
# - Sample correlation heatmap with hierarchical clustering
# - Annotation bars for condition and replicate
#
# USAGE:
# Rscript scripts/sample_correlation_heatmap.R
#===============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(DiffBind)
    library(pheatmap)
    library(RColorBrewer)
})

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Output directory
output_dir <- "results/07_analysis_narrow"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("Generating Sample Correlation Heatmap\n")
cat("========================================\n")

#-------------------------------------------------------------------------------
# Create DiffBind object and count reads
#-------------------------------------------------------------------------------

# Create sample sheet for DiffBind
samples <- data.frame(
    SampleID = c(
        "TES-1", "TES-2", "TES-3",
        "TEAD1-1", "TEAD1-2", "TEAD1-3"
    ),
    Tissue = rep("SNB19", 6),
    Factor = c(rep("TES", 3), rep("TEAD1", 3)),
    Condition = c(rep("TES", 3), rep("TEAD1", 3)),
    Replicate = rep(1:3, 2),
    bamReads = paste0(
        "results/04_filtered/",
        c("TES-1", "TES-2", "TES-3", "TEAD1-1", "TEAD1-2", "TEAD1-3"),
        "_filtered.bam"
    ),
    ControlID = c(rep("IggMs", 3), rep("IggRb", 3)),
    bamControl = c(
        rep("results/04_filtered/IggMs_filtered.bam", 3),
        rep("results/04_filtered/IggRb_filtered.bam", 3)
    ),
    Peaks = paste0(
        "results/05_peaks_narrow/",
        c("TES-1", "TES-2", "TES-3", "TEAD1-1", "TEAD1-2", "TEAD1-3"),
        "_peaks.narrowPeak"
    ),
    PeakCaller = rep("macs", 6)
)

# Check if all files exist
cat("\nChecking input files...\n")
bam_exist <- file.exists(samples$bamReads)
peak_exist <- file.exists(samples$Peaks)

cat("BAM files found:", sum(bam_exist), "/", length(bam_exist), "\n")
cat("Peak files found:", sum(peak_exist), "/", length(peak_exist), "\n")

if (!all(bam_exist) || !all(peak_exist)) {
    cat("\nMissing files:\n")
    if (!all(bam_exist)) {
        cat("BAM:", samples$bamReads[!bam_exist], "\n")
    }
    if (!all(peak_exist)) {
        cat("Peaks:", samples$Peaks[!peak_exist], "\n")
    }
    stop("Required files not found!")
}

# Try to load existing DiffBind object or create new one
dba_file <- file.path(output_dir, "diffbind_correlation_object.RData")

if (file.exists(dba_file)) {
    cat("\nLoading existing DiffBind object...\n")
    load(dba_file)
} else {
    cat("\nCreating DiffBind object...\n")

    # Save sample sheet
    samplesheet_file <- "config/diffbind_samplesheet_correlation.csv"
    write.csv(samples, samplesheet_file, row.names = FALSE, quote = FALSE)

    # Load data
    dba_obj <- dba(sampleSheet = samplesheet_file)

    # Count reads with simple parameters (avoid summits which can cause issues)
    cat("Counting reads in peaks (this may take a while)...\n")
    dba_obj <- dba.count(dba_obj, bParallel = TRUE)

    # Try to normalize if the function exists and works
    tryCatch({
        dba_obj <- dba.normalize(dba_obj)
        cat("Normalization applied.\n")
    }, error = function(e) {
        cat("Note: Skipping normalization (not critical for correlation analysis)\n")
    })

    # Save for future use
    save(dba_obj, file = dba_file)
    cat("DiffBind object saved to:", dba_file, "\n")
}

#-------------------------------------------------------------------------------
# Extract count matrix and compute correlations
#-------------------------------------------------------------------------------

cat("\nExtracting count matrix...\n")

# Get correlation matrix directly from DiffBind (most reliable method)
# This uses DiffBind's internal correlation calculation
cor_from_dba <- dba.plotHeatmap(dba_obj, correlations = TRUE, plot = FALSE)

# Sample names for display
sample_names <- c("TES1", "TES2", "TES3", "TEAD1_1", "TEAD1_2", "TEAD1_3")

# Create annotation data frame for samples
annotation_df <- data.frame(
    condition = c(rep("TES", 3), rep("TEAD1", 3)),
    samples = as.character(rep(1:3, 2)),
    row.names = sample_names
)

# Define colors matching the example figure style exactly
# TES = teal (#008B8B), TEAD1 = brown (#8B4513) - matching the fig1.png style
annotation_colors <- list(
    condition = c(TES = "#008B8B", TEAD1 = "#8B4513"),  # Teal and brown
    samples = c("1" = "#4D4D4D", "2" = "#B0B0B0", "3" = "#E34234")  # Dark gray, light gray, red
)

# Reorder annotation_df to put condition first (like in the reference figure)
annotation_df <- annotation_df[, c("condition", "samples")]

#-------------------------------------------------------------------------------
# Generate Euclidean Distance Heatmap
#-------------------------------------------------------------------------------

cat("\nGenerating distance heatmap...\n")

# Convert correlation to distance (1 - correlation gives a distance-like metric)
# Or use sqrt(2*(1-cor)) for proper Euclidean distance from correlation
dist_from_cor <- sqrt(2 * (1 - cor_from_dba))
rownames(dist_from_cor) <- sample_names
colnames(dist_from_cor) <- sample_names

# Create distance object for clustering
dist_obj <- as.dist(dist_from_cor)

# PDF version - Distance heatmap (matching fig1.png style)
pdf(file.path(output_dir, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(
    dist_from_cor,
    clustering_distance_rows = dist_obj,
    clustering_distance_cols = dist_obj,
    clustering_method = "complete",
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
    main = "",
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "white",
    cellwidth = 35,
    cellheight = 35,
    display_numbers = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 50,
    treeheight_col = 50
)
dev.off()

# PNG version
png(file.path(output_dir, "sample_distance_heatmap.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(
    dist_from_cor,
    clustering_distance_rows = dist_obj,
    clustering_distance_cols = dist_obj,
    clustering_method = "complete",
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
    main = "",
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "white",
    cellwidth = 35,
    cellheight = 35,
    display_numbers = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 50,
    treeheight_col = 50
)
dev.off()

#-------------------------------------------------------------------------------
# Generate Correlation Heatmap (alternative view)
#-------------------------------------------------------------------------------

cat("Generating correlation heatmap...\n")

# Use the correlation matrix from DiffBind
cor_matrix <- cor_from_dba
rownames(cor_matrix) <- sample_names
colnames(cor_matrix) <- sample_names

# PDF version
pdf(file.path(output_dir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(
    cor_matrix,
    clustering_method = "complete",
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                               "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D",
                               "#B2182B", "#67001F"))(100),
    main = "",
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "white",
    cellwidth = 35,
    cellheight = 35,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 50,
    treeheight_col = 50
)
dev.off()

# PNG version
png(file.path(output_dir, "sample_correlation_heatmap.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(
    cor_matrix,
    clustering_method = "complete",
    annotation_row = annotation_df,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE",
                               "#D1E5F0", "#FDDBC7", "#F4A582", "#D6604D",
                               "#B2182B", "#67001F"))(100),
    main = "",
    fontsize = 11,
    fontsize_row = 10,
    fontsize_col = 10,
    border_color = "white",
    cellwidth = 35,
    cellheight = 35,
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_names_row = TRUE,
    annotation_names_col = TRUE,
    annotation_legend = TRUE,
    legend = TRUE,
    treeheight_row = 50,
    treeheight_col = 50
)
dev.off()

#-------------------------------------------------------------------------------
# Print summary statistics
#-------------------------------------------------------------------------------

cat("\n========================================\n")
cat("Analysis Complete!\n")
cat("========================================\n")
cat("\nOutput files:\n")
cat("  - sample_distance_heatmap.pdf/png (Euclidean-like distance from correlation)\n")
cat("  - sample_correlation_heatmap.pdf/png (Pearson correlation)\n")
cat(paste("\nOutput directory:", output_dir, "\n"))

cat("\nCorrelation matrix:\n")
print(round(cor_matrix, 3))

cat("\nDistance matrix:\n")
print(round(dist_from_cor, 3))

# Print within-group vs between-group statistics
tes_idx <- 1:3
tead1_idx <- 4:6

within_tes_cor <- mean(cor_matrix[tes_idx, tes_idx][upper.tri(cor_matrix[tes_idx, tes_idx])])
within_tead1_cor <- mean(cor_matrix[tead1_idx, tead1_idx][upper.tri(cor_matrix[tead1_idx, tead1_idx])])
between_cor <- mean(cor_matrix[tes_idx, tead1_idx])

cat("\nWithin-group mean correlations:\n")
cat("  TES replicates:", round(within_tes_cor, 3), "\n")
cat("  TEAD1 replicates:", round(within_tead1_cor, 3), "\n")
cat("  Between TES-TEAD1:", round(between_cor, 3), "\n")

cat("\n========================================\n")
