#!/usr/bin/env Rscript
# Improved DiffBind Analysis with Enhanced Statistical Rigor and Biological Interpretation
# Focus: Narrow peaks only with proper statistical workflow

library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(GenomeInfoDb)
library(dplyr)
library(yaml)

# Load configuration
config <- yaml.load_file("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/config/pipeline_config.yaml")
setwd(config$directories$base)

# Create comprehensive output directory
output_dir <- file.path(config$directories$results, "diffbind_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to create sample sheet with error checking
create_sample_sheet <- function(config) {
  samples_list <- list()

  for (group_name in names(config$samples$groups)) {
    group_info <- config$samples$groups[[group_name]]
    for (i in seq_along(group_info$samples)) {
      sample_id <- group_info$samples[i]
      samples_list[[length(samples_list) + 1]] <- data.frame(
        SampleID = sample_id,
        Tissue = "SNB19",
        Factor = group_name,
        Condition = group_name,
        Replicate = i,
        bamReads = file.path(config$directories$results, "04_filtered", paste0(sample_id, "_filtered.bam")),
        ControlID = group_info$control,
        bamControl = file.path(config$directories$results, "04_filtered", paste0(group_info$control, "_filtered.bam")),
        Peaks = file.path(config$directories$results, "05_peaks", paste0(sample_id, "_peaks.narrowPeak")),
        PeakCaller = "macs",
        stringsAsFactors = FALSE
      )
    }
  }

  return(do.call(rbind, samples_list))
}

# Create and validate sample sheet
samples <- create_sample_sheet(config)
write.csv(samples, file.path(config$directories$config, "diffbind_samplesheet_improved.csv"),
  row.names = FALSE, quote = FALSE
)

# Validate input files exist
cat("Validating input files...\n")
missing_files <- c()
for (i in 1:nrow(samples)) {
  files_to_check <- c(samples$bamReads[i], samples$bamControl[i], samples$Peaks[i])
  for (file in files_to_check) {
    if (!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }
}

if (length(missing_files) > 0) {
  cat("WARNING: Missing files detected:\n")
  cat(paste(missing_files, collapse = "\n"))
  cat("\nPipeline may fail. Please check file paths.\n")
}

# Load DiffBind object with error handling
cat("Loading DiffBind data...\n")
tryCatch(
  {
    dba_obj <- dba(sampleSheet = file.path(config$directories$config, "diffbind_samplesheet_improved.csv"))

    # Initial exploration
    cat(paste("Loaded", nrow(dba_obj$samples), "samples\n"))
    cat("Sample overview:\n")
    print(dba_obj)

    # Generate initial correlation plot before counting
    pdf(file.path(output_dir, "01_initial_correlation.pdf"), width = 10, height = 8)
    dba.plotHeatmap(dba_obj, correlations = TRUE, main = "Initial Peak Overlap Correlations")
    dev.off()
  },
  error = function(e) {
    cat("Error loading DiffBind data:", e$message, "\n")
    quit(status = 1)
  }
)

# Count reads with quality control
cat("Counting reads in peaks...\n")
dba_obj <- dba.count(dba_obj, bParallel = TRUE, score = DBA_SCORE_TMM_MINUS_FULL)

# Quality control: Check read counts
count_summary <- dba.show(dba_obj, bCounts = TRUE)
write.csv(count_summary, file.path(output_dir, "read_count_summary.csv"))

# Flag samples with unusually low counts
median_counts <- median(count_summary$Reads)
low_count_samples <- count_summary[count_summary$Reads < median_counts * 0.5, ]
if (nrow(low_count_samples) > 0) {
  cat("WARNING: Samples with low read counts:\n")
  print(low_count_samples[, c("ID", "Reads")])
}

# Normalize data
cat("Normalizing data...\n")
dba_obj <- dba.normalize(dba_obj, method = DBA_DESEQ2)

# Enhanced correlation analysis after normalization
pdf(file.path(output_dir, "02_post_normalization_correlation.pdf"), width = 12, height = 10)
dba.plotHeatmap(dba_obj, correlations = TRUE, main = "Post-Normalization Correlations")
dev.off()

# PCA analysis with enhanced visualization
pdf(file.path(output_dir, "03_enhanced_PCA.pdf"), width = 12, height = 10)
dba.plotPCA(dba_obj, DBA_CONDITION,
  label = DBA_ID,
  main = "PCA Analysis - All Samples"
)
dev.off()

# Set up contrasts based on configuration
cat("Setting up statistical contrasts...\n")
for (i in seq_along(config$contrasts)) {
  contrast_info <- config$contrasts[[i]]

  tryCatch(
    {
      dba_obj <- dba.contrast(dba_obj,
        group1 = dba_obj$masks[[contrast_info$group1]],
        group2 = dba_obj$masks[[contrast_info$group2]],
        name1 = contrast_info$group1,
        name2 = contrast_info$group2
      )
      cat(paste("Added contrast:", contrast_info$name, "\n"))
    },
    error = function(e) {
      cat(paste("Error adding contrast", contrast_info$name, ":", e$message, "\n"))
    }
  )
}

# Perform differential analysis with enhanced parameters
cat("Performing differential binding analysis...\n")
dba_obj <- dba.analyze(dba_obj,
  method = DBA_DESEQ2,
  bFullLibrarySize = TRUE,
  bSubControl = TRUE
)

# Generate enhanced results for each contrast
for (i in 1:length(config$contrasts)) {
  contrast_info <- config$contrasts[[i]]
  contrast_name <- contrast_info$name

  cat(paste("Processing contrast:", contrast_name, "\n"))

  # Get results with multiple FDR thresholds
  res_strict <- dba.report(dba_obj, contrast = i, th = 0.01) # Strict
  res_standard <- dba.report(dba_obj, contrast = i, th = 0.05) # Standard
  res_all <- dba.report(dba_obj, contrast = i, th = 1) # All sites

  # Save results
  write.csv(
    as.data.frame(res_strict),
    file.path(output_dir, paste0(contrast_name, "_strict_FDR001.csv"))
  )
  write.csv(
    as.data.frame(res_standard),
    file.path(output_dir, paste0(contrast_name, "_standard_FDR005.csv"))
  )
  write.csv(
    as.data.frame(res_all),
    file.path(output_dir, paste0(contrast_name, "_all_sites.csv"))
  )

  # Enhanced MA plot
  pdf(file.path(output_dir, paste0(contrast_name, "_MA_plot.pdf")), width = 10, height = 8)
  dba.plotMA(dba_obj,
    contrast = i,
    main = paste("MA Plot:", contrast_info$description),
    sub = paste("FDR < 0.05:", length(res_standard), "sites")
  )
  dev.off()

  # Volcano plot
  if (length(res_all) > 0) {
    res_df <- as.data.frame(res_all)
    res_df$log2FoldChange <- log2(res_df$Fold)
    res_df$negLog10FDR <- -log10(res_df$FDR)
    res_df$significant <- res_df$FDR < 0.05 & abs(res_df$log2FoldChange) > log2(1.5)

    p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = negLog10FDR, color = significant)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c("grey70", "red")) +
      theme_classic() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "blue") +
      labs(
        title = paste("Volcano Plot:", contrast_info$description),
        x = "log2 Fold Change",
        y = "-log10(FDR)",
        subtitle = paste("Significant sites (FDR<0.05, |FC|>1.5):", sum(res_df$significant))
      ) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

    ggsave(file.path(output_dir, paste0(contrast_name, "_volcano_plot.pdf")),
      p_volcano,
      width = 10, height = 8
    )
  }

  # Statistics summary
  stats_summary <- data.frame(
    Contrast = contrast_name,
    Total_Sites = length(res_all),
    Significant_FDR005 = length(res_standard),
    Significant_FDR001 = length(res_strict),
    Upregulated_FDR005 = sum(res_standard$Fold > 1),
    Downregulated_FDR005 = sum(res_standard$Fold < 1)
  )

  if (i == 1) {
    all_stats <- stats_summary
  } else {
    all_stats <- rbind(all_stats, stats_summary)
  }
}

# Save comprehensive statistics
write.csv(all_stats, file.path(output_dir, "differential_binding_statistics.csv"), row.names = FALSE)

# Generate final summary plots
pdf(file.path(output_dir, "04_final_correlation_heatmap.pdf"), width = 12, height = 12)
dba.plotHeatmap(dba_obj,
  correlations = TRUE,
  main = "Final Sample Correlations\n(Post-Analysis)"
)
dev.off()

# Venn diagram if multiple contrasts
if (length(config$contrasts) > 1) {
  pdf(file.path(output_dir, "05_contrast_venn_diagram.pdf"), width = 12, height = 10)
  dba.plotVenn(dba_obj,
    contrast = 1:length(config$contrasts),
    main = "Overlap of Differential Binding Sites"
  )
  dev.off()
}

cat("\n=== IMPROVED DIFFBIND ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_dir, "\n")
cat("\nSummary of significant sites (FDR < 0.05):\n")
print(all_stats[, c("Contrast", "Significant_FDR005", "Upregulated_FDR005", "Downregulated_FDR005")])
