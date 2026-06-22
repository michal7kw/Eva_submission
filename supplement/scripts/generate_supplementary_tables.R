#!/usr/bin/env Rscript
# =============================================================================
# Generate Supplementary Tables for SRF_Eva Analysis
# =============================================================================
# Table S2: DEGs RNA-seq (significant only, padj < 0.05)
# Table S3: GO and GSEA categories (statistically significant)
# Table S5: CUT&Run peaks (TES and TEAD1, annotated)
# =============================================================================

# Set paths
base_dir <- "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission"
output_dir <- file.path(base_dir, "supplementary_tables")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("Generating Supplementary Tables\n")
cat("============================================================\n\n")

# =============================================================================
# TABLE S2: DEGs RNA-seq
# =============================================================================
cat(">>> Generating Table S2: DEGs RNA-seq\n")

deseq_file <- file.path(base_dir, "SRF_Eva_RNA/results/05_deseq2/significant_genes_TES_vs_GFP.txt")

if (file.exists(deseq_file)) {
  # Read DESeq2 results
  degs <- read.delim(deseq_file, stringsAsFactors = FALSE)

  # Clean up column names for publication
  colnames(degs) <- c("Gene_ID", "Gene_Symbol", "baseMean", "log2FoldChange",
                       "lfcSE", "stat", "pvalue", "padj")

  # Sort by adjusted p-value (most significant first)
  degs <- degs[order(degs$padj), ]

  # Add direction column
  degs$Direction <- ifelse(degs$log2FoldChange > 0, "UP", "DOWN")

  # Reorder columns
  degs <- degs[, c("Gene_ID", "Gene_Symbol", "Direction", "baseMean", "log2FoldChange",
                   "lfcSE", "stat", "pvalue", "padj")]

  # Export
  output_file <- file.path(output_dir, "Table_S2_DEGs_RNA_seq.csv")
  write.csv(degs, output_file, row.names = FALSE)

  # Summary
  cat(sprintf("   Total significant DEGs: %d\n", nrow(degs)))
  cat(sprintf("   Upregulated: %d\n", sum(degs$Direction == "UP")))
  cat(sprintf("   Downregulated: %d\n", sum(degs$Direction == "DOWN")))
  cat(sprintf("   Output: %s\n\n", output_file))
} else {
  cat("   ERROR: DESeq2 results file not found!\n\n")
}

# =============================================================================
# TABLE S3: GO and GSEA Categories
# =============================================================================
cat(">>> Generating Table S3: GO and GSEA Categories\n")

gsea_dir <- file.path(base_dir, "SRF_Eva_RNA/results/06_gsea")
go_dir <- file.path(base_dir, "SRF_Eva_RNA/results/06_go_enrichment")

# Initialize combined results
combined_results <- data.frame()

# --- GSEA Results ---
gsea_collections <- c("hallmark", "GO_BP", "KEGG", "reactome", "C6_oncogenic")
gsea_names <- c("Hallmark", "GO_BP", "KEGG", "Reactome", "Oncogenic")

for (i in seq_along(gsea_collections)) {
  gsea_file <- file.path(gsea_dir, paste0("gsea_results_", gsea_collections[i], ".csv"))

  if (file.exists(gsea_file)) {
    gsea_data <- read.csv(gsea_file, stringsAsFactors = FALSE)

    if (nrow(gsea_data) > 0) {
      # Filter by significance (FDR < 0.25 is standard for GSEA)
      gsea_sig <- gsea_data[gsea_data$p.adjust < 0.25, ]

      if (nrow(gsea_sig) > 0) {
        gsea_formatted <- data.frame(
          Analysis_Type = "GSEA",
          Collection = gsea_names[i],
          Direction = ifelse(gsea_sig$NES > 0, "UP", "DOWN"),
          Term_ID = gsea_sig$ID,
          Term_Name = gsea_sig$Description,
          Enrichment_Score = gsea_sig$NES,
          pvalue = gsea_sig$pvalue,
          FDR = gsea_sig$p.adjust,
          Gene_Count = gsea_sig$setSize,
          stringsAsFactors = FALSE
        )

        combined_results <- rbind(combined_results, gsea_formatted)
        cat(sprintf("   GSEA %s: %d significant terms\n", gsea_names[i], nrow(gsea_sig)))
      }
    }
  }
}

# --- GO Enrichment Results (Over-representation Analysis) ---
go_ontologies <- c("BP", "MF", "CC")
go_directions <- c("Upregulated", "Downregulated")

for (direction in go_directions) {
  for (ontology in go_ontologies) {
    go_file <- file.path(go_dir, tolower(direction),
                         paste0("GO_", ontology, "_", direction, "_results.csv"))

    if (file.exists(go_file)) {
      go_data <- read.csv(go_file, stringsAsFactors = FALSE)

      if (nrow(go_data) > 0 && "qvalue" %in% colnames(go_data)) {
        # Filter by significance (qvalue < 0.05)
        go_sig <- go_data[!is.na(go_data$qvalue) & go_data$qvalue < 0.05, ]

        if (nrow(go_sig) > 0) {
          go_formatted <- data.frame(
            Analysis_Type = "GO_ORA",
            Collection = paste0("GO:", ontology),
            Direction = ifelse(direction == "Upregulated", "UP", "DOWN"),
            Term_ID = go_sig$ID,
            Term_Name = go_sig$Description,
            Enrichment_Score = go_sig$FoldEnrichment,
            pvalue = go_sig$pvalue,
            FDR = go_sig$qvalue,
            Gene_Count = go_sig$Count,
            stringsAsFactors = FALSE
          )

          combined_results <- rbind(combined_results, go_formatted)
          cat(sprintf("   GO %s %s: %d significant terms\n", ontology, direction, nrow(go_sig)))
        }
      }
    }
  }
}

# Sort and export
if (nrow(combined_results) > 0) {
  combined_results <- combined_results[order(combined_results$Analysis_Type,
                                              combined_results$Collection,
                                              combined_results$Direction,
                                              combined_results$FDR), ]

  output_file <- file.path(output_dir, "Table_S3_GO_GSEA.csv")
  write.csv(combined_results, output_file, row.names = FALSE)

  cat(sprintf("\n   Total significant pathways: %d\n", nrow(combined_results)))
  cat(sprintf("   Output: %s\n\n", output_file))
} else {
  cat("   WARNING: No significant GO/GSEA results found!\n\n")
}

# =============================================================================
# TABLE S5: CUT&Run Peaks
# =============================================================================
cat(">>> Generating Table S5: CUT&Run Peaks\n")

peaks_dir <- file.path(base_dir, "SRF_Eva_Cut_and_Tag/results/07_analysis_narrow")

# Function to process peak file
process_peaks <- function(input_file, sample_name) {
  if (!file.exists(input_file)) {
    cat(sprintf("   ERROR: %s not found!\n", input_file))
    return(NULL)
  }

  peaks <- read.csv(input_file, stringsAsFactors = FALSE)

  # Clean up and rename columns
  peaks_clean <- data.frame(
    Peak_ID = peaks$V4,
    Chromosome = peaks$seqnames,
    Start = peaks$start,
    End = peaks$end,
    Width = peaks$width,
    Score = peaks$V5,
    Signal_Value = peaks$V7,
    pValue = peaks$V8,
    qValue = peaks$V9,
    Summit = peaks$V10,
    Annotation = peaks$annotation,
    Nearest_Gene_ID = peaks$geneId,
    Nearest_Transcript = peaks$transcriptId,
    Distance_to_TSS = peaks$distanceToTSS,
    stringsAsFactors = FALSE
  )

  # Sort by chromosome and position
  peaks_clean <- peaks_clean[order(peaks_clean$Chromosome, peaks_clean$Start), ]

  return(peaks_clean)
}

# Process TES peaks
tes_file <- file.path(peaks_dir, "TES_peaks_annotated.csv")
tes_peaks <- process_peaks(tes_file, "TES")

if (!is.null(tes_peaks)) {
  output_file <- file.path(output_dir, "Table_S5a_TES_peaks.csv")
  write.csv(tes_peaks, output_file, row.names = FALSE)
  cat(sprintf("   TES peaks: %d\n", nrow(tes_peaks)))
  cat(sprintf("   Output: %s\n", output_file))
}

# Process TEAD1 peaks
tead1_file <- file.path(peaks_dir, "TEAD1_peaks_annotated.csv")
tead1_peaks <- process_peaks(tead1_file, "TEAD1")

if (!is.null(tead1_peaks)) {
  output_file <- file.path(output_dir, "Table_S5b_TEAD1_peaks.csv")
  write.csv(tead1_peaks, output_file, row.names = FALSE)
  cat(sprintf("   TEAD1 peaks: %d\n", nrow(tead1_peaks)))
  cat(sprintf("   Output: %s\n", output_file))
}

# =============================================================================
# Summary
# =============================================================================
cat("\n============================================================\n")
cat("Summary of Generated Files\n")
cat("============================================================\n")

output_files <- list.files(output_dir, pattern = "Table_S.*\\.csv$", full.names = TRUE)
for (f in output_files) {
  file_info <- file.info(f)
  n_lines <- length(readLines(f)) - 1  # Subtract header
  cat(sprintf("   %s\n", basename(f)))
  cat(sprintf("      Rows: %d | Size: %.1f KB\n", n_lines, file_info$size / 1024))
}

cat("\n============================================================\n")
cat("Done!\n")
cat("============================================================\n")
