#!/usr/bin/env Rscript
# Integrated Peak Annotation and GO Analysis
# Improved biological interpretation with statistical rigor

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(GenomeInfoDb)
library(dplyr)
library(yaml)
library(ReactomePA)
library(DOSE)
library(enrichplot)

# Load configuration
config <- yaml.load_file("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/config/pipeline_config.yaml")
setwd(config$directories$base)

# Set up directories
output_dir <- file.path(config$directories$results, "integrated_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Function to read and process differential binding results
load_diffbind_results <- function(contrast_name, fdr_threshold = 0.05) {
  diffbind_dir <- file.path(config$directories$results, "diffbind_analysis")
  file_path <- file.path(diffbind_dir, paste0(contrast_name, "_standard_FDR005.csv"))

  if (!file.exists(file_path)) {
    cat("Warning: DiffBind results not found for", contrast_name, "\n")
    return(NULL)
  }

  results <- read.csv(file_path, row.names = 1)

  # Convert to GRanges for ChIPseeker
  gr <- GRanges(
    seqnames = results$seqnames,
    ranges = IRanges(start = results$start, end = results$end),
    strand = "*",
    fold_change = results$Fold,
    pvalue = results$p.value,
    padj = results$FDR
  )

  # Fix chromosome naming
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  GenomeInfoDb::genome(gr) <- "hg38"

  return(gr)
}

# Function for comprehensive peak annotation
annotate_peaks_comprehensive <- function(peaks, name) {
  if (is.null(peaks) || length(peaks) == 0) {
    return(NULL)
  }

  # Annotate peaks
  anno <- annotatePeak(peaks,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db"
  )

  # Save annotation results
  write.csv(
    as.data.frame(anno),
    file.path(output_dir, paste0(name, "_annotated_peaks.csv"))
  )

  return(anno)
}

# Function for enhanced GO and pathway analysis
perform_pathway_analysis <- function(gene_list, name, universe = NULL) {
  if (length(gene_list) < 10) {
    cat("Too few genes for", name, "(", length(gene_list), "). Skipping analysis.\n")
    return(NULL)
  }

  results_list <- list()

  # GO Biological Process
  ego_bp <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    universe = universe,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )

  if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
    results_list$GO_BP <- ego_bp
    write.csv(ego_bp@result, file.path(output_dir, paste0(name, "_GO_BP.csv")))
  }

  # GO Molecular Function
  ego_mf <- enrichGO(
    gene = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    universe = universe,
    readable = TRUE
  )

  if (!is.null(ego_mf) && nrow(ego_mf@result) > 0) {
    results_list$GO_MF <- ego_mf
    write.csv(ego_mf@result, file.path(output_dir, paste0(name, "_GO_MF.csv")))
  }

  # KEGG Pathway Analysis
  kegg_result <- enrichKEGG(
    gene = gene_list,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    results_list$KEGG <- kegg_result
    write.csv(kegg_result@result, file.path(output_dir, paste0(name, "_KEGG.csv")))
  }

  # Reactome Pathway Analysis
  reactome_result <- enrichPathway(
    gene = gene_list,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )

  if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
    results_list$Reactome <- reactome_result
    write.csv(reactome_result@result, file.path(output_dir, paste0(name, "_Reactome.csv")))
  }

  return(results_list)
}

# Function to create enhanced visualizations
create_pathway_plots <- function(pathway_results, name) {
  if (is.null(pathway_results) || length(pathway_results) == 0) {
    return()
  }

  plot_list <- list()

  for (analysis_type in names(pathway_results)) {
    result <- pathway_results[[analysis_type]]

    if (nrow(result@result) == 0) next

    # Dotplot
    p1 <- dotplot(result, showCategory = 15, font.size = 11) +
      ggtitle(paste(name, "-", analysis_type, "Enrichment")) +
      theme_classic()
    plot_list[[paste0(analysis_type, "_dot")]] <- p1

    # Barplot
    p2 <- barplot(result, showCategory = 10) +
      ggtitle(paste(name, "-", analysis_type, "Top Terms")) +
      theme_classic()
    plot_list[[paste0(analysis_type, "_bar")]] <- p2

    # Network plot for GO terms
    if (analysis_type %in% c("GO_BP", "GO_MF") && nrow(result@result) >= 5) {
      tryCatch(
        {
          p3 <- emapplot(pairwise_termsim(result), showCategory = 20) +
            ggtitle(paste(name, "-", analysis_type, "Network"))
          plot_list[[paste0(analysis_type, "_network")]] <- p3
        },
        error = function(e) {
          cat("Could not create network plot for", analysis_type, ":", e$message, "\n")
        }
      )
    }
  }

  # Save plots
  for (plot_name in names(plot_list)) {
    filename <- file.path(output_dir, paste0(name, "_", plot_name, ".pdf"))
    ggsave(filename, plot_list[[plot_name]], width = 12, height = 10)

    # Also save as PNG
    png_filename <- file.path(output_dir, paste0(name, "_", plot_name, ".png"))
    ggsave(png_filename, plot_list[[plot_name]], width = 12, height = 10, dpi = 300)
  }
}

# Main analysis workflow
cat("Starting integrated annotation and pathway analysis...\n")

# Process each contrast
for (contrast_info in config$contrasts) {
  contrast_name <- contrast_info$name
  cat("\nProcessing contrast:", contrast_name, "\n")

  # Load differential binding results
  diff_peaks <- load_diffbind_results(contrast_name)

  if (!is.null(diff_peaks)) {
    # Separate up and down-regulated peaks
    up_peaks <- diff_peaks[diff_peaks$fold_change > 1.5]
    down_peaks <- diff_peaks[diff_peaks$fold_change < 1 / 1.5]

    # Annotate peaks
    all_anno <- annotate_peaks_comprehensive(diff_peaks, paste0(contrast_name, "_all"))
    up_anno <- annotate_peaks_comprehensive(up_peaks, paste0(contrast_name, "_up"))
    down_anno <- annotate_peaks_comprehensive(down_peaks, paste0(contrast_name, "_down"))

    # Extract gene lists for pathway analysis
    if (!is.null(all_anno)) {
      all_genes <- as.data.frame(all_anno)$geneId
      all_genes <- all_genes[!is.na(all_genes)]

      # Get universe of all genes for background
      universe_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")

      # Pathway analysis for all differential sites
      all_pathways <- perform_pathway_analysis(all_genes, paste0(contrast_name, "_all"), universe_genes)
      create_pathway_plots(all_pathways, paste0(contrast_name, "_all"))

      # Pathway analysis for up-regulated sites
      if (!is.null(up_anno) && length(up_peaks) > 10) {
        up_genes <- as.data.frame(up_anno)$geneId
        up_genes <- up_genes[!is.na(up_genes)]
        up_pathways <- perform_pathway_analysis(up_genes, paste0(contrast_name, "_upregulated"), universe_genes)
        create_pathway_plots(up_pathways, paste0(contrast_name, "_upregulated"))
      }

      # Pathway analysis for down-regulated sites
      if (!is.null(down_anno) && length(down_peaks) > 10) {
        down_genes <- as.data.frame(down_anno)$geneId
        down_genes <- down_genes[!is.na(down_genes)]
        down_pathways <- perform_pathway_analysis(down_genes, paste0(contrast_name, "_downregulated"), universe_genes)
        create_pathway_plots(down_pathways, paste0(contrast_name, "_downregulated"))
      }
    }
  }
}

# Create comparative analysis for TES-specific vs TEAD1-specific
cat("\nPerforming TES vs TEAD1 specific analysis...\n")

# This requires the differential binding results to identify condition-specific sites
tes_vs_tead1_peaks <- load_diffbind_results("TES_vs_TEAD1")

if (!is.null(tes_vs_tead1_peaks)) {
  # TES-specific: significantly higher in TES
  tes_specific_peaks <- tes_vs_tead1_peaks[tes_vs_tead1_peaks$fold_change > 1.5 & tes_vs_tead1_peaks$padj < 0.05]

  # TEAD1-specific: significantly higher in TEAD1
  tead1_specific_peaks <- tes_vs_tead1_peaks[tes_vs_tead1_peaks$fold_change < 1 / 1.5 & tes_vs_tead1_peaks$padj < 0.05]

  if (length(tes_specific_peaks) > 0) {
    tes_specific_anno <- annotate_peaks_comprehensive(tes_specific_peaks, "TES_specific")
    if (!is.null(tes_specific_anno)) {
      tes_specific_genes <- as.data.frame(tes_specific_anno)$geneId
      tes_specific_genes <- tes_specific_genes[!is.na(tes_specific_genes)]

      # Save TES-specific genes
      tes_gene_symbols <- bitr(tes_specific_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
      write.table(tes_gene_symbols$SYMBOL, file.path(output_dir, "TES_specific_genes.txt"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )

      # Comprehensive pathway analysis
      tes_pathways <- perform_pathway_analysis(tes_specific_genes, "TES_specific")
      create_pathway_plots(tes_pathways, "TES_specific")
    }
  }

  if (length(tead1_specific_peaks) > 0) {
    tead1_specific_anno <- annotate_peaks_comprehensive(tead1_specific_peaks, "TEAD1_specific")
    if (!is.null(tead1_specific_anno)) {
      tead1_specific_genes <- as.data.frame(tead1_specific_anno)$geneId
      tead1_specific_genes <- tead1_specific_genes[!is.na(tead1_specific_genes)]

      # Save TEAD1-specific genes
      tead1_gene_symbols <- bitr(tead1_specific_genes, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
      write.table(tead1_gene_symbols$SYMBOL, file.path(output_dir, "TEAD1_specific_genes.txt"),
        quote = FALSE, row.names = FALSE, col.names = FALSE
      )

      # Comprehensive pathway analysis
      tead1_pathways <- perform_pathway_analysis(tead1_specific_genes, "TEAD1_specific")
      create_pathway_plots(tead1_pathways, "TEAD1_specific")
    }
  }
}

# Generate summary report
summary_file <- file.path(output_dir, "analysis_summary.txt")
cat("=== INTEGRATED ANNOTATION AND PATHWAY ANALYSIS SUMMARY ===\n", file = summary_file)
cat("Analysis completed on:", as.character(Sys.time()), "\n\n", file = summary_file, append = TRUE)

cat("Files generated in:", output_dir, "\n", file = summary_file, append = TRUE)
cat("- Peak annotations: *_annotated_peaks.csv\n", file = summary_file, append = TRUE)
cat("- GO Biological Process: *_GO_BP.csv\n", file = summary_file, append = TRUE)
cat("- GO Molecular Function: *_GO_MF.csv\n", file = summary_file, append = TRUE)
cat("- KEGG Pathways: *_KEGG.csv\n", file = summary_file, append = TRUE)
cat("- Reactome Pathways: *_Reactome.csv\n", file = summary_file, append = TRUE)
cat("- Visualizations: *.pdf and *.png files\n", file = summary_file, append = TRUE)

cat("\n=== INTEGRATED ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_dir, "\n")
