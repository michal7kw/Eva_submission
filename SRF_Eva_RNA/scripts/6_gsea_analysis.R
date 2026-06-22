#!/usr/bin/env Rscript

#===============================================================================
# GSEA Analysis Script
# Creates classic enrichment plots with NES and FDR annotations
#
# Usage:
#   Rscript gsea_analysis.R                # Default: plots only (use cached .rds results)
#   Rscript gsea_analysis.R --compute      # Run GSEA computation (reuse cache if available)
#   Rscript gsea_analysis.R --force        # Force recomputation of GSEA (ignore cache)
#===============================================================================

# Load required libraries
suppressPackageStartupMessages({
    library(clusterProfiler)
    library(enrichplot)
    library(org.Hs.eg.db)
    library(msigdbr)
    library(ggplot2)
    library(dplyr)
    library(tibble)
    library(gridExtra)
    library(patchwork)
})

#===============================================================================
# Parse command line arguments
#===============================================================================

args <- commandArgs(trailingOnly = TRUE)
# Default: PLOTS_ONLY=TRUE (skip GSEA computation, use cached results)
# Use --compute to run GSEA computation
# Use --force to force recomputation even if cache exists
PLOTS_ONLY <- !("--compute" %in% args)
FORCE_RECOMPUTE <- "--force" %in% args

cat("=== GSEA Analysis ===\n")
cat("Timestamp:", as.character(Sys.time()), "\n")
cat("Arguments:", if (length(args) > 0) paste(args, collapse = ", ") else "none", "\n")
cat("Force recompute:", FORCE_RECOMPUTE, "\n")
cat("Plots only:", PLOTS_ONLY, "\n\n")

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/06_gsea")

#===============================================================================
# Configuration
#===============================================================================

# GSEA parameters
MIN_GENE_SET_SIZE <- 15
MAX_GENE_SET_SIZE <- 500
N_PERMUTATIONS <- 10000
PVALUE_CUTOFF <- 0.05
FDR_CUTOFF <- 0.25

# Plot parameters - increased size to prevent cropping
PLOT_WIDTH <- 18
PLOT_HEIGHT <- 10
DPI <- 180

# Font sizes for publication
TITLE_SIZE <- 20
AXIS_TITLE_SIZE <- 18
AXIS_TEXT_SIZE <- 16
ANNOTATION_SIZE <- 8

# Number of top pathways to plot
N_TOP_PATHWAYS <- 20
N_INDIVIDUAL_PLOTS <- 5  # Top pathways for individual enrichment plots

# Comparison name for plot subtitles
COMPARISON_NAME <- "TES vs GFP"

# Color scheme matching heatmaps (from heatmap_publication.R)
# GFP = Brown (#8B4513), TES = Teal (#2E8B8B)
# Semi-saturated versions for gradient bar (more visible than fully faded):
COLOR_TES_GRADIENT <- "#5FBFBF"   # Medium teal for TES/upregulated (high rank) - more visible
COLOR_GFP_GRADIENT <- "#C9A86C"   # Medium tan/brown for GFP/downregulated (low rank) - more visible
COLOR_ENRICHMENT_LINE <- "#D73027"  # Red enrichment score line

# Alternative: fully faded versions if preferred
COLOR_TES_FADED <- "#A8D8D8"   # Light teal
COLOR_GFP_FADED <- "#DEB887"   # Light brown/burlywood

# Specific pathways to always plot (if significant)
PATHWAYS_OF_INTEREST <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "REACTOME_M_PHASE",
    "REACTOME_DNA_REPLICATION",
    "CORDENONSI_YAP_CONSERVED_SIGNATURE"
)

#===============================================================================
# Custom GSEA plot function with TES/GFP color gradient
#===============================================================================

create_gsea_plot_custom_colors <- function(gsea_result, pathway_id,
                                            title = NULL,
                                            add_stats = TRUE,
                                            color_high = COLOR_TES_GRADIENT,
                                            color_low = COLOR_GFP_GRADIENT,
                                            line_color = COLOR_ENRICHMENT_LINE) {
    # Get pathway data
    idx <- which(gsea_result@result$ID == pathway_id)
    if (length(idx) == 0) {
        return(NULL)
    }

    pathway_data <- gsea_result@result[idx, ]
    nes_val <- round(pathway_data$NES, 4)
    pval <- pathway_data$pvalue
    qval <- pathway_data$p.adjust

    if (is.null(title)) {
        title <- gsub("_", " ", pathway_data$Description)
        if (nchar(title) > 60) {
            title <- paste0(substr(title, 1, 57), "...")
        }
    }

    # Get gene list and calculate running enrichment score
    gene_list <- gsea_result@geneList
    gene_set <- gsea_result@geneSets[[pathway_id]]

    # Calculate running enrichment score
    n <- length(gene_list)
    gene_hits <- names(gene_list) %in% gene_set

    # Running sum statistics
    hit_indicator <- as.integer(gene_hits)
    no_hit_indicator <- 1 - hit_indicator

    # Calculate running ES
    Phit <- cumsum(abs(gene_list) * hit_indicator) / sum(abs(gene_list[gene_hits]))
    Pmiss <- cumsum(no_hit_indicator) / sum(no_hit_indicator)
    running_es <- Phit - Pmiss

    # Create data frames for plotting
    es_data <- data.frame(
        rank = 1:n,
        running_es = running_es
    )

    hit_data <- data.frame(
        rank = which(gene_hits),
        y = 0
    )

    rank_data <- data.frame(
        rank = 1:n,
        value = gene_list
    )

    # Panel 1: Running Enrichment Score
    p1 <- ggplot(es_data, aes(x = rank, y = running_es)) +
        geom_line(color = line_color, linewidth = 1.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        labs(y = "Enrichment Score", x = NULL, title = title) +
        theme_classic(base_size = AXIS_TEXT_SIZE) +
        theme(
            plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
            axis.title.y = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
            axis.text = element_text(size = AXIS_TEXT_SIZE, color = "black"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = margin(10, 15, 5, 15)
        )

    # Add NES/FDR annotation if requested
    if (add_stats) {
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)
        annotation_text <- sprintf("NES: %.3f\nFDR: %s", nes_val, qval_text)

        p1 <- p1 +
            annotate("text", x = n * 0.85, y = max(running_es) * 0.8,
                     label = annotation_text, hjust = 0.5, vjust = 1,
                     size = ANNOTATION_SIZE, fontface = "bold")
    }

    # Panel 2: Gene hit barcode
    p2 <- ggplot(hit_data, aes(x = rank, y = y)) +
        geom_segment(aes(xend = rank, yend = 1), color = "black", linewidth = 0.3) +
        labs(x = NULL, y = NULL) +
        theme_void() +
        theme(
            plot.margin = margin(0, 15, 0, 15)
        ) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(limits = c(1, n), expand = c(0, 0))

    # Panel 3: Ranked list gradient with custom colors
    # Create gradient based on POSITION (not value) for smooth color transition
    rank_data$position_color <- (rank_data$rank - 1) / (n - 1)  # 0 to 1

    p3 <- ggplot(rank_data, aes(x = rank, y = 1, fill = position_color)) +
        geom_raster(interpolate = TRUE) +  # Use geom_raster to avoid white line artifacts
        scale_fill_gradient(low = color_high, high = color_low, guide = "none") +  # TES (left) to GFP (right)
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(x = "Rank in Ordered Dataset", y = NULL) +
        theme_classic(base_size = AXIS_TEXT_SIZE) +
        theme(
            axis.title.x = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
            axis.text.x = element_text(size = AXIS_TEXT_SIZE, color = "black"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            plot.margin = margin(0, 15, 10, 15)
        )

    # Add labels for TES and GFP sides
    p3 <- p3 +
        annotate("text", x = n * 0.05, y = 1, label = "TES",
                 size = 4, color = "white", fontface = "bold") +
        annotate("text", x = n * 0.95, y = 1, label = "GFP",
                 size = 4, color = "white", fontface = "bold")

    # Combine panels
    combined_plot <- p1 / p2 / p3 +
        plot_layout(heights = c(3, 0.5, 1))

    return(combined_plot)
}

#===============================================================================
# Load and prepare data
#===============================================================================

cat("Loading DESeq2 results...\n")
res <- read.table("../05_deseq2/deseq2_results_TES_vs_GFP.txt",
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat(sprintf("Loaded %d genes\n", nrow(res)))

# Remove NA values
res <- res[!is.na(res$log2FoldChange) & !is.na(res$pvalue), ]
cat(sprintf("After removing NA: %d genes\n", nrow(res)))

#===============================================================================
# Create ranked gene list
#===============================================================================

cat("\nCreating ranked gene list...\n")

# Clean Ensembl IDs (remove version)
res$ensembl_clean <- gsub("\\..*", "", res$gene_id)

# Map to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = res$ensembl_clean,
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

res$entrez <- entrez_ids

# Remove genes without Entrez IDs
res_entrez <- res[!is.na(res$entrez), ]
cat(sprintf("Genes with Entrez IDs: %d\n", nrow(res_entrez)))

# Create ranked list for GSEA
# Use shrunken LFC (apeglm) if available — more robust ranking that reduces noise
# from low-count genes. Falls back to DESeq2 stat column if shrunken LFC not present.
if ("log2FoldChange_shrunk" %in% colnames(res_entrez) &&
    sum(!is.na(res_entrez$log2FoldChange_shrunk)) > 0) {
    cat("Using apeglm-shrunken log2FoldChange for GSEA ranking\n")
    res_entrez$rank_metric <- res_entrez$log2FoldChange_shrunk
} else if ("stat" %in% colnames(res_entrez)) {
    cat("Shrunken LFC not available, using DESeq2 Wald statistic for GSEA ranking\n")
    res_entrez$rank_metric <- res_entrez$stat
} else {
    cat("Using signed -log10(pvalue) for GSEA ranking\n")
    res_entrez$rank_metric <- -log10(res_entrez$pvalue) * sign(res_entrez$log2FoldChange)
}

# Handle infinite/NA values
res_entrez$rank_metric[is.infinite(res_entrez$rank_metric)] <-
    sign(res_entrez$rank_metric[is.infinite(res_entrez$rank_metric)]) * 300
res_entrez$rank_metric[is.na(res_entrez$rank_metric)] <- 0

# Create named vector for GSEA
gene_list <- res_entrez$rank_metric
names(gene_list) <- res_entrez$entrez

# Sort in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

cat(sprintf("Final ranked gene list: %d genes\n", length(gene_list)))
cat(sprintf("Range: %.2f to %.2f\n", min(gene_list), max(gene_list)))

#===============================================================================
# Get MSigDB gene sets (compatible with msigdbr v10.0.0+)
#===============================================================================

cat("\nLoading MSigDB gene sets...\n")

# Check msigdbr version and use appropriate syntax
msigdbr_version <- packageVersion("msigdbr")
cat(sprintf("msigdbr version: %s\n", msigdbr_version))

# First, check column names in msigdbr
cat("Checking msigdbr column names...\n")
test_df <- msigdbr(species = "Homo sapiens", collection = "H")
cat("Available columns:", paste(colnames(test_df), collapse = ", "), "\n")

# Determine the correct column name for Entrez gene IDs
entrez_col <- if ("entrez_gene" %in% colnames(test_df)) {
    "entrez_gene"
} else if ("entrez_id" %in% colnames(test_df)) {
    "entrez_id"
} else if ("ncbi_gene" %in% colnames(test_df)) {
    "ncbi_gene"
} else if ("human_entrez_gene" %in% colnames(test_df)) {
    "human_entrez_gene"
} else {
    # Find any column with "entrez" or "ncbi" in the name
    entrez_candidates <- grep("entrez|ncbi", colnames(test_df), value = TRUE, ignore.case = TRUE)
    if (length(entrez_candidates) > 0) entrez_candidates[1] else NULL
}

if (is.null(entrez_col)) {
    stop("Could not find Entrez gene ID column in msigdbr output")
}
cat(sprintf("Using Entrez column: %s\n", entrez_col))

# Function to get gene sets with version compatibility
get_msigdb_sets <- function(collection, subcollection = NULL) {
    tryCatch({
        if (is.null(subcollection)) {
            df <- msigdbr(species = "Homo sapiens", collection = collection)
        } else {
            df <- msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
        }

        # Select gene set name and entrez ID using the correct column
        result <- df %>%
            dplyr::select(gs_name, !!sym(entrez_col)) %>%
            dplyr::rename(entrez_gene = !!sym(entrez_col)) %>%
            as.data.frame()

        return(result)
    }, error = function(e) {
        cat(sprintf("  Error loading %s/%s: %s\n", collection, subcollection, e$message))
        return(NULL)
    })
}

# Hallmark gene sets (H)
hallmark_sets <- get_msigdb_sets("H")

# GO Biological Process - try different subcollection names
gobp_sets <- get_msigdb_sets("C5", "GO:BP")
if (is.null(gobp_sets)) {
    gobp_sets <- get_msigdb_sets("C5", "BP")
}
if (is.null(gobp_sets)) {
    # Fall back to all C5 and filter
    c5_all <- get_msigdb_sets("C5")
    if (!is.null(c5_all)) {
        gobp_sets <- c5_all[grepl("^GOBP_", c5_all$gs_name), ]
    }
}

# Reactome pathways - try different subcollection names
reactome_sets <- get_msigdb_sets("C2", "CP:REACTOME")
if (is.null(reactome_sets)) {
    reactome_sets <- get_msigdb_sets("C2", "REACTOME")
}
if (is.null(reactome_sets)) {
    # Fall back to all C2 and filter
    c2_all <- get_msigdb_sets("C2")
    if (!is.null(c2_all)) {
        reactome_sets <- c2_all[grepl("^REACTOME_", c2_all$gs_name), ]
    }
}

# KEGG pathways - try different subcollection names
kegg_sets <- get_msigdb_sets("C2", "CP:KEGG_MEDICUS")
if (is.null(kegg_sets) || nrow(kegg_sets) == 0) {
    kegg_sets <- get_msigdb_sets("C2", "KEGG")
}
if (is.null(kegg_sets) || nrow(kegg_sets) == 0) {
    # Fall back to all C2 and filter
    c2_all <- get_msigdb_sets("C2")
    if (!is.null(c2_all)) {
        kegg_sets <- c2_all[grepl("^KEGG_", c2_all$gs_name), ]
    }
}

# C6 Oncogenic signatures (for YAP, MYC, RAS, etc.)
c6_sets <- get_msigdb_sets("C6")

# Report loaded sets
if (!is.null(hallmark_sets)) cat(sprintf("Hallmark sets: %d\n", length(unique(hallmark_sets$gs_name))))
if (!is.null(gobp_sets)) cat(sprintf("GO:BP sets: %d\n", length(unique(gobp_sets$gs_name))))
if (!is.null(reactome_sets)) cat(sprintf("Reactome sets: %d\n", length(unique(reactome_sets$gs_name))))
if (!is.null(kegg_sets)) cat(sprintf("KEGG sets: %d\n", length(unique(kegg_sets$gs_name))))
if (!is.null(c6_sets)) cat(sprintf("C6 Oncogenic sets: %d\n", length(unique(c6_sets$gs_name))))

#===============================================================================
# Custom function to create publication-ready GSEA plot
#===============================================================================

create_gsea_plot <- function(gsea_result, pathway_id, title = NULL,
                              show_nes = TRUE, show_fdr = TRUE) {

    # Get pathway data
    idx <- which(gsea_result@result$ID == pathway_id)
    if (length(idx) == 0) {
        cat(sprintf("  Pathway not found: %s\n", pathway_id))
        return(NULL)
    }

    nes <- round(gsea_result@result$NES[idx], 4)
    fdr <- gsea_result@result$p.adjust[idx]
    fdr_text <- if (fdr < 0.001) "< 0.001" else sprintf("%.3f", fdr)

    if (is.null(title)) {
        title <- gsub("_", " ", gsea_result@result$Description[idx])
        # Truncate long titles
        if (nchar(title) > 60) {
            title <- paste0(substr(title, 1, 57), "...")
        }
    }

    # Create base plot
    p <- gseaplot2(gsea_result, geneSetID = pathway_id,
                   title = title,
                   color = "#D73027",
                   base_size = AXIS_TEXT_SIZE,
                   pvalue_table = FALSE)

    # Add NES and FDR annotations
    if (show_nes || show_fdr) {
        annotation_text <- ""
        if (show_nes) annotation_text <- paste0("NES: ", nes)
        if (show_fdr) annotation_text <- paste0(annotation_text, "\nFDR: ", fdr_text)

        p <- p + annotate("text",
                          x = Inf, y = Inf,
                          label = annotation_text,
                          hjust = 1.1, vjust = 1.5,
                          size = ANNOTATION_SIZE,
                          fontface = "bold")
    }

    return(p)
}

#===============================================================================
# Run GSEA for each gene set collection (with caching support)
#===============================================================================

# Function to check if cached results exist and are valid
check_cached_results <- function(name) {
    csv_file <- sprintf("gsea_results_%s.csv", name)
    rds_file <- sprintf("gsea_results_%s.rds", name)

    # Check if RDS file exists (contains full GSEA object for plotting)
    if (file.exists(rds_file)) {
        cat(sprintf("  Found cached RDS results: %s\n", rds_file))
        return(list(exists = TRUE, rds = TRUE, file = rds_file))
    }

    # Check if CSV file exists (contains results table only)
    if (file.exists(csv_file)) {
        cat(sprintf("  Found cached CSV results: %s\n", csv_file))
        return(list(exists = TRUE, rds = FALSE, file = csv_file))
    }

    return(list(exists = FALSE, rds = FALSE, file = NULL))
}

# Function to load cached GSEA results
load_cached_gsea <- function(name, gene_sets, gene_list) {
    rds_file <- sprintf("gsea_results_%s.rds", name)
    csv_file <- sprintf("gsea_results_%s.csv", name)

    # Try to load RDS file (full GSEA object - best for plotting)
    if (file.exists(rds_file)) {
        cat(sprintf("  Loading cached GSEA object from: %s\n", rds_file))
        tryCatch({
            gsea_result <- readRDS(rds_file)
            n_sig <- sum(gsea_result@result$p.adjust < FDR_CUTOFF, na.rm = TRUE)
            cat(sprintf("  Loaded %d pathways (%d significant at FDR < %.2f)\n",
                        nrow(gsea_result@result), n_sig, FDR_CUTOFF))
            return(gsea_result)
        }, error = function(e) {
            cat(sprintf("  Error loading RDS: %s\n", e$message))
            return(NULL)
        })
    }

    # If only CSV exists, we need to reconstruct the GSEA object
    # This is limited - plotting functions may not work fully
    if (file.exists(csv_file)) {
        cat(sprintf("  Note: Only CSV found. For full plotting support, rerun with --force\n"))
        cat(sprintf("  Loading results table from: %s\n", csv_file))
        tryCatch({
            results_df <- read.csv(csv_file, stringsAsFactors = FALSE)
            cat(sprintf("  Loaded %d pathways from CSV\n", nrow(results_df)))

            # Return NULL since we can't properly reconstruct the GSEA object
            # The plotting functions need the full object with geneList, etc.
            cat(sprintf("  Warning: Cannot create full GSEA object from CSV alone.\n"))
            cat(sprintf("  Run with --force to generate full results for plotting.\n"))
            return(NULL)
        }, error = function(e) {
            cat(sprintf("  Error loading CSV: %s\n", e$message))
            return(NULL)
        })
    }

    return(NULL)
}

run_gsea_analysis <- function(gene_sets, name, gene_list) {

    cat(sprintf("\n=== Running GSEA: %s ===\n", name))

    # Check for cached results unless force recompute is requested
    if (!FORCE_RECOMPUTE) {
        cached <- check_cached_results(name)
        if (cached$exists && cached$rds) {
            gsea_result <- load_cached_gsea(name, gene_sets, gene_list)
            if (!is.null(gsea_result)) {
                return(gsea_result)
            }
            cat("  Cached results invalid, recomputing...\n")
        } else if (cached$exists && !cached$rds) {
            cat("  CSV cache found but RDS needed for plotting. Recomputing...\n")
        }
    } else {
        cat("  Force recompute requested, ignoring cached results.\n")
    }

    # If plots-only mode and no valid cache, skip
    if (PLOTS_ONLY) {
        cat("  Plots-only mode: Skipping computation (no valid cache found)\n")
        return(NULL)
    }

    tryCatch({
        # Initialize RNG before GSEA to avoid '.Random.seed' not found error
        set.seed(42)

        gsea_result <- GSEA(gene_list,
                            TERM2GENE = gene_sets,
                            minGSSize = MIN_GENE_SET_SIZE,
                            maxGSSize = MAX_GENE_SET_SIZE,
                            pvalueCutoff = 1,  # Keep all for filtering later
                            pAdjustMethod = "BH",
                            nPermSimple = N_PERMUTATIONS,
                            verbose = FALSE,
                            seed = TRUE)  # Use TRUE instead of integer to use already set seed

        n_sig <- sum(gsea_result@result$p.adjust < FDR_CUTOFF, na.rm = TRUE)
        cat(sprintf("Significant pathways (FDR < %.2f): %d\n", FDR_CUTOFF, n_sig))

        # Save results as CSV (for easy viewing)
        results_df <- as.data.frame(gsea_result@result)
        write.csv(results_df, sprintf("gsea_results_%s.csv", name), row.names = FALSE)
        cat(sprintf("Saved: gsea_results_%s.csv\n", name))

        # Save full GSEA object as RDS (for plotting/reuse)
        saveRDS(gsea_result, sprintf("gsea_results_%s.rds", name))
        cat(sprintf("Saved: gsea_results_%s.rds (for reuse)\n", name))

        return(gsea_result)

    }, error = function(e) {
        cat(sprintf("Error in GSEA for %s: %s\n", name, e$message))
        return(NULL)
    })
}

#===============================================================================
# Run GSEA analyses
#===============================================================================

# Run GSEA only for successfully loaded gene sets
gsea_hallmark <- if (!is.null(hallmark_sets) && nrow(hallmark_sets) > 0) {
    run_gsea_analysis(hallmark_sets, "hallmark", gene_list)
} else { cat("Skipping Hallmark (not loaded)\n"); NULL }

gsea_gobp <- if (!is.null(gobp_sets) && nrow(gobp_sets) > 0) {
    run_gsea_analysis(gobp_sets, "GO_BP", gene_list)
} else { cat("Skipping GO:BP (not loaded)\n"); NULL }

gsea_reactome <- if (!is.null(reactome_sets) && nrow(reactome_sets) > 0) {
    run_gsea_analysis(reactome_sets, "reactome", gene_list)
} else { cat("Skipping Reactome (not loaded)\n"); NULL }

gsea_kegg <- if (!is.null(kegg_sets) && nrow(kegg_sets) > 0) {
    run_gsea_analysis(kegg_sets, "KEGG", gene_list)
} else { cat("Skipping KEGG (not loaded)\n"); NULL }

gsea_c6 <- if (!is.null(c6_sets) && nrow(c6_sets) > 0) {
    run_gsea_analysis(c6_sets, "C6_oncogenic", gene_list)
} else { cat("Skipping C6 Oncogenic (not loaded)\n"); NULL }

#===============================================================================
# Create publication-ready plots
#===============================================================================

cat("\n=== Creating GSEA Plots ===\n")

# Theme for publication with generous margins to prevent cropping
theme_gsea <- theme_classic(base_size = AXIS_TEXT_SIZE) +
    theme(
        plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
        axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
        axis.text = element_text(size = AXIS_TEXT_SIZE, color = "black"),
        legend.text = element_text(size = AXIS_TEXT_SIZE),
        legend.title = element_text(size = AXIS_TEXT_SIZE, face = "bold"),
        plot.margin = margin(20, 25, 20, 25)  # top, right, bottom, left - increased margins
    )

#-------------------------------------------------------------------------------
# Custom function to create GSEA plot with custom gradient colors
# Creates both annotated and clean versions
#-------------------------------------------------------------------------------

create_custom_gsea_plot <- function(gsea_result, pathway_id, title = NULL,
                                     add_stats = TRUE, nes_inside = FALSE) {

    # Get pathway data
    idx <- which(gsea_result@result$ID == pathway_id)
    if (length(idx) == 0) {
        cat(sprintf("  Pathway not found: %s\n", pathway_id))
        return(NULL)
    }

    pathway_data <- gsea_result@result[idx, ]
    nes_val <- round(pathway_data$NES, 4)
    pval <- pathway_data$pvalue
    qval <- pathway_data$p.adjust

    if (is.null(title)) {
        title <- gsub("_", " ", pathway_data$Description)
        if (nchar(title) > 60) {
            title <- paste0(substr(title, 1, 57), "...")
        }
    }

    # Get running enrichment score data
    gene_set_idx <- which(gsea_result@result$ID == pathway_id)

    # Extract data using gseaplot2's internal approach
    gsea_data <- gsea_result@result[gene_set_idx, ]
    gene_list_ranks <- gsea_result@geneList

    # Get the core enrichment genes
    core_genes <- strsplit(gsea_data$core_enrichment, "/")[[1]]
    leading_edge_idx <- which(names(gene_list_ranks) %in% core_genes)

    # Create base plot using gseaplot2 then modify
    base_plot <- gseaplot2(gsea_result, geneSetID = pathway_id,
                           title = title,
                           color = COLOR_ENRICHMENT_LINE,
                           base_size = AXIS_TEXT_SIZE,
                           pvalue_table = FALSE)

    # Modify the ranking gradient bar colors by accessing the ggplot layers
    # The gradient bar is typically in subplot 3
    # We need to modify the fill gradient

    # Create modified plot with custom gradient
    # Access the third panel (ranked list) and change its fill colors
    if (inherits(base_plot, "patchwork")) {
        # patchwork object - modify the third subplot
        n_panels <- length(base_plot$patches$plots) + 1

        if (n_panels >= 3) {
            # The third panel contains the rank gradient
            # We modify it by adding a new scale
            tryCatch({
                # Rebuild the third panel with custom colors
                base_plot[[3]] <- base_plot[[3]] +
                    scale_fill_gradient(
                        low = COLOR_GFP_FADED,
                        high = COLOR_TES_FADED,
                        guide = "none"
                    )
            }, error = function(e) {
                # If modification fails, continue with original
                cat(sprintf("    Note: Could not modify gradient colors: %s\n", e$message))
            })
        }
    }

    # Add statistics inside the plot if requested
    if (add_stats && nes_inside) {
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)

        # Add NES annotation inside the enrichment score panel
        annotation_text <- sprintf("NES: %.3f\nFDR: %s", nes_val, qval_text)

        base_plot <- base_plot +
            plot_annotation(
                caption = annotation_text,
                theme = theme(
                    plot.caption = element_text(
                        hjust = 0.95, vjust = 0.95,
                        size = ANNOTATION_SIZE * 1.5,
                        face = "bold"
                    )
                )
            )
    } else if (add_stats && !nes_inside) {
        # Add statistics as inset table
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)

        stats_table <- data.frame(
            Metric = c("NES", "p-value", "FDR"),
            Value = c(as.character(nes_val), pval_text, qval_text)
        )

        table_grob <- gridExtra::tableGrob(
            stats_table,
            rows = NULL,
            theme = gridExtra::ttheme_minimal(
                core = list(fg_params = list(fontsize = 10)),
                colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
            )
        )

        base_plot <- base_plot + patchwork::inset_element(
            table_grob,
            left = 0.7, bottom = 0.7, right = 0.98, top = 0.95
        )
    }

    return(base_plot)
}

#-------------------------------------------------------------------------------
# Function to create and save enrichment plots for top pathways
#-------------------------------------------------------------------------------

save_gsea_plots <- function(gsea_result, name, n_plots = N_INDIVIDUAL_PLOTS) {

    if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
        cat(sprintf("No results for %s, skipping plots\n", name))
        return()
    }

    cat(sprintf("\nCreating plots for %s...\n", name))

    # Get significant pathways sorted by NES
    sig_results <- gsea_result@result %>%
        filter(p.adjust < FDR_CUTOFF) %>%
        arrange(desc(abs(NES)))

    if (nrow(sig_results) == 0) {
        cat(sprintf("  No significant pathways for %s\n", name))
        # Still create plots for top pathways by p-value
        sig_results <- gsea_result@result %>%
            arrange(pvalue) %>%
            head(n_plots)
    }

    n_to_plot <- min(n_plots, nrow(sig_results))

    # Also check for pathways of interest
    poi_in_collection <- intersect(PATHWAYS_OF_INTEREST, gsea_result@result$ID)
    if (length(poi_in_collection) > 0) {
        cat(sprintf("  Found %d pathways of interest in %s\n", length(poi_in_collection), name))
    }

    # Create individual enrichment plots
    for (i in 1:n_to_plot) {
        pathway_id <- sig_results$ID[i]
        pathway_name <- gsub("_", " ", sig_results$Description[i])

        # Truncate long names for filename
        safe_name <- gsub("[^A-Za-z0-9]", "_", substr(pathway_name, 1, 50))

        cat(sprintf("  %d. %s\n", i, substr(pathway_name, 1, 60)))

        # Create plot with title including comparison name
        plot_title <- paste0(pathway_name, "\n(", COMPARISON_NAME, ")")

        # Get statistics for this pathway
        pathway_stats <- sig_results[i, ]
        nes_val <- round(pathway_stats$NES, 3)
        pval <- pathway_stats$pvalue
        qval <- pathway_stats$p.adjust
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)

        # VERSION 1: With statistics (custom colors with TES/GFP gradient)
        p_custom <- create_gsea_plot_custom_colors(
            gsea_result, pathway_id,
            title = plot_title,
            add_stats = TRUE
        )

        if (!is.null(p_custom)) {
            ggsave(sprintf("plots/gsea_%s_%02d_%s.png", name, i, safe_name),
                   p_custom, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
                   bg = "white")
        }

        # VERSION 2: Clean version (no statistics - for manual annotation)
        p_clean <- create_gsea_plot_custom_colors(
            gsea_result, pathway_id,
            title = plot_title,
            add_stats = FALSE
        )

        if (!is.null(p_clean)) {
            ggsave(sprintf("plots/gsea_%s_%02d_%s_clean.png", name, i, safe_name),
                   p_clean, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
                   bg = "white")
        }

        # VERSION 3: Original gseaplot2 style (with stats table inset)
        p_orig <- gseaplot2(gsea_result, geneSetID = pathway_id,
                            title = plot_title,
                            color = COLOR_ENRICHMENT_LINE,
                            base_size = AXIS_TEXT_SIZE,
                            pvalue_table = FALSE)

        stats_table <- data.frame(
            Metric = c("NES", "p-value", "FDR"),
            Value = c(as.character(nes_val), pval_text, qval_text)
        )

        table_grob <- gridExtra::tableGrob(
            stats_table,
            rows = NULL,
            theme = gridExtra::ttheme_minimal(
                core = list(fg_params = list(fontsize = 10)),
                colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
            )
        )

        p_orig_annotated <- p_orig + patchwork::inset_element(
            table_grob,
            left = 0.7, bottom = 0.7, right = 0.98, top = 0.95
        )

        ggsave(sprintf("plots/gsea_%s_%02d_%s_original.png", name, i, safe_name),
               p_orig_annotated, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
               bg = "white")
    }

    #---------------------------------------------------------------------------
    # Create plots for pathways of interest (if they exist in this collection)
    #---------------------------------------------------------------------------

    for (poi in poi_in_collection) {
        # Check if already plotted in top N
        if (poi %in% sig_results$ID[1:n_to_plot]) {
            next
        }

        poi_data <- gsea_result@result[gsea_result@result$ID == poi, ]
        if (nrow(poi_data) == 0) next

        pathway_name <- gsub("_", " ", poi_data$Description)
        safe_name <- gsub("[^A-Za-z0-9]", "_", substr(pathway_name, 1, 50))

        cat(sprintf("  POI: %s (NES=%.3f, FDR=%.3f)\n",
                    substr(pathway_name, 1, 50), poi_data$NES, poi_data$p.adjust))

        plot_title <- paste0(pathway_name, "\n(", COMPARISON_NAME, ")")

        nes_val <- round(poi_data$NES, 3)
        pval <- poi_data$pvalue
        qval <- poi_data$p.adjust
        pval_text <- if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)
        qval_text <- if (qval < 0.001) "< 0.001" else sprintf("%.3f", qval)

        # Custom colors version with stats
        p_custom <- create_gsea_plot_custom_colors(
            gsea_result, poi,
            title = plot_title,
            add_stats = TRUE
        )

        if (!is.null(p_custom)) {
            ggsave(sprintf("plots/gsea_%s_POI_%s.png", name, safe_name),
                   p_custom, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
                   bg = "white")
        }

        # Clean version (no stats)
        p_clean <- create_gsea_plot_custom_colors(
            gsea_result, poi,
            title = plot_title,
            add_stats = FALSE
        )

        if (!is.null(p_clean)) {
            ggsave(sprintf("plots/gsea_%s_POI_%s_clean.png", name, safe_name),
                   p_clean, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
                   bg = "white")
        }

        # Original gseaplot2 style version
        p_orig <- gseaplot2(gsea_result, geneSetID = poi,
                            title = plot_title,
                            color = COLOR_ENRICHMENT_LINE,
                            base_size = AXIS_TEXT_SIZE,
                            pvalue_table = FALSE)

        stats_table <- data.frame(
            Metric = c("NES", "p-value", "FDR"),
            Value = c(as.character(nes_val), pval_text, qval_text)
        )

        table_grob <- gridExtra::tableGrob(
            stats_table,
            rows = NULL,
            theme = gridExtra::ttheme_minimal(
                core = list(fg_params = list(fontsize = 10)),
                colhead = list(fg_params = list(fontsize = 10, fontface = "bold"))
            )
        )

        p_orig_annotated <- p_orig + patchwork::inset_element(
            table_grob,
            left = 0.7, bottom = 0.7, right = 0.98, top = 0.95
        )

        ggsave(sprintf("plots/gsea_%s_POI_%s_original.png", name, safe_name),
               p_orig_annotated, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = DPI,
               bg = "white")
    }

    #---------------------------------------------------------------------------
    # Create combined plot for top up and down regulated pathways
    # Limit to 3 each direction for cleaner visualization
    #---------------------------------------------------------------------------

    # Top upregulated (positive NES) - limit to 3 for readability
    up_pathways <- gsea_result@result %>%
        filter(NES > 0) %>%
        arrange(pvalue) %>%
        head(3)

    # Top downregulated (negative NES) - limit to 3 for readability
    down_pathways <- gsea_result@result %>%
        filter(NES < 0) %>%
        arrange(pvalue) %>%
        head(3)

    if (nrow(up_pathways) > 0 || nrow(down_pathways) > 0) {
        top_pathways <- bind_rows(up_pathways, down_pathways)

        if (nrow(top_pathways) > 1) {
            # Create combined plot with comparison name in title
            combined_title <- paste0("Top Enriched Pathways (", COMPARISON_NAME, ")")

            p_multi <- gseaplot2(gsea_result,
                                 geneSetID = top_pathways$ID,
                                 title = combined_title,
                                 pvalue_table = FALSE,
                                 base_size = 11,
                                 rel_heights = c(1.5, 0.3, 0.5),
                                 subplots = 1:3)

            # Create custom stats table for combined plot (without gene set names)
            combined_stats <- data.frame(
                NES = round(top_pathways$NES, 3),
                `p-value` = ifelse(top_pathways$pvalue < 0.001, "< 0.001",
                                   sprintf("%.3f", top_pathways$pvalue)),
                FDR = ifelse(top_pathways$p.adjust < 0.001, "< 0.001",
                                   sprintf("%.3f", top_pathways$p.adjust)),
                check.names = FALSE
            )

            combined_table_grob <- gridExtra::tableGrob(
                combined_stats,
                rows = NULL,
                theme = gridExtra::ttheme_minimal(
                    core = list(fg_params = list(fontsize = 9)),
                    colhead = list(fg_params = list(fontsize = 9, fontface = "bold"))
                )
            )

            # Add table to plot
            p_multi <- p_multi + patchwork::inset_element(
                combined_table_grob,
                left = 0.72, bottom = 0.75, right = 0.98, top = 0.95
            )

            ggsave(sprintf("plots/gsea_%s_top_combined.png", name),
                   p_multi, width = 16, height = 14, dpi = DPI,
                   bg = "white")

            cat(sprintf("  Saved combined plot: gsea_%s_top_combined.png\n", name))
        }
    }

    #---------------------------------------------------------------------------
    # Create dotplot
    #---------------------------------------------------------------------------

    if (nrow(gsea_result@result) > 0) {
        # Dotplot for top pathways
        n_dot <- min(N_TOP_PATHWAYS, nrow(gsea_result@result))

        p_dot <- dotplot(gsea_result, showCategory = n_dot,
                         split = ".sign", font.size = 10) +
            facet_grid(.~.sign) +
            ggtitle(paste0("GSEA Dotplot (", COMPARISON_NAME, ")")) +
            theme_gsea +
            theme(axis.text.y = element_text(size = 9),
                  strip.text = element_text(size = 12, face = "bold"))

        ggsave(sprintf("plots/gsea_%s_dotplot.png", name),
               p_dot, width = 16, height = 12, dpi = DPI,
               bg = "white")

        cat(sprintf("  Saved dotplot: gsea_%s_dotplot.png\n", name))
    }

    #---------------------------------------------------------------------------
    # Create ridge plot
    #---------------------------------------------------------------------------

    if (nrow(gsea_result@result) >= 5) {
        n_ridge <- min(15, nrow(gsea_result@result))

        tryCatch({
            p_ridge <- ridgeplot(gsea_result, showCategory = n_ridge) +
                theme_gsea +
                theme(axis.text.y = element_text(size = 9))

            ggsave(sprintf("plots/gsea_%s_ridgeplot.png", name),
                   p_ridge, width = 14, height = 12, dpi = DPI,
                   bg = "white")

            cat(sprintf("  Saved ridgeplot: gsea_%s_ridgeplot.png\n", name))
        }, error = function(e) {
            cat(sprintf("  Could not create ridgeplot: %s\n", e$message))
        })
    }
}

# Create plots for each gene set collection
# Increase n_plots to 10 for main collections to capture more interesting pathways
if (!is.null(gsea_hallmark)) save_gsea_plots(gsea_hallmark, "hallmark", n_plots = 10)
if (!is.null(gsea_gobp)) save_gsea_plots(gsea_gobp, "GO_BP", n_plots = 5)
if (!is.null(gsea_reactome)) save_gsea_plots(gsea_reactome, "reactome", n_plots = 10)
if (!is.null(gsea_kegg)) save_gsea_plots(gsea_kegg, "KEGG", n_plots = 5)
if (!is.null(gsea_c6)) save_gsea_plots(gsea_c6, "C6_oncogenic", n_plots = 10)

#===============================================================================
# Create summary table
#===============================================================================

cat("\n=== Creating Summary ===\n")

create_summary <- function(gsea_result, name) {
    if (is.null(gsea_result)) return(NULL)

    gsea_result@result %>%
        filter(p.adjust < FDR_CUTOFF) %>%
        mutate(Collection = name) %>%
        dplyr::select(Collection, ID, Description, NES, pvalue, p.adjust, setSize) %>%
        arrange(pvalue)
}

summary_list <- list(
    create_summary(gsea_hallmark, "Hallmark"),
    create_summary(gsea_gobp, "GO:BP"),
    create_summary(gsea_reactome, "Reactome"),
    create_summary(gsea_kegg, "KEGG"),
    create_summary(gsea_c6, "C6_Oncogenic")
)

summary_df <- bind_rows(summary_list)

if (nrow(summary_df) > 0) {
    write.csv(summary_df, "gsea_significant_summary.csv", row.names = FALSE)
    cat(sprintf("Total significant pathways across all collections: %d\n", nrow(summary_df)))

    # Print top pathways
    cat("\nTop 10 Significant Pathways:\n")
    print(head(summary_df[, c("Collection", "Description", "NES", "p.adjust")], 10))
} else {
    cat("No significant pathways found at FDR < 0.25\n")
}

#===============================================================================
# Final summary
#===============================================================================

cat("\n", strrep("=", 60), "\n")
cat("=== GSEA Analysis Complete ===\n")
cat(strrep("=", 60), "\n")

cat("\nOutput files:\n")
cat("  - gsea_results_*.csv: Full GSEA results for each collection\n")
cat("  - gsea_results_*.rds: Cached GSEA objects for reuse/plotting\n")
cat("  - gsea_significant_summary.csv: Summary of significant pathways\n")
cat("\nPlot types (for each pathway):\n")
cat("  - gsea_*_##_*.png: Custom colors (TES/GFP gradient) with NES/FDR inside\n")
cat("  - gsea_*_##_*_clean.png: Custom colors, NO stats (for manual annotation)\n")
cat("  - gsea_*_##_*_original.png: Standard gseaplot2 style with stats table\n")
cat("  - gsea_*_POI_*.png: Pathways of interest (custom colors with stats)\n")
cat("  - gsea_*_dotplot.png: Dotplots by collection\n")
cat("  - gsea_*_ridgeplot.png: Ridge plots by collection\n")
cat("  - gsea_*_top_combined.png: Combined top pathway plots\n")
cat("\nColor scheme:\n")
cat("  - Gradient bar: TES (teal, left/high) to GFP (brown, right/low)\n")
cat("  - Matching heatmap colors for consistency\n")
cat("\nCollections analyzed:\n")
cat("  - Hallmark (H): Curated cancer hallmarks\n")
cat("  - GO:BP (C5): Gene Ontology Biological Process\n")
cat("  - Reactome (C2): Curated pathway database\n")
cat("  - KEGG (C2): KEGG pathways\n")
cat("  - C6 Oncogenic: Cancer-related signatures (YAP, MYC, RAS, etc.)\n")

cat("\nAnalysis parameters:\n")
cat(sprintf("  Gene set size range: %d - %d\n", MIN_GENE_SET_SIZE, MAX_GENE_SET_SIZE))
cat(sprintf("  Permutations: %d\n", N_PERMUTATIONS))
cat(sprintf("  FDR cutoff: %.2f\n", FDR_CUTOFF))

cat("\nUsage:\n")
cat("  Rscript gsea_analysis.R                # Default: plots only (use cached .rds results)\n")
cat("  Rscript gsea_analysis.R --compute      # Run GSEA computation (reuse cache if available)\n")
cat("  Rscript gsea_analysis.R --force        # Force recomputation of GSEA (ignore cache)\n")

