#!/usr/bin/env Rscript
#
# ==============================================================================
# TRUE GSEA ANALYSIS: Rank-based Gene Set Enrichment Analysis
# ==============================================================================
#
# Gene Set Enrichment Analysis (GSEA) determines whether a set of genes (e.g.,
# a pathway) shows statistically significant differences between two biological
# states (e.g., TES vs GFP).
#
# DIFFERENCE FROM STANDARD ENRICHMENT (ORA):
# - Standard "Over-Representation Analysis" (ORA) takes a LIST of significant genes
#   (e.g., DEGs) and asks "Are apoptosis genes over-represented in this list?"
# - GSEA takes ALL genes, ranked by a metric (usually fold change). It asks:
#   "Do apoptosis genes cluster at the top or bottom of this sorted list?"
#   This allows detecting pathways where many genes change slightly in the same
#   direction, even if individual genes aren't significant.
#
# This script implements "True" GSEA using the entire transcriptome ranked by
# log2FoldChange.
#
# ==============================================================================

# %%
suppressPackageStartupMessages({
  # Core Analysis Packages
  library(fgsea) # Fast Gene Set Enrichment Analysis (main algorithm)
  library(dplyr) # specific data manipulation (filter, select, mutate)
  library(tidyr) # Data tidying
  library(readr) # Fast CSV/TSV reading used for loading large datasets
  library(stringr) # String manipulation

  # Visualization Packages
  library(ggplot2) # Main plotting engine
  library(ggrepel) # Handles non-overlapping text labels in plots

  # Annotation Databases
  # These provide the mapping between IDs (Ensembl, Entrez) and meaningful names,
  # as well as the Gene Ontology (GO) definitions.
  library(org.Hs.eg.db) # Human genome annotation database
  library(GO.db) # Gene Ontology database (definitions of biological processes)
  library(clusterProfiler) # Used here primarily for some helper functions
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")

cat("=== TRUE GSEA ANALYSIS: Rank-Based Enrichment ===\n")
cat("Analysis started:", as.character(Sys.time()), "\n\n")

# Create output directory
output_dir <- "output/01_true_gsea_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PHASE 1: LOAD RNA-SEQ DATA AND CREATE RANKED GENE LIST
# =============================================================================
# GSEA requires a *ranked list* of genes. The most common ranking metric is
# the "WALD statistic" or "log2FoldChange".
# We use log2FoldChange here:
# - Positive values = Upregulated in TES
# - Negative values = Downregulated in TES
# =============================================================================

cat("=== PHASE 1: Creating Ranked Gene List ===\n")

# %%
# Load complete RNA-seq results
# IMPORTANT: We load ALL genes, not just the significant ones.
# GSEA needs the "background" of unchanged genes to calculate enrichment scores.
rna_all <- read.delim("../../../SRF_Eva_RNA/results/05_deseq2/deseq2_results_TES_vs_GFP.txt",
  stringsAsFactors = FALSE
)

cat(sprintf("✓ Loaded %d genes from RNA-seq\n", nrow(rna_all)))

# Clean and prepare for ranking
# We need to map Ensembl IDs (ENSG...) to Gene Symbols (e.g., TP53) because
# most pathway databases use Symbols.
rna_all$ensembl_id <- gsub("\\..*", "", rna_all$gene_id)

head(rna_all)

# Get gene symbols
# mapIds is the standard function to translate between ID types.
rna_all$gene_symbol <- mapIds(org.Hs.eg.db,
  keys = rna_all$ensembl_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

head(rna_all)

# %%
# Remove genes without symbols or log2FC
# We cannot rank genes if they don't have a value or a name.
rna_ranked <- rna_all %>%
  filter(!is.na(gene_symbol) & !is.na(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

cat(sprintf("✓ Ranked %d genes with valid symbols and fold changes\n", nrow(rna_ranked)))

# Create ranked list (gene symbol → log2FC)
# The 'fgsea' function expects a named numeric vector:
#   Values = ranking metric (log2FoldChange)
#   Names = gene IDs (Symbols)
gene_ranks <- setNames(rna_ranked$log2FoldChange, rna_ranked$gene_symbol)

# Remove duplicates (keep first = highest ranking)
# Sometimes multiple Ensembl IDs map to the same Symbol. We keep the one
# with the strongest signal (or just the first one) to avoid duplicates.
gene_ranks <- gene_ranks[!duplicated(names(gene_ranks))]

cat(sprintf("✓ Final ranked list: %d unique genes\n", length(gene_ranks)))
cat(sprintf("  Range: %.2f to %.2f log2FC\n", min(gene_ranks), max(gene_ranks)))
cat(sprintf("  Most upregulated: %s (%.2f)\n", names(gene_ranks)[1], gene_ranks[1]))
cat(sprintf(
  "  Most downregulated: %s (%.2f)\n\n", names(gene_ranks)[length(gene_ranks)],
  gene_ranks[length(gene_ranks)]
))

# %%

# =============================================================================
# PHASE 2: PREPARE GENE SETS FOR GSEA
# =============================================================================
# A "Gene Set" is just a named list of gene symbols representing a biological concept.
# We will use two sources:
# 1. GO Biological Process (BP): Standard comprehensive database.
# 2. Custom Cancer Sets: Manually defined for this specific research question.
# =============================================================================

cat("=== PHASE 2: Preparing Gene Sets ===\n")

# Load GO gene sets from org.Hs.eg.db
cat("Loading GO Biological Process gene sets...\n")

# Use AnnotationDbi to extract all genes associated with "BP" (Biological Process)
go_bp <- AnnotationDbi::select(org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "GOALL"),
  columns = c("SYMBOL", "ONTOLOGYALL"),
  keytype = "GOALL"
) %>%
  filter(ONTOLOGYALL == "BP") %>%
  filter(!is.na(SYMBOL))

# Convert to list format for fgsea: List of vectors, where list name = Pathway ID
go_gene_sets <- split(go_bp$SYMBOL, go_bp$GOALL)

# Remove very small and very large gene sets
# - Too small (<10): Statistical noise, hard to be significant.
# - Too large (>500): Too broad to be biologically meaningful (e.g., "cellular process").
go_gene_sets <- go_gene_sets[sapply(go_gene_sets, length) >= 10 &
  sapply(go_gene_sets, length) <= 500]

cat(sprintf("✓ Loaded %d GO BP gene sets (10-500 genes each)\n", length(go_gene_sets)))

# Create cancer-focused gene sets manually
# These acts as "positive controls" or specific hypotheses we want to test directly.
cancer_gene_sets <- list(
  APOPTOSIS = c(
    "BAX", "BAK1", "BID", "BIM", "PUMA", "NOXA", "CASP3", "CASP8", "CASP9",
    "FAS", "FASL", "TNFRSF10A", "TNFRSF10B", "TP53", "BCL2", "BCL2L1",
    "MCL1", "APAF1", "CYCS", "DIABLO"
  ),
  ANTI_APOPTOSIS = c(
    "BCL2", "BCL2L1", "BCL2L2", "MCL1", "BCL2A1", "BCLW",
    "BIRC2", "BIRC3", "BIRC5", "XIAP", "NAIP", "CFLIP"
  ),
  CELL_CYCLE = c(
    "CDK1", "CDK2", "CDK4", "CDK6", "CCNA1", "CCNA2", "CCNB1", "CCNB2",
    "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2", "E2F1", "E2F2", "E2F3",
    "RB1", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B"
  ),
  MIGRATION = c(
    "ITGB1", "ITGB3", "ITGB5", "ITGA5", "ITGAV", "CDH1", "CDH2",
    "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "MMP2", "MMP9",
    "TIMP1", "TIMP2", "CXCR4", "CXCL12"
  ),
  HIPPO_YAP = c(
    "YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2",
    "STK3", "STK4", "MOB1A", "MOB1B", "SAV1", "NF2", "AMOT", "AMOTL1",
    "AMOTL2", "PTPN14", "FAT1", "FAT2", "FAT3", "FAT4"
  ),
  ANGIOGENESIS = c(
    "VEGFA", "VEGFB", "VEGFC", "VEGFD", "FLT1", "KDR", "FLT4",
    "ANGPT1", "ANGPT2", "TEK", "PDGFA", "PDGFB", "PDGFRA", "PDGFRB",
    "FGF2", "FGFR1", "FGFR2", "HIF1A", "EPAS1"
  ),
  EMT = c(
    "CDH1", "CDH2", "VIM", "FN1", "SNAI1", "SNAI2", "SLUG", "TWIST1", "TWIST2",
    "ZEB1", "ZEB2", "GSC", "FOXC2", "TCF3", "TCF4"
  ),
  GLIOBLASTOMA_CORE = c(
    "EGFR", "PTEN", "TP53", "CDKN2A", "CDKN2B", "NF1", "RB1",
    "PIK3CA", "PIK3R1", "PDGFRA", "MET", "BRAF", "IDH1", "IDH2",
    "ATRX", "H3F3A", "TERT", "MDM2", "MDM4"
  )
)

# Combine GO and custom gene sets into one large collection
all_gene_sets <- c(go_gene_sets, cancer_gene_sets)

cat(sprintf("✓ Total gene sets for GSEA: %d\n", length(all_gene_sets)))
cat(sprintf("  - GO BP: %d\n", length(go_gene_sets)))
cat(sprintf("  - Custom cancer sets: %d\n\n", length(cancer_gene_sets)))

# %%

# =============================================================================
# PHASE 3: RUN FGSEA
# =============================================================================
# The `fgsea` function calculates:
# 1. Enrichment Score (ES): How much the gene set is overrepresented at the top or bottom of the list.
# 2. Normalized Enrichment Score (NES): ES normalized for gene set size.
#    - NES > 0: Pathway is UPREGULATED (genes found at the top of the list)
#    - NES < 0: Pathway is DOWNREGULATED (genes found at the bottom)
# 3. p-value & padj: Statistical significance.
# =============================================================================

cat("=== PHASE 3: Running Gene Set Enrichment Analysis ===\n")
cat("This may take several minutes...\n\n")

set.seed(42) # For reproducibility of permutations

# Run fgsea
fgsea_results <- fgsea(
  pathways = all_gene_sets,
  stats = gene_ranks, # The named vector of log2FC
  minSize = 10, # Ignore sets smaller than this
  maxSize = 500, # Ignore sets larger than this
  nproc = 8, # Parallel processing
  nPermSimple = 10000 # Number of permutations for p-value calculation
)

# Filter significant results
# We use Adjusted P-value (padj) to control for False Discovery Rate.
fgsea_sig <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

cat(sprintf("✓ GSEA complete!\n"))
cat(sprintf("  Total gene sets tested: %d\n", nrow(fgsea_results)))
cat(sprintf("  Significant gene sets (FDR < 0.05): %d\n", nrow(fgsea_sig)))
cat(sprintf("  Upregulated (positive NES): %d\n", sum(fgsea_sig$NES > 0)))
cat(sprintf("  Downregulated (negative NES): %d\n\n", sum(fgsea_sig$NES < 0)))

# %%

# =============================================================================
# PHASE 4: MAP GO IDs TO DESCRIPTIONS AND FILTER FOR CANCER-RELEVANT PATHWAYS
# =============================================================================
# Raw GO IDs (e.g., GO:0006915) are hard to read. We maps them to descriptions
# (e.g., "apoptotic process").
# Then, we filter the huge list of results to find only those relevant to cancer,
# using a keyword search text-mining approach.
# =============================================================================

cat("=== PHASE 4: Mapping GO IDs to Descriptions ===\n")

# First, map GO IDs to their human-readable descriptions
go_ids <- fgsea_sig$pathway[grepl("^GO:", fgsea_sig$pathway)]
cat(sprintf("  Found %d GO terms to map\n", length(go_ids)))

if (length(go_ids) > 0) {
  # Get GO term descriptions from GO.db
  go_descriptions <- AnnotationDbi::select(GO.db,
    keys = go_ids,
    columns = "TERM",
    keytype = "GOID"
  )

  # Create a lookup table (ID -> Term)
  go_desc_lookup <- setNames(go_descriptions$TERM, go_descriptions$GOID)

  # Add description column to fgsea_sig
  # If it's a GO term, look it up; otherwise (custom sets), keep the original name.
  fgsea_sig$description <- ifelse(grepl("^GO:", fgsea_sig$pathway),
    go_desc_lookup[fgsea_sig$pathway],
    fgsea_sig$pathway
  )

  cat(sprintf(
    "  Successfully mapped %d GO terms to descriptions\n",
    sum(!is.na(fgsea_sig$description))
  ))
} else {
  fgsea_sig$description <- fgsea_sig$pathway
}

# Also add descriptions to all results (not just significant ones) for export
go_ids_all <- fgsea_results$pathway[grepl("^GO:", fgsea_results$pathway)]
if (length(go_ids_all) > 0) {
  go_descriptions_all <- AnnotationDbi::select(GO.db,
    keys = go_ids_all,
    columns = "TERM",
    keytype = "GOID"
  )
  go_desc_lookup_all <- setNames(go_descriptions_all$TERM, go_descriptions_all$GOID)
  fgsea_results$description <- ifelse(grepl("^GO:", fgsea_results$pathway),
    go_desc_lookup_all[fgsea_results$pathway],
    fgsea_results$pathway
  )
} else {
  fgsea_results$description <- fgsea_results$pathway
}

cat("\n=== Filtering Cancer-Relevant Pathways ===\n")

# Define keywords - EXPANDED for comprehensive cancer pathway capture
# We define a broad dictionary of terms related to cancer biology.
cancer_keywords <- c(
  # Cell death pathways
  "apoptosis", "apoptotic", "cell death", "programmed cell death",
  "necrosis", "necrotic", "ferroptosis", "pyroptosis", "anoikis",
  "autophagy", "autophagic", "survival", "viability", "senescence",

  # Cell cycle and proliferation
  "proliferation", "proliferative", "cell cycle", "mitosis", "mitotic",
  "cell division", "growth", "G1/S", "G2/M", "S phase", "M phase",
  "DNA replication", "chromosome", "spindle", "cytokinesis",
  "cyclin", "checkpoint", "DNA repair", "DNA damage",

  # Migration and invasion
  "migration", "migratory", "motility", "invasion", "invasive",
  "chemotaxis", "chemotactic", "cell movement", "locomotion",
  "adhesion", "cell adhesion", "focal adhesion", "cell junction",
  "cytoskeleton", "actin", "tubulin", "microtubule",

  # Angiogenesis and vasculature
  "angiogenesis", "angiogenic", "blood vessel", "vasculature",
  "endothelial", "VEGF", "vascular",

  # Signaling pathways
  "signaling", "signal transduction", "kinase", "phosphorylation",
  "growth factor", "receptor", "activation", "cascade",

  # Metabolism
  "glycolysis", "metabolism", "metabolic", "glucose", "ATP",
  "oxidative", "respiration", "biosynthesis", "catabolic",

  # Transcription and chromatin
  "transcription", "gene expression", "chromatin", "histone",
  "RNA processing", "splicing", "translation",

  # Cancer-specific terms
  "tumor", "cancer", "oncogenic", "transformation",
  "EMT", "epithelial", "mesenchymal", "stemness",
  "Hippo", "YAP", "TEAD", "Wnt", "Notch", "Hedgehog",

  # Stress response
  "stress", "oxidative stress", "hypoxia", "ER stress",
  "unfolded protein", "heat shock", "inflammatory"
)

# Convert keywords to a regex pattern "term1|term2|term3"
pattern <- paste(cancer_keywords, collapse = "|")

# Filter using the DESCRIPTION column (not pathway ID) for proper matching
fgsea_cancer <- fgsea_sig %>%
  filter(grepl(pattern, description, ignore.case = TRUE))

cat(sprintf("✓ Cancer-relevant pathways: %d\n", nrow(fgsea_cancer)))
cat(sprintf("  Upregulated: %d\n", sum(fgsea_cancer$NES > 0)))
cat(sprintf("  Downregulated: %d\n\n", sum(fgsea_cancer$NES < 0)))

# Categorize pathways into functional groups (priority order matters - first match wins)
# This assigns a high-level label ("Cell Death") to detailed terms ("Regulation of apoptosis")
# to simplify the visualization.
fgsea_cancer$category <- NA

# Cell death (highest priority for apoptosis-related terms)
fgsea_cancer$category[grepl("apoptosis|apoptotic|death|necrosis|ferroptosis|pyroptosis|anoikis|autophagy|survival|senescence",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Cell Death"

# Cell cycle and proliferation
fgsea_cancer$category[grepl("proliferation|cell cycle|mitosis|mitotic|division|growth|G1|G2|S phase|M phase|replication|chromosome|spindle|cyclin|checkpoint",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Proliferation"

# Migration and invasion
fgsea_cancer$category[grepl("migration|invasion|motility|chemotaxis|locomotion|adhesion|cytoskeleton|actin|tubulin",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Migration"

# Angiogenesis
fgsea_cancer$category[grepl("angiogenesis|angiogenic|blood vessel|vascular|vasculature|endothelial|VEGF",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Angiogenesis"

# Metabolism
fgsea_cancer$category[grepl("glycolysis|metabolism|metabolic|glucose|ATP|oxidative|respiration|biosynthesis",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Metabolism"

# Signaling pathways
fgsea_cancer$category[grepl("signaling|signal transduction|kinase|phosphorylation|growth factor|receptor|Hippo|YAP|TEAD|Wnt|Notch",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Signaling"

# Transcription and chromatin regulation
fgsea_cancer$category[grepl("transcription|gene expression|chromatin|histone|RNA processing|splicing",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Transcription"

# Stress response
fgsea_cancer$category[grepl("stress|hypoxia|ER stress|unfolded protein|heat shock|inflammatory",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "Stress Response"

# EMT and transformation
fgsea_cancer$category[grepl("EMT|epithelial|mesenchymal|transformation|tumor|cancer|oncogenic",
  fgsea_cancer$description,
  ignore.case = TRUE
)] <- "EMT/Transformation"

# %%

# =============================================================================
# PHASE 5: EXPORT RESULTS
# =============================================================================

cat("=== PHASE 5: Exporting Results ===\n")

# Convert list columns to character strings for CSV export
# 'leadingEdge' is a list of genes driving the enrichment (the "core enrichment" genes).
# We unlist it into a semicolon-separated string to save it in a CSV cell.
fgsea_results_export <- fgsea_results %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  select(pathway, description, everything()) # Put description after pathway

fgsea_sig_export <- fgsea_sig %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  select(pathway, description, everything())

fgsea_cancer_export <- fgsea_cancer %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
  select(pathway, description, everything())

# Export all results
write.csv(fgsea_results_export, file.path(output_dir, "fgsea_all_pathways.csv"), row.names = FALSE)
write.csv(fgsea_sig_export, file.path(output_dir, "fgsea_significant_pathways.csv"), row.names = FALSE)
write.csv(fgsea_cancer_export, file.path(output_dir, "fgsea_cancer_pathways.csv"), row.names = FALSE)

cat(sprintf("✓ Results exported:\n"))
cat(sprintf("  - All pathways: %d\n", nrow(fgsea_results_export)))
cat(sprintf("  - Significant pathways: %d\n", nrow(fgsea_sig_export)))
cat(sprintf("  - Cancer-relevant pathways: %d\n\n", nrow(fgsea_cancer_export)))

# %%

# =============================================================================
# PHASE 6: VISUALIZATIONS
# =============================================================================

cat("=== PHASE 6: Creating Visualizations ===\n")

# =============================================================================
# PLOT 1: TOP SIGNIFICANT PATHWAYS (ALL - no keyword filtering)
# =============================================================================
cat("Creating Plot 1: Top significant pathways (all pathways)...\n")

# Get top pathways from ALL significant results (no cancer keyword filtering)
top_all_up <- fgsea_sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(15)

top_all_down <- fgsea_sig %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(15)

top_all_pathways <- rbind(top_all_up, top_all_down)

if (nrow(top_all_pathways) > 0) {
  # Truncate long descriptions for display
  top_all_pathways$display_name <- ifelse(
    nchar(top_all_pathways$description) > 60,
    paste0(substr(top_all_pathways$description, 1, 57), "..."),
    top_all_pathways$description
  )

  pdf(file.path(output_dir, "01_top_pathways_barplot.pdf"), width = 14, height = 12)
  p1 <- ggplot(top_all_pathways, aes(x = reorder(display_name, NES), y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = "Top 30 Significant Pathways from GSEA",
      subtitle = "All significant pathways (FDR < 0.05) - Normalized Enrichment Score",
      x = "Pathway",
      y = "NES",
      fill = "Direction"
    ) +
    scale_fill_manual(
      values = c("FALSE" = "#1F78B4", "TRUE" = "#E31A1C"),
      labels = c("Downregulated in TES", "Upregulated in TES")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
  print(p1)
  dev.off()
  cat(sprintf("  ✓ Saved: 01_top_pathways_barplot.pdf (%d pathways)\n", nrow(top_all_pathways)))
}

# =============================================================================
# PLOT 1b: TOP CANCER-RELEVANT PATHWAYS (keyword filtered)
# =============================================================================
cat("Creating Plot 1b: Top cancer-relevant pathways...\n")

top_cancer_up <- fgsea_cancer %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  head(15)

top_cancer_down <- fgsea_cancer %>%
  filter(NES < 0) %>%
  arrange(NES) %>%
  head(15)

top_cancer_pathways <- rbind(top_cancer_up, top_cancer_down)

if (nrow(top_cancer_pathways) > 0) {
  # Truncate long descriptions for display
  top_cancer_pathways$display_name <- ifelse(
    nchar(top_cancer_pathways$description) > 60,
    paste0(substr(top_cancer_pathways$description, 1, 57), "..."),
    top_cancer_pathways$description
  )

  pdf(file.path(output_dir, "01b_cancer_pathways_barplot.pdf"), width = 14, height = 12)
  p1b <- ggplot(top_cancer_pathways, aes(x = reorder(display_name, NES), y = NES, fill = NES > 0)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    coord_flip() +
    labs(
      title = "Top Cancer-Relevant Pathways from GSEA",
      subtitle = "Filtered by cancer keywords - Normalized Enrichment Score",
      x = "Pathway",
      y = "NES",
      fill = "Direction"
    ) +
    scale_fill_manual(
      values = c("FALSE" = "#1F78B4", "TRUE" = "#E31A1C"),
      labels = c("Downregulated in TES", "Upregulated in TES")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
  print(p1b)
  dev.off()
  cat(sprintf("  ✓ Saved: 01b_cancer_pathways_barplot.pdf (%d pathways)\n", nrow(top_cancer_pathways)))
} else {
  cat("  ⚠ No cancer-relevant pathways found for barplot\n")
}

# Plot 2: Enrichment plots for specific pathways
# These plots show the "Mountain Plot" or "Random Walk", visualizing exactly
# WHERE the genes in the pathway fall in our ranked list.
cat("Creating detailed enrichment plots for key pathways...\n")

# Select interesting pathways
key_pathway_names <- c()
if ("APOPTOSIS" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "APOPTOSIS")
if ("CELL_CYCLE" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "CELL_CYCLE")
if ("MIGRATION" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "MIGRATION")
if ("HIPPO_YAP" %in% names(all_gene_sets)) key_pathway_names <- c(key_pathway_names, "HIPPO_YAP")

if (length(key_pathway_names) > 0) {
  pdf(file.path(output_dir, "02_detailed_enrichment_plots.pdf"), width = 12, height = 8)
  for (pathway_name in key_pathway_names) {
    p <- plotEnrichment(all_gene_sets[[pathway_name]], gene_ranks) +
      labs(title = paste("GSEA Enrichment Plot:", pathway_name)) +
      theme_bw(base_size = 14)
    print(p)
  }
  dev.off()
}

# Plot 3: Bubble plot of cancer pathways
# A multi-dimensional plot:
# - X axis: NES (Direction and strength)
# - Y axis: Significance (-log10 FDR)
# - Size: Pathway size (gene count)
# - Color: Biological Category
cat("Creating bubble plot...\n")

if (nrow(fgsea_cancer) > 0 && !all(is.na(fgsea_cancer$category))) {
  # Add text labels for top pathways
  fgsea_cancer_plot <- fgsea_cancer %>%
    filter(!is.na(category)) %>%
    mutate(short_desc = ifelse(nchar(description) > 40,
      paste0(substr(description, 1, 37), "..."),
      description
    ))

  pdf(file.path(output_dir, "03_cancer_pathways_bubble.pdf"), width = 14, height = 10)
  p3 <- ggplot(
    fgsea_cancer_plot,
    aes(x = NES, y = -log10(padj), color = category, size = size)
  ) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    # Add labels for top pathways
    ggrepel::geom_text_repel(
      data = fgsea_cancer_plot %>%
        group_by(category) %>%
        slice_max(abs(NES), n = 2) %>%
        ungroup(),
      aes(label = short_desc),
      size = 2.5,
      max.overlaps = 20,
      show.legend = FALSE
    ) +
    labs(
      title = "Cancer Pathway Enrichment (GSEA)",
      subtitle = "Bubble size = number of genes in pathway",
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(Adjusted p-value)",
      color = "Category",
      size = "Gene Set Size"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    ) +
    scale_color_brewer(palette = "Set1")
  print(p3)
  dev.off()
  cat(sprintf("  ✓ Saved: 03_cancer_pathways_bubble.pdf (%d pathways)\n", nrow(fgsea_cancer_plot)))
} else {
  cat("  ⚠ No cancer pathways with categories found for bubble plot\n")
}

# Plot 4: Category summary
# A simple bar chart showing how many pathways in each category are Up vs Down.
cat("Creating category summary plot...\n")

if (!all(is.na(fgsea_cancer$category))) {
  category_summary <- fgsea_cancer %>%
    filter(!is.na(category)) %>%
    mutate(direction = ifelse(NES > 0, "Up", "Down")) %>%
    group_by(category, direction) %>%
    summarise(count = n(), .groups = "drop")

  pdf(file.path(output_dir, "04_category_summary.pdf"), width = 10, height = 7)
  p4 <- ggplot(category_summary, aes(x = category, y = count, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    labs(
      title = "Cancer Pathway Categories (GSEA)",
      subtitle = "Number of significantly enriched pathways per category",
      x = "Category",
      y = "Number of Pathways",
      fill = "Direction"
    ) +
    scale_fill_manual(values = c("Up" = "#E31A1C", "Down" = "#1F78B4")) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  print(p4)
  dev.off()
}

cat("✓ All visualizations created\n\n")

# %%

# =============================================================================
# PHASE 7: SUMMARY REPORT
# =============================================================================

cat("=== PHASE 7: Generating Summary Report ===\n")

summary_file <- file.path(output_dir, "GSEA_SUMMARY.txt")
cat("TRUE GSEA ANALYSIS SUMMARY\n", file = summary_file)
cat("==========================\n\n", file = summary_file, append = TRUE)
cat(paste("Generated:", Sys.time(), "\n\n"), file = summary_file, append = TRUE)

cat("METHOD:\n", file = summary_file, append = TRUE)
cat("  - Algorithm: fgsea (Fast Gene Set Enrichment Analysis)\n", file = summary_file, append = TRUE)
cat("  - Ranking metric: log2FoldChange (TES vs GFP)\n", file = summary_file, append = TRUE)
cat(sprintf("  - Total genes ranked: %d\n", length(gene_ranks)), file = summary_file, append = TRUE)
cat(sprintf("  - Gene sets tested: %d\n", nrow(fgsea_results)), file = summary_file, append = TRUE)
cat(sprintf("  - Permutations: 10,000\n\n"), file = summary_file, append = TRUE)

cat("RESULTS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Significant pathways (FDR < 0.05): %d\n", nrow(fgsea_sig)), file = summary_file, append = TRUE)
cat(sprintf("    - Positively enriched (NES > 0): %d\n", sum(fgsea_sig$NES > 0)), file = summary_file, append = TRUE)
cat(sprintf("    - Negatively enriched (NES < 0): %d\n\n", sum(fgsea_sig$NES < 0)), file = summary_file, append = TRUE)

cat("CANCER-RELEVANT PATHWAYS:\n", file = summary_file, append = TRUE)
cat(sprintf("  Total cancer pathways: %d\n", nrow(fgsea_cancer)), file = summary_file, append = TRUE)

if (!all(is.na(fgsea_cancer$category))) {
  for (cat_name in unique(fgsea_cancer$category[!is.na(fgsea_cancer$category)])) {
    cat_pathways <- fgsea_cancer %>% filter(category == cat_name)
    cat(
      sprintf(
        "    %s: %d (%d up, %d down)\n",
        cat_name,
        nrow(cat_pathways),
        sum(cat_pathways$NES > 0),
        sum(cat_pathways$NES < 0)
      ),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n\nTOP 10 UPREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_up <- fgsea_sig %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(10)
  for (i in 1:min(nrow(top_up), 10)) {
    cat(sprintf("  %d. %s (NES=%.2f, FDR=%.2e)\n", i, top_up$pathway[i], top_up$NES[i], top_up$padj[i]),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n\nTOP 10 DOWNREGULATED PATHWAYS:\n", file = summary_file, append = TRUE)
if (nrow(fgsea_sig) > 0) {
  top_down <- fgsea_sig %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(10)
  for (i in 1:min(nrow(top_down), 10)) {
    cat(sprintf("  %d. %s (NES=%.2f, FDR=%.2e)\n", i, top_down$pathway[i], top_down$NES[i], top_down$padj[i]),
      file = summary_file, append = TRUE
    )
  }
}

cat("\n✓ Summary report saved\n\n")

cat("========================================\n")
cat("TRUE GSEA ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Completed:", as.character(Sys.time()), "\n")
cat(sprintf("Output directory: %s\n\n", output_dir))
cat("Key files:\n")
cat("  CSV Results (with GO term descriptions):\n")
cat("  - fgsea_all_pathways.csv (all tested pathways with descriptions)\n")
cat("  - fgsea_significant_pathways.csv (FDR < 0.05 with descriptions)\n")
cat("  - fgsea_cancer_pathways.csv (cancer-relevant subset with descriptions)\n")
cat("\n  Visualizations:\n")
cat("  - 01_top_pathways_barplot.pdf (TOP 30 significant pathways - ALL)\n")
cat("  - 01b_cancer_pathways_barplot.pdf (TOP cancer-relevant pathways only)\n")
cat("  - 02_detailed_enrichment_plots.pdf (enrichment curves)\n")
cat("  - 03_cancer_pathways_bubble.pdf (cancer pathways by category)\n")
cat("  - 04_category_summary.pdf (category counts)\n")
cat("  - GSEA_SUMMARY.txt\n\n")
cat("This is TRUE GSEA using ranked gene lists!\n")
cat("All genes contribute to enrichment score calculation.\n")
cat("GO terms now mapped to human-readable descriptions.\n")
