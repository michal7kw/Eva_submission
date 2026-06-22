library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(GenomeInfoDb)

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Load TxDb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read broad peak files (TESmut removed - failed sample)
peak_files <- c(
    TES = "results/05_peaks_broad/TES_peaks.broadPeak",
    TEAD1 = "results/05_peaks_broad/TEAD1_peaks.broadPeak"
)

peak_list <- lapply(peak_files, readPeakFile)

# Filter out empty peak files
peak_list <- peak_list[sapply(peak_list, function(x) length(x) > 0)]

# Fix chromosome naming - convert to UCSC format (add 'chr' prefix)
peak_list <- lapply(peak_list, function(peaks) {
    if (length(peaks) > 0) {
        # Harmonise chromosome naming and genome metadata with UCSC hg38
        GenomeInfoDb::seqlevelsStyle(peaks) <- "UCSC"
        GenomeInfoDb::genome(peaks) <- "hg38"
    }
    return(peaks)
})

# Annotate broad peaks
anno_list <- lapply(peak_list, annotatePeak, tssRegion = c(-3000, 3000), TxDb = txdb)

# Create output directory
dir.create("results/07_analysis_broad", recursive = TRUE, showWarnings = FALSE)

# Plot annotation statistics
pdf("results/07_analysis_broad/peak_annotation_barplot.pdf", width = 12, height = 8)
plotAnnoBar(anno_list)
dev.off()

# Convert to PNG
png("results/07_analysis_broad/peak_annotation_barplot.png", width = 12, height = 8, units = "in", res = 300)
plotAnnoBar(anno_list)
dev.off()

# Distance to TSS plot
pdf("results/07_analysis_broad/peak_distance_to_tss.pdf", width = 12, height = 8)
plotDistToTSS(anno_list)
dev.off()

# Convert to PNG
png("results/07_analysis_broad/peak_distance_to_tss.png", width = 12, height = 8, units = "in", res = 300)
plotDistToTSS(anno_list)
dev.off()

# Perform downstream analysis only if both TES and TEAD1 peaks are present
if ("TES" %in% names(anno_list) && "TEAD1" %in% names(anno_list)) {
    # Get gene lists
    tes_anno_df <- as.data.frame(anno_list[["TES"]])
    if ("geneId" %in% colnames(tes_anno_df) && nrow(tes_anno_df) > 0) {
        tes_gene_ids <- tes_anno_df$geneId
    } else {
        tes_gene_ids <- character(0)
    }

    tead1_anno_df <- as.data.frame(anno_list[["TEAD1"]])
    if ("geneId" %in% colnames(tead1_anno_df) && nrow(tead1_anno_df) > 0) {
        tead1_gene_ids <- tead1_anno_df$geneId
    } else {
        tead1_gene_ids <- character(0)
    }

    # Convert ENTREZ IDs to SYMBOLs
    if (length(tes_gene_ids) > 0) {
        tes_symbols <- bitr(tes_gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        tes_genes <- tes_symbols$SYMBOL
    } else {
        tes_genes <- character(0)
    }

    if (length(tead1_gene_ids) > 0) {
        tead1_symbols <- bitr(tead1_gene_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
        tead1_genes <- tead1_symbols$SYMBOL
    } else {
        tead1_genes <- character(0)
    }

    # Remove NAs from gene lists
    tes_genes <- tes_genes[!is.na(tes_genes)]
    tead1_genes <- tead1_genes[!is.na(tead1_genes)]

    # Find overlapping genes
    common_genes <- intersect(tes_genes, tead1_genes)
    tes_specific <- setdiff(tes_genes, tead1_genes)
    tead1_specific <- setdiff(tead1_genes, tes_genes)

    # Save gene lists
    write.table(common_genes, "results/07_analysis_broad/TES_TEAD1_common_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(tes_specific, "results/07_analysis_broad/TES_specific_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(tead1_specific, "results/07_analysis_broad/TEAD1_specific_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # GO enrichment for TES target genes
    ego <- enrichGO(
        gene = na.omit(tes_gene_ids),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    if (!is.null(ego)) {
        pdf("results/07_analysis_broad/TES_GO_enrichment.pdf", width = 10, height = 8)
        dotplot(ego, showCategory = 20)
        dev.off()

        png("results/07_analysis_broad/TES_GO_enrichment.png", width = 10, height = 8, units = "in", res = 300)
        dotplot(ego, showCategory = 20)
        dev.off()
    }
}

# Save results - standardized naming
# Save annotated peaks
for (name in names(anno_list)) {
    write.csv(
        as.data.frame(anno_list[[name]]),
        paste0("results/07_analysis_broad/", name, "_peaks_annotated.csv")
    )
}

cat("Broad peak annotation complete!\n")
cat("Results saved to results/07_analysis_broad/\n")
