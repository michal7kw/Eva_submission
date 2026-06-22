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

# Read narrow peak files (TES and TEAD1 only - TESmut removed)
peak_files <- c(
    TES = "results/05_peaks_narrow/TES_peaks.narrowPeak",
    TEAD1 = "results/05_peaks_narrow/TEAD1_peaks.narrowPeak"
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

# Annotate narrow peaks
anno_list <- lapply(peak_list, annotatePeak, tssRegion = c(-3000, 3000), TxDb = txdb)

# Create output directory
dir.create("results/07_analysis_narrow", recursive = TRUE, showWarnings = FALSE)

# Plot annotation statistics
pdf("results/07_analysis_narrow/peak_annotation_barplot.pdf", width = 12, height = 8)
plotAnnoBar(anno_list)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/peak_annotation_barplot.png", width = 12, height = 8, units = "in", res = 300)
plotAnnoBar(anno_list)
dev.off()

# Distance to TSS plot
pdf("results/07_analysis_narrow/peak_distance_to_tss.pdf", width = 12, height = 8)
plotDistToTSS(anno_list)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/peak_distance_to_tss.png", width = 12, height = 8, units = "in", res = 300)
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
    write.table(common_genes, "results/07_analysis_narrow/TES_TEAD1_common_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(tes_specific, "results/07_analysis_narrow/TES_specific_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(tead1_specific, "results/07_analysis_narrow/TEAD1_specific_genes.txt",
        quote = FALSE, row.names = FALSE, col.names = FALSE
    )

    # GO enrichment for TES target genes
    if (length(tes_gene_ids) > 10) {
        ego_tes <- enrichGO(
            gene = na.omit(tes_gene_ids),
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05
        )

        if (!is.null(ego_tes) && nrow(ego_tes@result) > 0) {
            pdf("results/07_analysis_narrow/TES_GO_enrichment.pdf", width = 12, height = 10)
            print(dotplot(ego_tes, showCategory = 20) + ggtitle("TES GO Enrichment - Biological Process"))
            dev.off()

            png("results/07_analysis_narrow/TES_GO_enrichment.png", width = 12, height = 10, units = "in", res = 300)
            print(dotplot(ego_tes, showCategory = 20) + ggtitle("TES GO Enrichment - Biological Process"))
            dev.off()

            # Save GO results
            write.csv(ego_tes@result, "results/07_analysis_narrow/TES_GO_results.csv")
        } else {
            cat("No significant GO terms found for TES\n")
        }
    }

    # GO enrichment for TEAD1 target genes
    if (length(tead1_gene_ids) > 10) {
        ego_tead1 <- enrichGO(
            gene = na.omit(tead1_gene_ids),
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05
        )

        if (!is.null(ego_tead1) && nrow(ego_tead1@result) > 0) {
            pdf("results/07_analysis_narrow/TEAD1_GO_enrichment.pdf", width = 12, height = 10)
            print(dotplot(ego_tead1, showCategory = 20) + ggtitle("TEAD1 GO Enrichment - Biological Process"))
            dev.off()

            png("results/07_analysis_narrow/TEAD1_GO_enrichment.png", width = 12, height = 10, units = "in", res = 300)
            print(dotplot(ego_tead1, showCategory = 20) + ggtitle("TEAD1 GO Enrichment - Biological Process"))
            dev.off()

            # Save GO results
            write.csv(ego_tead1@result, "results/07_analysis_narrow/TEAD1_GO_results.csv")
        } else {
            cat("No significant GO terms found for TEAD1\n")
        }
    }

    # Comparative GO analysis between TES and TEAD1
    if (length(tes_gene_ids) > 10 && length(tead1_gene_ids) > 10) {
        gene_list <- list(TES = na.omit(tes_gene_ids), TEAD1 = na.omit(tead1_gene_ids))
        compGO <- compareCluster(gene_list,
            fun = "enrichGO",
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05
        )

        if (!is.null(compGO) && nrow(compGO@compareClusterResult) > 0) {
            pdf("results/07_analysis_narrow/TES_TEAD1_comparative_GO.pdf", width = 14, height = 10)
            print(dotplot(compGO, showCategory = 10) + ggtitle("TES vs TEAD1 GO Enrichment Comparison"))
            dev.off()

            png("results/07_analysis_narrow/TES_TEAD1_comparative_GO.png", width = 14, height = 10, units = "in", res = 300)
            print(dotplot(compGO, showCategory = 10) + ggtitle("TES vs TEAD1 GO Enrichment Comparison"))
            dev.off()

            # Save comparative results
            write.csv(compGO@compareClusterResult, "results/07_analysis_narrow/TES_TEAD1_comparative_GO.csv")
        }
    }
}

# Save results
# Save annotated peaks
for (name in names(anno_list)) {
    write.csv(
        as.data.frame(anno_list[[name]]),
        paste0("results/07_analysis_narrow/", name, "_peaks_annotated.csv")
    )
}

cat("Narrow peak annotation complete!\n")
cat("Results saved to results/07_analysis_narrow/\n")
