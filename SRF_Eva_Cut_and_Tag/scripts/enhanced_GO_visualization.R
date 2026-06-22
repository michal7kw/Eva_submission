#!/usr/bin/env Rscript
# Enhanced GO visualization for TES-specific vs TEAD1 differential binding analysis
# This script creates more informative and interpretable GO enrichment plots

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(forcats)
library(gridExtra)
library(enrichplot)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Read the TES-specific genes (enriched in TES but not TEAD1)
tes_specific_genes <- read.table("results/07_analysis_narrow/TES_specific_genes.txt",
    stringsAsFactors = FALSE
)$V1

# Convert gene symbols to ENTREZ IDs
tes_entrez <- bitr(tes_specific_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

print(paste("TES-specific genes:", length(tes_specific_genes)))
print(paste("Converted to ENTREZ:", nrow(tes_entrez)))

# Enhanced GO enrichment analysis for TES-specific genes
ego_tes_specific <- enrichGO(
    gene = tes_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
)

# Create output directory
dir.create("results/07_analysis_narrow/enhanced_GO", recursive = TRUE, showWarnings = FALSE)

if (!is.null(ego_tes_specific) && nrow(ego_tes_specific@result) > 0) {
    # 1. Enhanced dotplot with better formatting
    p1 <- dotplot(ego_tes_specific,
        showCategory = 15,
        font.size = 12
    ) +
        ggtitle("TES-Specific Pathways\n(Enriched in TES but not TEAD1)",
            subtitle = paste("Based on", length(tes_specific_genes), "TES-specific genes")
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.y = element_text(size = 11),
            axis.text.x = element_text(size = 10),
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        ) +
        scale_color_gradient(low = "blue", high = "red", name = "p.adjust") +
        xlab("Gene Ratio") + ylab("GO Terms")

    # 2. Barplot showing enrichment scores
    p2 <- barplot(ego_tes_specific,
        showCategory = 15
    ) +
        ggtitle("TES-Specific GO Term Enrichment Scores") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 9)
        ) +
        scale_fill_gradient(low = "lightblue", high = "darkred", name = "p.adjust") +
        xlab("Enrichment Score") + ylab("GO Terms")

    # 3. Enhanced comparison plot showing TES vs TEAD1 differences
    # First, load TEAD1 data for comparison
    tead1_specific_genes <- read.table("results/07_analysis_narrow/TEAD1_specific_genes.txt",
        stringsAsFactors = FALSE
    )$V1
    tead1_entrez <- bitr(tead1_specific_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    # GO enrichment for TEAD1-specific genes
    ego_tead1_specific <- enrichGO(
        gene = tead1_entrez$ENTREZID,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.01,
        qvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
    )

    # Comparative analysis
    if (!is.null(ego_tead1_specific) && nrow(ego_tead1_specific@result) > 0) {
        gene_list_specific <- list(
            "TES-specific" = tes_entrez$ENTREZID,
            "TEAD1-specific" = tead1_entrez$ENTREZID
        )

        comp_specific <- compareCluster(gene_list_specific,
            fun = "enrichGO",
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.01,
            qvalueCutoff = 0.05
        )

        # Enhanced comparative plot
        if (!is.null(comp_specific) && nrow(comp_specific@compareClusterResult) > 0) {
            p3 <- dotplot(comp_specific,
                showCategory = 12,
                font.size = 11
            ) +
                ggtitle("TES-specific vs TEAD1-specific Pathways",
                    subtitle = "Comparing unique binding targets of each transcription factor"
                ) +
                theme_classic() +
                theme(
                    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, size = 11),
                    axis.text.y = element_text(size = 10),
                    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
                    legend.title = element_text(size = 10)
                ) +
                scale_color_gradient(low = "#3B82F6", high = "#DC2626", name = "p.adjust") +
                xlab("Condition") + ylab("GO Terms")
        }
    }

    # 4. Create summary statistics plot
    tes_top_terms <- head(ego_tes_specific@result, 10)
    tes_top_terms$NegLog10P <- -log10(tes_top_terms$p.adjust)
    tes_top_terms$Description <- factor(tes_top_terms$Description,
        levels = rev(tes_top_terms$Description)
    )

    p4 <- ggplot(tes_top_terms, aes(
        x = NegLog10P, y = Description,
        fill = Count, size = as.numeric(sub("/.*", "", GeneRatio))
    )) +
        geom_point(shape = 21, color = "black", stroke = 0.3) +
        scale_fill_gradient(low = "#FEF3C7", high = "#DC2626", name = "Gene\nCount") +
        scale_size_continuous(range = c(3, 10), name = "Gene\nRatio") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.title = element_text(size = 10)
        ) +
        ggtitle("TES-Specific Pathway Significance and Gene Counts") +
        xlab("-log10(adjusted p-value)") +
        ylab("GO Terms")

    # Save all plots
    pdf("results/07_analysis_narrow/enhanced_GO/TES_specific_enhanced_dotplot.pdf", width = 14, height = 10)
    print(p1)
    dev.off()

    png("results/07_analysis_narrow/enhanced_GO/TES_specific_enhanced_dotplot.png",
        width = 14, height = 10, units = "in", res = 300
    )
    print(p1)
    dev.off()

    pdf("results/07_analysis_narrow/enhanced_GO/TES_specific_barplot.pdf", width = 12, height = 10)
    print(p2)
    dev.off()

    png("results/07_analysis_narrow/enhanced_GO/TES_specific_barplot.png",
        width = 12, height = 10, units = "in", res = 300
    )
    print(p2)
    dev.off()

    if (exists("p3")) {
        pdf("results/07_analysis_narrow/enhanced_GO/TES_vs_TEAD1_specific_comparison.pdf", width = 16, height = 12)
        print(p3)
        dev.off()

        png("results/07_analysis_narrow/enhanced_GO/TES_vs_TEAD1_specific_comparison.png",
            width = 16, height = 12, units = "in", res = 300
        )
        print(p3)
        dev.off()
    }

    pdf("results/07_analysis_narrow/enhanced_GO/TES_specific_significance_plot.pdf", width = 12, height = 8)
    print(p4)
    dev.off()

    png("results/07_analysis_narrow/enhanced_GO/TES_specific_significance_plot.png",
        width = 12, height = 8, units = "in", res = 300
    )
    print(p4)
    dev.off()

    # 5. Network visualization of GO terms
    if (nrow(ego_tes_specific@result) >= 5) {
        # Create GO term network
        p5 <- emapplot(pairwise_termsim(ego_tes_specific),
            showCategory = 20,
            cex_label_category = 0.8
        ) +
            ggtitle("TES-Specific GO Term Network\nClusters show related biological processes") +
            theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

        pdf("results/07_analysis_narrow/enhanced_GO/TES_specific_network.pdf", width = 14, height = 14)
        print(p5)
        dev.off()

        png("results/07_analysis_narrow/enhanced_GO/TES_specific_network.png",
            width = 14, height = 14, units = "in", res = 300
        )
        print(p5)
        dev.off()
    }

    # Save detailed results
    write.csv(ego_tes_specific@result,
        "results/07_analysis_narrow/enhanced_GO/TES_specific_detailed_GO_results.csv",
        row.names = FALSE
    )

    # Create summary report
    summary_stats <- data.frame(
        Metric = c(
            "Total TES-specific genes", "Genes with ENTREZ ID", "Significant GO terms (p.adj < 0.05)",
            "Most significant term", "Top enriched process"
        ),
        Value = c(
            length(tes_specific_genes), nrow(tes_entrez), nrow(ego_tes_specific@result),
            ego_tes_specific@result$Description[1], ego_tes_specific@result$Description[1]
        ),
        Details = c(
            "Genes found in TES but not TEAD1 peaks", "Successfully converted to ENTREZ",
            "After multiple testing correction", paste("p.adj =", format(ego_tes_specific@result$p.adjust[1], scientific = TRUE)),
            paste("Gene ratio =", ego_tes_specific@result$GeneRatio[1])
        )
    )

    write.csv(summary_stats,
        "results/07_analysis_narrow/enhanced_GO/TES_specific_analysis_summary.csv",
        row.names = FALSE
    )

    print("Enhanced GO analysis complete!")
    print(paste("Generated", nrow(ego_tes_specific@result), "significant GO terms for TES-specific genes"))
} else {
    print("No significant GO terms found for TES-specific genes")
}
