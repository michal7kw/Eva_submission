library(VennDiagram)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Peak counts from shell script (TESmut removed - failed sample)
tes_total=79841
tead1_total=61624

# Direct overlap counts (peaks with any overlap)
tes_overlap_tead1=28243
tead1_overlap_tes=28887

# For Venn diagrams, use overlap counts (peaks that have any overlap)
tes_tead1_overlap=28243

tes_unique_tead1=51598
tead1_unique_tes=32737

# Summit-based overlaps
summit_tes_tead1=26663

# Reciprocal overlaps (>=50%)
reciprocal_tes_tead1=15852

# Jaccard indices
jaccard_tes_tead1=0.249

cat(paste("Narrow Peak counts (TESmut removed - failed sample):\n"))
cat(paste("TES:", tes_total, "\n"))
cat(paste("TEAD1:", tead1_total, "\n\n"))

cat(paste("Direct Overlaps:\n"))
cat(paste("TES overlapping TEAD1:", tes_overlap_tead1, "\n\n"))

# Create Venn diagrams for pairwise comparisons (TESmut removed - failed sample)
output_dir <- "results/12_pairwise_overlap_narrow/venn"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# TES vs TEAD1 (only comparison now - TESmut removed)
cat("Creating TES vs TEAD1 Venn diagram...\n")
png(file.path(output_dir, "TES_vs_TEAD1_narrow_venn.png"),
    width = 3000, height = 3000, res = 300)
venn.plot <- draw.pairwise.venn(
    area1 = tes_total,
    area2 = tead1_total,
    cross.area = tes_tead1_overlap,
    category = c("TES", "TEAD1"),
    col = "transparent",
    fill = c("lightblue", "lightgreen"),
    alpha = 0.6,
    cex = 1.5,
    fontfamily = "serif",
    main = "TES vs TEAD1 Narrow Peak Overlap",
    main.cex = 1.8
)
grid.draw(venn.plot)
dev.off()

# Create comprehensive summary statistics table (TESmut removed)
stats_table <- data.frame(
    Comparison = c("TES vs TEAD1"),
    Direct_Overlap = c(tes_tead1_overlap),
    Summit_Overlap = c(summit_tes_tead1),
    Reciprocal_50pct = c(reciprocal_tes_tead1),
    TES_Only = c(tes_unique_tead1),
    TEAD1_Only = c(tead1_unique_tes),
    TES_Total = c(tes_total),
    TEAD1_Total = c(tead1_total),
    Overlap_Pct_TES = round(c(tes_tead1_overlap/tes_total*100), 1),
    Overlap_Pct_TEAD1 = round(c(tes_tead1_overlap/tead1_total*100), 1),
    Jaccard_Index = c(jaccard_tes_tead1)
)

# Write statistics to file
write.table(stats_table, file.path(output_dir, "narrow_peak_overlap_statistics.txt"),
            sep="\t", row.names=FALSE, quote=FALSE)

# Create bar plot comparing overlap methods
overlap_comparison <- data.frame(
    Comparison = rep(c("TES vs TEAD1"), 3),
    Method = c("Direct Overlap", "Summit-based (200bp)", "Reciprocal (50%)"),
    Count = c(tes_tead1_overlap, summit_tes_tead1, reciprocal_tes_tead1)
)

p_methods <- ggplot(overlap_comparison, aes(x = Comparison, y = Count, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparison of Overlap Methods for Narrow Peaks",
         y = "Number of Overlapping Peaks", x = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set2")

ggsave(file.path(output_dir, "overlap_methods_comparison.pdf"), p_methods, width = 10, height = 6)
ggsave(file.path(output_dir, "overlap_methods_comparison.png"), p_methods, width = 10, height = 6, dpi = 300)

cat("Narrow peak Venn diagram analysis completed (TESmut excluded - failed sample).\n")
cat("\nSummary Statistics:\n")
print(stats_table)

cat("\nInterpretation Guide:\n")
cat("- Direct Overlap: Any genomic overlap (most permissive)\n")
cat("- Summit-based: Peak centers within 200bp (appropriate for TF binding)\n")
cat("- Reciprocal 50%: At least 50% mutual overlap (most stringent)\n")
cat("- Jaccard Index: Overlap similarity score (0-1, higher = more similar)\n")
cat("- NOTE: TESmut samples removed from analysis (failed sample)\n")
