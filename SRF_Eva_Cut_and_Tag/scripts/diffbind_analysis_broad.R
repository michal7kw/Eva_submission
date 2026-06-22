library(DiffBind)
library(ggplot2)
library(pheatmap)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Create sample sheet for broad peaks (TESmut removed - failed sample)
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
    c(
      "TES-1", "TES-2", "TES-3",
      "TEAD1-1", "TEAD1-2", "TEAD1-3"
    ),
    "_filtered.bam"
  ),
  ControlID = c(rep("IggMs", 3), rep("IggRb", 3)),
  bamControl = c(
    rep("results/04_filtered/IggMs_filtered.bam", 3),
    rep("results/04_filtered/IggRb_filtered.bam", 3)
  ),
  Peaks = paste0(
    "results/05_peaks_broad/",
    c(
      "TES-1", "TES-2", "TES-3",
      "TEAD1-1", "TEAD1-2", "TEAD1-3"
    ),
    "_peaks.broadPeak"
  ),
  PeakCaller = rep("macs", 6)
)

# Create output directory
dir.create("results/07_analysis_broad", recursive = TRUE, showWarnings = FALSE)

# Save sample sheet
write.csv(samples, "config/diffbind_samplesheet_broad.csv", row.names = FALSE, quote = FALSE)

# Load data
dba_obj <- dba(sampleSheet = "config/diffbind_samplesheet_broad.csv")

# Count reads
dba_obj <- dba.count(dba_obj, bParallel = TRUE)

# Normalize
dba_obj <- dba.normalize(dba_obj)

# TESmut removed - failed sample (no TES vs TESmut contrast)

# Establish contrast: TES vs TEAD1
dba_obj <- dba.contrast(dba_obj,
  group1 = dba_obj$masks$TES,
  group2 = dba_obj$masks$TEAD1,
  name1 = "TES", name2 = "TEAD1"
)

# Perform differential analysis
dba_obj <- dba.analyze(dba_obj)

# Get results (TESmut removed - only TES vs TEAD1 contrast)
res_tes_vs_tead1 <- dba.report(dba_obj, contrast = 1, th = 1)

# Save results - standardized naming to match narrow peaks
write.csv(
  as.data.frame(res_tes_vs_tead1),
  "results/07_analysis_broad/DiffBind_TES_vs_TEAD1.csv"
)

# Plot MA - standardized naming (TESmut removed - only 1 contrast)
pdf("results/07_analysis_broad/DiffBind_MA_plots.pdf", width = 8, height = 6)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

# Convert to PNG
png("results/07_analysis_broad/DiffBind_MA_plots.png", width = 8, height = 6, units = "in", res = 300)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

# Plot PCA - standardized naming
pdf("results/07_analysis_broad/DiffBind_PCA.pdf", width = 8, height = 8)
dba.plotPCA(dba_obj, label = DBA_ID)
dev.off()

# Convert to PNG
png("results/07_analysis_broad/DiffBind_PCA.png", width = 8, height = 8, units = "in", res = 300)
dba.plotPCA(dba_obj, label = DBA_ID)
dev.off()

# Plot correlation heatmap - standardized naming
pdf("results/07_analysis_broad/DiffBind_correlation.pdf", width = 10, height = 10)
dba.plotHeatmap(dba_obj, correlations = TRUE)
dev.off()

# Convert to PNG
png("results/07_analysis_broad/DiffBind_correlation.png", width = 10, height = 10, units = "in", res = 300)
dba.plotHeatmap(dba_obj, correlations = TRUE)
dev.off()

# Venn diagram removed - only 1 contrast now (TESmut removed)

print("Broad peak DiffBind analysis complete!")
cat("Results saved to results/07_analysis_broad/\n")
