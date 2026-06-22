library(DiffBind)
library(ggplot2)
library(pheatmap)

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")

# Create sample sheet for narrow peaks (TES and TEAD1 only - TESmut removed)
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
    "results/05_peaks_narrow/",
    c(
      "TES-1", "TES-2", "TES-3",
      "TEAD1-1", "TEAD1-2", "TEAD1-3"
    ),
    "_peaks.narrowPeak"
  ),
  PeakCaller = rep("macs", 6)
)

# Create output directory
dir.create("results/07_analysis_narrow", recursive = TRUE, showWarnings = FALSE)

# Save sample sheet
write.csv(samples, "config/diffbind_samplesheet.csv", row.names = FALSE, quote = FALSE)

# Load data
dba_obj <- dba(sampleSheet = "config/diffbind_samplesheet.csv")

# ============================================================================
# Apply ENCODE blacklist and generate greylists from IgG controls
# ============================================================================
# Greylists identify sample-specific artifact regions from IgG control data,
# complementing the ENCODE blacklist which covers known problematic regions.
# DiffBind v3+ supports this via dba.blacklist() with greylist=TRUE.

cat("Applying blacklist and greylist filtering...\n")

# Try to apply blacklist + greylist; fall back gracefully if GreyListChIP
# is not installed or if the BSgenome for hg38 is not available
tryCatch({
    # Apply ENCODE blacklist for hg38 and generate greylists from IgG controls
    dba_obj <- dba.blacklist(dba_obj, blacklist = DBA_BLACKLIST_HG38, greylist = TRUE)

    # Report how many peaks were removed
    cat("Blacklist/greylist filtering applied successfully\n")

    # Save greylist regions for reference
    greylist_dir <- "results/07_analysis_narrow/greylists"
    dir.create(greylist_dir, recursive = TRUE, showWarnings = FALSE)

    # Export greylist if available
    if (!is.null(dba_obj$greylist)) {
        greylist_gr <- dba_obj$greylist$master
        if (!is.null(greylist_gr) && length(greylist_gr) > 0) {
            rtracklayer::export(greylist_gr,
                                file.path(greylist_dir, "combined_greylist.bed"))
            cat(sprintf("  Greylist regions saved: %d regions\n", length(greylist_gr)))
        }
    }
}, error = function(e) {
    cat("WARNING: Greylist generation failed:", conditionMessage(e), "\n")
    cat("Falling back to ENCODE blacklist only...\n")
    tryCatch({
        dba_obj <<- dba.blacklist(dba_obj, blacklist = DBA_BLACKLIST_HG38,
                                   greylist = FALSE)
        cat("ENCODE blacklist applied successfully\n")
    }, error = function(e2) {
        cat("WARNING: Blacklist filtering also failed:", conditionMessage(e2), "\n")
        cat("Proceeding without blacklist/greylist filtering\n")
    })
})

# Count reads
dba_obj <- dba.count(dba_obj, bParallel = TRUE)

# Normalize
dba_obj <- dba.normalize(dba_obj)

# Establish contrast: TES vs TEAD1 (only comparison now)
dba_obj <- dba.contrast(dba_obj,
  group1 = dba_obj$masks$TES,
  group2 = dba_obj$masks$TEAD1,
  name1 = "TES", name2 = "TEAD1"
)

# Perform differential analysis
dba_obj <- dba.analyze(dba_obj)

# Get results
res_tes_vs_tead1 <- dba.report(dba_obj, contrast = 1, th = 1)

# Save results
write.csv(
  as.data.frame(res_tes_vs_tead1),
  "results/07_analysis_narrow/DiffBind_TES_vs_TEAD1.csv"
)

# Plot MA
pdf("results/07_analysis_narrow/DiffBind_MA_plots.pdf", width = 8, height = 6)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/DiffBind_MA_plots.png", width = 8, height = 6, units = "in", res = 300)
dba.plotMA(dba_obj, contrast = 1)
dev.off()

# Plot PCA
pdf("results/07_analysis_narrow/DiffBind_PCA.pdf", width = 8, height = 8)
dba.plotPCA(dba_obj, label = DBA_ID)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/DiffBind_PCA.png", width = 8, height = 8, units = "in", res = 300)
dba.plotPCA(dba_obj, label = DBA_ID)
dev.off()

# Plot correlation heatmap
pdf("results/07_analysis_narrow/DiffBind_correlation.pdf", width = 10, height = 10)
dba.plotHeatmap(dba_obj, correlations = TRUE)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/DiffBind_correlation.png", width = 10, height = 10, units = "in", res = 300)
dba.plotHeatmap(dba_obj, correlations = TRUE)
dev.off()

# Venn diagram showing peak overlaps between TES and TEAD1
pdf("results/07_analysis_narrow/DiffBind_venn.pdf", width = 8, height = 8)
dba.plotVenn(dba_obj, dba_obj$masks$TES | dba_obj$masks$TEAD1)
dev.off()

# Convert to PNG
png("results/07_analysis_narrow/DiffBind_venn.png", width = 8, height = 8, units = "in", res = 300)
dba.plotVenn(dba_obj, dba_obj$masks$TES | dba_obj$masks$TEAD1)
dev.off()

print("Narrow peak DiffBind analysis complete!")
cat("Results saved to results/07_analysis_narrow/\n")
cat("Note: TESmut samples have been removed from analysis\n")
