#!/bin/bash
# meDIP NOTE: reads DNA-methylation (meDIP) data from the ARCHIVED tree SRF_Eva_top/Archive/meDIP, which is NOT self-contained in Eva_submission. If the Archive is moved or removed, repoint or copy the meDIP inputs before running.
#SBATCH --job-name=a1_32_encode_enhancer
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/32_encode_enhancer_methylation.out
#SBATCH --error=logs/32_encode_enhancer_methylation.err

# =============================================================================
# ENCODE ENHANCER METHYLATION ANALYSIS
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE-defined enhancers to determine if
#          TES binding causes methylation at enhancer regions.
#
# Key Question: Do TES-bound enhancers show HYPERMETHYLATION (TES > GFP)?
#               Or do they show HYPOMETHYLATION like gene bodies?
#
# =============================================================================

echo "=========================================="
echo "ENCODE ENHANCER METHYLATION ANALYSIS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/32_encode_enhancer"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: DOWNLOAD ENCODE cCRE ANNOTATIONS
# =============================================================================

echo "=== STEP 1: Downloading ENCODE cCRE Annotations ==="
echo ""

CCRE_FILE="${OUTDIR}/GRCh38-cCREs.bed"

if [ ! -f "${CCRE_FILE}" ]; then
    echo "Downloading ENCODE cCRE annotations..."

    # Download BigBed from UCSC
    wget -q -O "${OUTDIR}/encodeCcreCombined.bb" \
        "http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb"

    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to download ENCODE cCRE annotations"
        exit 1
    fi

    # Download bigBedToBed if not present
    if [ ! -f "${OUTDIR}/bigBedToBed" ]; then
        wget -q -O "${OUTDIR}/bigBedToBed" \
            "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed"
        chmod +x "${OUTDIR}/bigBedToBed"
    fi

    # Convert BigBed to BED
    "${OUTDIR}/bigBedToBed" "${OUTDIR}/encodeCcreCombined.bb" "${CCRE_FILE}"

    echo "Downloaded and converted: ${CCRE_FILE}"
else
    echo "ENCODE cCRE file already exists: ${CCRE_FILE}"
fi

# Check file
N_CCRE=$(wc -l < "${CCRE_FILE}")
echo "Total cCREs: ${N_CCRE}"
echo ""

# =============================================================================
# STEP 2: FILTER FOR ENHANCER TYPES
# =============================================================================

echo "=== STEP 2: Filtering for Enhancer Types ==="
echo ""

# ENCODE cCRE types (column 11 in UCSC format):
# - PLS: Promoter-like signature
# - pELS: Proximal enhancer-like signature (within 2kb of TSS)
# - dELS: Distal enhancer-like signature (true enhancers, >2kb from TSS)
# - CTCF-only: CTCF binding only
# - DNase-H3K4me3: DNase + H3K4me3 (promoter-like)

# Extract distal enhancers (dELS) - these are the true enhancers
awk -F'\t' '$11 == "dELS"' "${CCRE_FILE}" > "${OUTDIR}/ENCODE_distal_enhancers.bed"
N_DELS=$(wc -l < "${OUTDIR}/ENCODE_distal_enhancers.bed")
echo "Distal enhancers (dELS): ${N_DELS}"

# Extract proximal enhancers (pELS) - near promoters
awk -F'\t' '$11 == "pELS"' "${CCRE_FILE}" > "${OUTDIR}/ENCODE_proximal_enhancers.bed"
N_PELS=$(wc -l < "${OUTDIR}/ENCODE_proximal_enhancers.bed")
echo "Proximal enhancers (pELS): ${N_PELS}"

# Also get promoters for comparison
awk -F'\t' '$11 == "PLS"' "${CCRE_FILE}" > "${OUTDIR}/ENCODE_promoters.bed"
N_PLS=$(wc -l < "${OUTDIR}/ENCODE_promoters.bed")
echo "Promoters (PLS): ${N_PLS}"
echo ""

# =============================================================================
# STEP 3: PREPARE ENHANCER SETS BY TES/TEAD1 BINDING
# =============================================================================

echo "=== STEP 3: Preparing Enhancer Sets by Binding Status ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 32_encode_enhancer_methylation.R

# Check if BED files were created
if [ ! -f "${OUTDIR}/TES_bound_enhancers.bed" ]; then
    echo "ERROR: R script failed to create BED files"
    exit 1
fi

# Get counts
N_TES_BOUND=$(wc -l < "${OUTDIR}/TES_bound_enhancers.bed")
N_TEAD1_BOUND=$(wc -l < "${OUTDIR}/TEAD1_bound_enhancers.bed")
N_UNBOUND=$(wc -l < "${OUTDIR}/Unbound_enhancers.bed")
N_TES_ONLY=$(wc -l < "${OUTDIR}/TES_only_enhancers.bed")
N_TEAD1_ONLY=$(wc -l < "${OUTDIR}/TEAD1_only_enhancers.bed")
N_BOTH=$(wc -l < "${OUTDIR}/Both_bound_enhancers.bed")

echo ""
echo "Enhancer counts by binding status:"
echo "  TES-bound enhancers:        ${N_TES_BOUND}"
echo "  TEAD1-bound enhancers:      ${N_TEAD1_BOUND}"
echo "  TES-only enhancers:         ${N_TES_ONLY}"
echo "  TEAD1-only enhancers:       ${N_TEAD1_ONLY}"
echo "  Both TES+TEAD1 enhancers:   ${N_BOTH}"
echo "  Unbound enhancers:          ${N_UNBOUND}"
echo ""

# =============================================================================
# STEP 4: COMPUTE METHYLATION MATRIX
# =============================================================================

echo "=== STEP 4: Computing Methylation Matrix ==="
echo ""

conda activate tg

# BigWig files
TES_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/05_bigwig/TES_average_IP.bw"
GFP_METH="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/05_bigwig/GFP_average.bw"
TES_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/06_bigwig/TES_comb.bw"
TEAD1_BIND="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/06_bigwig/TEAD1_comb.bw"

# Verify inputs
for f in "$TES_METH" "$GFP_METH" "$TES_BIND" "$TEAD1_BIND"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

# Subsample unbound enhancers to match TES-bound count (for faster computation)
N_TES_BOUND=$(wc -l < "${OUTDIR}/TES_bound_enhancers.bed")
shuf -n ${N_TES_BOUND} ${OUTDIR}/Unbound_enhancers.bed | sort -k1,1 -k2,2n > ${OUTDIR}/Unbound_enhancers_subsampled.bed
N_SUBSAMPLE=$(wc -l < "${OUTDIR}/Unbound_enhancers_subsampled.bed")
echo "Subsampled unbound enhancers: ${N_SUBSAMPLE} (from ${N_UNBOUND})"

# Compute matrix with 3 groups: TES-bound, TEAD1-only, Unbound (subsampled)
computeMatrix reference-point \
    --referencePoint center \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${OUTDIR}/TES_bound_enhancers.bed \
       ${OUTDIR}/TEAD1_only_enhancers.bed \
       ${OUTDIR}/Unbound_enhancers_subsampled.bed \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/encode_enhancer_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

echo "Matrix created: ${OUTDIR}/encode_enhancer_matrix.gz"
echo ""

# =============================================================================
# STEP 5: GENERATE PROFILE PLOTS
# =============================================================================

echo "=== STEP 5: Generating Profile Plots ==="
echo ""

# Main comparison plot (all signals)
plotProfile -m ${OUTDIR}/encode_enhancer_matrix.gz \
    -out ${OUTDIR}/MAIN_ENCODE_Enhancer_Profile.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --refPointLabel "Enhancer Center" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "TES-bound (n=${N_TES_BOUND})" \
                   "TEAD1-only (n=${N_TEAD1_ONLY})" \
                   "Unbound subsampled (n=${N_SUBSAMPLE})" \
    --plotTitle "ENCODE Enhancers: Methylation by TES Binding Status" \
    --plotHeight 14 \
    --plotWidth 16 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Created: MAIN_ENCODE_Enhancer_Profile.png"

# Methylation only plot
plotProfile -m ${OUTDIR}/encode_enhancer_matrix.gz \
    -out ${OUTDIR}/METHYLATION_ENCODE_Enhancer_Profile.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "TES-bound" "TEAD1-only" "Unbound" \
    --plotTitle "Methylation at ENCODE Enhancers" \
    --yMin 0 \
    --dpi 300

echo "  Created: METHYLATION_ENCODE_Enhancer_Profile.png"

# Heatmap
plotHeatmap -m ${OUTDIR}/encode_enhancer_matrix.gz \
    -out ${OUTDIR}/ENCODE_Enhancer_Heatmap.png \
    --colorMap RdBu_r \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "TES-bound" "TEAD1-only" "Unbound" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --zMin 0 \
    --dpi 200

echo "  Created: ENCODE_Enhancer_Heatmap.png"
echo ""

# =============================================================================
# STEP 6: STATISTICAL QUANTIFICATION
# =============================================================================

echo "=== STEP 6: Statistical Quantification ==="
echo ""

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")
OUTPUT_DIR <- "output/32_encode_enhancer"

# Read Matrix
matrix_file <- file.path(OUTPUT_DIR, "encode_enhancer_matrix.gz")
con <- gzfile(matrix_file, "rt")
header_line <- readLines(con, n = 1)
close(con)
header_json <- fromJSON(gsub("^@", "", header_line))

cat("Matrix info:\n")
cat(sprintf("  Samples: %s\n", paste(header_json$sample_labels, collapse=", ")))
cat(sprintf("  Groups: %s\n", paste(header_json$group_labels, collapse=", ")))
cat(sprintf("  Group boundaries: %s\n", paste(header_json$group_boundaries, collapse=", ")))

mat <- fread(cmd = paste("zcat", matrix_file, "| tail -n +2"), header = FALSE)

# Groups
group_bounds <- header_json$group_boundaries
n_g1 <- group_bounds[2]
n_g2 <- group_bounds[3] - group_bounds[2]
n_g3 <- group_bounds[4] - group_bounds[3]

g1_rows <- 1:n_g1  # TES-bound
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)  # TEAD1-only
g3_rows <- (n_g1 + n_g2 + 1):(n_g1 + n_g2 + n_g3)  # Unbound

cat(sprintf("\n  Group 1 (TES-bound): %d enhancers\n", n_g1))
cat(sprintf("  Group 2 (TEAD1-only): %d enhancers\n", n_g2))
cat(sprintf("  Group 3 (Unbound): %d enhancers\n", n_g3))

# Bins: Center is at 5000bp. Total 10000bp. 50bp bins -> 200 bins.
# Center bins: roughly 90-110 (±500bp from center)
n_bins_per_sample <- as.integer(header_json$upstream[1] + header_json$downstream[1]) / header_json$`bin size`[1]
center_bins <- 90:110

tes_meth_cols <- 6 + center_bins
gfp_meth_cols <- 6 + n_bins_per_sample + center_bins

# Calculate methylation difference (TES - GFP) at center
calc_diff <- function(rows) {
    tes <- rowMeans(as.matrix(mat[rows, ..tes_meth_cols]), na.rm=TRUE)
    gfp <- rowMeans(as.matrix(mat[rows, ..gfp_meth_cols]), na.rm=TRUE)
    return(tes - gfp)
}

g1_diff <- calc_diff(g1_rows)
g2_diff <- calc_diff(g2_rows)
g3_diff <- calc_diff(g3_rows)

# Statistics
cat("\n========================================\n")
cat("METHYLATION DIFFERENCE (TES - GFP) at ENCODE Enhancers\n")
cat("(Center ±500bp)\n")
cat("========================================\n\n")

cat(sprintf("TES-bound enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% positive (TES > GFP): %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% negative (TES < GFP): %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("TEAD1-only enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% positive (TES > GFP): %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% negative (TES < GFP): %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("Unbound enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% positive (TES > GFP): %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% negative (TES < GFP): %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

# Statistical tests
cat("========================================\n")
cat("STATISTICAL TESTS\n")
cat("========================================\n\n")

# TES-bound vs Unbound
res_1v3 <- wilcox.test(g1_diff, g3_diff)
cat("TES-bound vs Unbound:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v3$p.value))
if (mean(g1_diff, na.rm=TRUE) > mean(g3_diff, na.rm=TRUE)) {
    cat("  Direction: TES-bound shows MORE methylation\n\n")
} else {
    cat("  Direction: TES-bound shows LESS methylation\n\n")
}

# TEAD1-only vs Unbound
res_2v3 <- wilcox.test(g2_diff, g3_diff)
cat("TEAD1-only vs Unbound:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_2v3$p.value))
if (mean(g2_diff, na.rm=TRUE) > mean(g3_diff, na.rm=TRUE)) {
    cat("  Direction: TEAD1-only shows MORE methylation\n\n")
} else {
    cat("  Direction: TEAD1-only shows LESS methylation\n\n")
}

# TES-bound vs TEAD1-only
res_1v2 <- wilcox.test(g1_diff, g2_diff)
cat("TES-bound vs TEAD1-only:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v2$p.value))
if (mean(g1_diff, na.rm=TRUE) > mean(g2_diff, na.rm=TRUE)) {
    cat("  Direction: TES-bound shows MORE methylation than TEAD1-only\n\n")
} else {
    cat("  Direction: TES-bound shows LESS methylation than TEAD1-only\n\n")
}

# KEY BIOLOGICAL QUESTION
cat("========================================\n")
cat("KEY BIOLOGICAL QUESTION\n")
cat("========================================\n\n")

if (mean(g1_diff, na.rm=TRUE) > 0 && res_1v3$p.value < 0.05) {
    cat("RESULT: TES binding CAUSES HYPERMETHYLATION at ENCODE enhancers!\n")
    cat("        (TES-bound enhancers show TES > GFP, higher than unbound)\n")
} else if (mean(g1_diff, na.rm=TRUE) < 0 && res_1v3$p.value < 0.05) {
    cat("RESULT: TES binding is associated with HYPOMETHYLATION at ENCODE enhancers\n")
    cat("        (TES-bound enhancers show TES < GFP, lower than unbound)\n")
} else {
    cat("RESULT: No significant difference in methylation between TES-bound and unbound enhancers\n")
}

# Save results
results_df <- data.frame(
    Group = c("TES_bound", "TEAD1_only", "Unbound"),
    N = c(n_g1, n_g2, n_g3),
    Mean_diff = c(mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE), mean(g3_diff, na.rm=TRUE)),
    Median_diff = c(median(g1_diff, na.rm=TRUE), median(g2_diff, na.rm=TRUE), median(g3_diff, na.rm=TRUE)),
    Pct_hypermethylated = c(100*mean(g1_diff > 0, na.rm=TRUE), 100*mean(g2_diff > 0, na.rm=TRUE), 100*mean(g3_diff > 0, na.rm=TRUE)),
    Pct_hypomethylated = c(100*mean(g1_diff < 0, na.rm=TRUE), 100*mean(g2_diff < 0, na.rm=TRUE), 100*mean(g3_diff < 0, na.rm=TRUE))
)

write.csv(results_df, file.path(OUTPUT_DIR, "methylation_statistics.csv"), row.names=FALSE)
cat("\nResults saved to: methylation_statistics.csv\n")

# Save detailed statistics
sink(file.path(OUTPUT_DIR, "statistical_results.txt"))
cat("ENCODE Enhancer Methylation Analysis\n")
cat("====================================\n\n")
cat("Methylation Difference (TES - GFP)\n\n")
print(results_df)
cat("\n\nStatistical Tests:\n")
cat("\nTES-bound vs Unbound:\n")
print(res_1v3)
cat("\nTEAD1-only vs Unbound:\n")
print(res_2v3)
cat("\nTES-bound vs TEAD1-only:\n")
print(res_1v2)
sink()

RSCRIPT_QUANTIFY

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -la ${OUTDIR}/*.bed ${OUTDIR}/*.png ${OUTDIR}/*.csv ${OUTDIR}/*.txt 2>/dev/null | awk '{print "  " $NF}'
