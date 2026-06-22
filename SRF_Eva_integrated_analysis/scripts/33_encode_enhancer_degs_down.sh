#!/bin/bash
# meDIP NOTE: reads DNA-methylation (meDIP) data from the ARCHIVED tree SRF_Eva_top/Archive/meDIP, which is NOT self-contained in Eva_submission. If the Archive is moved or removed, repoint or copy the meDIP inputs before running.
#SBATCH --job-name=a1_33_encode_degs_down
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/33_encode_enhancer_degs_down.out
#SBATCH --error=logs/33_encode_enhancer_degs_down.err

# =============================================================================
# ENCODE ENHANCERS OF DEGs DOWN - METHYLATION ANALYSIS
# =============================================================================
#
# Purpose: Analyze methylation at ENCODE enhancers specifically associated
#          with downregulated DEGs, stratified by TES/TEAD1 binding status.
#
# Key Question: Do TES-bound enhancers of DEGs DOWN show a different
#               methylation pattern than all TES-bound enhancers?
#
# =============================================================================

echo "=========================================="
echo "ENCODE ENHANCERS OF DEGs DOWN"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1

OUTDIR="output/33_encode_enhancer_degs_down"
mkdir -p ${OUTDIR}
mkdir -p logs

# =============================================================================
# STEP 1: PREPARE BED FILES (R SCRIPT)
# =============================================================================

echo "=== STEP 1: Preparing Enhancer BED Files ==="
echo ""

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate r_chipseq_env

Rscript 33_encode_enhancer_degs_down.R

# Check if BED files were created
BED_TES="${OUTDIR}/TES_bound_enhancers_DEGs_DOWN.bed"
BED_TEAD1="${OUTDIR}/TEAD1_only_enhancers_DEGs_DOWN.bed"
BED_UNBOUND="${OUTDIR}/Unbound_enhancers_DEGs_DOWN.bed"
BED_CONTROL="${OUTDIR}/Control_enhancers.bed"

if [[ ! -f "$BED_TES" ]]; then
    echo "ERROR: BED files not created. Check R script output."
    exit 1
fi

echo ""
echo "BED files created successfully."
echo ""

# Get counts
N_TES=$(wc -l < ${BED_TES})
N_TEAD1=$(wc -l < ${BED_TEAD1})
N_UNBOUND=$(wc -l < ${BED_UNBOUND})
N_CONTROL=$(wc -l < ${BED_CONTROL})

echo "Enhancer Counts (DEGs DOWN):"
echo "  TES-bound enhancers:    ${N_TES}"
echo "  TEAD1-only enhancers:   ${N_TEAD1}"
echo "  Unbound enhancers:      ${N_UNBOUND}"
echo "  Control enhancers:      ${N_CONTROL}"
echo ""

# =============================================================================
# STEP 2: SETUP DEEPTOOLS
# =============================================================================

echo "=== STEP 2: Setting up deepTools ==="

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

# =============================================================================
# STEP 3: COMPUTE METHYLATION MATRIX
# =============================================================================

echo "=== STEP 3: Computing Methylation Matrix ==="
echo ""

# Compute matrix with 4 groups
computeMatrix reference-point \
    --referencePoint center \
    -S $TES_METH $GFP_METH $TES_BIND $TEAD1_BIND \
    -R ${BED_TES} \
       ${BED_TEAD1} \
       ${BED_UNBOUND} \
       ${BED_CONTROL} \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --binSize 50 \
    --skipZeros \
    --missingDataAsZero \
    -o ${OUTDIR}/encode_degs_down_matrix.gz \
    -p 16 \
    2>&1 | grep -v "Skipping\|did not match"

echo "Matrix created: ${OUTDIR}/encode_degs_down_matrix.gz"
echo ""

# =============================================================================
# STEP 4: GENERATE PROFILE PLOTS
# =============================================================================

echo "=== STEP 4: Generating Profile Plots ==="
echo ""

# Main comparison plot (all signals)
plotProfile -m ${OUTDIR}/encode_degs_down_matrix.gz \
    -out ${OUTDIR}/MAIN_DEGs_DOWN_Enhancer_Profile.png \
    --perGroup \
    --colors "#7B3294" "#636363" "#E31A1C" "#377EB8" \
    --refPointLabel "Enhancer Center" \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "TES-bound DEGs DOWN (n=${N_TES})" \
                   "TEAD1-only DEGs DOWN (n=${N_TEAD1})" \
                   "Unbound DEGs DOWN (n=${N_UNBOUND})" \
                   "Control (n=${N_CONTROL})" \
    --plotTitle "ENCODE Enhancers of DEGs DOWN: Methylation by TES Binding" \
    --plotHeight 14 \
    --plotWidth 18 \
    --legendLocation "upper-left" \
    --yMin 0 \
    --dpi 300

echo "  Created: MAIN_DEGs_DOWN_Enhancer_Profile.png"

# Methylation only plot
plotProfile -m ${OUTDIR}/encode_degs_down_matrix.gz \
    -out ${OUTDIR}/METHYLATION_DEGs_DOWN_Profile.png \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --samplesLabel "TES meth" "GFP meth" \
    --regionsLabel "TES-bound DEGs DOWN" "TEAD1-only DEGs DOWN" "Unbound DEGs DOWN" "Control" \
    --plotTitle "Methylation at Enhancers of DEGs DOWN" \
    --yMin 0 \
    --dpi 300

echo "  Created: METHYLATION_DEGs_DOWN_Profile.png"

# Heatmap
plotHeatmap -m ${OUTDIR}/encode_degs_down_matrix.gz \
    -out ${OUTDIR}/DEGs_DOWN_Enhancer_Heatmap.png \
    --colorMap RdBu_r \
    --samplesLabel "TES meth" "GFP meth" "TES bind" "TEAD1 bind" \
    --regionsLabel "TES-bound" "TEAD1-only" "Unbound" "Control" \
    --sortUsing mean \
    --sortUsingSamples 1 \
    --zMin 0 \
    --dpi 200

echo "  Created: DEGs_DOWN_Enhancer_Heatmap.png"
echo ""

# =============================================================================
# STEP 5: STATISTICAL QUANTIFICATION
# =============================================================================

echo "=== STEP 5: Statistical Quantification ==="
echo ""

conda activate r_chipseq_env

Rscript - << 'RSCRIPT_QUANTIFY'
suppressPackageStartupMessages({
    library(data.table)
    library(jsonlite)
    library(dplyr)
})

setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1")
OUTPUT_DIR <- "output/33_encode_enhancer_degs_down"

# Read Matrix
matrix_file <- file.path(OUTPUT_DIR, "encode_degs_down_matrix.gz")
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
n_g4 <- group_bounds[5] - group_bounds[4]

g1_rows <- 1:n_g1  # TES-bound DEGs DOWN
g2_rows <- (n_g1 + 1):(n_g1 + n_g2)  # TEAD1-only DEGs DOWN
g3_rows <- (n_g1 + n_g2 + 1):(n_g1 + n_g2 + n_g3)  # Unbound DEGs DOWN
g4_rows <- (n_g1 + n_g2 + n_g3 + 1):(n_g1 + n_g2 + n_g3 + n_g4)  # Control

cat(sprintf("\n  Group 1 (TES-bound DEGs DOWN): %d enhancers\n", n_g1))
cat(sprintf("  Group 2 (TEAD1-only DEGs DOWN): %d enhancers\n", n_g2))
cat(sprintf("  Group 3 (Unbound DEGs DOWN): %d enhancers\n", n_g3))
cat(sprintf("  Group 4 (Control): %d enhancers\n", n_g4))

# Bins for center region
n_bins_per_sample <- as.integer(header_json$upstream[1] + header_json$downstream[1]) / header_json$`bin size`[1]
center_bins <- 90:110  # ±500bp from center

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
g4_diff <- calc_diff(g4_rows)

# Statistics
cat("\n========================================\n")
cat("METHYLATION DIFFERENCE (TES - GFP) at ENCODE Enhancers of DEGs DOWN\n")
cat("(Center ±500bp)\n")
cat("========================================\n\n")

cat(sprintf("TES-bound DEGs DOWN enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g1_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g1_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated (TES > GFP): %.1f%%\n", 100*mean(g1_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated (TES < GFP): %.1f%%\n\n", 100*mean(g1_diff < 0, na.rm=TRUE)))

cat(sprintf("TEAD1-only DEGs DOWN enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g2_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g2_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated (TES > GFP): %.1f%%\n", 100*mean(g2_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated (TES < GFP): %.1f%%\n\n", 100*mean(g2_diff < 0, na.rm=TRUE)))

cat(sprintf("Unbound DEGs DOWN enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g3_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g3_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated (TES > GFP): %.1f%%\n", 100*mean(g3_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated (TES < GFP): %.1f%%\n\n", 100*mean(g3_diff < 0, na.rm=TRUE)))

cat(sprintf("Control enhancers:\n"))
cat(sprintf("  Mean diff:   %.4f\n", mean(g4_diff, na.rm=TRUE)))
cat(sprintf("  Median diff: %.4f\n", median(g4_diff, na.rm=TRUE)))
cat(sprintf("  %% hypermethylated (TES > GFP): %.1f%%\n", 100*mean(g4_diff > 0, na.rm=TRUE)))
cat(sprintf("  %% hypomethylated (TES < GFP): %.1f%%\n\n", 100*mean(g4_diff < 0, na.rm=TRUE)))

# Statistical tests
cat("========================================\n")
cat("STATISTICAL TESTS\n")
cat("========================================\n\n")

# TES-bound DEGs DOWN vs Unbound DEGs DOWN
res_1v3 <- wilcox.test(g1_diff, g3_diff)
cat("TES-bound DEGs DOWN vs Unbound DEGs DOWN:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v3$p.value))
if (mean(g1_diff, na.rm=TRUE) > mean(g3_diff, na.rm=TRUE)) {
    cat("  Direction: TES-bound shows MORE methylation\n\n")
} else {
    cat("  Direction: TES-bound shows LESS methylation\n\n")
}

# TES-bound DEGs DOWN vs Control
res_1v4 <- wilcox.test(g1_diff, g4_diff)
cat("TES-bound DEGs DOWN vs Control:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_1v4$p.value))
if (mean(g1_diff, na.rm=TRUE) > mean(g4_diff, na.rm=TRUE)) {
    cat("  Direction: TES-bound DEGs DOWN shows MORE methylation than Control\n\n")
} else {
    cat("  Direction: TES-bound DEGs DOWN shows LESS methylation than Control\n\n")
}

# Unbound DEGs DOWN vs Control
res_3v4 <- wilcox.test(g3_diff, g4_diff)
cat("Unbound DEGs DOWN vs Control:\n")
cat(sprintf("  Wilcoxon p-value: %.4e\n", res_3v4$p.value))
if (mean(g3_diff, na.rm=TRUE) > mean(g4_diff, na.rm=TRUE)) {
    cat("  Direction: Unbound DEGs DOWN shows MORE methylation than Control\n\n")
} else {
    cat("  Direction: Unbound DEGs DOWN shows LESS methylation than Control\n\n")
}

# KEY BIOLOGICAL QUESTION
cat("========================================\n")
cat("KEY BIOLOGICAL QUESTION\n")
cat("========================================\n\n")

if (mean(g1_diff, na.rm=TRUE) > 0 && res_1v3$p.value < 0.05) {
    cat("RESULT: TES binding CAUSES HYPERMETHYLATION at enhancers of DEGs DOWN!\n")
} else if (mean(g1_diff, na.rm=TRUE) < 0 && res_1v3$p.value < 0.05) {
    cat("RESULT: TES binding is associated with HYPOMETHYLATION at enhancers of DEGs DOWN\n")
    cat("        (Same pattern as all enhancers - TES binding does NOT cause methylation)\n")
} else {
    cat("RESULT: No significant difference in methylation between TES-bound and unbound enhancers of DEGs DOWN\n")
}

# Save results
results_df <- data.frame(
    Group = c("TES_bound_DEGs_DOWN", "TEAD1_only_DEGs_DOWN", "Unbound_DEGs_DOWN", "Control"),
    N = c(n_g1, n_g2, n_g3, n_g4),
    Mean_diff = c(mean(g1_diff, na.rm=TRUE), mean(g2_diff, na.rm=TRUE), mean(g3_diff, na.rm=TRUE), mean(g4_diff, na.rm=TRUE)),
    Median_diff = c(median(g1_diff, na.rm=TRUE), median(g2_diff, na.rm=TRUE), median(g3_diff, na.rm=TRUE), median(g4_diff, na.rm=TRUE)),
    Pct_hypermethylated = c(100*mean(g1_diff > 0, na.rm=TRUE), 100*mean(g2_diff > 0, na.rm=TRUE), 100*mean(g3_diff > 0, na.rm=TRUE), 100*mean(g4_diff > 0, na.rm=TRUE)),
    Pct_hypomethylated = c(100*mean(g1_diff < 0, na.rm=TRUE), 100*mean(g2_diff < 0, na.rm=TRUE), 100*mean(g3_diff < 0, na.rm=TRUE), 100*mean(g4_diff < 0, na.rm=TRUE))
)

write.csv(results_df, file.path(OUTPUT_DIR, "methylation_statistics.csv"), row.names=FALSE)
cat("\nResults saved to: methylation_statistics.csv\n")

RSCRIPT_QUANTIFY

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -la ${OUTDIR}/*.bed ${OUTDIR}/*.png ${OUTDIR}/*.csv 2>/dev/null | awk '{print "  " $NF}'
