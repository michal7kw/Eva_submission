#!/bin/bash
# meDIP NOTE: reads DNA-methylation (meDIP) data from the ARCHIVED tree SRF_Eva_top/Archive/meDIP, which is NOT self-contained in Eva_submission. If the Archive is moved or removed, repoint or copy the meDIP inputs before running.
#SBATCH --job-name=34_verify_norm_bias
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/34_verify_normalization_bias.out
#SBATCH --error=logs/34_verify_normalization_bias.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# =============================================================================
# VERIFY NORMALIZATION BIAS IN meDIP ANALYSIS
# =============================================================================
#
# Purpose: Check whether the observed hypomethylation at TES-bound enhancers
#          could be an artifact of INPUT normalization (or lack thereof)
#
# Three verification tasks:
# 1. Check INPUT signal at TES-bound vs Unbound enhancers
# 2. Create INPUT-subtracted BigWig files and repeat analysis
# 3. Compare MEDIPS IP-only vs INPUT-normalized results
#
# =============================================================================

echo "=========================================="
echo "VERIFY NORMALIZATION BIAS"
echo "=========================================="
echo "Started: $(date)"
echo ""

cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts/analysis_1

# =============================================================================
# PATH CONFIGURATION
# =============================================================================

OUTDIR="output/34_verify_normalization_bias"
mkdir -p ${OUTDIR}

# BigWig files
TES_IP="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/05_bigwig/TES_average_IP.bw"
GFP_IP="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/05_bigwig/GFP_average.bw"
INPUT_BW="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/05_bigwig/TESmut-1-INPUT_RPKM.bw"

# Enhancer BED files from script 32
TES_BOUND_ENH="output/32_encode_enhancer/TES_bound_enhancers.bed"
UNBOUND_ENH="output/32_encode_enhancer/Unbound_enhancers_subsampled.bed"

# MEDIPS directories
MEDIPS_IP_ONLY="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/07_differential_MEDIPS"
MEDIPS_INPUT_NORM="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Archive/meDIP/results/07_differential_MEDIPS_INPUT_normalized"

# Verify files exist
echo "=== Verifying Input Files ==="
for f in "$TES_IP" "$GFP_IP" "$INPUT_BW" "$TES_BOUND_ENH" "$UNBOUND_ENH"; do
    if [ -f "$f" ]; then
        echo "  Found: $f"
    else
        echo "  MISSING: $f"
        exit 1
    fi
done
echo ""

# =============================================================================
# TASK 1: CHECK INPUT SIGNAL AT ENHANCERS
# =============================================================================

echo "=== TASK 1: Checking INPUT Signal at Enhancers ==="
echo ""

# Activate deepTools environment (same as script 32)
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate tg

# Compute matrix for INPUT signal at TES-bound vs Unbound enhancers
echo "Computing INPUT signal matrix..."
computeMatrix reference-point \
    --referencePoint center \
    -S ${INPUT_BW} \
    -R ${TES_BOUND_ENH} ${UNBOUND_ENH} \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --binSize 50 \
    -o ${OUTDIR}/input_signal_matrix.gz \
    -p 8

# Generate profile plot for INPUT signal
echo "Generating INPUT signal profile..."
plotProfile \
    -m ${OUTDIR}/input_signal_matrix.gz \
    --perGroup \
    --colors "#E31A1C" "#636363" \
    --samplesLabel "INPUT" \
    --regionsLabel "TES-bound enhancers" "Unbound enhancers" \
    -out ${OUTDIR}/INPUT_signal_at_enhancers_profile.png \
    --plotTitle "INPUT Signal at ENCODE Enhancers"

echo "  Created: INPUT_signal_at_enhancers_profile.png"

# =============================================================================
# TASK 2: CREATE INPUT-SUBTRACTED BIGWIG FILES
# =============================================================================

echo ""
echo "=== TASK 2: Creating INPUT-Subtracted BigWig Files ==="
echo ""

# Subtract INPUT from TES IP
echo "Creating TES INPUT-subtracted BigWig..."
bigwigCompare \
    -b1 ${TES_IP} \
    -b2 ${INPUT_BW} \
    --operation subtract \
    -o ${OUTDIR}/TES_INPUT_subtracted.bw \
    -p 8

# Subtract INPUT from GFP IP
echo "Creating GFP INPUT-subtracted BigWig..."
bigwigCompare \
    -b1 ${GFP_IP} \
    -b2 ${INPUT_BW} \
    --operation subtract \
    -o ${OUTDIR}/GFP_INPUT_subtracted.bw \
    -p 8

echo "  Created: TES_INPUT_subtracted.bw"
echo "  Created: GFP_INPUT_subtracted.bw"

# =============================================================================
# TASK 2B: REPEAT ENHANCER ANALYSIS WITH INPUT-SUBTRACTED FILES
# =============================================================================

echo ""
echo "Computing methylation matrix with INPUT-subtracted files..."

computeMatrix reference-point \
    --referencePoint center \
    -S ${OUTDIR}/TES_INPUT_subtracted.bw ${OUTDIR}/GFP_INPUT_subtracted.bw \
    -R ${TES_BOUND_ENH} ${UNBOUND_ENH} \
    --beforeRegionStartLength 5000 \
    --afterRegionStartLength 5000 \
    --binSize 50 \
    -o ${OUTDIR}/input_subtracted_enhancer_matrix.gz \
    -p 8

# Profile plot with INPUT-subtracted data
echo "Generating INPUT-subtracted profile..."
plotProfile \
    -m ${OUTDIR}/input_subtracted_enhancer_matrix.gz \
    --perGroup \
    --colors "#7B3294" "#636363" \
    --samplesLabel "TES-INPUT" "GFP-INPUT" \
    --regionsLabel "TES-bound enhancers" "Unbound enhancers" \
    -out ${OUTDIR}/INPUT_subtracted_enhancer_profile.png \
    --plotTitle "INPUT-Subtracted Methylation at ENCODE Enhancers"

echo "  Created: INPUT_subtracted_enhancer_profile.png"

# =============================================================================
# TASK 3: RUN R SCRIPT FOR STATISTICAL ANALYSIS
# =============================================================================

echo ""
echo "=== TASK 3: Statistical Analysis (R) ==="
echo ""

# Switch to R environment
conda deactivate
conda activate r_chipseq_env

Rscript 34_verify_normalization_bias.R

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=========================================="
echo "VERIFICATION COMPLETE"
echo "=========================================="
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -la ${OUTDIR}/
echo ""
