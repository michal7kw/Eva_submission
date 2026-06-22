#!/bin/bash
#SBATCH --job-name=homer_peaks
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --output=logs/homer_peaks_%j.out
#SBATCH --error=logs/homer_peaks_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kubacki.michal@hsr.it

# ============================================================================
# HOMER Motif Analysis for TES and TEAD1 Cut&Tag Peaks
# ============================================================================
# Purpose: Identify TF binding motifs in Cut&Tag peaks to confirm TEAD consensus
#          Expected: TEAD1/TEAD4 consensus (GGAATG / MCAT motif)
#
# Input:  Narrow peaks from TES and TEAD1 Cut&Tag
# Output: HOMER motif enrichment results for both samples
# ============================================================================

set -e
set -o pipefail

# Navigate to Cut&Tag directory
cd /beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag

# Create output directories
mkdir -p results/08_homer_motifs/TES_peaks
mkdir -p results/08_homer_motifs/TEAD1_peaks
mkdir -p logs

echo "=============================================="
echo "HOMER Motif Analysis - Cut&Tag Peaks"
echo "Started: $(date)"
echo "=============================================="

# Activate conda environment with HOMER
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate homer_env

# Check HOMER is available
if ! command -v findMotifsGenome.pl &> /dev/null; then
    echo "Trying genomics_env..."
    conda activate genomics_env
    if ! command -v findMotifsGenome.pl &> /dev/null; then
        echo "ERROR: HOMER not available"
        exit 1
    fi
fi

echo "HOMER version: $(findMotifsGenome.pl 2>&1 | head -1 || echo 'Version check failed')"

# Define peak files
# First check for consensus peaks, then individual peaks
if [[ -f "results/11_combined_replicates_narrow/TES_consensus_peaks.bed" ]]; then
    TES_PEAKS="results/11_combined_replicates_narrow/TES_consensus_peaks.bed"
elif [[ -f "results/05_peaks_narrow/TES_peaks.narrowPeak" ]]; then
    TES_PEAKS="results/05_peaks_narrow/TES_peaks.narrowPeak"
else
    # Find any TES peak file
    TES_PEAKS=$(find results -name "*TES*peak*" -type f | head -1)
fi

if [[ -f "results/11_combined_replicates_narrow/TEAD1_consensus_peaks.bed" ]]; then
    TEAD1_PEAKS="results/11_combined_replicates_narrow/TEAD1_consensus_peaks.bed"
elif [[ -f "results/05_peaks_narrow/TEAD1_peaks.narrowPeak" ]]; then
    TEAD1_PEAKS="results/05_peaks_narrow/TEAD1_peaks.narrowPeak"
else
    TEAD1_PEAKS=$(find results -name "*TEAD1*peak*" -type f | head -1)
fi

echo "Using TES peaks: $TES_PEAKS"
echo "Using TEAD1 peaks: $TEAD1_PEAKS"

# Verify files exist
for PEAKS in "$TES_PEAKS" "$TEAD1_PEAKS"; do
    if [[ -z "$PEAKS" ]] || [[ ! -f "$PEAKS" ]]; then
        echo "ERROR: Peak file not found: $PEAKS"
        exit 1
    fi
done

# Count peaks
TES_COUNT=$(wc -l < "$TES_PEAKS")
TEAD1_COUNT=$(wc -l < "$TEAD1_PEAKS")
echo "TES peaks: $TES_COUNT"
echo "TEAD1 peaks: $TEAD1_COUNT"

# ============================================
# Run HOMER on TES peaks
# ============================================
echo ""
echo "Running HOMER on TES peaks..."
echo "=============================================="

findMotifsGenome.pl \
    "$TES_PEAKS" \
    hg38 \
    results/08_homer_motifs/TES_peaks/ \
    -size 200 \
    -p 8 \
    -len 8,10,12 \
    -S 25 \
    2>&1 | tee logs/homer_TES_peaks.log

# ============================================
# Run HOMER on TEAD1 peaks
# ============================================
echo ""
echo "Running HOMER on TEAD1 peaks..."
echo "=============================================="

findMotifsGenome.pl \
    "$TEAD1_PEAKS" \
    hg38 \
    results/08_homer_motifs/TEAD1_peaks/ \
    -size 200 \
    -p 8 \
    -len 8,10,12 \
    -S 25 \
    2>&1 | tee logs/homer_TEAD1_peaks.log

# ============================================
# Summarize results
# ============================================
echo ""
echo "=============================================="
echo "Analysis Complete"
echo "=============================================="

for SAMPLE in TES TEAD1; do
    OUTDIR="results/08_homer_motifs/${SAMPLE}_peaks"
    echo ""
    echo "$SAMPLE peak motifs:"
    echo "  Results: $OUTDIR/homerResults.html"

    if [[ -f "$OUTDIR/knownResults.txt" ]]; then
        echo "  Top 5 known motifs:"
        head -6 "$OUTDIR/knownResults.txt" | tail -5 | sed 's/^/    /'

        echo ""
        echo "  TEAD motifs found:"
        grep -i "TEAD" "$OUTDIR/knownResults.txt" | head -5 | sed 's/^/    /' || echo "    (None in results)"
    fi
done

echo ""
echo "=============================================="
echo "EXPECTED: TEAD1/TEAD4 consensus (GGAATG / MCAT)"
echo "=============================================="
echo ""
echo "Finished: $(date)"
