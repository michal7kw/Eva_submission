#!/bin/bash

#===============================================================================
# SCRIPT: 8b_homer_motifs.sh
# PURPOSE: HOMER de novo motif discovery for Cut&Tag narrow peaks
#
# DESCRIPTION:
# Performs de novo motif discovery using HOMER findMotifsGenome.pl
# on consensus narrow peaks from TES and TEAD1 conditions.
#
# KEY OPERATIONS:
# 1. Prepares peak files with chr prefix for HOMER compatibility
# 2. Runs de novo motif discovery for each condition
# 3. Identifies enriched known motifs from HOMER database
# 4. Generates motif logos and HTML reports
#
# INPUTS:
# - Consensus narrow peaks from results/11_combined_replicates_narrow/peaks/
#
# OUTPUTS:
# - results/08_homer_motifs/{condition}_peaks/
#   - homerResults.html - Main results report
#   - knownResults.html - Known motif enrichment
#   - homerMotifs.all.motifs - All discovered motifs
#   - Motif logos and statistics
#
# USAGE:
# sbatch 8b_homer_motifs.sh
#===============================================================================

#SBATCH --job-name=8b_homer_motifs
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8b_homer_motifs.err"
#SBATCH --output="./logs/8b_homer_motifs.out"

set -e

# Activate HOMER environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate homer_env

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
PEAKS_DIR="${BASE_DIR}/results/11_combined_replicates_narrow/peaks"
OUTPUT_DIR="${BASE_DIR}/results/08_homer_motifs"
TMP_DIR="${BASE_DIR}/results/08_homer_motifs/tmp"

# Create output directories
mkdir -p ${OUTPUT_DIR}
mkdir -p ${TMP_DIR}

echo "=========================================="
echo "HOMER De Novo Motif Analysis"
echo "Date: $(date)"
echo "=========================================="

# Function to run HOMER motif analysis
run_homer_motifs() {
    local condition=$1
    local input_peaks="${PEAKS_DIR}/${condition}_consensus_peaks.bed"
    local output_dir="${OUTPUT_DIR}/${condition}_peaks"
    local tmp_file="${TMP_DIR}/${condition}_peaks_homer.bed"

    echo ""
    echo "Processing ${condition}..."

    # Check if input file exists
    if [ ! -s "${input_peaks}" ]; then
        echo "  WARNING: Peak file not found or empty: ${input_peaks}"
        return 1
    fi

    # Count peaks
    peak_count=$(wc -l < "${input_peaks}")
    echo "  Input peaks: ${peak_count}"

    # Prepare peaks for HOMER (add chr prefix if needed)
    echo "  Preparing peak file..."
    awk 'BEGIN{OFS="\t"} {
        if ($1 !~ /^chr/) {
            print "chr"$1, $2, $3, $4, $5
        } else {
            print $1, $2, $3, $4, $5
        }
    }' "${input_peaks}" > "${tmp_file}"

    # Run HOMER findMotifsGenome.pl
    echo "  Running HOMER motif discovery..."
    findMotifsGenome.pl \
        "${tmp_file}" \
        hg38 \
        "${output_dir}" \
        -size 200 \
        -mask \
        -p 8

    echo "  Completed ${condition} motif analysis"

    # Clean up temp file
    rm -f "${tmp_file}"
}

# Run motif analysis for each condition (TESmut removed - failed sample)
for condition in TES TEAD1; do
    run_homer_motifs "${condition}" || echo "  Failed or skipped ${condition}"
done

# Clean up temp directory
rmdir ${TMP_DIR} 2>/dev/null || true

echo ""
echo "=========================================="
echo "HOMER motif analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}/"
echo "Date: $(date)"
echo "=========================================="
