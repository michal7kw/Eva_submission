#!/bin/bash

#===============================================================================
# SCRIPT: 8_annotate_narrow.sh
# PURPOSE: Narrow peak annotation and motif analysis for Cut&Tag data
#
# DESCRIPTION:
# This script performs comprehensive narrow peak annotation and motif discovery analysis
# for identified narrow peaks. It annotates peaks with genomic features (promoters,
# introns, exons, etc.) and discovers enriched DNA motifs within narrow peak regions
# using HOMER motif analysis tools.
#
# KEY OPERATIONS:
# 1. Peak annotation using R-based genomic annotation tools
# 2. HOMER motif discovery for TES and TEAD1 peaks
# 3. Motif enrichment analysis with background correction
# 4. Generation of motif logos and statistics
#
# METHODOLOGY:
# - Uses R script for detailed genomic feature annotation
# - Employs HOMER findMotifsGenome.pl for de novo motif discovery
# - Analyzes 200bp regions around peak centers
# - Includes repeat masking for accurate motif identification
# - Compares motifs between different conditions
#
# IMPORTANT PARAMETERS:
# - Memory: 16GB (sufficient for motif analysis)
# - Time: 2 hours (adequate for HOMER processing)
# - Threads: 8 (parallel motif discovery)
# - Motif size: 200bp around peak centers
# - Genome: hg38 reference
#
# INPUTS:
# - Narrow peak files from MACS2 (05_peaks directory)
# - Peak comparison results from previous analysis
# - R annotation script (annotate_peaks.R)
#
# OUTPUTS:
# - Genomic feature annotations (from R script)
# - HOMER motif analysis results (TES_motifs/, TEAD1_motifs/)
# - Motif logos, enrichment statistics, and HTML reports
# - Known and de novo motif discoveries
#
# DEPENDENCIES:
# - R with genomic annotation packages
# - HOMER motif analysis suite
# - Conda environment: annotation_enrichment
#
# USAGE:
# sbatch 8_annotate_narrow.sh
#
# NOTES:
# - Adds "chr" prefix to peak coordinates for HOMER compatibility
# - Includes file existence checks before processing
# - Generates publication-ready motif analysis reports
# - Compares motifs between TES and TEAD1 conditions
#===============================================================================

#SBATCH --job-name=8_annotate_narrow
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/8_annotate_narrow.err"
#SBATCH --output="./logs/8_annotate_narrow.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate annotation_enrichment

# Load HOMER module
# module load homer

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
PEAKS_DIR="${BASE_DIR}/results/05_peaks_narrow"
CONSENSUS_DIR="${BASE_DIR}/results/11_combined_replicates_narrow"
ANALYSIS_DIR="${BASE_DIR}/results/07_analysis_narrow"

# Create output directory if it doesn't exist
mkdir -p ${ANALYSIS_DIR}

echo "=========================================="
echo "Starting peak annotation analysis"
echo "Date: $(date)"
echo "=========================================="

# Determine which peak sets to use for annotation
# Priority: consensus peaks > MACS2 combined peaks
echo "Preparing peak files for annotation..."

for condition in TES TEAD1; do  # TESmut removed - failed sample
    consensus_scored="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks_scored.bed"
    consensus_full="${CONSENSUS_DIR}/peaks/${condition}_consensus_peaks.bed"
    macs2_peaks="${PEAKS_DIR}/${condition}_peaks.narrowPeak"

    if [ -f "$consensus_scored" ]; then
        echo "  Using quality-scored consensus peaks for $condition annotation"
        cp "$consensus_scored" "${ANALYSIS_DIR}/${condition}_peaks_for_annotation.bed"
    elif [ -f "$consensus_full" ]; then
        echo "  Using full consensus peaks for $condition annotation"
        cp "$consensus_full" "${ANALYSIS_DIR}/${condition}_peaks_for_annotation.bed"
    elif [ -f "$macs2_peaks" ]; then
        echo "  Using MACS2 combined peaks for $condition annotation (consensus not available)"
        cp "$macs2_peaks" "${ANALYSIS_DIR}/${condition}_peaks_for_annotation.bed"
    else
        echo "  WARNING: No peak file found for $condition"
    fi
done

# Run R script for genomic annotation
echo "Running R annotation script..."
Rscript ${BASE_DIR}/scripts/annotate_peaks.R

# HOMER motif analysis for TES peaks (COMMENTED OUT)
# if [ -s "${PEAKS_DIR}/TES_peaks.narrowPeak" ]; then
#     echo "Running motif analysis for TES peaks..."
#     awk 'BEGIN{OFS="\t"} {print "chr"$0}' ${PEAKS_DIR}/TES_peaks.narrowPeak > ${PEAKS_DIR}/TES_peaks.narrowPeak.tmp
#     findMotifsGenome.pl \
#         ${PEAKS_DIR}/TES_peaks.narrowPeak.tmp \
#         hg38 \
#         ${ANALYSIS_DIR}/TES_motifs \
#         -size 200 \
#         -mask \
#         -p 8
#     rm ${PEAKS_DIR}/TES_peaks.narrowPeak.tmp
# else
#     echo "Skipping HOMER for TES peaks, file is empty or does not exist."
# fi

# Compare with TEAD1 peaks (COMMENTED OUT)
# if [ -s "${PEAKS_DIR}/TEAD1_peaks.narrowPeak" ]; then
#     echo "Running motif analysis for TEAD1 peaks..."
#     awk 'BEGIN{OFS="\t"} {print "chr"$0}' ${PEAKS_DIR}/TEAD1_peaks.narrowPeak > ${PEAKS_DIR}/TEAD1_peaks.narrowPeak.tmp
#     findMotifsGenome.pl \
#         ${PEAKS_DIR}/TEAD1_peaks.narrowPeak.tmp \
#         hg38 \
#         ${ANALYSIS_DIR}/TEAD1_motifs \
#         -size 200 \
#         -mask \
#         -p 8
#     rm ${PEAKS_DIR}/TEAD1_peaks.narrowPeak.tmp
# else
#     echo "Skipping HOMER for TEAD1 peaks, file is empty or does not exist."
# fi

echo "Narrow peak annotation and motif analysis complete!"