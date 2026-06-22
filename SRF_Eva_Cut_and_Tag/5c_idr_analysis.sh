#!/bin/bash

#===============================================================================
# SCRIPT: 5c_idr_analysis.sh
# PURPOSE: Replicate reproducibility assessment (IDR-style overlap analysis)
#
# DESCRIPTION:
# Assesses peak reproducibility between biological replicates using pairwise
# overlap rates. This is the ENCODE-recommended approach for validating that
# peaks are consistently detected across replicates. Reports overlap rates
# and correlation matrices per experimental group.
#
# USAGE:
# sbatch 5c_idr_analysis.sh
#
# INPUTS:
# - Individual replicate narrowPeak files from step 5
#
# OUTPUTS:
# - results/idr_analysis/replicate_overlap_summary.csv
# - results/idr_analysis/correlation_matrix.pdf
# - results/idr_analysis/overlap_rates_by_group.pdf
#===============================================================================

#SBATCH --job-name=5c_idr
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/5c_idr.err"
#SBATCH --output="./logs/5c_idr.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate diffbind_analysis

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"

echo "=========================================="
echo "IDR / Replicate Reproducibility Analysis"
echo "Date: $(date)"
echo "=========================================="

mkdir -p ${BASE_DIR}/results/idr_analysis
mkdir -p ${BASE_DIR}/logs

# Verify peak files exist
MISSING=0
for SAMPLE in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3; do
    PEAK_FILE="${BASE_DIR}/results/05_peaks_narrow/${SAMPLE}_peaks.narrowPeak"
    if [[ ! -f "${PEAK_FILE}" ]]; then
        echo "WARNING: Missing peak file for ${SAMPLE}: ${PEAK_FILE}"
        MISSING=$((MISSING + 1))
    fi
done

if [[ ${MISSING} -gt 0 ]]; then
    echo "WARNING: ${MISSING} peak files missing. Analysis may be incomplete."
fi

echo "Running IDR analysis R script..."
Rscript ${BASE_DIR}/scripts/idr_analysis.R

if [[ $? -eq 0 ]]; then
    echo ""
    echo "IDR analysis complete!"
    echo "Results: ${BASE_DIR}/results/idr_analysis/"
    ls -lh ${BASE_DIR}/results/idr_analysis/ 2>/dev/null
else
    echo "ERROR: IDR analysis failed. Check R script and logs."
    exit 1
fi

echo ""
echo "=========================================="
echo "IDR analysis completed"
echo "Date: $(date)"
echo "=========================================="
