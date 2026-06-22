#!/bin/bash

#===============================================================================
# SCRIPT: rerun_integration.sh
# PURPOSE: Re-run integrative analysis after Cut&Tag and RNA-seq updates
#
# DESCRIPTION:
# Submits the updated 10_final_integrative_analysis.R which now includes:
# - LFC threshold for direct target classification
# - Distance threshold sensitivity analysis (5kb, 10kb, 25kb, 50kb, 100kb)
# - Sensitivity plot showing robustness of default 50kb threshold
#
# PREREQUISITES:
# Both Cut&Tag and RNA-seq pipelines must have completed first:
# - Cut&Tag peaks: results/05_peaks_narrow/ (re-generated with new filtering)
# - RNA-seq DEGs: results/05_deseq2/deseq2_results_TES_vs_GFP.txt (with apeglm)
#
# NEW OUTPUTS:
# - output/10_final_integrative_analysis/direct_targets/distance_threshold_sensitivity.csv
# - output/10_final_integrative_analysis/plots/06_distance_sensitivity.pdf
#
# USAGE:
#   ./rerun_integration.sh              # Run integration
#   ./rerun_integration.sh --dry-run    # Show what would be submitted
#===============================================================================

set -e

# Configuration
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_integrated_analysis/scripts"
cd "$SCRIPT_DIR"

# Parse arguments
DRY_RUN=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run]"
            exit 1
            ;;
    esac
done

# Create logs directory
mkdir -p logs

echo "=============================================="
echo "Integration Analysis RERUN (Revised Code)"
echo "=============================================="
echo "Script directory: $SCRIPT_DIR"
echo "Dry run: $DRY_RUN"
echo ""
echo "New features in 10_final_integrative_analysis.R:"
echo "  - LFC threshold for direct target classification"
echo "  - Distance threshold sensitivity analysis (5-100kb)"
echo "  - Sensitivity plot validating 50kb default threshold"
echo ""

#===============================================================================
# Submit integration analysis
#===============================================================================
echo "=== Submitting Integrative Analysis ==="

if [[ ! -f "10_final_integrative_analysis.sh" ]]; then
    echo "ERROR: 10_final_integrative_analysis.sh not found in $SCRIPT_DIR"
    exit 1
fi

if [[ "$DRY_RUN" == "true" ]]; then
    echo "  [DRY-RUN] sbatch 10_final_integrative_analysis.sh"
else
    RESULT=$(sbatch 10_final_integrative_analysis.sh 2>&1)
    JOB_ID=$(echo "$RESULT" | grep -oP '\d+$')
    if [[ -n "$JOB_ID" ]]; then
        echo "  Submitted: 10_final_integrative_analysis.sh -> Job ID: $JOB_ID"
    else
        echo "  ERROR: $RESULT"
        exit 1
    fi
fi

echo ""
echo "=============================================="
echo "Integration Rerun Submitted!"
echo "=============================================="
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f ${SCRIPT_DIR}/logs/10_final_integrative_analysis.out"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
