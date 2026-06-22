#!/bin/bash

#===============================================================================
# SCRIPT: rerun_all_revised.sh
# PURPOSE: Master orchestrator to re-run all pipelines after code revisions
#
# DESCRIPTION:
# Coordinates the re-execution of Cut&Tag, RNA-seq, and integration pipelines
# after the code revisions in commits 9230d72 and 7a6d72f.
#
# EXECUTION ORDER:
#   1. Cut&Tag (from 4_filter.sh) and RNA-seq (from 4_quantify.sh) in PARALLEL
#   2. Integration analysis AFTER both pipelines complete
#
# CHANGES SUMMARY:
# Cut&Tag:
#   - 4_filter.sh: Mark-only duplicates (bio-duplicates preserved for CUT&Tag)
#   - New: 3b_spikein, 4b_fragment, 5b_seacr, 5c_idr, 5d_frip
#   - All downstream re-generated (peaks, annotation, DiffBind, etc.)
#
# RNA-seq:
#   - 4_quantify.sh: Strandedness auto-detection from STAR counts
#   - deseq2_analysis.R: Improved filtering + apeglm LFC shrinkage
#   - 6_go/gsea: Text scaling and parameter improvements
#   - New: 3b_rseqc.sh (comprehensive QC)
#
# Integration:
#   - 10_final_integrative_analysis.R: LFC threshold + distance sensitivity analysis
#
# USAGE:
#   ./rerun_all_revised.sh              # Run everything
#   ./rerun_all_revised.sh --dry-run    # Preview without submitting
#===============================================================================

set -e

PROJ_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission"
cd "$PROJ_DIR"

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

DRY_FLAG=""
if [[ "$DRY_RUN" == "true" ]]; then
    DRY_FLAG="--dry-run"
fi

echo "=============================================="
echo "MASTER RERUN: All Pipelines (Revised Code)"
echo "=============================================="
echo "Project directory: $PROJ_DIR"
echo "Dry run: $DRY_RUN"
echo "Timestamp: $(date)"
echo ""
echo "Execution plan:"
echo "  Phase 1 (parallel): Cut&Tag + RNA-seq pipelines"
echo "  Phase 2 (sequential): Integration (after both complete)"
echo ""

#===============================================================================
# PHASE 1: Cut&Tag and RNA-seq in parallel
#===============================================================================
echo "###############################################"
echo "# PHASE 1: Cut&Tag Pipeline (from filtering) #"
echo "###############################################"
echo ""

ALL_JOB_IDS=""

if [[ "$DRY_RUN" == "true" ]]; then
    bash rerun_cutandtag.sh --dry-run
else
    # Capture output and extract ALL job IDs
    CUTANDTAG_OUTPUT=$(bash rerun_cutandtag.sh 2>&1)
    echo "$CUTANDTAG_OUTPUT"
    # Extract every Job ID from the output
    CUTANDTAG_IDS=$(echo "$CUTANDTAG_OUTPUT" | grep -oP 'Job ID: \K\d+' | tr '\n' ':' | sed 's/:$//')
    if [[ -n "$CUTANDTAG_IDS" ]]; then
        ALL_JOB_IDS="$CUTANDTAG_IDS"
        echo ""
        echo "  Cut&Tag job IDs: $CUTANDTAG_IDS"
    fi
fi

echo ""
echo "###############################################"
echo "# PHASE 1: RNA-seq Pipeline (from quantify)  #"
echo "###############################################"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    bash rerun_rnaseq.sh --dry-run
else
    RNASEQ_OUTPUT=$(bash rerun_rnaseq.sh 2>&1)
    echo "$RNASEQ_OUTPUT"
    RNASEQ_IDS=$(echo "$RNASEQ_OUTPUT" | grep -oP 'Job ID: \K\d+' | tr '\n' ':' | sed 's/:$//')
    if [[ -n "$RNASEQ_IDS" ]]; then
        ALL_JOB_IDS="${ALL_JOB_IDS:+${ALL_JOB_IDS}:}${RNASEQ_IDS}"
        echo ""
        echo "  RNA-seq job IDs: $RNASEQ_IDS"
    fi
fi

#===============================================================================
# PHASE 2: Integration (after both pipelines complete)
#===============================================================================
echo ""
echo "###############################################"
echo "# PHASE 2: Integration Analysis              #"
echo "###############################################"
echo ""

INTEGRATION_DIR="$PROJ_DIR/SRF_Eva_integrated_analysis/scripts"
mkdir -p "$INTEGRATION_DIR/logs"

if [[ "$DRY_RUN" == "true" ]]; then
    echo "  [DRY-RUN] Would submit integration after Cut&Tag + RNA-seq complete"
    echo "  [DRY-RUN] sbatch --dependency=afterok:<all_pipeline_jobs> $INTEGRATION_DIR/10_final_integrative_analysis.sh"
elif [[ -n "$ALL_JOB_IDS" ]]; then
    echo "  Waiting for all upstream jobs: $ALL_JOB_IDS"
    RESULT=$(sbatch --dependency=afterok:$ALL_JOB_IDS "$INTEGRATION_DIR/10_final_integrative_analysis.sh" 2>&1)
    JOB_ID=$(echo "$RESULT" | grep -oP '\d+$')
    if [[ -n "$JOB_ID" ]]; then
        echo "  Submitted: 10_final_integrative_analysis.sh -> Job ID: $JOB_ID"
    else
        echo "  ERROR: $RESULT"
    fi
else
    echo "  ERROR: No upstream job IDs captured. Cannot submit integration with dependencies."
    echo "  Please run the pipelines individually first, then run ./rerun_integration.sh"
    exit 1
fi

#===============================================================================
# Summary
#===============================================================================
echo ""
echo "=============================================="
echo "Master Rerun Submission Complete!"
echo "=============================================="
echo "Timestamp: $(date)"
echo ""
echo "Monitor all jobs:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs:"
echo "  tail -f $PROJ_DIR/SRF_Eva_Cut_and_Tag/logs/*.out"
echo "  tail -f $PROJ_DIR/SRF_Eva_RNA/logs/*.out"
echo "  tail -f $PROJ_DIR/SRF_Eva_integrated_analysis/scripts/logs/*.out"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
