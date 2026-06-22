#!/bin/bash

#===============================================================================
# SCRIPT: master_rerun_all.sh
# PURPOSE: Master orchestration script to run all Cut&Tag pipeline steps
#
# DESCRIPTION:
# This script submits all pipeline steps with proper SLURM dependencies,
# allowing maximum parallelism where possible while respecting data flow.
#
# PIPELINE STRUCTURE:
# Stage 1 (Parallel): QC + Trimming
#   - 1_qc.sh (raw FASTQ QC)
#   - 2_adapter_trim.sh (adapter trimming)
#
# Stage 2: Alignment (depends on trimming)
#   - 3_align.sh
#
# Stage 3: Filtering (depends on alignment)
#   - 4_filter.sh
#
# Stage 4 (Parallel): Peak calling + BigWig (depends on filtering)
#   - 5_peak_calling_narrow.sh
#   - 6_bigwig.sh
#
# Stage 5 (Parallel): Analysis steps (depends on peaks + bigwig)
#   - 7_compare_narrow.sh
#   - 8_annotate_narrow.sh
#   - 8b_homer_motifs.sh
#   - 9_diff_bind_narrow.sh
#   - 11_combine_replicates_narrow.sh
#   - 13_combine_bigwig.sh
#
# Stage 6 (Parallel): Overlap analyses (depends on combine steps)
#   - 12_pairwise_overlap_narrow.sh
#   - 13_ses_overlap.sh
#
# Stage 7: SES plots (depends on SES overlap)
#   - 14_ses_overlap_plots.sh
#
# Stage 8 (Parallel): Final reports (depends on most steps)
#   - 10_summary.sh
#   - 16_run_additional_analyses.sh
#   - multiqc.sh
#
# USAGE:
#   ./master_rerun_all.sh           # Run full pipeline
#   ./master_rerun_all.sh --dry-run # Show what would be submitted
#   ./master_rerun_all.sh --from 5  # Start from stage 5
#===============================================================================

set -e

# Configuration
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
cd "$BASE_DIR"

# Parse arguments
DRY_RUN=false
START_STAGE=1

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --from)
            START_STAGE=$2
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--dry-run] [--from STAGE]"
            exit 1
            ;;
    esac
done

# Create logs directory if needed
mkdir -p logs

echo "=============================================="
echo "Cut&Tag Pipeline Master Orchestration Script"
echo "=============================================="
echo "Base directory: $BASE_DIR"
echo "Starting from stage: $START_STAGE"
echo "Dry run: $DRY_RUN"
echo ""

# Function to submit job and capture job ID
submit_job() {
    local script=$1
    local dependency=$2
    local job_id=""

    if [[ ! -f "$script" ]]; then
        echo "WARNING: Script not found: $script - skipping"
        return 1
    fi

    local dep_arg=""
    if [[ -n "$dependency" ]]; then
        dep_arg="--dependency=afterok:$dependency"
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
        echo "  [DRY-RUN] sbatch $dep_arg $script"
        echo "DRYRUN_$(basename $script .sh)"
    else
        local result=$(sbatch $dep_arg "$script" 2>&1)
        job_id=$(echo "$result" | grep -oP '\d+$')
        if [[ -n "$job_id" ]]; then
            echo "  Submitted: $script -> Job ID: $job_id"
            echo "$job_id"
        else
            echo "  ERROR submitting $script: $result"
            return 1
        fi
    fi
}

# Function to join job IDs with colon for dependencies
join_ids() {
    local IFS=':'
    echo "$*"
}

# Track job IDs for each stage
declare -A STAGE_JOBS

#===============================================================================
# STAGE 1: QC and Trimming (Parallel)
#===============================================================================
if [[ $START_STAGE -le 1 ]]; then
    echo ""
    echo "=== STAGE 1: Quality Control and Adapter Trimming ==="

    JOB_QC=$(submit_job "1_qc.sh" "")
    JOB_TRIM=$(submit_job "2_adapter_trim.sh" "")

    STAGE_JOBS[1]="$JOB_QC:$JOB_TRIM"

    # Stage 2 only depends on trimming (QC is informational)
    STAGE1_DEP="$JOB_TRIM"
else
    echo ""
    echo "=== STAGE 1: Skipped (starting from stage $START_STAGE) ==="
    STAGE1_DEP=""
fi

#===============================================================================
# STAGE 2: Alignment
#===============================================================================
if [[ $START_STAGE -le 2 ]]; then
    echo ""
    echo "=== STAGE 2: Alignment ==="

    JOB_ALIGN=$(submit_job "3_align.sh" "$STAGE1_DEP")
    STAGE_JOBS[2]="$JOB_ALIGN"
    STAGE2_DEP="$JOB_ALIGN"
else
    echo ""
    echo "=== STAGE 2: Skipped (starting from stage $START_STAGE) ==="
    STAGE2_DEP=""
fi

#===============================================================================
# STAGE 3: Filtering
#===============================================================================
if [[ $START_STAGE -le 3 ]]; then
    echo ""
    echo "=== STAGE 3: BAM Filtering ==="

    JOB_FILTER=$(submit_job "4_filter.sh" "$STAGE2_DEP")
    STAGE_JOBS[3]="$JOB_FILTER"
    STAGE3_DEP="$JOB_FILTER"
else
    echo ""
    echo "=== STAGE 3: Skipped (starting from stage $START_STAGE) ==="
    STAGE3_DEP=""
fi

#===============================================================================
# STAGE 4: Peak Calling and BigWig Generation (Parallel)
#===============================================================================
if [[ $START_STAGE -le 4 ]]; then
    echo ""
    echo "=== STAGE 4: Peak Calling and BigWig Generation ==="

    JOB_PEAKS=$(submit_job "5_peak_calling_narrow.sh" "$STAGE3_DEP")
    JOB_BIGWIG=$(submit_job "6_bigwig.sh" "$STAGE3_DEP")

    STAGE_JOBS[4]="$JOB_PEAKS:$JOB_BIGWIG"
    STAGE4_PEAKS_DEP="$JOB_PEAKS"
    STAGE4_BIGWIG_DEP="$JOB_BIGWIG"
    STAGE4_DEP="$JOB_PEAKS:$JOB_BIGWIG"
else
    echo ""
    echo "=== STAGE 4: Skipped (starting from stage $START_STAGE) ==="
    STAGE4_PEAKS_DEP=""
    STAGE4_BIGWIG_DEP=""
    STAGE4_DEP=""
fi

#===============================================================================
# STAGE 5: Analysis Steps (Parallel)
#===============================================================================
if [[ $START_STAGE -le 5 ]]; then
    echo ""
    echo "=== STAGE 5: Peak Analysis ==="

    # Jobs that depend on peaks only
    JOB_ANNOTATE=$(submit_job "8_annotate_narrow.sh" "$STAGE4_PEAKS_DEP")
    JOB_HOMER=$(submit_job "8b_homer_motifs.sh" "$STAGE4_PEAKS_DEP")
    JOB_DIFFBIND=$(submit_job "9_diff_bind_narrow.sh" "$STAGE4_PEAKS_DEP")

    # Jobs that depend on both peaks and bigwig
    JOB_COMPARE=$(submit_job "7_compare_narrow.sh" "$STAGE4_DEP")
    JOB_COMBINE_REP=$(submit_job "11_combine_replicates_narrow.sh" "$STAGE4_DEP")

    # Job that depends on bigwig only
    JOB_COMBINE_BW=$(submit_job "13_combine_bigwig.sh" "$STAGE4_BIGWIG_DEP")

    STAGE_JOBS[5]="$JOB_ANNOTATE:$JOB_HOMER:$JOB_DIFFBIND:$JOB_COMPARE:$JOB_COMBINE_REP:$JOB_COMBINE_BW"
    STAGE5_COMBINE_REP_DEP="$JOB_COMBINE_REP"
else
    echo ""
    echo "=== STAGE 5: Skipped (starting from stage $START_STAGE) ==="
    STAGE5_COMBINE_REP_DEP=""
fi

#===============================================================================
# STAGE 6: Overlap Analyses (Parallel)
#===============================================================================
if [[ $START_STAGE -le 6 ]]; then
    echo ""
    echo "=== STAGE 6: Overlap Analyses ==="

    JOB_PAIRWISE=$(submit_job "12_pairwise_overlap_narrow.sh" "$STAGE5_COMBINE_REP_DEP")
    JOB_SES=$(submit_job "13_ses_overlap.sh" "$STAGE5_COMBINE_REP_DEP")

    STAGE_JOBS[6]="$JOB_PAIRWISE:$JOB_SES"
    STAGE6_SES_DEP="$JOB_SES"
else
    echo ""
    echo "=== STAGE 6: Skipped (starting from stage $START_STAGE) ==="
    STAGE6_SES_DEP=""
fi

#===============================================================================
# STAGE 7: SES Overlap Plots
#===============================================================================
if [[ $START_STAGE -le 7 ]]; then
    echo ""
    echo "=== STAGE 7: SES Overlap Visualization ==="

    JOB_SES_PLOTS=$(submit_job "14_ses_overlap_plots.sh" "$STAGE6_SES_DEP")

    STAGE_JOBS[7]="$JOB_SES_PLOTS"
    STAGE7_DEP="$JOB_SES_PLOTS"
else
    echo ""
    echo "=== STAGE 7: Skipped (starting from stage $START_STAGE) ==="
    STAGE7_DEP=""
fi

#===============================================================================
# STAGE 8: Final Reports (Parallel)
#===============================================================================
if [[ $START_STAGE -le 8 ]]; then
    echo ""
    echo "=== STAGE 8: Final Reports and Summary ==="

    # Collect all previous job IDs for final dependency
    ALL_PRIOR_JOBS=""
    for stage in "${!STAGE_JOBS[@]}"; do
        if [[ -n "${STAGE_JOBS[$stage]}" ]]; then
            if [[ -n "$ALL_PRIOR_JOBS" ]]; then
                ALL_PRIOR_JOBS="$ALL_PRIOR_JOBS:${STAGE_JOBS[$stage]}"
            else
                ALL_PRIOR_JOBS="${STAGE_JOBS[$stage]}"
            fi
        fi
    done

    # Remove empty entries from dependency list
    ALL_PRIOR_JOBS=$(echo "$ALL_PRIOR_JOBS" | sed 's/::/:/g' | sed 's/^://' | sed 's/:$//')

    JOB_SUMMARY=$(submit_job "10_summary.sh" "$ALL_PRIOR_JOBS")
    JOB_ADDITIONAL=$(submit_job "16_run_additional_analyses.sh" "$ALL_PRIOR_JOBS")
    JOB_MULTIQC=$(submit_job "multiqc.sh" "$ALL_PRIOR_JOBS")

    STAGE_JOBS[8]="$JOB_SUMMARY:$JOB_ADDITIONAL:$JOB_MULTIQC"
fi

#===============================================================================
# Summary
#===============================================================================
echo ""
echo "=============================================="
echo "Pipeline Submission Complete!"
echo "=============================================="
echo ""
echo "Jobs submitted per stage:"
for stage in $(echo "${!STAGE_JOBS[@]}" | tr ' ' '\n' | sort -n); do
    echo "  Stage $stage: ${STAGE_JOBS[$stage]}"
done

echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f logs/*.out"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
