#!/bin/bash

#===============================================================================
# SCRIPT: rerun_cutandtag.sh
# PURPOSE: Re-run Cut&Tag pipeline from filtering onwards after code revisions
#
# DESCRIPTION:
# Submits Cut&Tag pipeline jobs with SLURM dependency chains, starting from
# 4_filter.sh (skipping 3_align.sh which uses existing aligned BAMs).
#
# CHANGES TRIGGERING THIS RERUN:
# 1. 4_filter.sh: Mark-only duplicates (REMOVE_DUPLICATES=false) instead of
#    removing them; changed -F 1804 to -F 2828 to keep CUT&Tag biological
#    duplicates (Tn5 insertion at same position ≠ PCR artifact)
# 2. New QC/analysis scripts: 3b_spikein_calibration, 4b_fragment_stratification,
#    5b_peak_calling_seacr, 5c_idr_analysis, 5d_frip_calculation
# 3. New scripts/diffbind_analysis.R for DiffBind
#
# DOWNSTREAM CASCADE:
# Since filtered BAMs change → peaks change → all downstream analyses change
#
# USAGE:
#   ./rerun_cutandtag.sh              # Run pipeline from filtering onwards
#   ./rerun_cutandtag.sh --dry-run    # Show what would be submitted
#===============================================================================

set -e

# Configuration
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
ORIG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_CUTandTAG"
cd "$SCRIPT_DIR"

#===============================================================================
# SETUP: Create symlinks for input data that only exists in original directory
# (Aligned BAMs and trimmed FASTQs were not copied to Eva_submission)
#===============================================================================
echo "=== Setting up symlinks for input data ==="

# Symlink aligned BAMs from original (needed by 4_filter.sh)
if [[ ! -L "results/03_aligned" && -d "results/03_aligned" ]]; then
    # Remove empty directory and replace with symlink
    if [[ -z "$(ls -A results/03_aligned/*.bam 2>/dev/null)" ]]; then
        rm -rf results/03_aligned
        ln -s "${ORIG_DIR}/results/03_aligned" results/03_aligned
        echo "  Symlinked: results/03_aligned -> ${ORIG_DIR}/results/03_aligned"
    else
        echo "  results/03_aligned already has BAM files, skipping symlink"
    fi
elif [[ -L "results/03_aligned" ]]; then
    echo "  results/03_aligned symlink already exists"
else
    ln -s "${ORIG_DIR}/results/03_aligned" results/03_aligned
    echo "  Symlinked: results/03_aligned -> ${ORIG_DIR}/results/03_aligned"
fi

# Symlink trimmed FASTQs from original (needed by 3b_spikein_calibration.sh)
# NOTE: results/02_trimmed/ may contain QC reports (trimming_report.txt, fastqc.html)
# that should NOT be deleted. We symlink individual .fq.gz files instead.
if [[ -z "$(ls -A results/02_trimmed/*.fq.gz 2>/dev/null)" ]]; then
    mkdir -p results/02_trimmed
    for fq in "${ORIG_DIR}"/results/02_trimmed/*.fq.gz; do
        if [[ -f "$fq" ]]; then
            fname=$(basename "$fq")
            if [[ ! -e "results/02_trimmed/$fname" ]]; then
                ln -s "$fq" "results/02_trimmed/$fname"
            fi
        fi
    done
    echo "  Symlinked trimmed FASTQ files from ${ORIG_DIR}/results/02_trimmed/"
else
    echo "  results/02_trimmed already has FASTQ files, skipping symlinks"
fi

echo ""

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
echo "Cut&Tag Pipeline RERUN (Revised Code)"
echo "=============================================="
echo "Script directory: $SCRIPT_DIR"
echo "Original data:    $ORIG_DIR"
echo "Dry run: $DRY_RUN"
echo "Starting from: 4_filter.sh (skipping alignment)"
echo ""

# Function to submit job and capture job ID
submit_job() {
    local script=$1
    local dependency=$2

    if [[ ! -f "$script" ]]; then
        echo "  WARNING: Script not found: $script - skipping"
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
        local job_id=$(echo "$result" | grep -oP '\d+$')
        if [[ -n "$job_id" ]]; then
            echo "  Submitted: $script -> Job ID: $job_id"
            echo "$job_id"
        else
            echo "  ERROR submitting $script: $result"
            return 1
        fi
    fi
}

#===============================================================================
# STAGE 1: BAM Filtering (modified - mark-only duplicates, keep bio-duplicates)
#===============================================================================
echo ""
echo "=== STAGE 1: BAM Filtering (MODIFIED) ==="
echo "  Change: REMOVE_DUPLICATES=false, -F 2828 (keep CUT&Tag bio-duplicates)"

JOB_FILTER=$(submit_job "4_filter.sh" "")

#===============================================================================
# STAGE 1b: Spike-in Calibration (NEW - runs on trimmed FASTQs, parallel)
#===============================================================================
echo ""
echo "=== STAGE 1b: Spike-in Calibration (NEW) ==="

JOB_SPIKEIN=$(submit_job "3b_spikein_calibration.sh" "")

#===============================================================================
# STAGE 2: Peak Calling + BigWig + New QC (depends on filtering)
#===============================================================================
echo ""
echo "=== STAGE 2: Peak Calling, BigWig, and Fragment Analysis ==="

JOB_PEAKS=$(submit_job "5_peak_calling_narrow.sh" "$JOB_FILTER")
JOB_BIGWIG=$(submit_job "6_bigwig.sh" "$JOB_FILTER")

# NEW: SEACR peak calling (alternative to MACS2)
JOB_SEACR=$(submit_job "5b_peak_calling_seacr.sh" "$JOB_FILTER")

# NEW: Fragment size stratification analysis
JOB_FRAG=$(submit_job "4b_fragment_stratification.sh" "$JOB_FILTER")

#===============================================================================
# STAGE 3: Peak Analysis (depends on peaks)
#===============================================================================
echo ""
echo "=== STAGE 3: Peak Analysis and Annotation ==="

PEAKS_DEP="$JOB_PEAKS"
BOTH_DEP="$JOB_PEAKS:$JOB_BIGWIG"

# Downstream of peak calling
JOB_ANNOTATE=$(submit_job "8_annotate_narrow.sh" "$PEAKS_DEP")
JOB_HOMER=$(submit_job "8b_homer_motifs.sh" "$PEAKS_DEP")
JOB_DIFFBIND=$(submit_job "9_diff_bind_narrow.sh" "$PEAKS_DEP")
JOB_COMPARE=$(submit_job "7_compare_narrow.sh" "$BOTH_DEP")
JOB_COMBINE_REP=$(submit_job "11_combine_replicates_narrow.sh" "$BOTH_DEP")
JOB_COMBINE_BW=$(submit_job "13_combine_bigwig.sh" "$JOB_BIGWIG")

# NEW: IDR analysis (depends on peak calling)
JOB_IDR=$(submit_job "5c_idr_analysis.sh" "$PEAKS_DEP")

# NEW: FRiP calculation (depends on filtering + peaks)
JOB_FRIP=$(submit_job "5d_frip_calculation.sh" "$PEAKS_DEP")

#===============================================================================
# STAGE 4: Overlap Analyses (depends on combined replicates)
#===============================================================================
echo ""
echo "=== STAGE 4: Overlap Analyses ==="

JOB_PAIRWISE=$(submit_job "12_pairwise_overlap_narrow.sh" "$JOB_COMBINE_REP")
JOB_SES=$(submit_job "13_ses_overlap.sh" "$JOB_COMBINE_REP")

#===============================================================================
# STAGE 5: SES Overlap Plots (depends on SES overlap)
#===============================================================================
echo ""
echo "=== STAGE 5: SES Overlap Visualization ==="

JOB_SES_PLOTS=$(submit_job "14_ses_overlap_plots.sh" "$JOB_SES")

#===============================================================================
# STAGE 6: Sample Correlation (depends on bigwig)
#===============================================================================
echo ""
echo "=== STAGE 6: Sample Correlation ==="

JOB_CORR=$(submit_job "15_sample_correlation.sh" "$JOB_BIGWIG")

#===============================================================================
# STAGE 7: Final Reports (depends on most steps)
#===============================================================================
echo ""
echo "=== STAGE 7: Summary and Reports ==="

# Build dependency list from all major analysis jobs
ALL_DEPS=""
for jid in $JOB_ANNOTATE $JOB_HOMER $JOB_DIFFBIND $JOB_COMPARE $JOB_COMBINE_REP $JOB_COMBINE_BW $JOB_PAIRWISE $JOB_SES_PLOTS $JOB_IDR $JOB_FRIP $JOB_CORR $JOB_SPIKEIN $JOB_FRAG; do
    if [[ -n "$jid" && "$jid" != DRYRUN_* ]]; then
        ALL_DEPS="${ALL_DEPS:+${ALL_DEPS}:}${jid}"
    fi
done

JOB_SUMMARY=$(submit_job "10_summary.sh" "$ALL_DEPS")
JOB_ADDITIONAL=$(submit_job "16_run_additional_analyses.sh" "$ALL_DEPS")
JOB_MULTIQC=$(submit_job "multiqc.sh" "$ALL_DEPS")

#===============================================================================
# Summary
#===============================================================================
echo ""
echo "=============================================="
echo "Cut&Tag Rerun Submission Complete!"
echo "=============================================="
echo ""
echo "Pipeline starts from: 4_filter.sh (uses existing aligned BAMs)"
echo "New scripts included: 3b_spikein, 4b_fragment, 5b_seacr, 5c_idr, 5d_frip"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f ${SCRIPT_DIR}/logs/*.out"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
