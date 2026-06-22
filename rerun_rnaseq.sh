#!/bin/bash

#===============================================================================
# SCRIPT: rerun_rnaseq.sh
# PURPOSE: Re-run RNA-seq pipeline from quantification onwards after code revisions
#
# DESCRIPTION:
# Submits RNA-seq pipeline jobs with SLURM dependency chains:
#   1. 3b_rseqc.sh (NEW - comprehensive QC, standalone)
#   2. 4_quantify.sh (MODIFIED - strandedness auto-detection)
#   3. rerun_deseq2.sh (runs UPDATED R script directly, not 5_deseq2.sh)
#   4. 6_go_enrichment.sh + 6_gsea_analysis.sh (MODIFIED - text scaling, params)
#   5. volcano_plot_publication.sh + heatmap_publication.sh
#
# CHANGES TRIGGERING THIS RERUN:
# 1. 4_quantify.sh: Auto-detects library strandedness from STAR counts
#    (selects column 2/3/4 based on >70% threshold)
# 2. deseq2_analysis.R: Improved filtering (>=10 counts in >=3 samples),
#    apeglm LFC shrinkage, alpha=0.05
# 3. 6_go_enrichment.R: TEXT_SCALE parameter for publication readability
# 4. 6_gsea_analysis.R: Parameter adjustments
# 5. 3b_rseqc.sh: NEW comprehensive QC (RSeQC metrics)
#
# IMPORTANT: Uses rerun_deseq2.sh instead of 5_deseq2.sh to avoid the
# heredoc overwriting the updated deseq2_analysis.R with the old version.
#
# USAGE:
#   ./rerun_rnaseq.sh              # Run pipeline
#   ./rerun_rnaseq.sh --dry-run    # Show what would be submitted
#===============================================================================

set -e

# Configuration
SCRIPT_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
ORIG_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/SRF_Eva_RNA"
cd "$SCRIPT_DIR"

#===============================================================================
# SETUP: Create symlinks for input data that only exists in original directory
# (Aligned BAMs and STAR counts were not copied to Eva_submission)
#===============================================================================
echo "=== Setting up symlinks for input data ==="

# Symlink aligned BAMs + STAR counts from original (needed by 3b_rseqc, 4_quantify)
# NOTE: results/03_aligned/ may contain subdirectories with QC stats that should be kept.
# We symlink individual BAM and counts.tab files instead of replacing the entire directory.
if [[ -z "$(ls -A results/03_aligned/*_sorted.bam 2>/dev/null)" ]]; then
    mkdir -p results/03_aligned
    for f in "${ORIG_DIR}"/results/03_aligned/*_sorted.bam "${ORIG_DIR}"/results/03_aligned/*_counts.tab; do
        if [[ -f "$f" ]]; then
            fname=$(basename "$f")
            if [[ ! -e "results/03_aligned/$fname" ]]; then
                ln -s "$f" "results/03_aligned/$fname"
            fi
        fi
    done
    echo "  Symlinked BAM and counts.tab files from ${ORIG_DIR}/results/03_aligned/"
else
    echo "  results/03_aligned already has BAM files, skipping symlinks"
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
echo "RNA-seq Pipeline RERUN (Revised Code)"
echo "=============================================="
echo "Script directory: $SCRIPT_DIR"
echo "Original data:    $ORIG_DIR"
echo "Dry run: $DRY_RUN"
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
# STAGE 1: New QC + Quantification (parallel - independent)
#===============================================================================
echo ""
echo "=== STAGE 1a: RSeQC Quality Control (NEW) ==="
echo "  Strandedness validation, gene body coverage, read distribution"

JOB_RSEQC=$(submit_job "scripts/3b_rseqc.sh" "")

echo ""
echo "=== STAGE 1b: Quantification (MODIFIED) ==="
echo "  Change: Auto-detects library strandedness from STAR counts"

JOB_QUANTIFY=$(submit_job "scripts/4_quantify.sh" "")

#===============================================================================
# STAGE 2: DESeq2 with updated R script
#===============================================================================
echo ""
echo "=== STAGE 2: DESeq2 Differential Expression (UPDATED R script) ==="
echo "  Changes: Improved filtering, apeglm LFC shrinkage, alpha=0.05"
echo "  NOTE: Using rerun_deseq2.sh (not 5_deseq2.sh) to preserve updated R code"

JOB_DESEQ2=$(submit_job "scripts/rerun_deseq2.sh" "$JOB_QUANTIFY")

#===============================================================================
# STAGE 3: Functional Analysis (depends on DESeq2)
#===============================================================================
echo ""
echo "=== STAGE 3: Functional Enrichment (MODIFIED) ==="

JOB_GO=$(submit_job "scripts/6_go_enrichment.sh" "$JOB_DESEQ2")
JOB_GSEA=$(submit_job "scripts/6_gsea_analysis.sh" "$JOB_DESEQ2")

#===============================================================================
# STAGE 4: Publication Visualizations (depends on DESeq2)
#===============================================================================
echo ""
echo "=== STAGE 4: Publication Visualizations ==="

JOB_VOLCANO=$(submit_job "scripts/volcano_plot_publication.sh" "$JOB_DESEQ2")
JOB_HEATMAP=$(submit_job "scripts/heatmap_publication.sh" "$JOB_DESEQ2")

#===============================================================================
# Summary
#===============================================================================
echo ""
echo "=============================================="
echo "RNA-seq Rerun Submission Complete!"
echo "=============================================="
echo ""
echo "Dependency chain:"
echo "  3b_rseqc.sh (standalone QC, no downstream)"
echo "  4_quantify.sh → rerun_deseq2.sh → 6_go + 6_gsea + volcano + heatmap"
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo "  tail -f ${SCRIPT_DIR}/logs/*.out"
echo ""

if [[ "$DRY_RUN" == "true" ]]; then
    echo "NOTE: This was a dry run. No jobs were actually submitted."
fi
