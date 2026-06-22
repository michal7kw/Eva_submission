#!/bin/bash

#===============================================================================
# SCRIPT: 10_summary.sh
# PURPOSE: Generate comprehensive MultiQC report and dual-track analysis summary
#
# DESCRIPTION:
# This script creates a comprehensive summary report of the entire Cut&Tag
# analysis pipeline using MultiQC for quality control aggregation and
# generates a detailed analysis summary document with key findings from both
# narrow and broad peak analysis tracks, conclusions, and next steps.
#
# KEY OPERATIONS:
# 1. MultiQC report generation from all pipeline outputs
# 2. Aggregation of QC metrics across all analysis steps
# 3. Creation of comprehensive dual-track analysis summary document
# 4. Integration of peak statistics from both narrow and broad peak analyses
# 5. Comparison between narrow and broad peak methodologies
# 6. Documentation of experimental design and conclusions
#
# METHODOLOGY:
# - Uses MultiQC to parse and aggregate QC reports
# - Combines results from FastQC, alignment, peak calling, etc.
# - Creates publication-ready HTML report
# - Generates markdown summary with biological interpretation
# - Compares narrow vs broad peak analysis results
# - Includes recommendations for follow-up experiments
#
# IMPORTANT PARAMETERS:
# - Memory: 8GB (sufficient for report generation)
# - Time: 1 hour (adequate for MultiQC processing)
# - Threads: 4 (parallel processing for report generation)
# - Output format: HTML report with interactive plots
#
# INPUTS:
# - All QC outputs from previous pipeline steps
# - FastQC reports, alignment statistics, peak counts
# - Analysis results from both narrow and broad peak tracks
# - Peak statistics from comparative analysis
#
# OUTPUTS:
# - MultiQC HTML report (TES_TEAD1_CutTag_Report.html)
# - Comprehensive dual-track analysis summary (ANALYSIS_SUMMARY.md)
# - Integrated QC metrics and visualizations
# - Biological interpretation and track comparisons
# - Recommendations for follow-up studies
#
# DEPENDENCIES:
# - MultiQC for report aggregation
# - Conda environment: multiqc
# - Both narrow and broad pipeline outputs
#
# USAGE:
# sbatch 10_summary.sh
#
# NOTES:
# - Final step of both analysis pipelines
# - Provides comprehensive dual-track project overview
# - Includes biological interpretation of results
# - Compares narrow vs broad peak methodologies
# - Serves as basis for manuscript preparation
# - Documents experimental design and methodology
#===============================================================================

#SBATCH --job-name=10_report
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/10_report.err"
#SBATCH --output="./logs/10_report.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate multiqc

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"

echo "=========================================="
echo "Starting comprehensive dual-track summary"
echo "Date: $(date)"
echo "=========================================="

# Create output directory
mkdir -p ${BASE_DIR}/results/report

# Run MultiQC to aggregate all QC reports
echo "Generating MultiQC report..."
multiqc \
    ${BASE_DIR}/results \
    -o ${BASE_DIR}/results/report \
    -n TES_TEAD1_CutTag_Report

# Function to safely append file content or placeholder
append_file_or_placeholder() {
    local file_path=$1
    local placeholder_msg=$2
    local output_file=$3

    if [ -f "$file_path" ]; then
        cat "$file_path" >> "$output_file"
    else
        echo "$placeholder_msg" >> "$output_file"
    fi
}

# Create comprehensive dual-track summary report
echo "Creating comprehensive analysis summary..."
cat > ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'
# Cut&Tag Dual-Track Analysis Summary: TES/TEAD1 Epigenetic Study

## Experimental Design
- **TES**: Epigenetic silencer with TEAD1 binding domain + KRAB + DNMT3A/3L
- **TEAD1**: Endogenous TEAD1 for comparison
- **Controls**: IggMs (for TES), IggRb (for TEAD1)
- **NOTE**: TESmut samples excluded from analysis (failed sample)

## Analysis Methodology
This study employed dual-track peak calling methodology:
- **Track A (Narrow Peaks)**: Precise transcription factor binding site identification
- **Track B (Broad Peaks)**: Enhanced overlap detection and chromatin domain analysis

## Track A: Narrow Peak Analysis Results

### Peak Statistics (Narrow)
EOF

# Add narrow peak statistics
append_file_or_placeholder "${BASE_DIR}/results/07_analysis_narrow/peak_stats.txt" \
    "*Narrow peak statistics not available - run narrow peak pipeline*" \
    "${BASE_DIR}/results/ANALYSIS_SUMMARY.md"

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

### Narrow Peak Key Findings
EOF

# Add narrow peak analysis results if available
if [ -f "${BASE_DIR}/results/07_analysis_narrow/TES_TEAD1_common_genes.txt" ]; then
    narrow_common_genes=$(wc -l < "${BASE_DIR}/results/07_analysis_narrow/TES_TEAD1_common_genes.txt")
    echo "- Common target genes (TES/TEAD1): $narrow_common_genes" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/07_analysis_narrow/TES_specific_genes.txt" ]; then
    narrow_tes_specific=$(wc -l < "${BASE_DIR}/results/07_analysis_narrow/TES_specific_genes.txt")
    echo "- TES-specific target genes: $narrow_tes_specific" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/07_analysis_narrow/TEAD1_specific_genes.txt" ]; then
    narrow_tead1_specific=$(wc -l < "${BASE_DIR}/results/07_analysis_narrow/TEAD1_specific_genes.txt")
    echo "- TEAD1-specific target genes: $narrow_tead1_specific" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

## Track B: Broad Peak Analysis Results

### Peak Statistics (Broad)
EOF

# Add broad peak statistics
append_file_or_placeholder "${BASE_DIR}/results/07_analysis_broad/peak_stats.txt" \
    "*Broad peak statistics not available - run broad peak pipeline*" \
    "${BASE_DIR}/results/ANALYSIS_SUMMARY.md"

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

### Broad Peak Key Findings
EOF

# Add broad peak analysis results if available
if [ -f "${BASE_DIR}/results/07_analysis_broad/TES_TEAD1_common_genes.txt" ]; then
    broad_common_genes=$(wc -l < "${BASE_DIR}/results/07_analysis_broad/TES_TEAD1_common_genes.txt")
    echo "- Common target genes (TES/TEAD1): $broad_common_genes" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/07_analysis_broad/TES_specific_genes.txt" ]; then
    broad_tes_specific=$(wc -l < "${BASE_DIR}/results/07_analysis_broad/TES_specific_genes.txt")
    echo "- TES-specific target genes: $broad_tes_specific" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/07_analysis_broad/TEAD1_specific_genes.txt" ]; then
    broad_tead1_specific=$(wc -l < "${BASE_DIR}/results/07_analysis_broad/TEAD1_specific_genes.txt")
    echo "- TEAD1-specific target genes: $broad_tead1_specific" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

## Methodology Comparison

### Peak Calling Strategy Comparison
| Analysis Type | Narrow Peaks | Broad Peaks |
|---------------|--------------|-------------|
| **Purpose** | Precise TF binding sites | Chromatin domain detection |
| **Sensitivity** | High specificity | High sensitivity |
| **Overlap Detection** | Moderate | Enhanced |
| **Best For** | Motif analysis | Chromatin accessibility |

### Track Performance Summary
EOF

# Generate track comparison
echo "Comparing analysis tracks..." >&2

# Narrow peak counts
if [ -f "${BASE_DIR}/results/05_peaks_narrow/TES_peaks.narrowPeak" ]; then
    narrow_tes=$(wc -l < "${BASE_DIR}/results/05_peaks_narrow/TES_peaks.narrowPeak")
    echo "- Narrow TES peaks: $narrow_tes" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/05_peaks_narrow/TEAD1_peaks.narrowPeak" ]; then
    narrow_tead1=$(wc -l < "${BASE_DIR}/results/05_peaks_narrow/TEAD1_peaks.narrowPeak")
    echo "- Narrow TEAD1 peaks: $narrow_tead1" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

# Broad peak counts
if [ -f "${BASE_DIR}/results/05_peaks_broad/TES_peaks.broadPeak" ]; then
    broad_tes=$(wc -l < "${BASE_DIR}/results/05_peaks_broad/TES_peaks.broadPeak")
    echo "- Broad TES peaks: $broad_tes" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

if [ -f "${BASE_DIR}/results/05_peaks_broad/TEAD1_peaks.broadPeak" ]; then
    broad_tead1=$(wc -l < "${BASE_DIR}/results/05_peaks_broad/TEAD1_peaks.broadPeak")
    echo "- Broad TEAD1 peaks: $broad_tead1" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

## Consensus Analysis Results (Improved Pipeline)

### Key Features of Improved Consensus Analysis:
- **Summit-based consensus**: Precise TF binding site identification
- **Quality-weighted peaks**: Signal strength prioritization
- **IDR reproducibility**: ENCODE-standard statistical framework
- **Multiple peak sets**: Different stringency levels for various analyses

EOF

# Add consensus results if available
if [ -f "${BASE_DIR}/results/11_combined_replicates_narrow/REPLICATE_COMBINATION_SUMMARY_NARROW.txt" ]; then
    echo "### Narrow Peak Consensus (Improved Pipeline)" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    echo '```' >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    append_file_or_placeholder "${BASE_DIR}/results/11_combined_replicates_narrow/REPLICATE_COMBINATION_SUMMARY_NARROW.txt" \
        "*Narrow consensus analysis not available*" \
        "${BASE_DIR}/results/ANALYSIS_SUMMARY.md"
    echo '```' >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    echo "" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

# Add IDR statistics if available
CONSENSUS_NARROW="${BASE_DIR}/results/11_combined_replicates_narrow"
if [ -d "${CONSENSUS_NARROW}/idr" ]; then
    echo "### IDR Reproducibility Metrics (Narrow Peaks)" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    echo "" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md

    for condition in TES TEAD1; do  # TESmut removed - failed sample
        if [ -f "${CONSENSUS_NARROW}/idr/${condition}_rep1_vs_rep2_idr.txt" ]; then
            idr_12=$(wc -l < "${CONSENSUS_NARROW}/idr/${condition}_rep1_vs_rep2_idr.txt" 2>/dev/null || echo 0)
            idr_13=$(wc -l < "${CONSENSUS_NARROW}/idr/${condition}_rep1_vs_rep3_idr.txt" 2>/dev/null || echo 0)
            idr_23=$(wc -l < "${CONSENSUS_NARROW}/idr/${condition}_rep2_vs_rep3_idr.txt" 2>/dev/null || echo 0)
            idr_union=$(wc -l < "${CONSENSUS_NARROW}/idr/${condition}_idr_union.bed" 2>/dev/null || echo 0)
            idr_intersect=$(wc -l < "${CONSENSUS_NARROW}/idr/${condition}_idr_intersection.bed" 2>/dev/null || echo 0)

            echo "**${condition} IDR Results:**" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "- Rep1 vs Rep2: $idr_12 peaks" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "- Rep1 vs Rep3: $idr_13 peaks" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "- Rep2 vs Rep3: $idr_23 peaks" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "- IDR union: $idr_union peaks" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "- IDR intersection (most stringent): $idr_intersect peaks" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
            echo "" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
        fi
    done
fi

if [ -f "${BASE_DIR}/results/11_combined_replicates_broad/REPLICATE_COMBINATION_SUMMARY_BROAD.txt" ]; then
    echo "### Broad Peak Consensus" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    echo '```' >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    append_file_or_placeholder "${BASE_DIR}/results/11_combined_replicates_broad/REPLICATE_COMBINATION_SUMMARY_BROAD.txt" \
        "*Broad consensus analysis not available*" \
        "${BASE_DIR}/results/ANALYSIS_SUMMARY.md"
    echo '```' >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    echo "" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
fi

cat >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md << 'EOF'

## Biological Conclusions
1. **TES binding specificity**: Compare overlap between TES and TEAD1 peaks across both methodologies
2. **TESmut validation**: Confirm loss of binding in mutant using both narrow and broad peak approaches
3. **Target genes**: Analyze common targets between TES and TEAD1 from both analysis tracks
4. **Methodological insights**: Broad peaks provide enhanced sensitivity for overlap detection
5. **Epigenetic modifications**: Validate KRAB/DNMT activity at target sites identified by both methods

## Recommendations

### Peak Set Selection for Downstream Analysis

#### For Motif Analysis:
- **Primary**: Summit-based consensus peaks (`*_consensus_summits.bed`)
- **Alternative**: IDR intersection peaks (most stringent)

#### For Functional Studies (RNA-seq correlation, etc.):
- **Primary**: Quality-scored consensus peaks (`*_consensus_peaks_scored.bed`) ⭐ RECOMMENDED
- **Alternative**: IDR union peaks (permissive but reproducible)

#### For Publication:
- **Primary**: IDR-based peaks (intersection for stringent, union for permissive)
- **Include**: Reproducibility statistics from IDR analysis

### Analysis Track Selection
- **Use Narrow Peaks for**: Precise motif analysis, specific binding site identification
- **Use Broad Peaks for**: Overlap analysis, chromatin domain studies, consensus calling
- **Use Both Tracks for**: Comprehensive analysis and methodological validation

### Quality Thresholds:
- **Excellent data**: IDR reproducibility >70%
- **Good data**: IDR reproducibility 50-70%
- **Acceptable data**: IDR reproducibility 30-50%
- **Poor data**: IDR reproducibility <30% (consider repeating)

### Next Steps
1. Validate key targets by ChIP-qPCR using IDR intersection peaks (highest confidence)
2. RNA-seq to measure transcriptional repression at quality-scored consensus targets
3. Bisulfite sequencing for DNA methylation analysis at high-confidence IDR sites
4. ATAC-seq to assess chromatin accessibility changes using consensus peaks
5. Motif analysis using summit-based consensus peaks
6. Compare results between different peak stringency levels

## Output Files Directory Structure

### Shared Components
- FastQC reports: `results/01_fastqc/`
- Trimmed reads: `results/02_trimmed/`
- Aligned BAMs: `results/03_aligned/`
- Filtered BAMs: `results/04_filtered/`
- BigWig tracks: `results/06_bigwig/`
- MultiQC report: `results/report/`

### Track A: Narrow Peak Analysis
- Peak files: `results/05_peaks_narrow/`
- Analysis results: `results/07_analysis_narrow/`
- Consensus results: `results/11_combined_replicates_narrow/`

### Track B: Broad Peak Analysis
- Peak files: `results/05_peaks_broad/`
- Analysis results: `results/07_analysis_broad/`
- Consensus results: `results/11_combined_replicates_broad/`

### Comparative Analysis
- Pairwise overlaps: `results/12_pairwise_overlap/`
- SES overlaps: `results/13_ses_overlap/`

## Analysis Quality Metrics
EOF

# Add quality metrics summary
echo "### Pipeline Completion Status" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md

# Check for key output files
check_file_status() {
    local file_path=$1
    local description=$2
    if [ -f "$file_path" ]; then
        echo "✅ $description: Completed" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    else
        echo "❌ $description: Not completed" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
    fi
}

check_file_status "${BASE_DIR}/results/05_peaks_narrow/TES_peaks.narrowPeak" "Narrow peak calling"
check_file_status "${BASE_DIR}/results/05_peaks_broad/TES_peaks.broadPeak" "Broad peak calling"
check_file_status "${BASE_DIR}/results/07_analysis_narrow/peak_stats.txt" "Narrow peak analysis"
check_file_status "${BASE_DIR}/results/07_analysis_broad/peak_stats.txt" "Broad peak analysis"
check_file_status "${BASE_DIR}/results/12_pairwise_overlap/PAIRWISE_OVERLAP_SUMMARY.md" "Pairwise overlap analysis"

echo "" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
echo "*Analysis completed on: $(date)*" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md
echo "*Pipeline version: Dual-track Cut&Tag analysis*" >> ${BASE_DIR}/results/ANALYSIS_SUMMARY.md

echo "Analysis pipeline complete! Check results/ANALYSIS_SUMMARY.md for summary."