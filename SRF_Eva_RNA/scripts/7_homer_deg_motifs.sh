#!/bin/bash

#===============================================================================
# SCRIPT: 7_homer_deg_motifs.sh
# PURPOSE: HOMER de novo motif discovery in promoters of differentially expressed genes
#
# DESCRIPTION:
# Performs de novo motif discovery using HOMER findMotifsGenome.pl on promoter
# regions (TSS ± 2kb) of significantly differentially expressed genes.
# Analyzes upregulated and downregulated genes separately to identify distinct
# regulatory motifs associated with each direction of expression change.
#
# KEY OPERATIONS:
# 1. Extracts significant DEGs from DESeq2 results (padj < 0.05)
# 2. Separates genes into upregulated (log2FC > 0) and downregulated (log2FC < 0)
# 3. Generates promoter BED files using GTF annotation
# 4. Runs HOMER de novo motif discovery for each direction
# 5. Identifies enriched known motifs from HOMER database
# 6. Generates comparative analysis summary
#
# INPUTS:
# - DESeq2 results: results/05_deseq2/deseq2_results_TES_vs_GFP.txt
# - GTF annotation: gencode.v44.annotation.gtf
#
# OUTPUTS:
# - results/07_homer_motifs/upregulated_promoters/
#   - homerResults.html - Main results report
#   - knownResults.txt - Known motif enrichment
#   - homerMotifs.all.motifs - All discovered motifs
# - results/07_homer_motifs/downregulated_promoters/
#   - Same structure as above
# - results/07_homer_motifs/all_degs_promoters/
#   - Combined analysis (optional)
# - results/07_homer_motifs/summary/
#   - Comparative analysis and gene lists
#
# USAGE:
# sbatch scripts/7_homer_deg_motifs.sh
#
# PARAMETERS (can be modified below):
# - PADJ_THRESHOLD: Adjusted p-value threshold (default: 0.05)
# - LOG2FC_THRESHOLD: Optional fold change threshold (default: 0, no filter)
# - PROMOTER_SIZE: Distance upstream/downstream of TSS (default: 2000)
#===============================================================================

#SBATCH --job-name=7_homer_deg_motifs
#SBATCH --account=kubacki.michal
#SBATCH --partition=workq
#SBATCH --mem=32GB
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/7_homer_deg_motifs.err"
#SBATCH --output="logs/7_homer_deg_motifs.out"

set -e
set -o pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

# Base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
DESEQ_RESULTS="${BASE_DIR}/results/05_deseq2/deseq2_results_TES_vs_GFP.txt"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"
OUTPUT_DIR="${BASE_DIR}/results/07_homer_motifs"
TMP_DIR="${OUTPUT_DIR}/tmp"

# Analysis parameters
PADJ_THRESHOLD=0.05
LOG2FC_THRESHOLD=0  # Set >0 to filter by fold change (e.g., 1 for 2-fold)
PROMOTER_SIZE=2000  # TSS ± this value
NCPUS=8

# HOMER parameters
MOTIF_LENGTHS="8,10,12"  # De novo motif lengths
TOP_MOTIFS=25            # Number of top motifs to report
GENOME="hg38"

#===============================================================================
# SETUP
#===============================================================================

cd "${BASE_DIR}"

# Create output directories
mkdir -p "${OUTPUT_DIR}/upregulated_promoters"
mkdir -p "${OUTPUT_DIR}/downregulated_promoters"
mkdir -p "${OUTPUT_DIR}/all_degs_promoters"
mkdir -p "${OUTPUT_DIR}/summary"
mkdir -p "${TMP_DIR}"
mkdir -p logs

echo "=============================================="
echo "HOMER Motif Analysis - DEG Promoters"
echo "=============================================="
echo "Date: $(date)"
echo "Working directory: ${BASE_DIR}"
echo ""
echo "Parameters:"
echo "  Adjusted p-value threshold: ${PADJ_THRESHOLD}"
echo "  Log2FC threshold: ${LOG2FC_THRESHOLD}"
echo "  Promoter size: TSS ± ${PROMOTER_SIZE} bp"
echo "  Motif lengths: ${MOTIF_LENGTHS}"
echo "=============================================="

#===============================================================================
# ACTIVATE ENVIRONMENT
#===============================================================================

echo ""
echo "Activating conda environment..."
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate homer_env

# Verify HOMER is available
if ! command -v findMotifsGenome.pl &> /dev/null; then
    echo "ERROR: HOMER not found in homer_env. Trying genomics_env..."
    conda activate genomics_env
    if ! command -v findMotifsGenome.pl &> /dev/null; then
        echo "ERROR: HOMER not available. Please install HOMER."
        exit 1
    fi
fi

echo "HOMER location: $(which findMotifsGenome.pl)"

#===============================================================================
# STEP 1: EXTRACT GENE TSS COORDINATES FROM GTF
#===============================================================================

echo ""
echo "=============================================="
echo "STEP 1: Extracting gene TSS from GTF annotation"
echo "=============================================="

# Create gene TSS BED file from GTF
# Format: chr, start, end, gene_id, score, strand
TSS_BED="${TMP_DIR}/all_genes_tss.bed"

echo "Parsing GTF file..."
awk -F'\t' '
    $3 == "gene" {
        # Extract gene_id from attributes
        match($9, /gene_id "([^"]+)"/, arr)
        gene_id = arr[1]

        # Remove version number from gene_id
        gsub(/\.[0-9]+$/, "", gene_id)

        chr = $1
        strand = $7

        # TSS is at start for + strand, end for - strand
        if (strand == "+") {
            tss = $4
        } else {
            tss = $5
        }

        # Only keep standard chromosomes
        if (chr ~ /^chr[0-9XY]+$/) {
            print chr"\t"tss"\t"tss"\t"gene_id"\t0\t"strand
        }
    }
' "${GTF_FILE}" | sort -k1,1 -k2,2n -u > "${TSS_BED}"

TSS_COUNT=$(wc -l < "${TSS_BED}")
echo "Extracted ${TSS_COUNT} unique gene TSS positions"

#===============================================================================
# STEP 2: EXTRACT AND CLASSIFY DEGS
#===============================================================================

echo ""
echo "=============================================="
echo "STEP 2: Extracting differentially expressed genes"
echo "=============================================="

# Check input file
if [[ ! -f "${DESEQ_RESULTS}" ]]; then
    echo "ERROR: DESeq2 results file not found: ${DESEQ_RESULTS}"
    exit 1
fi

# Count total genes in DESeq2 results
TOTAL_GENES=$(tail -n +2 "${DESEQ_RESULTS}" | wc -l)
echo "Total genes in DESeq2 results: ${TOTAL_GENES}"

# Extract significant upregulated genes
UP_GENES="${OUTPUT_DIR}/summary/upregulated_genes.txt"
awk -F'\t' -v padj="${PADJ_THRESHOLD}" -v lfc="${LOG2FC_THRESHOLD}" '
    NR > 1 && $8 != "NA" && $8 < padj && $4 > lfc {
        # Remove version from gene_id
        gene_id = $1
        gsub(/\.[0-9]+$/, "", gene_id)
        print gene_id"\t"$2"\t"$4"\t"$8
    }
' "${DESEQ_RESULTS}" > "${UP_GENES}"

UP_COUNT=$(wc -l < "${UP_GENES}")
echo "Upregulated genes (padj < ${PADJ_THRESHOLD}, log2FC > ${LOG2FC_THRESHOLD}): ${UP_COUNT}"

# Extract significant downregulated genes
DOWN_GENES="${OUTPUT_DIR}/summary/downregulated_genes.txt"
awk -F'\t' -v padj="${PADJ_THRESHOLD}" -v lfc="${LOG2FC_THRESHOLD}" '
    NR > 1 && $8 != "NA" && $8 < padj && $4 < -lfc {
        # Remove version from gene_id
        gene_id = $1
        gsub(/\.[0-9]+$/, "", gene_id)
        print gene_id"\t"$2"\t"$4"\t"$8
    }
' "${DESEQ_RESULTS}" > "${DOWN_GENES}"

DOWN_COUNT=$(wc -l < "${DOWN_GENES}")
echo "Downregulated genes (padj < ${PADJ_THRESHOLD}, log2FC < -${LOG2FC_THRESHOLD}): ${DOWN_COUNT}"

# All DEGs
ALL_GENES="${OUTPUT_DIR}/summary/all_degs.txt"
cat "${UP_GENES}" "${DOWN_GENES}" > "${ALL_GENES}"
ALL_COUNT=$(wc -l < "${ALL_GENES}")
echo "Total DEGs: ${ALL_COUNT}"

#===============================================================================
# STEP 3: CREATE PROMOTER BED FILES
#===============================================================================

echo ""
echo "=============================================="
echo "STEP 3: Creating promoter BED files"
echo "=============================================="

create_promoter_bed() {
    local gene_list=$1
    local output_bed=$2
    local description=$3

    echo "  Creating ${description} promoter BED..."

    # Join gene list with TSS coordinates and extend to promoter regions
    awk -F'\t' 'NR==FNR {genes[$1]=1; next} $4 in genes' \
        "${gene_list}" "${TSS_BED}" | \
    awk -v size="${PROMOTER_SIZE}" 'BEGIN{OFS="\t"} {
        start = $2 - size
        if (start < 0) start = 0
        end = $2 + size
        print $1, start, end, $4, $5, $6
    }' | sort -k1,1 -k2,2n > "${output_bed}"

    local count=$(wc -l < "${output_bed}")
    echo "    Generated ${count} promoter regions"
}

# Create promoter BED for each gene set
UP_PROMOTERS="${TMP_DIR}/upregulated_promoters.bed"
DOWN_PROMOTERS="${TMP_DIR}/downregulated_promoters.bed"
ALL_PROMOTERS="${TMP_DIR}/all_degs_promoters.bed"
BG_PROMOTERS="${TMP_DIR}/background_promoters.bed"

# Extract just gene IDs for matching
cut -f1 "${UP_GENES}" > "${TMP_DIR}/up_gene_ids.txt"
cut -f1 "${DOWN_GENES}" > "${TMP_DIR}/down_gene_ids.txt"
cut -f1 "${ALL_GENES}" > "${TMP_DIR}/all_gene_ids.txt"

create_promoter_bed "${TMP_DIR}/up_gene_ids.txt" "${UP_PROMOTERS}" "upregulated"
create_promoter_bed "${TMP_DIR}/down_gene_ids.txt" "${DOWN_PROMOTERS}" "downregulated"
create_promoter_bed "${TMP_DIR}/all_gene_ids.txt" "${ALL_PROMOTERS}" "all DEGs"

# Create background: all gene promoters
echo "  Creating background (all genes) promoter BED..."
awk -v size="${PROMOTER_SIZE}" 'BEGIN{OFS="\t"} {
    start = $2 - size
    if (start < 0) start = 0
    end = $2 + size
    print $1, start, end, $4, $5, $6
}' "${TSS_BED}" | sort -k1,1 -k2,2n > "${BG_PROMOTERS}"

BG_COUNT=$(wc -l < "${BG_PROMOTERS}")
echo "    Background promoters: ${BG_COUNT}"

#===============================================================================
# STEP 4: RUN HOMER MOTIF ANALYSIS
#===============================================================================

echo ""
echo "=============================================="
echo "STEP 4: Running HOMER motif analysis"
echo "=============================================="

run_homer_analysis() {
    local input_bed=$1
    local output_dir=$2
    local description=$3
    local bg_bed=$4

    local region_count=$(wc -l < "${input_bed}")

    echo ""
    echo "  Processing ${description}..."
    echo "    Input regions: ${region_count}"
    echo "    Output directory: ${output_dir}"

    if [[ ${region_count} -lt 50 ]]; then
        echo "    WARNING: Only ${region_count} regions. Minimum recommended: 50"
        echo "    Proceeding anyway, but results may be unreliable."
    fi

    # Run HOMER findMotifsGenome
    # -size given: Use exact BED coordinates
    # -bg: Use all gene promoters as background
    # -p: Number of processors
    # -len: De novo motif lengths
    # -S: Number of top motifs
    findMotifsGenome.pl \
        "${input_bed}" \
        "${GENOME}" \
        "${output_dir}" \
        -size given \
        -bg "${bg_bed}" \
        -p ${NCPUS} \
        -len ${MOTIF_LENGTHS} \
        -S ${TOP_MOTIFS} \
        2>&1 | tee "${output_dir}/homer.log"

    echo "    Completed ${description} motif analysis"

    # Extract summary if knownResults.txt exists
    if [[ -f "${output_dir}/knownResults.txt" ]]; then
        echo ""
        echo "    Top 10 known motifs for ${description}:"
        head -11 "${output_dir}/knownResults.txt" | tail -10 | while read line; do
            echo "      ${line}" | cut -f1,3
        done
    fi
}

# Run analysis for each gene set
echo ""
echo "--- Upregulated Gene Promoters ---"
run_homer_analysis \
    "${UP_PROMOTERS}" \
    "${OUTPUT_DIR}/upregulated_promoters" \
    "upregulated genes" \
    "${BG_PROMOTERS}"

echo ""
echo "--- Downregulated Gene Promoters ---"
run_homer_analysis \
    "${DOWN_PROMOTERS}" \
    "${OUTPUT_DIR}/downregulated_promoters" \
    "downregulated genes" \
    "${BG_PROMOTERS}"

echo ""
echo "--- All DEG Promoters ---"
run_homer_analysis \
    "${ALL_PROMOTERS}" \
    "${OUTPUT_DIR}/all_degs_promoters" \
    "all DEGs" \
    "${BG_PROMOTERS}"

#===============================================================================
# STEP 5: GENERATE SUMMARY REPORT
#===============================================================================

echo ""
echo "=============================================="
echo "STEP 5: Generating summary report"
echo "=============================================="

SUMMARY_FILE="${OUTPUT_DIR}/summary/HOMER_MOTIF_SUMMARY.txt"

cat > "${SUMMARY_FILE}" << EOF
===============================================================================
HOMER MOTIF ANALYSIS SUMMARY
TES vs GFP Differential Expression - Promoter Motifs
===============================================================================

Analysis Date: $(date)
Working Directory: ${BASE_DIR}

PARAMETERS
----------
Adjusted p-value threshold: ${PADJ_THRESHOLD}
Log2FC threshold: ${LOG2FC_THRESHOLD}
Promoter region: TSS ± ${PROMOTER_SIZE} bp
De novo motif lengths: ${MOTIF_LENGTHS}
Genome: ${GENOME}

INPUT DATA
----------
DESeq2 results: ${DESEQ_RESULTS}
Total genes tested: ${TOTAL_GENES}

GENE COUNTS
-----------
Upregulated genes (TES > GFP): ${UP_COUNT}
Downregulated genes (TES < GFP): ${DOWN_COUNT}
Total DEGs: ${ALL_COUNT}

PROMOTER REGIONS
----------------
Upregulated promoters: $(wc -l < "${UP_PROMOTERS}")
Downregulated promoters: $(wc -l < "${DOWN_PROMOTERS}")
All DEG promoters: $(wc -l < "${ALL_PROMOTERS}")
Background (all genes): ${BG_COUNT}

RESULTS LOCATIONS
-----------------
Upregulated: ${OUTPUT_DIR}/upregulated_promoters/
  - homerResults.html (de novo motifs)
  - knownResults.txt (known motif enrichment)

Downregulated: ${OUTPUT_DIR}/downregulated_promoters/
  - homerResults.html (de novo motifs)
  - knownResults.txt (known motif enrichment)

All DEGs: ${OUTPUT_DIR}/all_degs_promoters/
  - homerResults.html (de novo motifs)
  - knownResults.txt (known motif enrichment)

Gene Lists: ${OUTPUT_DIR}/summary/
  - upregulated_genes.txt
  - downregulated_genes.txt
  - all_degs.txt

===============================================================================
TOP KNOWN MOTIFS BY CATEGORY
===============================================================================

EOF

# Add top motifs from each analysis to summary
for direction in upregulated downregulated all_degs; do
    known_file="${OUTPUT_DIR}/${direction}_promoters/knownResults.txt"
    if [[ -f "${known_file}" ]]; then
        echo "" >> "${SUMMARY_FILE}"
        echo "${direction^^} PROMOTERS - Top 15 Known Motifs:" >> "${SUMMARY_FILE}"
        echo "----------------------------------------" >> "${SUMMARY_FILE}"
        head -16 "${known_file}" | tail -15 | \
            awk -F'\t' 'BEGIN{OFS="\t"} {print NR". "$1, "p="$3}' >> "${SUMMARY_FILE}"
    fi
done

# Check for TEAD motifs specifically
echo "" >> "${SUMMARY_FILE}"
echo "===============================================================================" >> "${SUMMARY_FILE}"
echo "TEAD-RELATED MOTIFS (Key TES Targets)" >> "${SUMMARY_FILE}"
echo "===============================================================================" >> "${SUMMARY_FILE}"

for direction in upregulated downregulated all_degs; do
    known_file="${OUTPUT_DIR}/${direction}_promoters/knownResults.txt"
    if [[ -f "${known_file}" ]]; then
        echo "" >> "${SUMMARY_FILE}"
        echo "${direction^^}:" >> "${SUMMARY_FILE}"
        grep -i "TEAD" "${known_file}" | head -5 | \
            awk -F'\t' '{print "  "$1, "p="$3}' >> "${SUMMARY_FILE}" || \
        echo "  (No TEAD motifs found in enriched results)" >> "${SUMMARY_FILE}"
    fi
done

echo "" >> "${SUMMARY_FILE}"
echo "===============================================================================" >> "${SUMMARY_FILE}"
echo "Analysis completed: $(date)" >> "${SUMMARY_FILE}"
echo "===============================================================================" >> "${SUMMARY_FILE}"

# Display summary
cat "${SUMMARY_FILE}"

#===============================================================================
# CLEANUP
#===============================================================================

echo ""
echo "=============================================="
echo "Cleanup"
echo "=============================================="

# Keep promoter BED files in summary for reference
cp "${UP_PROMOTERS}" "${OUTPUT_DIR}/summary/upregulated_promoters.bed"
cp "${DOWN_PROMOTERS}" "${OUTPUT_DIR}/summary/downregulated_promoters.bed"
cp "${ALL_PROMOTERS}" "${OUTPUT_DIR}/summary/all_degs_promoters.bed"

# Remove temporary files
rm -rf "${TMP_DIR}"
echo "Temporary files cleaned up"

#===============================================================================
# COMPLETION
#===============================================================================

echo ""
echo "=============================================="
echo "HOMER Motif Analysis Complete!"
echo "=============================================="
echo ""
echo "Results saved to: ${OUTPUT_DIR}/"
echo ""
echo "Key output files:"
echo "  - ${OUTPUT_DIR}/upregulated_promoters/homerResults.html"
echo "  - ${OUTPUT_DIR}/downregulated_promoters/homerResults.html"
echo "  - ${OUTPUT_DIR}/all_degs_promoters/homerResults.html"
echo "  - ${OUTPUT_DIR}/summary/HOMER_MOTIF_SUMMARY.txt"
echo ""
echo "To view HTML reports, copy to local machine:"
echo "  scp -r user@server:${OUTPUT_DIR}/*.html /local/path/"
echo ""
echo "Completed: $(date)"
