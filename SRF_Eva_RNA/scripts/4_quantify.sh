#!/bin/bash

#===============================================================================
# SCRIPT: 4_quantify.sh
# PURPOSE: Gene Expression Quantification for RNA-seq Data
#
# DESCRIPTION:
# This script performs gene expression quantification using featureCounts
# or combines STAR gene counts into a unified count matrix. It is the fourth
# step in the RNA-seq analysis pipeline and prepares count data for
# differential expression analysis.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with Subread (featureCounts)
# 2. Either runs featureCounts on BAM files OR combines STAR gene counts
# 3. Generates a unified count matrix for all samples
# 4. Creates sample metadata file
# 5. Performs basic count statistics and quality checks
#
# METHODOLOGY:
# - Option A: Use featureCounts for precise gene quantification
# - Option B: Combine STAR gene counts (faster, already generated)
# - Generates count matrix suitable for DESeq2 analysis
# - Includes gene annotation and sample metadata
#
# EXPECTED INPUTS:
# - BAM files: ${ALIGNED_DIR}/${SAMPLE}_sorted.bam
# - OR STAR counts: ${ALIGNED_DIR}/${SAMPLE}_counts.tab
# - GTF annotation: ${GTF_FILE}
# - Sample design: ${BASE_DIR}/config/design.txt
#
# EXPECTED OUTPUTS:
# - Count matrix: ${OUTPUT_DIR}/count_matrix.txt
# - Sample metadata: ${OUTPUT_DIR}/sample_metadata.txt
# - Summary statistics: ${OUTPUT_DIR}/quantification_summary.txt
# - SLURM logs: ./logs/4_quantify.out and ./logs/4_quantify.err
#
# DEPENDENCIES:
# - Subread (featureCounts) via conda environment
# - R (for matrix processing)
# - GTF annotation file
# - SLURM workload manager
#
# USAGE:
# sbatch 4_quantify.sh
#===============================================================================

#SBATCH --job-name=4_quantify
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/4_quantify.err"
#SBATCH --output="./logs/4_quantify.out"

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
ALIGNED_DIR="${BASE_DIR}/results/03_aligned"
OUTPUT_DIR="${BASE_DIR}/results/04_quantified"

# Annotation file
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

echo "=== Starting gene expression quantification ==="
echo "Timestamp: $(date)"
echo "Input directory: ${ALIGNED_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "GTF file: ${GTF_FILE}"
echo ""

# Read samples
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))

#===============================================================================
# STRANDEDNESS CHECK
# Determine library strandedness using the first sample's BAM file.
# STAR ReadsPerGene.out.tab columns: 1=gene_id, 2=unstranded, 3=sense, 4=antisense
# Most modern RNA-seq kits are stranded (column 4 for dUTP-based reverse-strand).
# We infer strandedness by comparing STAR's N_noFeature counts across columns.
#===============================================================================
echo "=== Checking library strandedness ==="
FIRST_SAMPLE=${SAMPLES[0]}
FIRST_COUNT_FILE="${ALIGNED_DIR}/${FIRST_SAMPLE}_counts.tab"

if [[ -f "${FIRST_COUNT_FILE}" ]]; then
    # Extract summary rows from STAR ReadsPerGene.out.tab
    # Column 2: unstranded, Column 3: sense (++,--), Column 4: antisense (+-,-+)
    UNSTRANDED=$(awk '$1=="N_noFeature" {print $2}' "${FIRST_COUNT_FILE}")
    SENSE_NOFEAT=$(awk '$1=="N_noFeature" {print $3}' "${FIRST_COUNT_FILE}")
    ANTISENSE_NOFEAT=$(awk '$1=="N_noFeature" {print $4}' "${FIRST_COUNT_FILE}")

    # Compare total gene counts across columns to infer strandedness
    TOTAL_UNSTRANDED=$(awk '$1!~/^N_/ {sum+=$2} END {print sum}' "${FIRST_COUNT_FILE}")
    TOTAL_SENSE=$(awk '$1!~/^N_/ {sum+=$3} END {print sum}' "${FIRST_COUNT_FILE}")
    TOTAL_ANTISENSE=$(awk '$1!~/^N_/ {sum+=$4} END {print sum}' "${FIRST_COUNT_FILE}")

    echo "  Sample: ${FIRST_SAMPLE}"
    echo "  Total gene counts - Unstranded: ${TOTAL_UNSTRANDED}, Sense: ${TOTAL_SENSE}, Antisense: ${TOTAL_ANTISENSE}"

    # If sense or antisense captures >70% of unstranded counts, library is stranded
    if [[ ${TOTAL_UNSTRANDED} -eq 0 ]]; then
        echo "  WARNING: No unstranded counts detected. Defaulting to unstranded."
        STRAND_COL=2
    else
        SENSE_PCT=$((TOTAL_SENSE * 100 / TOTAL_UNSTRANDED))
        ANTISENSE_PCT=$((TOTAL_ANTISENSE * 100 / TOTAL_UNSTRANDED))

        echo "  Sense/Unstranded: ${SENSE_PCT}%, Antisense/Unstranded: ${ANTISENSE_PCT}%"

        if [[ ${ANTISENSE_PCT} -gt 70 ]]; then
            STRAND_COL=4  # Reverse-stranded (dUTP, most common)
            echo "  RESULT: Library is REVERSE-STRANDED (using column 4 = antisense)"
        elif [[ ${SENSE_PCT} -gt 70 ]]; then
            STRAND_COL=3  # Forward-stranded
            echo "  RESULT: Library is FORWARD-STRANDED (using column 3 = sense)"
        else
            STRAND_COL=2  # Unstranded
            echo "  RESULT: Library is UNSTRANDED (using column 2 = unstranded)"
        fi
    fi
else
    echo "  WARNING: Cannot check strandedness, count file not found. Defaulting to unstranded."
    STRAND_COL=2
fi

echo "  Using STAR count column: ${STRAND_COL}"
echo ""

# Method selection: use STAR counts (faster) or featureCounts (more precise)
USE_STAR_COUNTS=true

if [[ "$USE_STAR_COUNTS" == "true" ]]; then
    echo "=== Using STAR gene counts (Method A) ==="

    # Check if STAR count files exist
    echo "Checking STAR count files..."
    for SAMPLE in "${SAMPLES[@]}"; do
        COUNT_FILE="${ALIGNED_DIR}/${SAMPLE}_counts.tab"
        if [[ ! -f "${COUNT_FILE}" ]]; then
            echo "ERROR: STAR count file not found: ${COUNT_FILE}"
            exit 1
        fi
        echo "  Found: ${COUNT_FILE}"
    done

    echo "=== Combining STAR gene counts ==="

    # Write strandedness column to file for R script to read
    echo "${STRAND_COL}" > ${OUTPUT_DIR}/strand_column.txt

    # Create R script for combining counts
    cat > ${OUTPUT_DIR}/combine_star_counts.R << 'EOF'
#!/usr/bin/env Rscript

# Load required libraries
if (!require("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "https://cran.r-project.org/")
    library(tidyverse)
}

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/04_quantified")

# Read sample names
samples <- readLines("../../config/samples.txt")
aligned_dir <- "../03_aligned"

# Initialize count matrix
count_matrix <- NULL
gene_info <- NULL

# Read the strandedness column determined by the shell script
# 2 = unstranded, 3 = sense (forward-stranded), 4 = antisense (reverse-stranded)
strand_col <- as.integer(readLines("strand_column.txt"))
col_names <- c("gene_id", "unstranded", "sense", "antisense")
strand_label <- col_names[strand_col]
cat(sprintf("Using STAR count column %d (%s) based on strandedness check\n", strand_col, strand_label))

cat("Processing STAR count files...\n")

for (i in seq_along(samples)) {
    sample <- samples[i]
    count_file <- file.path(aligned_dir, paste0(sample, "_counts.tab"))

    cat(sprintf("Processing %s (%d/%d)...\n", sample, i, length(samples)))

    # Read STAR count file
    # Columns: gene_id, unstranded, sense(++,--), antisense(+-,-+)
    counts <- read.table(count_file, header = FALSE, stringsAsFactors = FALSE,
                        col.names = col_names)

    # Remove summary lines at the end
    counts <- counts[!grepl("^N_", counts$gene_id), ]

    # For first sample, initialize matrix with gene info
    if (is.null(count_matrix)) {
        gene_info <- counts[, "gene_id", drop = FALSE]
        count_matrix <- data.frame(gene_id = counts$gene_id)
    }

    # Add sample counts using the strand-appropriate column
    count_matrix[[sample]] <- counts[[strand_label]]
}

# Write count matrix
cat("Writing count matrix...\n")
write.table(count_matrix, "count_matrix.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("Count matrix saved to count_matrix.txt\n")
cat(sprintf("Dimensions: %d genes x %d samples\n", nrow(count_matrix), ncol(count_matrix)-1))

# Basic statistics
total_counts <- colSums(count_matrix[, -1])
cat("\nTotal counts per sample:\n")
for (i in 1:length(total_counts)) {
    cat(sprintf("  %s: %s\n", names(total_counts)[i], format(total_counts[i], big.mark = ",")))
}

# Save summary statistics
summary_stats <- data.frame(
    sample = names(total_counts),
    total_counts = total_counts,
    genes_detected = colSums(count_matrix[, -1] > 0)
)

write.table(summary_stats, "quantification_summary.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nQuantification complete!\n")
EOF

    # Run R script
    echo "Running R script to combine counts..."
    Rscript ${OUTPUT_DIR}/combine_star_counts.R

else
    echo "=== Using featureCounts (Method B) ==="

    # Check if BAM files exist
    echo "Checking BAM files..."
    BAM_FILES=()
    for SAMPLE in "${SAMPLES[@]}"; do
        BAM_FILE="${ALIGNED_DIR}/${SAMPLE}_sorted.bam"
        if [[ ! -f "${BAM_FILE}" ]]; then
            echo "ERROR: BAM file not found: ${BAM_FILE}"
            exit 1
        fi
        BAM_FILES+=("${BAM_FILE}")
        echo "  Found: ${BAM_FILE}"
    done

    echo "=== Running featureCounts ==="

    # Run featureCounts
    featureCounts \
        -T 16 \
        -p \
        -t exon \
        -g gene_id \
        -a ${GTF_FILE} \
        -o ${OUTPUT_DIR}/featureCounts_output.txt \
        "${BAM_FILES[@]}"

    echo "=== Processing featureCounts output ==="

    # Create R script for processing featureCounts output
    cat > ${OUTPUT_DIR}/process_featurecounts.R << 'EOF'
#!/usr/bin/env Rscript

# Load required libraries
if (!require("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "https://cran.r-project.org/")
    library(tidyverse)
}

# Set working directory
setwd("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/04_quantified")

# Read featureCounts output
cat("Reading featureCounts output...\n")
fc_data <- read.table("featureCounts_output.txt", header = TRUE,
                     sep = "\t", stringsAsFactors = FALSE, comment.char = "#")

# Extract count matrix (remove first 6 annotation columns)
count_matrix <- fc_data[, -c(1:6)]

# Add gene_id as first column
count_matrix <- cbind(gene_id = fc_data$Geneid, count_matrix)

# Clean sample names (remove path and .bam extension)
colnames(count_matrix)[-1] <- gsub(".*\\/([^/]+)_sorted\\.bam", "\\1", colnames(count_matrix)[-1])

# Write count matrix
cat("Writing count matrix...\n")
write.table(count_matrix, "count_matrix.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("Count matrix saved to count_matrix.txt\n")
cat(sprintf("Dimensions: %d genes x %d samples\n", nrow(count_matrix), ncol(count_matrix)-1))

# Basic statistics
total_counts <- colSums(count_matrix[, -1])
cat("\nTotal counts per sample:\n")
for (i in 1:length(total_counts)) {
    cat(sprintf("  %s: %s\n", names(total_counts)[i], format(total_counts[i], big.mark = ",")))
}

# Save summary statistics
summary_stats <- data.frame(
    sample = names(total_counts),
    total_counts = total_counts,
    genes_detected = colSums(count_matrix[, -1] > 0)
)

write.table(summary_stats, "quantification_summary.txt",
           sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nQuantification complete!\n")
EOF

    # Run R script
    echo "Running R script to process featureCounts output..."
    Rscript ${OUTPUT_DIR}/process_featurecounts.R

fi

echo "=== Creating sample metadata ==="

# Copy and format design file
cp ${BASE_DIR}/config/design.txt ${OUTPUT_DIR}/sample_metadata.txt

echo "=== Checking output files ==="

# Check if count matrix was created
if [[ -f "${OUTPUT_DIR}/count_matrix.txt" ]]; then
    echo "SUCCESS: Count matrix created"
    MATRIX_SIZE=$(ls -lh "${OUTPUT_DIR}/count_matrix.txt" | awk '{print $5}')
    echo "  File: ${OUTPUT_DIR}/count_matrix.txt (${MATRIX_SIZE})"

    # Count dimensions
    GENE_COUNT=$(tail -n +2 "${OUTPUT_DIR}/count_matrix.txt" | wc -l)
    SAMPLE_COUNT=$(($(head -n 1 "${OUTPUT_DIR}/count_matrix.txt" | tr '\t' '\n' | wc -l) - 1))
    echo "  Dimensions: ${GENE_COUNT} genes × ${SAMPLE_COUNT} samples"
else
    echo "ERROR: Count matrix not created!"
    exit 1
fi

# Check summary statistics
if [[ -f "${OUTPUT_DIR}/quantification_summary.txt" ]]; then
    echo "SUCCESS: Summary statistics created"
    echo "  File: ${OUTPUT_DIR}/quantification_summary.txt"
    echo ""
    echo "Summary statistics:"
    cat "${OUTPUT_DIR}/quantification_summary.txt"
else
    echo "WARNING: Summary statistics not created!"
fi

echo ""
echo "=== Gene expression quantification complete ==="
echo "Timestamp: $(date)"
echo ""
echo "Next step: Run differential expression analysis (5_deseq2.sh)"