#!/bin/bash

#===============================================================================
# SCRIPT: 3_align.sh
# PURPOSE: Read Alignment for RNA-seq Sequencing Data
#
# DESCRIPTION:
# This script performs alignment of trimmed paired-end FASTQ files to the
# reference genome using STAR aligner. It is the third step in the RNA-seq
# analysis pipeline and generates sorted, indexed BAM files with comprehensive
# alignment statistics optimized for RNA-seq data.
#
# KEY OPERATIONS:
# 1. Activates the conda environment with STAR and SAMtools
# 2. Validates input trimmed FASTQ files exist
# 3. Performs STAR alignment with RNA-seq optimized parameters
# 4. Converts SAM output to sorted BAM format
# 5. Indexes BAM files for downstream analysis
# 6. Generates alignment statistics and read counts
#
# METHODOLOGY:
# - Uses STAR aligner optimized for RNA-seq data
# - 2-pass alignment mode for improved splice junction detection
# - Outputs both genomic and transcriptomic alignments
# - Generates gene counts table during alignment
# - Multi-threading for optimal performance
#
# IMPORTANT PARAMETERS:
# - Alignment mode: 2-pass for splice junction discovery
# - Output format: BAM, sorted by coordinates
# - Gene counting: during alignment using GeneCounts
# - Splice junction filtering: canonical and semi-canonical
# - Multi-mapping: up to 20 locations
#
# EXPECTED INPUTS:
# - Trimmed FASTQ files: ${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz
# - Trimmed FASTQ files: ${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz
# - STAR genome index: ${GENOME_INDEX}
# - GTF annotation: ${GTF_FILE}
# - Sample list: ${BASE_DIR}/config/samples.txt
#
# EXPECTED OUTPUTS:
# - Sorted BAM files: ${OUTPUT_DIR}/${SAMPLE}Aligned.sortedByCoord.out.bam
# - Gene counts: ${OUTPUT_DIR}/${SAMPLE}ReadsPerGene.out.tab
# - Alignment logs: ${OUTPUT_DIR}/${SAMPLE}Log.final.out
# - Junction files: ${OUTPUT_DIR}/${SAMPLE}SJ.out.tab
# - SLURM logs: ./logs/3_align_${ARRAY_ID}.out and ./logs/3_align_${ARRAY_ID}.err
#
# DEPENDENCIES:
# - STAR aligner (via conda environment)
# - SAMtools (for BAM processing and indexing)
# - Reference genome STAR index
# - GTF annotation file
# - SLURM workload manager
#
# USAGE:
# sbatch 3_align.sh
#===============================================================================

#SBATCH --job-name=3_align
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --array=0-5
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3_align_%a.err"
#SBATCH --output="./logs/3_align_%a.out"

# Set up conda environment with required tools
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
TRIMMED_DIR="${BASE_DIR}/results/02_trimmed"
OUTPUT_DIR="${BASE_DIR}/results/03_aligned"

# Reference genome and annotation - adjust paths as needed
GENOME_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCh38"
GTF_FILE="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=== Starting alignment for sample: ${SAMPLE} ==="
echo "Timestamp: $(date)"
echo "Input files:"
echo "  R1: ${TRIMMED_DIR}/${SAMPLE}_val_1.fq.gz"
echo "  R2: ${TRIMMED_DIR}/${SAMPLE}_val_2.fq.gz"
echo "Output directory: ${OUTPUT_DIR}"
echo "Genome index: ${GENOME_INDEX}"
echo "GTF file: ${GTF_FILE}"
echo ""

# Check if input files exist
if [[ ! -f "${TRIMMED_DIR}/${SAMPLE}_val_1.fq.gz" ]]; then
    echo "ERROR: R1 file not found: ${TRIMMED_DIR}/${SAMPLE}_val_1.fq.gz"
    exit 1
fi
if [[ ! -f "${TRIMMED_DIR}/${SAMPLE}_val_2.fq.gz" ]]; then
    echo "ERROR: R2 file not found: ${TRIMMED_DIR}/${SAMPLE}_val_2.fq.gz"
    exit 1
fi

# Check if genome index exists
if [[ ! -d "${GENOME_INDEX}" ]]; then
    echo "ERROR: STAR genome index not found: ${GENOME_INDEX}"
    echo "Please create STAR index first or adjust the path"
    exit 1
fi

# Create sample-specific output directory
SAMPLE_OUTPUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p ${SAMPLE_OUTPUT}

echo "=== Step 1: Running STAR alignment ==="
echo "Timestamp: $(date)"

# STAR alignment with RNA-seq optimized parameters
STAR \
    --runThreadN 32 \
    --genomeDir ${GENOME_INDEX} \
    --sjdbGTFfile ${GTF_FILE} \
    --readFilesIn ${TRIMMED_DIR}/${SAMPLE}_val_1.fq.gz ${TRIMMED_DIR}/${SAMPLE}_val_2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix ${SAMPLE_OUTPUT}/${SAMPLE} \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --quantMode GeneCounts \
    --twopassMode Basic

echo "=== Step 2: Indexing BAM file ==="
echo "Timestamp: $(date)"

# Index the BAM file
samtools index -@ 32 ${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam

echo "=== Step 3: Generating alignment statistics ==="
echo "Timestamp: $(date)"

# Generate alignment statistics
samtools flagstat ${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT}/${SAMPLE}_flagstat.txt

# Generate additional statistics
samtools stats ${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam > ${SAMPLE_OUTPUT}/${SAMPLE}_stats.txt

echo "=== Step 4: Checking output files ==="
BAM_FILE="${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam"
COUNTS_FILE="${SAMPLE_OUTPUT}/${SAMPLE}ReadsPerGene.out.tab"
LOG_FILE="${SAMPLE_OUTPUT}/${SAMPLE}Log.final.out"

if [[ -f "${BAM_FILE}" ]]; then
    BAM_SIZE=$(ls -lh "${BAM_FILE}" | awk '{print $5}')
    echo "SUCCESS: BAM file created - ${BAM_FILE} (${BAM_SIZE})"

    # Quick read count check
    READ_COUNT=$(samtools view -c "${BAM_FILE}")
    echo "Total reads in BAM: ${READ_COUNT}"
else
    echo "ERROR: BAM file not created!"
    exit 1
fi

if [[ -f "${COUNTS_FILE}" ]]; then
    GENE_COUNT=$(tail -n +5 "${COUNTS_FILE}" | wc -l)
    echo "SUCCESS: Gene counts file created - ${COUNTS_FILE}"
    echo "Genes with counts: ${GENE_COUNT}"
else
    echo "WARNING: Gene counts file not created!"
fi

if [[ -f "${LOG_FILE}" ]]; then
    echo "SUCCESS: STAR log file created - ${LOG_FILE}"
    echo "Key alignment metrics:"
    grep -E "(Uniquely mapped reads|Number of reads mapped to multiple loci|% of reads mapped to multiple loci)" "${LOG_FILE}"
else
    echo "WARNING: STAR log file not found!"
fi

# Copy files to main output directory for easier access
cp ${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam ${OUTPUT_DIR}/${SAMPLE}_sorted.bam
cp ${SAMPLE_OUTPUT}/${SAMPLE}Aligned.sortedByCoord.out.bam.bai ${OUTPUT_DIR}/${SAMPLE}_sorted.bam.bai
cp ${SAMPLE_OUTPUT}/${SAMPLE}ReadsPerGene.out.tab ${OUTPUT_DIR}/${SAMPLE}_counts.tab

echo "=== Alignment complete for ${SAMPLE} ==="
echo "Timestamp: $(date)"