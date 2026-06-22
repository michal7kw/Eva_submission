#!/bin/bash

#===============================================================================
# SCRIPT: 3b_spikein_calibration.sh
# PURPOSE: E. coli spike-in calibration for CUT&Tag normalization
#
# DESCRIPTION:
# CUT&Tag uses pA-Tn5 produced in E. coli. Residual E. coli DNA serves as a
# natural spike-in for calibration. This script aligns reads to E. coli K12
# genome and computes calibration scale factors for each sample.
#
# The scale factor is: scale_i = C / ecoli_reads_i
# where C = min(ecoli_reads) across all samples, ensuring the sample with
# fewest E. coli reads gets scale factor 1.0 (no downsampling).
#
# USAGE:
# sbatch 3b_spikein_calibration.sh
#
# INPUTS:
# - Trimmed FASTQ files from step 2
# - E. coli K12 Bowtie2 index
#
# OUTPUTS:
# - results/spikein/spikein_counts.txt (E. coli read counts per sample)
# - results/spikein/spikein_scale_factors.txt (calibration factors)
# - results/06_bigwig/*_spikein.bw (spike-in calibrated BigWig files)
#===============================================================================

#SBATCH --job-name=3b_spikein
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/3b_spikein.err"
#SBATCH --output="./logs/3b_spikein.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/alignment

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
TRIMMED_DIR="${BASE_DIR}/results/02_trimmed"
OUTPUT_DIR="${BASE_DIR}/results/spikein"
BIGWIG_DIR="${BASE_DIR}/results/06_bigwig"
BAM_DIR="${BASE_DIR}/results/04_filtered"

# E. coli K12 genome index — build if not present
ECOLI_INDEX="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/ecoli_K12"

mkdir -p ${OUTPUT_DIR}

echo "=== E. coli Spike-in Calibration ==="
echo "Timestamp: $(date)"

#===============================================================================
# Step 1: Build E. coli index if needed
#===============================================================================
if [[ ! -f "${ECOLI_INDEX}.1.bt2" ]]; then
    echo "E. coli K12 Bowtie2 index not found at ${ECOLI_INDEX}"
    echo "Please download E. coli K12 genome (GCF_000005845.2) and build index:"
    echo "  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    echo "  gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz"
    echo "  bowtie2-build GCF_000005845.2_ASM584v2_genomic.fna ${ECOLI_INDEX}"
    echo ""
    echo "Skipping spike-in calibration (index not available)."
    echo "This is a known limitation — document in the Discussion section."
    exit 0
fi

#===============================================================================
# Step 2: Align each sample to E. coli and count reads
#===============================================================================
SAMPLES=($(cat ${BASE_DIR}/config/samples.txt))
COUNT_FILE="${OUTPUT_DIR}/spikein_counts.txt"

echo "sample	ecoli_reads	human_reads" > ${COUNT_FILE}

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE}..."

    R1="${TRIMMED_DIR}/${SAMPLE}_R1_001_val_1.fq.gz"
    R2="${TRIMMED_DIR}/${SAMPLE}_R2_001_val_2.fq.gz"

    if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
        echo "  WARNING: Trimmed files not found for ${SAMPLE}, skipping"
        continue
    fi

    # Align to E. coli (use --very-sensitive-local, same as human alignment)
    ECOLI_READS=$(bowtie2 \
        -p 16 \
        --very-sensitive-local \
        --dovetail \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -I 10 \
        -X 700 \
        -x ${ECOLI_INDEX} \
        -1 ${R1} \
        -2 ${R2} \
        2> ${OUTPUT_DIR}/${SAMPLE}_ecoli_bowtie2.log | \
        samtools view -@ 16 -bS -f 2 - | \
        samtools view -c -)

    # Get human-aligned read count from filtered BAM
    HUMAN_READS=0
    if [[ -f "${BAM_DIR}/${SAMPLE}_filtered.bam" ]]; then
        HUMAN_READS=$(samtools view -c ${BAM_DIR}/${SAMPLE}_filtered.bam)
    fi

    echo "${SAMPLE}	${ECOLI_READS}	${HUMAN_READS}" >> ${COUNT_FILE}
    echo "  E. coli reads: ${ECOLI_READS}, Human reads: ${HUMAN_READS}"
done

echo ""
echo "=== Spike-in counts ==="
cat ${COUNT_FILE}

#===============================================================================
# Step 3: Compute scale factors
#===============================================================================
echo ""
echo "=== Computing scale factors ==="

SCALE_FILE="${OUTPUT_DIR}/spikein_scale_factors.txt"

# Use awk to compute scale factors: scale_i = min(ecoli) / ecoli_i
awk 'BEGIN {FS="\t"; OFS="\t"; min_ecoli=999999999}
     NR==1 {print $0, "scale_factor"; next}
     {
         if ($2+0 > 0 && $2+0 < min_ecoli) min_ecoli = $2+0;
         samples[NR] = $0;
         ecoli[NR] = $2+0;
     }
     END {
         for (i in samples) {
             if (ecoli[i] > 0) {
                 sf = min_ecoli / ecoli[i];
             } else {
                 sf = 1.0;
             }
             printf "%s\t%.6f\n", samples[i], sf;
         }
     }' ${COUNT_FILE} > ${SCALE_FILE}

echo "Scale factors:"
cat ${SCALE_FILE}

#===============================================================================
# Step 4: Generate spike-in calibrated BigWig files
#===============================================================================
echo ""
echo "=== Generating spike-in calibrated BigWig files ==="

# Read scale factors into associative array
declare -A SCALE_FACTORS
while IFS=$'\t' read -r sample ecoli human sf; do
    [[ "$sample" == "sample" ]] && continue
    SCALE_FACTORS[$sample]=$sf
done < ${SCALE_FILE}

# Activate bigwig conda env for bamCoverage
conda activate bigwig

for SAMPLE in "${SAMPLES[@]}"; do
    SF=${SCALE_FACTORS[$SAMPLE]:-1.0}
    BAM_FILE="${BAM_DIR}/${SAMPLE}_filtered.bam"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "  Skipping ${SAMPLE}: BAM not found"
        continue
    fi

    echo "  ${SAMPLE}: scale factor = ${SF}"

    bamCoverage \
        -b ${BAM_FILE} \
        -o ${BIGWIG_DIR}/${SAMPLE}_spikein.bw \
        --scaleFactor ${SF} \
        --binSize 10 \
        --numberOfProcessors 16 \
        --extendReads
done

echo ""
echo "=== Spike-in calibration complete ==="
echo "Counts: ${COUNT_FILE}"
echo "Scale factors: ${SCALE_FILE}"
echo "BigWig files: ${BIGWIG_DIR}/*_spikein.bw"
echo "Timestamp: $(date)"
