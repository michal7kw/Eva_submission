#!/bin/bash

#===============================================================================
# SCRIPT: 10_tead1_exon_collapsed.sh
# PURPOSE: Exon-collapsed RNA-seq coverage snapshot of the TEAD1 locus,
#          Mock (red) vs TES (blue), 3 replicates each. Signal is drawn ONLY
#          over exons; introns are rendered as narrow fixed-width skipped gaps,
#          so the readable look matches the (intronless) SOX2 reference panel.
#
#          This is an ALTERNATIVE view to 9_tead1_locus_snapshot.sh (which keeps
#          true genomic coordinates via pyGenomeTracks). Both are kept.
#
# DESCRIPTION:
# 1. Regenerates a TEAD1-only GTF subset from the gencode.v44 GTF (self-contained;
#    same awk filter as script 9), so this step depends only on the step-8 bigWigs.
# 2. Calls 10_tead1_exon_collapsed.py, which takes the UNION of TEAD1 exons,
#    lays them side-by-side at true bp width separated by narrow gaps, pulls
#    per-base coverage from each bigWig (pyBigWig) and renders 6 stacked panels
#    + a gene-model strip + an in-exon scale bar to PDF + PNG.
#
# PREREQUISITE: run 8_rna_coverage_bigwig.sh first (produces results/08_bigwig/*.bw)
#
# DEPENDENCIES / SETUP (one-time, on the cluster):
#   conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks
#   (pyBigWig, numpy and matplotlib ship with pyGenomeTracks)
#
# USAGE:
# sbatch 10_tead1_exon_collapsed.sh        # or: bash 10_tead1_exon_collapsed.sh
#===============================================================================

#SBATCH --job-name=10_tead1_exon_collapsed
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/10_tead1_exon_collapsed.err"
#SBATCH --output="./logs/10_tead1_exon_collapsed.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate pygenometracks

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
BW_DIR="${BASE_DIR}/results/08_bigwig"
OUT_DIR="${BASE_DIR}/results/08_bigwig"
SCRIPT_DIR="${BASE_DIR}/scripts"
GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

GENE="TEAD1"
GENE_GTF="${OUT_DIR}/${GENE}_gene.gtf"
OUT_PREFIX="${OUT_DIR}/${GENE}_exon_collapsed_Mock_vs_TES"

# Mock = GFP control; sample IDs from config/samples.txt
MOCK_BW=("${BW_DIR}/ASE01_GFP1.bw" "${BW_DIR}/ASE01_GFP2.bw" "${BW_DIR}/ASE01_GFP3.bw")
TES_BW=("${BW_DIR}/ASE01_TES1.bw" "${BW_DIR}/ASE01_TES2.bw" "${BW_DIR}/ASE01_TES3.bw")

mkdir -p "${OUT_DIR}"

# --- Sanity checks -----------------------------------------------------------
[[ -f "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}"; exit 1; }
for bw in "${MOCK_BW[@]}" "${TES_BW[@]}"; do
    [[ -f "${bw}" ]] || { echo "ERROR: BigWig not found: ${bw} (run 8_rna_coverage_bigwig.sh first)"; exit 1; }
done

# --- 1. TEAD1-only GTF subset (self-contained; same filter as script 9) ------
echo "Writing ${GENE} GTF subset..."
awk -v g="\"${GENE}\"" '$0 ~ ("gene_name " g)' "${GTF}" > "${GENE_GTF}"
echo "  ${GENE_GTF} ($(wc -l < "${GENE_GTF}") lines)"

# --- 2. Render exon-collapsed snapshot ---------------------------------------
echo "Rendering exon-collapsed snapshot..."
python "${SCRIPT_DIR}/10_tead1_exon_collapsed.py" \
    --gtf "${GENE_GTF}" \
    --gene "${GENE}" \
    --bw-mock "${MOCK_BW[@]}" \
    --bw-tes "${TES_BW[@]}" \
    --out-prefix "${OUT_PREFIX}" \
    --dpi 300

echo "=== Done ==="
echo "Outputs: ${OUT_PREFIX}.pdf / .png"
