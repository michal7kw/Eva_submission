#!/bin/bash

#===============================================================================
# SCRIPT: 9_tead1_locus_snapshot.sh
# PURPOSE: IGV-style genome-browser snapshot of the TEAD1 locus showing
#          per-replicate RNA-seq coverage for Mock (red) vs TES (blue),
#          with the gene model and a scale bar below (pyGenomeTracks).
#          Mirrors the SES-paper SOX2 track figure.
#
# DESCRIPTION:
# 1. Reads the exact TEAD1 coordinates from the gencode.v44 GTF (TEAD1 is on
#    chr11p15.3) and adds a flank for the view region.
# 2. Builds a TEAD1-only GTF subset for a clean gene model track.
# 3. Computes a single shared y-axis max across all 6 bigWigs over the region
#    (via pyBigWig, which ships with pyGenomeTracks) so Mock and TES are
#    directly comparable -- the reference figure uses one fixed scale.
# 4. Writes tracks.ini (3 Mock + 3 TES coverage tracks + gene track + x-axis)
#    and renders to PDF + PNG.
#
# PREREQUISITE: run 8_rna_coverage_bigwig.sh first (produces results/08_bigwig/*.bw)
#
# DEPENDENCIES / SETUP (one-time, on the cluster):
#   conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks
#   (pyBigWig is pulled in as a pyGenomeTracks dependency)
#
# USAGE:
# sbatch 9_tead1_locus_snapshot.sh        # or: bash 9_tead1_locus_snapshot.sh
#===============================================================================

#SBATCH --job-name=9_tead1_snapshot
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/9_tead1_snapshot.err"
#SBATCH --output="./logs/9_tead1_snapshot.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate pygenometracks

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
BW_DIR="${BASE_DIR}/results/08_bigwig"
OUT_DIR="${BASE_DIR}/results/08_bigwig"
GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

GENE="TEAD1"
FLANK=5000          # bp added on each side of the gene for the view region
INI="${OUT_DIR}/tracks_${GENE}.ini"
GENE_GTF="${OUT_DIR}/${GENE}_gene.gtf"

# Mock = GFP control; sample IDs from config/samples.txt
MOCK_BW=("${BW_DIR}/ASE01_GFP1.bw" "${BW_DIR}/ASE01_GFP2.bw" "${BW_DIR}/ASE01_GFP3.bw")
TES_BW=("${BW_DIR}/ASE01_TES1.bw" "${BW_DIR}/ASE01_TES2.bw" "${BW_DIR}/ASE01_TES3.bw")
MOCK_COLOR="#B22222"   # red  (Mock)
TES_COLOR="#1F78B4"    # blue (TES)

mkdir -p "${OUT_DIR}"

# --- Sanity checks -----------------------------------------------------------
[[ -f "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}"; exit 1; }
for bw in "${MOCK_BW[@]}" "${TES_BW[@]}"; do
    [[ -f "${bw}" ]] || { echo "ERROR: BigWig not found: ${bw} (run 8_rna_coverage_bigwig.sh first)"; exit 1; }
done

# --- 1. Extract TEAD1 gene coordinates from the GTF --------------------------
echo "Extracting ${GENE} coordinates from GTF..."
read -r CHROM GSTART GEND <<< "$(awk -v g="\"${GENE}\"" '
    $3=="gene" && $0 ~ ("gene_name " g) { print $1, $4, $5; exit }
' "${GTF}")" || true

if [[ -z "${CHROM:-}" ]]; then
    echo "ERROR: Could not find gene ${GENE} in ${GTF}"
    exit 1
fi

START=$(( GSTART - FLANK )); (( START < 1 )) && START=1
END=$(( GEND + FLANK ))
REGION="${CHROM}:${START}-${END}"
echo "  ${GENE} = ${CHROM}:${GSTART}-${GEND}  ->  view region ${REGION}"

# --- 2. TEAD1-only GTF subset for a clean gene-model track -------------------
awk -v g="\"${GENE}\"" '$0 ~ ("gene_name " g)' "${GTF}" > "${GENE_GTF}"
echo "  Wrote gene-model subset: ${GENE_GTF} ($(wc -l < "${GENE_GTF}") lines)"

# --- 3. Shared y-axis max across all 6 bigWigs over the region ---------------
echo "Computing shared y-axis max over the region..."
YMAX=$(python - "$REGION" "${MOCK_BW[@]}" "${TES_BW[@]}" <<'PY'
import sys, pyBigWig
region = sys.argv[1]
bws = sys.argv[2:]
chrom, rest = region.split(":")
start, end = (int(x) for x in rest.split("-"))
m = 0.0
for f in bws:
    bw = pyBigWig.open(f)
    # per-bin max over the region (binned coverage is what gets plotted)
    vals = bw.stats(chrom, start, end, type="max", nBins=max(1, (end - start) // 10))
    bw.close()
    for v in vals:
        if v is not None and v > m:
            m = v
# add ~8% headroom; guard against zero
print(round(m * 1.08, 4) if m > 0 else 1.0)
PY
)
echo "  Shared max_value = ${YMAX}"

# --- 4. Build tracks.ini -----------------------------------------------------
echo "Writing ${INI}..."
{
echo "[spacer]"
echo "height = 0.3"
echo ""

idx=1
for bw in "${MOCK_BW[@]}"; do
    echo "[Mock ${idx}]"
    echo "file = ${bw}"
    echo "title = Mock ${idx}"
    echo "color = ${MOCK_COLOR}"
    echo "min_value = 0"
    echo "max_value = ${YMAX}"
    echo "height = 1.4"
    echo "number_of_bins = 700"
    echo "nans_to_zeros = true"
    echo "type = fill"
    echo "show_data_range = true"
    echo "file_type = bigwig"
    echo ""
    idx=$((idx+1))
done

idx=1
for bw in "${TES_BW[@]}"; do
    echo "[TES ${idx}]"
    echo "file = ${bw}"
    echo "title = TES ${idx}"
    echo "color = ${TES_COLOR}"
    echo "min_value = 0"
    echo "max_value = ${YMAX}"
    echo "height = 1.4"
    echo "number_of_bins = 700"
    echo "nans_to_zeros = true"
    echo "type = fill"
    echo "show_data_range = true"
    echo "file_type = bigwig"
    echo ""
    idx=$((idx+1))
done

echo "[spacer]"
echo "height = 0.3"
echo ""
echo "[genes]"
echo "file = ${GENE_GTF}"
echo "title = ${GENE}"
echo "prefered_name = gene_name"
echo "merge_transcripts = true"
echo "color = darkblue"
echo "height = 1.2"
echo "fontsize = 10"
echo "file_type = gtf"
echo ""
echo "[x-axis]"
echo "where = bottom"
echo "fontsize = 10"
} > "${INI}"

# --- 5. Render PDF + PNG -----------------------------------------------------
echo "Rendering snapshot for ${REGION}..."
for ext in pdf png; do
    pyGenomeTracks \
        --tracks "${INI}" \
        --region "${REGION}" \
        --dpi 300 \
        --trackLabelFraction 0.12 \
        --width 24 \
        --title "${GENE} locus  (RNA-seq coverage, CPM)" \
        -o "${OUT_DIR}/${GENE}_locus_Mock_vs_TES.${ext}"
    echo "  Wrote ${OUT_DIR}/${GENE}_locus_Mock_vs_TES.${ext}"
done

echo "=== Done ==="
echo "Region: ${REGION}   shared max_value: ${YMAX}"
