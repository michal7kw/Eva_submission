#!/bin/bash

#===============================================================================
# SCRIPT: tead1_cutandtag_snapshot.sh
# PURPOSE: IGV-style genome-browser snapshot of the TEAD1 locus for the Cut&Tag
#          data: per-replicate binding coverage for TES (3), TEAD1 (3) and the
#          IgG / "Mock" background controls (IggMs = TES control, IggRb = TEAD1
#          control), optionally with the MACS-called peak boxes below. Mirrors
#          the RNA-seq TEAD1 snapshot (SRF_Eva_RNA/scripts/9_tead1_locus_snapshot.sh).
#
# RENDERS 4 figures = {CPM, RPGC normalization} x {with peaks, coverage only}.
#
# CHROM NAMING: the Cut&Tag bigWigs and narrowPeaks use Ensembl names ("11"),
#   while the gencode GTF uses "chr11". We strip the "chr" prefix from the TEAD1
#   gene-model subset and query the region as "11:..." so bigWigs, peaks and the
#   gene model all share one chrom convention (same GRCh38 coordinates).
#
# PREREQUISITE: results/06_bigwig/*.bw (CPM + RPGC) and the combined narrowPeaks
#   under results/11_combined_replicates_narrow/peaks/ (already present).
#
# DEPENDENCIES: conda env `pygenometracks` (pyGenomeTracks + pyBigWig).
#
# USAGE: sbatch tead1_cutandtag_snapshot.sh   # or: bash tead1_cutandtag_snapshot.sh
#===============================================================================

#SBATCH --job-name=tead1_cnt_snapshot
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/tead1_cnt_snapshot.err"
#SBATCH --output="./logs/tead1_cnt_snapshot.out"

set -euo pipefail

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate pygenometracks

BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
BW_DIR="${BASE_DIR}/results/06_bigwig"
PEAK_DIR="${BASE_DIR}/results/11_combined_replicates_narrow/peaks"
OUT_DIR="${BASE_DIR}/results/06_bigwig"
GTF="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

GENE="TEAD1"
FLANK=5000
GENE_GTF="${OUT_DIR}/${GENE}_gene.cnt.gtf"   # chr-stripped gene-model subset

# Colors: TES blue, TEAD1 green, IgG/Mock grey (peaks match their group color)
TES_COLOR="#1F78B4"
TEAD1_COLOR="#33A02C"
IGG_COLOR="#777777"

# Combined per-group MACS peaks (Ensembl chrom naming, already "11")
TES_PEAKS="${PEAK_DIR}/TES_combined_peaks.narrowPeak"
TEAD1_PEAKS="${PEAK_DIR}/TEAD1_combined_peaks.narrowPeak"

# Each (normalization x peak-mode) is rendered twice:
#   - "full"   : auto shared y-max, standard track height
#   - "_zoomed": y-max divided by ZOOM_DIVISOR (top peaks saturate, mid signal
#                fills the track) and taller coverage tracks -> taller figure
ZOOM_DIVISOR=4          # zoomed y-max = full y-max / this
FULL_TRACK_H=1.2        # coverage track height, full view
ZOOM_TRACK_H=2.6        # coverage track height, zoomed view (taller = more visible)

mkdir -p "${OUT_DIR}" logs

# --- Sanity checks -----------------------------------------------------------
[[ -f "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}"; exit 1; }
for s in TES-1 TES-2 TES-3 TEAD1-1 TEAD1-2 TEAD1-3 IggMs IggRb; do
    for n in CPM RPGC; do
        [[ -f "${BW_DIR}/${s}_${n}.bw" ]] || { echo "ERROR: missing bigWig ${BW_DIR}/${s}_${n}.bw"; exit 1; }
    done
done
[[ -f "${TES_PEAKS}" ]]   || { echo "ERROR: missing ${TES_PEAKS}"; exit 1; }
[[ -f "${TEAD1_PEAKS}" ]] || { echo "ERROR: missing ${TEAD1_PEAKS}"; exit 1; }

# --- 1. TEAD1 coordinates (strip "chr" so chrom matches bigWigs/peaks) -------
echo "Extracting ${GENE} coordinates from GTF..."
read -r CHROM GSTART GEND <<< "$(awk -v g="\"${GENE}\"" '
    $3=="gene" && $0 ~ ("gene_name " g) { sub(/^chr/,"",$1); print $1, $4, $5; exit }
' "${GTF}")" || true
[[ -n "${CHROM:-}" ]] || { echo "ERROR: gene ${GENE} not found in ${GTF}"; exit 1; }

START=$(( GSTART - FLANK )); (( START < 1 )) && START=1
END=$(( GEND + FLANK ))
REGION="${CHROM}:${START}-${END}"
echo "  ${GENE} = ${CHROM}:${GSTART}-${GEND}  ->  view region ${REGION}"

# --- 2. Gene-model GTF subset, chr-stripped ----------------------------------
# NB: FS=OFS="\t" is required -- a GTF's 9th column is space-filled free text, so
# splitting on default whitespace would shred it when sub() rebuilds the record.
awk -v g="\"${GENE}\"" 'BEGIN{FS=OFS="\t"} $0 ~ ("gene_name " g) { sub(/^chr/,"",$1); print }' \
    "${GTF}" > "${GENE_GTF}"
echo "  Wrote gene-model subset: ${GENE_GTF} ($(wc -l < "${GENE_GTF}") lines)"

# --- helper: shared y-axis max over a list of bigWigs in the region ----------
compute_ymax () {
    python - "$REGION" "$@" <<'PY'
import sys, pyBigWig
region = sys.argv[1]; bws = sys.argv[2:]
chrom, rest = region.split(":"); start, end = (int(x) for x in rest.split("-"))
m = 0.0
for f in bws:
    bw = pyBigWig.open(f)
    vals = bw.stats(chrom, start, end, type="max", nBins=max(1, (end - start)//10))
    bw.close()
    for v in vals:
        if v is not None and v > m: m = v
print(round(m * 1.08, 4) if m > 0 else 1.0)
PY
}

# --- helper: append a bigWig coverage track to the ini -----------------------
# Uses the global ${TRACK_H} (set per view) for height.
bw_track () {  # $1 title  $2 file  $3 color  $4 ymax
    cat <<EOF
[${1}]
file = ${2}
title = ${1}
color = ${3}
min_value = 0
max_value = ${4}
height = ${TRACK_H}
number_of_bins = 700
nans_to_zeros = true
type = fill
show_data_range = true
file_type = bigwig

EOF
}

# --- 3. Build + render each (normalization x peak-mode) ----------------------
for NORM in CPM RPGC; do
    TES_BW=("${BW_DIR}/TES-1_${NORM}.bw" "${BW_DIR}/TES-2_${NORM}.bw" "${BW_DIR}/TES-3_${NORM}.bw")
    TEAD1_BW=("${BW_DIR}/TEAD1-1_${NORM}.bw" "${BW_DIR}/TEAD1-2_${NORM}.bw" "${BW_DIR}/TEAD1-3_${NORM}.bw")
    IGG_BW=("${BW_DIR}/IggMs_${NORM}.bw" "${BW_DIR}/IggRb_${NORM}.bw")

    echo "Computing shared y-axis max (${NORM})..."
    YMAX=$(compute_ymax "${TES_BW[@]}" "${TEAD1_BW[@]}" "${IGG_BW[@]}")
    echo "  ${NORM} shared max_value = ${YMAX}"

    for VIEW in full zoomed; do
        if [[ "${VIEW}" == "zoomed" ]]; then
            VM=$(awk -v y="${YMAX}" -v d="${ZOOM_DIVISOR}" 'BEGIN{printf "%.4f", y/d}')
            TRACK_H="${ZOOM_TRACK_H}"; vsuffix="_zoomed"; vtitle="  [zoomed: y/${ZOOM_DIVISOR}]"
        else
            VM="${YMAX}"; TRACK_H="${FULL_TRACK_H}"; vsuffix=""; vtitle=""
        fi
        echo "  ${VIEW} view: max_value=${VM}, track height=${TRACK_H}"

        for PEAKMODE in peaks nopeaks; do
            INI="${OUT_DIR}/tracks_${GENE}_cnt_${NORM}_${PEAKMODE}${vsuffix}.ini"
            {
                echo "[spacer]"; echo "height = 0.3"; echo ""
                i=1; for bw in "${TES_BW[@]}";   do bw_track "TES ${i}"   "${bw}" "${TES_COLOR}"   "${VM}"; i=$((i+1)); done
                i=1; for bw in "${TEAD1_BW[@]}"; do bw_track "TEAD1 ${i}" "${bw}" "${TEAD1_COLOR}" "${VM}"; i=$((i+1)); done
                bw_track "IgG (Ms)" "${IGG_BW[0]}" "${IGG_COLOR}" "${VM}"
                bw_track "IgG (Rb)" "${IGG_BW[1]}" "${IGG_COLOR}" "${VM}"

                if [[ "${PEAKMODE}" == "peaks" ]]; then
                    echo "[spacer]"; echo "height = 0.2"; echo ""
                    echo "[TES peaks]";   echo "file = ${TES_PEAKS}";   echo "title = TES peaks"
                    echo "color = ${TES_COLOR}";   echo "height = 0.5"; echo "display = collapsed"
                    echo "labels = false"; echo "border_color = none"; echo "file_type = bed"; echo ""
                    echo "[TEAD1 peaks]"; echo "file = ${TEAD1_PEAKS}"; echo "title = TEAD1 peaks"
                    echo "color = ${TEAD1_COLOR}"; echo "height = 0.5"; echo "display = collapsed"
                    echo "labels = false"; echo "border_color = none"; echo "file_type = bed"; echo ""
                fi

                echo "[spacer]"; echo "height = 0.3"; echo ""
                echo "[genes]"; echo "file = ${GENE_GTF}"; echo "title = ${GENE}"
                echo "prefered_name = gene_name"; echo "merge_transcripts = true"
                echo "color = darkblue"; echo "height = 1.0"; echo "fontsize = 10"; echo "file_type = gtf"; echo ""
                echo "[x-axis]"; echo "where = bottom"; echo "fontsize = 10"
            } > "${INI}"

            suffix=""; [[ "${PEAKMODE}" == "peaks" ]] && suffix="_peaks"
            TITLE="${GENE} locus  (Cut&Tag, ${NORM})${vtitle}"
            echo "Rendering ${NORM} ${VIEW} (${PEAKMODE})..."
            for ext in pdf png; do
                pyGenomeTracks --tracks "${INI}" --region "${REGION}" \
                    --dpi 300 --trackLabelFraction 0.12 --width 24 --title "${TITLE}" \
                    -o "${OUT_DIR}/${GENE}_cutandtag_${NORM}${suffix}${vsuffix}.${ext}"
            done
            echo "  Wrote ${OUT_DIR}/${GENE}_cutandtag_${NORM}${suffix}${vsuffix}.png / .pdf"
        done
    done
done

echo "=== Done ==="
echo "Region: ${REGION}"
echo "Outputs in ${OUT_DIR}: ${GENE}_cutandtag_{CPM,RPGC}{,_peaks}{,_zoomed}.png/.pdf"
