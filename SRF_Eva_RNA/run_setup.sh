#!/bin/bash

#SBATCH --job-name=setup
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/setup.err"
#SBATCH --output="./logs/setup.out"


THREADS=32
LOCAL_SCRATCH_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/scratch"

# One-shot: download GENCODE v44 (GRCh38), verify checksums (SHA256 if present, else MD5),
# decompress, place GTF at your exact path, and build a STAR index.
# Multicore: uses THREADS, aria2c (if available) for downloads, pigz (if available) for decompress.
set -euo pipefail

# ========== USER PATHS (as requested) ==========
STAR_INDEX_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/genome/STAR_GRCh38"
GTF_PATH="/beegfs/scratch/ric.sessa/kubacki.michal/COMMONS/annotation/gencode.v44.annotation.gtf"

# Optional: build index on fast local scratch then copy to STAR_INDEX_DIR (empty to skip)
# Example: LOCAL_SCRATCH_DIR="/scratch/$USER/star_build_tmp"
LOCAL_SCRATCH_DIR="${LOCAL_SCRATCH_DIR:-}"

# ========== SETTINGS ==========
# For STAR sjdbOverhang: typical read length is 101 → overhang=100
READ_LENGTH="${READ_LENGTH:-101}"
SJDB_OVERHANG=$((READ_LENGTH - 1))

# Threads for STAR, pigz, aria2c; override with: THREADS=16 ./setup...
THREADS="${THREADS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)}"

# Use a downloads working area on the same filesystem for speed
WORKDIR="$(dirname "$STAR_INDEX_DIR")/downloads_grch38_gencode_v44"
mkdir -p "$WORKDIR" "$(dirname "$GTF_PATH")" "$STAR_INDEX_DIR"

# Data sources (GENCODE v44, GRCh38 primary assembly)
BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44"
FA_URL="$BASE/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="$BASE/gencode.v44.annotation.gtf.gz"

FA_GZ="$WORKDIR/GRCh38.primary_assembly.genome.fa.gz"
GTF_GZ="$WORKDIR/gencode.v44.annotation.gtf.gz"

# Redownload control: 0 = keep existing files; 1 = force re-download
FORCE_REDOWNLOAD="${FORCE_REDOWNLOAD:-0}"

# ========== HELPERS ==========
log(){ printf "[%s] %s\n" "$(date +'%F %T')" "$*"; }
fail(){ echo "ERROR: $*" >&2; exit 1; }

need_cmd(){ command -v "$1" >/dev/null 2>&1 || fail "Missing required command: $1"; }

# ========== CHECKS ==========
need_cmd awk
need_cmd grep
need_cmd sed
need_cmd sort
need_cmd tr
need_cmd cut
need_cmd mkdir
need_cmd cp

# Set up conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate /beegfs/scratch/ric.sessa/kubacki.michal/conda/envs/rnaseq-quant

# Pick downloader: aria2c (fast) → wget (fallback)
DOWNLOADER="wget"
if command -v aria2c >/dev/null 2>&1; then
  DOWNLOADER="aria2c"
fi

download_file () {
  local url="$1" out="$2"
  if [[ "$FORCE_REDOWNLOAD" -eq 1 ]]; then
    rm -f "$out"
  fi
  if [[ -s "$out" ]]; then
    log "Exists: $(basename "$out") — skipping download."
    return 0
  fi
  log "Downloading: $(basename "$out")"
  if [[ "$DOWNLOADER" == "aria2c" ]]; then
    aria2c -x 8 -s 8 -o "$out" "$url"
  else
    wget -O "$out" "$url"
  fi
}

# ========== DOWNLOAD FASTA & GTF ==========
download_file "$FA_URL"  "$FA_GZ"
download_file "$GTF_URL" "$GTF_GZ"

# ========== CHECKSUMS (prefer SHA256, fallback to MD5) ==========
SUMS_FILE=""
SUMS_TOOL=""

log "Fetching checksums (trying SHA256, then MD5)..."
if wget -q -O "$WORKDIR/SHA256SUMS" "$BASE/SHA256SUMS"; then
  SUMS_FILE="$WORKDIR/SHA256SUMS"
  SUMS_TOOL="sha256sum -c"
elif wget -q -O "$WORKDIR/SHA256SUMS.gz" "$BASE/SHA256SUMS.gz"; then
  gunzip -f "$WORKDIR/SHA256SUMS.gz"
  SUMS_FILE="$WORKDIR/SHA256SUMS"
  SUMS_TOOL="sha256sum -c"
elif wget -q -O "$WORKDIR/MD5SUMS" "$BASE/MD5SUMS"; then
  SUMS_FILE="$WORKDIR/MD5SUMS"
  SUMS_TOOL="md5sum -c"
else
  log "WARNING: No checksum file found. Skipping verification."
fi

if [[ -n "$SUMS_FILE" ]]; then
  log "Verifying checksums from $(basename "$SUMS_FILE")"
  awk '
    /GRCh38\.primary_assembly\.genome\.fa\.gz/ {print}
    /gencode\.v44\.annotation\.gtf\.gz/ {print}
  ' "$SUMS_FILE" > "$WORKDIR/SUMS_subset"
  ( cd "$WORKDIR" && $SUMS_TOOL SUMS_subset )
  log "Checksum verification: OK"
fi

# ========== DECOMPRESS (pigz if available) ==========
FA="$WORKDIR/GRCh38.primary_assembly.genome.fa"
GTF="$WORKDIR/gencode.v44.annotation.gtf"

log "Decompressing FASTA & GTF..."
if command -v pigz >/dev/null 2>&1; then
  pigz -df -p "$THREADS" "$FA_GZ"
  pigz -df -p "$THREADS" "$GTF_GZ"
else
  gunzip -fk "$FA_GZ"
  gunzip -fk "$GTF_GZ"
fi
[[ -s "$FA" ]]  || fail "Missing decompressed FASTA: $FA"
[[ -s "$GTF" ]] || fail "Missing decompressed GTF: $GTF"

# ========== PLACE GTF EXACTLY WHERE YOU WANT IT ==========
log "Copying GTF to: $GTF_PATH"
cp -f "$GTF" "$GTF_PATH"

# ========== BUILD STAR INDEX ==========
BUILD_DIR="$STAR_INDEX_DIR"
if [[ -n "${LOCAL_SCRATCH_DIR}" ]]; then
  mkdir -p "$LOCAL_SCRATCH_DIR"
  BUILD_DIR="$LOCAL_SCRATCH_DIR/STAR_GRCh38_build"
  mkdir -p "$BUILD_DIR"
  log "Using local scratch for indexing: $BUILD_DIR"
fi

log "Running STAR genomeGenerate with THREADS=$THREADS, sjdbOverhang=$SJDB_OVERHANG"
STAR \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$BUILD_DIR" \
  --genomeFastaFiles "$FA" \
  --sjdbGTFfile "$GTF_PATH" \
  --sjdbOverhang "$SJDB_OVERHANG"

# Optionally copy from local scratch to final destination
if [[ "$BUILD_DIR" != "$STAR_INDEX_DIR" ]]; then
  log "Copying STAR index from scratch to final dir: $STAR_INDEX_DIR"
  rsync -a --info=progress2 "$BUILD_DIR"/ "$STAR_INDEX_DIR"/
fi

# ========== QUICK SANITY CHECK ==========
if [[ ! -s "$STAR_INDEX_DIR/Genome" && ! -s "$STAR_INDEX_DIR/SA" ]]; then
  fail "STAR index appears incomplete in $STAR_INDEX_DIR"
fi

log "All done."
log "STAR genome index: $STAR_INDEX_DIR"
log "GTF annotation:    $GTF_PATH"
log "You can verify with: ./verify_paths.sh"

