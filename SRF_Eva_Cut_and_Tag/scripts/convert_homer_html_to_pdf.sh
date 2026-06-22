#!/bin/bash

#===============================================================================
# SCRIPT: convert_homer_html_to_pdf.sh
# PURPOSE: Convert HTML reports to PDF format
#
# DESCRIPTION:
# Converts HTML reports to PDF format for easier viewing and sharing.
#
# REQUIREMENTS:
# - Python with weasyprint package (installed automatically if needed)
#===============================================================================

#SBATCH --job-name=homer_to_pdf
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="./logs/homer_to_pdf.err"
#SBATCH --output="./logs/homer_to_pdf.out"

set -e

# Base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag"
HOMER_DIR="${BASE_DIR}/results/08_homer_motifs"

# Activate conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if weasyprint environment exists, create if not
ENV_NAME="weasyprint_env"
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment with weasyprint..."
    conda create -n ${ENV_NAME} python=3.10 -y
    conda activate ${ENV_NAME}
    pip install weasyprint
else
    conda activate ${ENV_NAME}
fi

echo "=========================================="
echo "HTML to PDF Converter"
echo "Date: $(date)"
echo "=========================================="

# Conditions to process (can be overridden by command line argument)
if [ $# -gt 0 ]; then
    # Pass all command line arguments as conditions
    CONDITIONS_ARGS="--conditions $@"
else
    # Let Python script use its defaults (TES, TEAD1, TESmut)
    CONDITIONS_ARGS=""
fi

# Run the Python converter
SCRIPT_PATH="${BASE_DIR}/scripts/convert_homer_html_to_pdf.py"

if [ ! -f "$SCRIPT_PATH" ]; then
    echo "ERROR: Python script not found at $SCRIPT_PATH"
    exit 1
fi

echo "Running Python conversion script..."
echo "Command: python3 $SCRIPT_PATH --base-dir $HOMER_DIR $CONDITIONS_ARGS"

python3 "$SCRIPT_PATH" \
    --base-dir "$HOMER_DIR" \
    $CONDITIONS_ARGS

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "PDF conversion complete!"
else
    echo ""
    echo "PDF conversion failed with exit code ${EXIT_CODE}"
fi

echo "Date: $(date)"
exit $EXIT_CODE
