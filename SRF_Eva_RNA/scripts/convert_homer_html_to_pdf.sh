#!/bin/bash

#===============================================================================
# SCRIPT: convert_homer_html_to_pdf.sh
# PURPOSE: Convert HOMER HTML reports to PDF format
#
# DESCRIPTION:
# Converts HOMER HTML reports (knownResults.html, homerResults.html) to PDF format
# for easier viewing and sharing.
#
# REQUIREMENTS:
# - Python with pdfkit package
# - wkhtmltopdf binary (installed via conda)
#
# USAGE:
# sbatch scripts/convert_homer_html_to_pdf.sh
# sbatch scripts/convert_homer_html_to_pdf.sh upregulated_promoters downregulated_promoters
#===============================================================================

#SBATCH --job-name=homer_to_pdf
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/homer_to_pdf.err"
#SBATCH --output="logs/homer_to_pdf.out"

set -e

# Base directories
BASE_DIR="/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA"
HOMER_DIR="${BASE_DIR}/results/07_homer_motifs"

cd "${BASE_DIR}"
mkdir -p logs

# Activate conda
source /opt/common/tools/ric.cosr/miniconda3/bin/activate

# Check if pdfkit environment exists, create if not
ENV_NAME="pdfkit_env"
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "Creating conda environment with pdfkit and wkhtmltopdf..."
    conda create -n ${ENV_NAME} python=3.10 -y
    conda activate ${ENV_NAME}
    # Install wkhtmltopdf from conda-forge (includes the binary)
    conda install -c conda-forge wkhtmltopdf -y
    # Install pdfkit python package
    pip install pdfkit
else
    conda activate ${ENV_NAME}
fi

# Verify wkhtmltopdf is available
if ! command -v wkhtmltopdf &> /dev/null; then
    echo "ERROR: wkhtmltopdf not found. Installing..."
    conda install -c conda-forge wkhtmltopdf -y
fi

echo "=========================================="
echo "HOMER HTML to PDF Converter"
echo "Date: $(date)"
echo "=========================================="
echo ""
echo "Python: $(which python3)"
echo "wkhtmltopdf: $(which wkhtmltopdf)"
echo "HOMER directory: ${HOMER_DIR}"
echo ""

# Conditions to process (can be overridden by command line argument)
if [ $# -gt 0 ]; then
    # Pass all command line arguments as conditions
    CONDITIONS_ARGS="--conditions $@"
else
    # Let Python script use its defaults (all_degs_promoters, downregulated_promoters, upregulated_promoters)
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
echo ""

python3 "$SCRIPT_PATH" \
    --base-dir "$HOMER_DIR" \
    $CONDITIONS_ARGS

EXIT_CODE=$?

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "PDF conversion complete!"
    echo ""
    echo "Output PDFs can be found in:"
    echo "  ${HOMER_DIR}/upregulated_promoters/*.pdf"
    echo "  ${HOMER_DIR}/downregulated_promoters/*.pdf"
    echo "  ${HOMER_DIR}/all_degs_promoters/*.pdf"
else
    echo ""
    echo "PDF conversion failed with exit code ${EXIT_CODE}"
fi

echo ""
echo "Completed: $(date)"
exit $EXIT_CODE
