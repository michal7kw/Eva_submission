#!/usr/bin/env python3
"""
Convert HOMER HTML reports to PDF format.

This script converts HOMER HTML reports to PDF format for easier viewing and sharing.
Uses pdfkit (wkhtmltopdf) which has fewer system dependencies than weasyprint.

Requirements:
    pip install pdfkit
    System: wkhtmltopdf (can be installed via conda: conda install -c conda-forge wkhtmltopdf)
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Default base directory
DEFAULT_HOMER_DIR = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_RNA/results/07_homer_motifs")

# HTML files to convert
HTML_FILES = ['knownResults.html', 'homerResults.html']

# Conditions (subdirectory names - these are the full directory names, not prefixes)
CONDITIONS = ['all_degs_promoters', 'downregulated_promoters', 'upregulated_promoters']


def check_dependencies():
    """Check if required packages are installed."""
    # Check for pdfkit
    try:
        import pdfkit
    except ImportError:
        print("ERROR: pdfkit is not installed.")
        print("\nTo install pdfkit, run:")
        print("  pip install pdfkit")
        return False

    # Check for wkhtmltopdf binary
    try:
        result = subprocess.run(['wkhtmltopdf', '--version'],
                                capture_output=True, text=True, timeout=10)
        if result.returncode != 0:
            raise FileNotFoundError()
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("ERROR: wkhtmltopdf is not installed or not in PATH.")
        print("\nTo install wkhtmltopdf:")
        print("  conda install -c conda-forge wkhtmltopdf")
        print("  # or on Debian/Ubuntu: sudo apt-get install wkhtmltopdf")
        return False


def get_pdfkit_options():
    """Return pdfkit options for better PDF rendering."""
    return {
        'page-size': 'A4',
        'orientation': 'Landscape',
        'margin-top': '10mm',
        'margin-right': '10mm',
        'margin-bottom': '10mm',
        'margin-left': '10mm',
        'encoding': 'UTF-8',
        'no-outline': None,
        'enable-local-file-access': None,  # Required for local images
        'quiet': None,
    }


def convert_html_to_pdf(html_path: Path, pdf_path: Path = None, verbose: bool = True) -> bool:
    """
    Convert a single HTML file to PDF.

    Args:
        html_path: Path to input HTML file
        pdf_path: Path to output PDF file (default: same location as HTML)
        verbose: Print progress messages

    Returns:
        True if conversion succeeded, False otherwise
    """
    import pdfkit

    html_path = Path(html_path)

    if not html_path.exists():
        if verbose:
            print(f"  ERROR: File not found: {html_path}")
        return False

    if pdf_path is None:
        pdf_path = html_path.with_suffix('.pdf')
    else:
        pdf_path = Path(pdf_path)

    try:
        if verbose:
            print(f"  Converting {html_path.name}...", end=' ', flush=True)

        pdfkit.from_file(str(html_path), str(pdf_path), options=get_pdfkit_options())

        if verbose:
            print(f"OK -> {pdf_path.name}")
        return True

    except Exception as e:
        if verbose:
            print(f"FAILED: {e}")
        return False


def convert_condition(condition: str, base_dir: Path = DEFAULT_HOMER_DIR, verbose: bool = True) -> tuple:
    """
    Convert all HTML files for a specific condition directory.

    Args:
        condition: Condition directory name (e.g., 'upregulated_promoters', 'downregulated_promoters')
        base_dir: Base HOMER results directory
        verbose: Print progress messages

    Returns:
        Tuple of (converted_count, failed_count)
    """
    # Use the condition name directly as the subdirectory name
    # (CONDITIONS already contains full directory names like 'upregulated_promoters')
    condition_dir = base_dir / condition

    if not condition_dir.exists():
        if verbose:
            print(f"\nSkipping {condition}: directory not found at {condition_dir}")
        return (0, 0)

    if verbose:
        print(f"\nProcessing {condition}...")
        print(f"  Directory: {condition_dir}")

    converted = 0
    failed = 0

    for html_file in HTML_FILES:
        html_path = condition_dir / html_file

        if not html_path.exists():
            if verbose:
                print(f"  Skipping {html_file}: file not found")
            continue

        if convert_html_to_pdf(html_path, verbose=verbose):
            converted += 1
        else:
            failed += 1

    return (converted, failed)


def main():
    parser = argparse.ArgumentParser(
        description="Convert HTML reports to PDF format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                           # Convert all default files
  %(prog)s --input file.html         # Convert specific file
  %(prog)s --base-dir /path/to/dir   # Use custom base directory
        """
    )

    parser.add_argument(
        '--input', '-i',
        type=Path,
        help='Convert a specific HTML file'
    )

    parser.add_argument(
        '--output', '-o',
        type=Path,
        help='Output PDF path (only with --input)'
    )

    parser.add_argument(
        '--base-dir', '-d',
        type=Path,
        default=DEFAULT_HOMER_DIR,
        help=f'Base project directory (default: {DEFAULT_HOMER_DIR})'
    )

    parser.add_argument(
        '--conditions', '-c',
        nargs='+',
        default=CONDITIONS,
        help=f'Conditions to process (default: {" ".join(CONDITIONS)})'
    )

    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Suppress progress messages'
    )

    args = parser.parse_args()
    verbose = not args.quiet

    # Check dependencies first
    if not check_dependencies():
        sys.exit(1)

    if verbose:
        print("=" * 50)
        print("HTML to PDF Converter")
        print("=" * 50)

    # Handle single file conversion
    if args.input:
        success = convert_html_to_pdf(args.input, args.output, verbose=verbose)
        sys.exit(0 if success else 1)

    # Convert all files in the list
    total_converted = 0
    total_failed = 0

    for condition in args.conditions:
        c, f = convert_condition(condition, args.base_dir, verbose=verbose)
        total_converted += c
        total_failed += f

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"Conversion complete!")
        print(f"  Converted: {total_converted}")
        print(f"  Failed: {total_failed}")
        print(f"{'=' * 50}")

    sys.exit(0 if total_failed == 0 else 1)


if __name__ == '__main__':
    main()
