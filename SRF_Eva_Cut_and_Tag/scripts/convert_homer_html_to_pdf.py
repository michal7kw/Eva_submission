#!/usr/bin/env python3
"""
Convert HTML reports to PDF format.

This script converts HTML reports to PDF format for easier viewing and sharing.

Requirements:
    pip install weasyprint

Usage:
    python convert_homer_html_to_pdf.py                     # Convert all default files
    python convert_homer_html_to_pdf.py --input file.html   # Specific file

Examples:
    # Convert all default HTML reports
    python scripts/convert_homer_html_to_pdf.py

    # Convert a specific file
    python scripts/convert_homer_html_to_pdf.py --input results/multiqc/summary_report.html
"""

import argparse
import sys
from pathlib import Path

# Default base directory
DEFAULT_HOMER_DIR = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag/results/08_homer_motifs")

# HTML files to convert
HTML_FILES = ['knownResults.html', 'homerResults.html']

# Conditions (subdirectory prefixes)
CONDITIONS = ['TES', 'TEAD1']


def check_dependencies():
    """Check if required packages are installed."""
    try:
        from weasyprint import HTML, CSS
        return True
    except ImportError:
        print("ERROR: weasyprint is not installed.")
        print("\nTo install weasyprint, run:")
        print("  pip install weasyprint")
        print("\nOr create a conda environment:")
        print("  conda create -n weasyprint_env python=3.10")
        print("  conda activate weasyprint_env")
        print("  pip install weasyprint")
        return False


def get_custom_css():
    """Return custom CSS for better PDF rendering."""
    from weasyprint import CSS
    return CSS(string='''
        @page {
            size: A4 landscape;
            margin: 1cm;
        }
        body {
            font-family: Arial, sans-serif;
            font-size: 10pt;
        }
        table {
            font-size: 8pt;
            border-collapse: collapse;
            width: 100%;
        }
        td, th {
            padding: 2px 4px;
            border: 1px solid #ccc;
            vertical-align: top;
        }
        h1 {
            font-size: 14pt;
            margin-bottom: 10px;
        }
        svg {
            max-width: 100%;
            height: auto;
        }
        a {
            color: #0066cc;
            text-decoration: none;
        }
    ''')


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
    from weasyprint import HTML

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

        html = HTML(filename=str(html_path))
        html.write_pdf(str(pdf_path), stylesheets=[get_custom_css()])

        if verbose:
            print(f"OK -> {pdf_path.name}")
        return True

    except Exception as e:
        if verbose:
            print(f"FAILED: {e}")
        return False

def convert_condition(condition: str, base_dir: Path = DEFAULT_HOMER_DIR, verbose: bool = True) -> tuple:
    """
    Convert all HTML files for a specific condition.

    Args:
        condition: Condition name (e.g., 'TES', 'TEAD1')
        base_dir: Base HOMER results directory
        verbose: Print progress messages

    Returns:
        Tuple of (converted_count, failed_count)
    """
    condition_dir = base_dir / f"{condition}_peaks"

    if not condition_dir.exists():
        if verbose:
            print(f"\nSkipping {condition}: directory not found")
        return (0, 0)

    if verbose:
        print(f"\nProcessing {condition}...")

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
