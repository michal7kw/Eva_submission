#!/usr/bin/env python3
"""
Summarize HOMER motif analysis results for stakeholder reports.

Extracts the most meaningful information from HOMER outputs:
- Top enriched known motifs (by p-value)
- De novo discovered motifs with best matches
- Cross-condition comparison
- TEAD-family motif enrichment (project-specific)

Outputs:
- CSV tables for each condition
- Combined comparison table
- Executive summary report
"""

import os
import re
import sys
import pandas as pd
from pathlib import Path

def parse_known_motifs(filepath, top_n=25):
    """Parse knownResults.txt and extract top motifs."""
    df = pd.read_csv(filepath, sep='\t')
    df.columns = [
        'Motif_Name', 'Consensus', 'P_value', 'Log_P_value', 'q_value',
        'Target_Seqs_with_Motif', 'Pct_Target', 'Background_Seqs', 'Pct_Background'
    ]

    # Clean up percentage columns
    df['Pct_Target'] = df['Pct_Target'].str.rstrip('%').astype(float)
    df['Pct_Background'] = df['Pct_Background'].str.rstrip('%').astype(float)

    # Calculate fold enrichment
    df['Fold_Enrichment'] = df['Pct_Target'] / df['Pct_Background'].replace(0, 0.01)

    # Extract TF name and family from motif name
    df['TF_Name'] = df['Motif_Name'].apply(lambda x: x.split('(')[0])
    df['TF_Family'] = df['Motif_Name'].apply(lambda x: x.split('(')[1].split(')')[0] if '(' in x else 'Unknown')

    return df.head(top_n)


def parse_denovo_motifs(filepath):
    """Parse homerMotifs.all.motifs for de novo discovered motifs."""
    motifs = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    consensus = parts[0][1:]  # Remove '>'
                    name = parts[1]
                    log_p = float(parts[3])
                    p_value = f"1e{int(log_p)}" if log_p < 0 else str(10**log_p)

                    # Parse target/background stats
                    stats = parts[5] if len(parts) > 5 else ''
                    match = re.search(r'T:([\d.]+)\(([\d.]+)%\),B:([\d.]+)\(([\d.]+)%\)', stats)
                    if match:
                        target_n = float(match.group(1))
                        target_pct = float(match.group(2))
                        bg_n = float(match.group(3))
                        bg_pct = float(match.group(4))
                        fold_enrich = target_pct / bg_pct if bg_pct > 0 else 0
                    else:
                        target_pct, bg_pct, fold_enrich = 0, 0, 0

                    motifs.append({
                        'Rank': len(motifs) + 1,
                        'Consensus': consensus,
                        'P_value': p_value,
                        'Log_P_value': log_p,
                        'Pct_Target': target_pct,
                        'Pct_Background': bg_pct,
                        'Fold_Enrichment': round(fold_enrich, 2)
                    })

    return pd.DataFrame(motifs)


def find_tead_motifs(df):
    """Extract TEAD/TEA family motifs from results."""
    tead_mask = df['Motif_Name'].str.contains(r'TEAD|TEA\)', case=False, regex=True)
    return df[tead_mask].copy()


def generate_summary_report(output_dir, conditions=['TES', 'TEAD1']):
    """Generate comprehensive summary report."""

    results = {}
    all_known = []
    all_denovo = []
    tead_summary = []

    for condition in conditions:
        cond_dir = output_dir / f"{condition}_peaks"
        if not cond_dir.exists():
            print(f"Warning: {cond_dir} not found, skipping")
            continue

        known_file = cond_dir / "knownResults.txt"
        denovo_file = cond_dir / "homerMotifs.all.motifs"

        if known_file.exists():
            known_df = parse_known_motifs(known_file, top_n=25)
            known_df['Condition'] = condition
            all_known.append(known_df)

            # Extract TEAD motifs
            tead_df = find_tead_motifs(known_df)
            if not tead_df.empty:
                tead_df['Condition'] = condition
                tead_summary.append(tead_df)

        if denovo_file.exists():
            denovo_df = parse_denovo_motifs(denovo_file)
            denovo_df['Condition'] = condition
            all_denovo.append(denovo_df.head(10))  # Top 10 de novo per condition

    return {
        'known': pd.concat(all_known, ignore_index=True) if all_known else pd.DataFrame(),
        'denovo': pd.concat(all_denovo, ignore_index=True) if all_denovo else pd.DataFrame(),
        'tead': pd.concat(tead_summary, ignore_index=True) if tead_summary else pd.DataFrame()
    }


def create_comparison_table(known_df, metric='Pct_Target'):
    """Create cross-condition comparison for top TF families."""
    if known_df.empty:
        return pd.DataFrame()

    # Get most enriched TF families across all conditions
    top_families = known_df.groupby('TF_Family')['Fold_Enrichment'].mean().nlargest(15).index

    # Create pivot table
    subset = known_df[known_df['TF_Family'].isin(top_families)]
    pivot = subset.pivot_table(
        index='TF_Family',
        columns='Condition',
        values=[metric, 'Fold_Enrichment'],
        aggfunc='first'
    )

    return pivot


def format_executive_summary(results, output_path):
    """Create human-readable executive summary."""
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("HOMER MOTIF ANALYSIS - EXECUTIVE SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        # TEAD family summary (project-relevant)
        f.write("## TEAD FAMILY MOTIF ENRICHMENT\n")
        f.write("-" * 40 + "\n")
        if not results['tead'].empty:
            tead_df = results['tead'][['Condition', 'TF_Name', 'Pct_Target', 'Fold_Enrichment', 'P_value']]
            tead_df = tead_df.sort_values(['Condition', 'Fold_Enrichment'], ascending=[True, False])
            f.write(tead_df.to_string(index=False) + "\n\n")
        else:
            f.write("No TEAD family motifs found in top results.\n\n")

        # Top motifs per condition
        f.write("\n## TOP 10 ENRICHED MOTIFS PER CONDITION\n")
        f.write("-" * 40 + "\n")
        if not results['known'].empty:
            for condition in results['known']['Condition'].unique():
                f.write(f"\n### {condition}\n")
                cond_df = results['known'][results['known']['Condition'] == condition].head(10)
                display_cols = ['TF_Name', 'TF_Family', 'Consensus', 'Pct_Target', 'Fold_Enrichment']
                f.write(cond_df[display_cols].to_string(index=False) + "\n")

        # De novo motifs
        f.write("\n\n## TOP DE NOVO DISCOVERED MOTIFS\n")
        f.write("-" * 40 + "\n")
        if not results['denovo'].empty:
            for condition in results['denovo']['Condition'].unique():
                f.write(f"\n### {condition}\n")
                cond_df = results['denovo'][results['denovo']['Condition'] == condition].head(5)
                display_cols = ['Rank', 'Consensus', 'Pct_Target', 'Fold_Enrichment', 'P_value']
                f.write(cond_df[display_cols].to_string(index=False) + "\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("Key: Pct_Target = % of peaks containing motif\n")
        f.write("     Fold_Enrichment = Target% / Background%\n")
        f.write("=" * 80 + "\n")


def main():
    base_dir = Path("/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Eva_top/Eva_submission/SRF_Eva_Cut_and_Tag")
    homer_dir = base_dir / "results" / "08_homer_motifs"
    summary_dir = homer_dir / "summary"
    summary_dir.mkdir(exist_ok=True)

    print("Generating HOMER motif summary...")

    # Generate summary
    results = generate_summary_report(homer_dir)

    # Save CSV outputs
    if not results['known'].empty:
        # Top known motifs
        results['known'].to_csv(summary_dir / "top_known_motifs_all_conditions.csv", index=False)

        # Simplified version for stakeholders
        stakeholder_cols = ['Condition', 'TF_Name', 'TF_Family', 'Consensus',
                           'Pct_Target', 'Fold_Enrichment', 'P_value']
        results['known'][stakeholder_cols].to_csv(
            summary_dir / "motif_summary_simple.csv", index=False
        )
        print(f"  Saved: top_known_motifs_all_conditions.csv")

    if not results['denovo'].empty:
        results['denovo'].to_csv(summary_dir / "denovo_motifs_all_conditions.csv", index=False)
        print(f"  Saved: denovo_motifs_all_conditions.csv")

    if not results['tead'].empty:
        results['tead'].to_csv(summary_dir / "tead_family_motifs.csv", index=False)
        print(f"  Saved: tead_family_motifs.csv")

    # Create comparison table
    comparison = create_comparison_table(results['known'])
    if not comparison.empty:
        comparison.to_csv(summary_dir / "motif_comparison_across_conditions.csv")
        print(f"  Saved: motif_comparison_across_conditions.csv")

    # Generate executive summary
    format_executive_summary(results, summary_dir / "executive_summary.txt")
    print(f"  Saved: executive_summary.txt")

    print(f"\nAll summaries saved to: {summary_dir}/")
    print("\nKey files for stakeholders:")
    print("  - motif_summary_simple.csv     (spreadsheet-friendly)")
    print("  - tead_family_motifs.csv       (TEAD-specific results)")
    print("  - executive_summary.txt        (human-readable report)")


if __name__ == "__main__":
    main()
