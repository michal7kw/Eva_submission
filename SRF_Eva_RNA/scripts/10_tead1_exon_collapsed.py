#!/usr/bin/env python3
# =============================================================================
# 10_tead1_exon_collapsed.py
#
# Exon-collapsed RNA-seq coverage snapshot of the TEAD1 locus, Mock (red) vs
# TES (blue), 3 replicates each. Unlike the genomic-coordinate pyGenomeTracks
# view (9_tead1_locus_snapshot.sh), this draws signal ONLY over exons and
# renders introns as narrow fixed-width skipped gaps -- so the exons (where the
# RNA-seq reads actually pile up) get the horizontal space and the figure reads
# like the intronless SOX2 reference panel.
#
# Approach:
#   - take the UNION of all TEAD1 exons (merged meta-gene, matching the
#     merge_transcripts=true gene model used by script 9);
#   - lay each exon side-by-side at its true bp width (exons are to-scale),
#     separated by a small constant gap (introns, NOT to scale);
#   - pull per-base coverage per exon from each bigWig and fill_between only
#     over exon segments, leaving the gaps blank.
#
# Dependencies (all already in the `pygenometracks` conda env):
#   pyBigWig, numpy, matplotlib
#
# Invoked by 10_tead1_exon_collapsed.sh; can also be run standalone (see --help).
# =============================================================================

import argparse
import sys

import numpy as np
import pyBigWig

import matplotlib
matplotlib.use("Agg")  # headless cluster rendering
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# -----------------------------------------------------------------------------
# Defaults / tunables (CLI overrides everything below)
# -----------------------------------------------------------------------------
MOCK_COLOR = "#B22222"   # red  (Mock = GFP control)   -- matches script 9
TES_COLOR = "#1F78B4"    # blue (TES)                  -- matches script 9


def parse_exons(gtf_path, gene):
    """Return (chrom, strand, [(start,end), ...]) of UNION exons (1-based incl)."""
    raw = []
    chrom = strand = None
    needle = 'gene_name "%s"' % gene
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            if needle not in f[8]:
                continue
            chrom = f[0]
            strand = f[6]
            raw.append((int(f[3]), int(f[4])))
    if not raw:
        sys.exit("ERROR: no exons for %s found in %s" % (gene, gtf_path))

    # Merge overlapping / adjacent exons into a union set
    raw.sort()
    merged = [list(raw[0])]
    for s, e in raw[1:]:
        if s <= merged[-1][1] + 1:          # overlap or directly adjacent
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return chrom, strand, [(s, e) for s, e in merged]


def exon_coverage(bw, chrom, s, e):
    """Per-base coverage over [s, e] (1-based incl); NaN -> 0."""
    arr = np.asarray(bw.values(chrom, s - 1, e), dtype=float)  # bigWig is 0-based half-open
    return np.nan_to_num(arr, nan=0.0)


def smooth(y, win):
    """Centered rolling mean of width `win` bp (no-op if win <= 1)."""
    if win is None or win <= 1 or y.size < win:
        return y
    kernel = np.ones(win) / win
    return np.convolve(y, kernel, mode="same")


def build_layout(exons, strand, gap_bp, flank_bp):
    """Order exons 5'->3' and assign collapsed x-coordinates.

    Returns a list of dicts: {gstart, gend, width, x0, x1} in plot order, where
    x0..x1 is the collapsed (drawing) coordinate span for that exon. A flank is
    added to the outermost (5' and 3') exons only.
    """
    order = list(range(len(exons)))
    if strand == "-":
        order = order[::-1]  # so the leftmost drawn exon is the 5' end

    layout = []
    x = 0.0
    for rank, i in enumerate(order):
        gs, ge = exons[i]
        # extend the first/last *drawn* exon outward by flank_bp (visual only)
        left_pad = flank_bp if rank == 0 else 0
        right_pad = flank_bp if rank == len(order) - 1 else 0
        gs_p, ge_p = gs - left_pad, ge + right_pad
        width = ge_p - gs_p + 1
        layout.append(dict(gstart=gs_p, gend=ge_p, width=width, x0=x, x1=x + width))
        x += width + gap_bp
    return layout, order


def main():
    ap = argparse.ArgumentParser(description="Exon-collapsed coverage snapshot.")
    ap.add_argument("--gtf", required=True, help="GTF (subset) containing the gene's exons")
    ap.add_argument("--gene", default="TEAD1")
    ap.add_argument("--bw-mock", nargs="+", required=True, help="Mock bigWigs (red)")
    ap.add_argument("--bw-tes", nargs="+", required=True, help="TES bigWigs (blue)")
    ap.add_argument("--out-prefix", required=True, help="output path prefix (no extension)")
    ap.add_argument("--gap-bp", type=int, default=None, help="intron gap width (bp); default ~1%% of exonic length")
    ap.add_argument("--flank-bp", type=int, default=0, help="bp added to the outer 5'/3' exons")
    ap.add_argument("--smooth-bp", type=int, default=0, help="rolling-mean window (bp); 0 = off")
    ap.add_argument("--ymax", type=float, default=None, help="fixed shared y-axis max (override)")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    chrom, strand, exons = parse_exons(args.gtf, args.gene)
    exonic_len = sum(e - s + 1 for s, e in exons)
    gap_bp = args.gap_bp if args.gap_bp is not None else max(50, round(0.01 * exonic_len))
    print("%s: %d union exons, %d exonic bp, strand %s, gap=%d bp"
          % (args.gene, len(exons), exonic_len, strand, gap_bp))

    layout, order = build_layout(exons, strand, gap_bp, args.flank_bp)
    total_x = layout[-1]["x1"]

    mock_files = args.bw_mock
    tes_files = args.bw_tes
    samples = ([("Mock %d" % (i + 1), f, MOCK_COLOR) for i, f in enumerate(mock_files)]
               + [("TES %d" % (i + 1), f, TES_COLOR) for i, f in enumerate(tes_files)])

    # --- gather per-exon coverage for every sample, and the shared y-max -------
    cov = {}  # sample label -> list of per-exon y arrays (in plot order)
    ymax_data = 0.0
    for label, fpath, _ in samples:
        bw = pyBigWig.open(fpath)
        per_exon = []
        for blk in layout:
            y = exon_coverage(bw, chrom, blk["gstart"], blk["gend"])
            y = smooth(y, args.smooth_bp)
            per_exon.append(y)
            if y.size:
                ymax_data = max(ymax_data, float(y.max()))
        bw.close()
        cov[label] = per_exon
    ymax = args.ymax if args.ymax is not None else (ymax_data * 1.08 if ymax_data > 0 else 1.0)
    print("  shared y-axis max = %.3f%s" % (ymax, " (override)" if args.ymax else ""))

    # --- figure: 6 coverage panels + 1 gene-model strip -----------------------
    n = len(samples)
    fig, axes = plt.subplots(
        n + 1, 1, figsize=(16, 1.05 * n + 2.2), sharex=True,
        gridspec_kw=dict(height_ratios=[1.0] * n + [0.55], hspace=0.12),
    )

    def shade_exons(ax):
        for k, blk in enumerate(layout):
            if k % 2 == 0:
                ax.axvspan(blk["x0"], blk["x1"], color="0.5", alpha=0.06, lw=0)

    # coverage panels
    for ax, (label, _f, color) in zip(axes[:n], samples):
        shade_exons(ax)
        for blk, y in zip(layout, cov[label]):
            x = blk["x0"] + np.arange(y.size)
            # minus-strand: reverse within-exon so the axis reads 5'->3' L->R
            yy = y[::-1] if strand == "-" else y
            ax.fill_between(x, 0, yy, color=color, linewidth=0)
        ax.set_xlim(0, total_x)
        ax.set_ylim(0, ymax)
        ax.set_yticks([])
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(bottom=False)
        # pyGenomeTracks-style data-range tag (top-left) + sample label (right)
        ax.text(0.002, 0.92, "[0 - %.0f]" % ymax, transform=ax.transAxes,
                va="top", ha="left", fontsize=8, color="0.35")
        ax.text(1.004, 0.5, label, transform=ax.transAxes, va="center",
                ha="left", fontsize=11, color=color, fontweight="bold")

    # gene-model strip (exon boxes joined by a thin intron line)
    gax = axes[n]
    gax.set_ylim(0, 1)
    gax.set_xlim(0, total_x)
    gax.axis("off")
    ymid, box_h = 0.5, 0.5
    # intron line across the whole span
    gax.plot([layout[0]["x0"], layout[-1]["x1"]], [ymid, ymid],
             color="darkblue", lw=1.0, zorder=1)
    for blk in layout:
        gax.add_patch(Rectangle((blk["x0"], ymid - box_h / 2), blk["width"], box_h,
                                facecolor="darkblue", edgecolor="darkblue", zorder=2))
    # axis is oriented 5'->3' left-to-right for either strand
    gax.text(0.0, -0.15, "%s  5'→  union exons (%s)" % (args.gene, strand),
             transform=gax.transAxes, va="top", ha="left", fontsize=10,
             color="darkblue", fontweight="bold")

    # in-exon scale bar (valid because exons are to-scale); 1 kb if it fits
    bar_bp = 1000 if total_x >= 3000 else max(100, int(total_x / 5))
    bar_label = "%d kb" % (bar_bp // 1000) if bar_bp >= 1000 else "%d bp" % bar_bp
    x_bar0 = total_x - bar_bp
    gax.plot([x_bar0, x_bar0 + bar_bp], [-0.18, -0.18], transform=gax.get_xaxis_transform(),
             color="black", lw=2.2, clip_on=False)
    gax.text(x_bar0 + bar_bp / 2, -0.30, bar_label, transform=gax.get_xaxis_transform(),
             va="top", ha="center", fontsize=9, clip_on=False)

    fig.suptitle("%s locus — RNA-seq coverage (CPM), exon-collapsed" % args.gene,
                 fontsize=14, fontweight="bold", y=0.985)
    fig.text(0.5, 0.945, "exons to scale; introns shown as narrow gaps (not to scale)",
             ha="center", fontsize=9, color="0.4")

    fig.subplots_adjust(left=0.04, right=0.93, top=0.93, bottom=0.07)
    for ext in ("pdf", "png"):
        out = "%s.%s" % (args.out_prefix, ext)
        fig.savefig(out, dpi=args.dpi)
        print("  wrote %s" % out)


if __name__ == "__main__":
    main()
