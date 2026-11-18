#!/usr/bin/env python
"""Step 6: depth along a taxon's mitogenome, with the barcode gene marked.

Answers two questions a breadth percentage hides: is the coverage spread over the genome or
piled on a few spots, and does it happen to fall on the barcode gene. The masking threshold
is drawn, because a base only enters the consensus above it.
"""

import argparse
import csv
import gzip

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, SURFACE, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
DEPTH, GENE, THRESHOLD = "#2a78d6", "#eb6834", "#52514e"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--site-counts", nargs="+", required=True)
    parser.add_argument("--taxa", nargs="+", required=True)
    parser.add_argument("--accessions", nargs="+", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--min-depth", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.cds, newline="") as fh:
        genes = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}

    fig, axes = plt.subplots(len(args.taxa), 1, figsize=(11, 3.1 * len(args.taxa)), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for ax, taxon, accession, path in zip(axes[:, 0], args.taxa, args.accessions, args.site_counts):
        gene = genes[accession]
        start, end = int(gene["start"]), int(gene["end"])
        depth = {}
        with gzip.open(path, "rt", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["id"] == "samples pooled":
                    depth[int(row["position"])] = sum(int(row[b]) for b in "ACGT")
        length = max(max(depth, default=1), end)
        values = [depth.get(i, 0) for i in range(1, length + 1)]
        in_gene = values[start - 1:end]

        ax.axvspan(start / 1000, end / 1000, color=GENE, alpha=0.13, zorder=0)
        ax.fill_between([i / 1000 for i in range(1, length + 1)], values, color=DEPTH,
                        linewidth=0, zorder=2)
        ax.axhline(args.min_depth, color=THRESHOLD, linewidth=1, linestyle=(0, (4, 3)), zorder=3)
        covered = sum(v > 0 for v in values) / length
        deep = sum(v >= args.min_depth for v in values) / length
        gene_covered = sum(v > 0 for v in in_gene) / len(in_gene)
        gene_deep = sum(v >= args.min_depth for v in in_gene) / len(in_gene)
        ax.set_title(
            f"{taxon.replace('_', ' ')} ({accession}) — whole genome: {covered:.0%} has any read, "
            f"{deep:.0%} reaches {args.min_depth}x\n"
            f"{gene['gene']} (shaded, {start:,}-{end:,}): {gene_covered:.0%} has any read, "
            f"{gene_deep:.0%} reaches {args.min_depth}x",
            color=INK, fontsize=10, loc="left")
        ax.set_xlim(0, length / 1000)
        ax.set_ylim(0, max(max(values), args.min_depth + 1) * 1.1)
        ax.set_xlabel("position in the mitogenome (kb)", color=MUTED, fontsize=9)
        ax.set_ylabel("reads covering", color=MUTED, fontsize=9)
        ax.annotate(f"{args.min_depth}x, the masking threshold", (length / 1000, args.min_depth),
                    xytext=(-4, 4), textcoords="offset points", ha="right", fontsize=7,
                    color=THRESHOLD)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150, facecolor=SURFACE)


if __name__ == "__main__":
    main()
