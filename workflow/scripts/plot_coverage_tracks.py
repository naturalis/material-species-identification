#!/usr/bin/env python
"""Step 4: depth along the genome for the panel species with the most deduplicated reads.

Samples are pooled for the picture only, after per-library deduplication; the blank is drawn
on the same axes. Each bar shows the deepest position in its bin, so short piles stay visible.
"""

import argparse
import csv
import gzip

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLE, BLANK = "#2a78d6", "#eb6834"
SURFACE, INK, MUTED = "#fcfcfb", "#0b0b0b", "#52514e"
BINS = 400


def binned(intervals, length):
    change = [0] * (length + 1)
    for start, end in intervals:
        change[start - 1] += 1
        change[end] -= 1
    depth, running = [], 0
    for step in change[:-1]:
        running += step
        depth.append(running)
    width = max(1, -(-length // BINS))
    edges = [i / 1000 for i in range(0, length, width)]
    peaks = [max(depth[i:i + width]) for i in range(0, length, width)]
    return edges + [length / 1000], peaks + [peaks[-1]]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--blank", required=True)
    parser.add_argument("--assigned", nargs="+", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--top", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.sheet, newline="") as fh:
        sheet = {row["library"]: row for row in csv.DictReader(fh, delimiter="\t")}
    totals, lengths, names = {}, {}, {}
    with open(args.metrics, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lengths[row["accession"]] = int(row["genome_length"])
            names[row["accession"]] = row["species"]
            if row["role"] == "sample":
                totals[row["accession"]] = totals.get(row["accession"], 0) + int(row["dedup_reads"])
    chosen = [a for a in sorted(totals, key=lambda a: -totals[a]) if totals[a] > 0][: args.top]

    intervals = {(a, group): [] for a in chosen for group in ("samples", "blank")}
    for path in args.assigned:
        library = path.split("/")[-1].removesuffix(".assigned.tsv.gz")
        info = sheet[library]
        group = "samples" if info["role"] == "sample" else "blank" if info["id"] == args.blank else None
        if group is None:
            continue
        with gzip.open(path, "rt", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row["accession"] in chosen and row["duplicate"] == "0":
                    intervals[row["accession"], group].append((int(row["start"]), int(row["end"])))

    columns = 3
    rows = max(1, -(-len(chosen) // columns))
    figure, axes = plt.subplots(rows, columns, figsize=(12, 2.9 * rows + 0.6),
                                facecolor=SURFACE, squeeze=False)
    for ax in axes.flat[len(chosen):]:
        ax.set_visible(False)
    for ax, accession in zip(axes.flat, chosen):
        length = lengths[accession]
        x, sample = binned(intervals[accession, "samples"], length)
        _, blank = binned(intervals[accession, "blank"], length)
        ax.fill_between(x, sample, step="post", color=SAMPLE, linewidth=0, label="samples pooled")
        ax.step(x, blank, where="post", color=BLANK, linewidth=1.2, label=args.blank)
        ax.set_facecolor(SURFACE)
        ax.set_xlim(0, length / 1000)
        ax.set_ylim(0, max(max(sample), max(blank), 1) * 1.1)
        ax.yaxis.get_major_locator().set_params(integer=True)
        ax.set_title(
            f"{names[accession]} ({accession})\n"
            f"{len(intervals[accession, 'samples']):,} deduplicated reads in the samples · "
            f"{len(intervals[accession, 'blank']):,} in {args.blank}",
            fontsize=9.5, color=INK, loc="left",
        )
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(axis="y", color=MUTED, alpha=0.15, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlabel("position (kb)", fontsize=8.5, color=MUTED)
        ax.set_ylabel("reads covering position", fontsize=8.5, color=MUTED)
    if chosen:
        axes.flat[0].legend(frameon=False, fontsize=8.5, loc="upper right")
    figure.suptitle("Step 4: reads assigned to one panel species, after position deduplication",
                    fontsize=11.5, color=INK, x=0.01, ha="left")
    figure.text(0.01, 0.005, f"Each panel has its own scales; each of {BINS} bars shows its "
                "deepest position. Samples pooled for display only.", fontsize=8, color=MUTED)
    figure.tight_layout(rect=(0, 0.02, 1, 0.97))
    figure.savefig(args.out, dpi=150, facecolor=SURFACE)


if __name__ == "__main__":
    main()
