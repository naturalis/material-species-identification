#!/usr/bin/env python
"""Step 5: summarise damage, fragment length and library complexity, with the blank alongside.

Damage and lengths are given per library and for the sample libraries pooled. Pooling is for
statistical power only; it comes after per-library deduplication. A raised first position
against the interior is the deamination signal. With full USER treatment none is expected.
"""

import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

POOLED = "samples pooled"
FIVE, THREE, OUTLINE = "#2a78d6", "#eb6834", "#52514e"
SURFACE, INK, MUTED = "#fcfcfb", "#0b0b0b", "#52514e"
MULTIPLIERS = (2, 5, 10)
MIN_PLOTTED = 50  # reads; damage_summary.tsv still lists every taxon


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def wilson(k, n, z=1.96):
    if n == 0:
        return "", ""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return round(max(0, centre - half), 4), round(min(1, centre + half), 4)


def rate(k, n):
    return round(k / n, 4) if n else ""


def interpolate(points, x):
    """Linear interpolation on preseq's (total reads, expected distinct) curve."""
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            return int(y0 + (y1 - y0) * (x - x0) / (x1 - x0)) if x1 > x0 else int(y0)
    return ""


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color=MUTED, alpha=0.15, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--blank", required=True)
    parser.add_argument("--taxa", nargs="+", required=True)
    parser.add_argument("--interior", required=True, help="first-last, e.g. 10-25")
    parser.add_argument("--damage", nargs="+", required=True)
    parser.add_argument("--lengths", nargs="+", required=True)
    parser.add_argument("--curves", nargs="+", required=True)
    parser.add_argument("--histograms", nargs="+", required=True)
    parser.add_argument("--library-lengths", required=True)
    parser.add_argument("--damage-out", required=True)
    parser.add_argument("--damage-plot", required=True)
    parser.add_argument("--lengths-out", required=True)
    parser.add_argument("--lengths-plot", required=True)
    parser.add_argument("--complexity-out", required=True)
    args = parser.parse_args()

    with open(args.sheet, newline="") as fh:
        sheet = {row["library"]: row for row in csv.DictReader(fh, delimiter="\t")}
    first, last = map(int, args.interior.split("-"))
    order = [row["id"] for row in sheet.values() if row["role"] == "sample"]
    order += [POOLED] + [row["id"] for row in sheet.values() if row["role"] == "control"]

    def groups_of(library):
        info = sheet[library]
        return [info["id"], POOLED] if info["role"] == "sample" else [info["id"]]

    # Damage: (group, species, end, position) -> [reference bases, substitutions]
    damage, reads = defaultdict(lambda: [0, 0]), Counter()
    for row in read_rows(args.damage):
        for group in groups_of(row["library"]):
            key = (group, row["species"], row["end"], int(row["position"]))
            damage[key][0] += int(row["reference_bases"])
            damage[key][1] += int(row["substitutions"])
            if row["end"] == "other":
                reads[group, row["species"]] += int(row["reads"])

    def window(group, species, end, positions):
        n = sum(damage[group, species, end, p][0] for p in positions)
        k = sum(damage[group, species, end, p][1] for p in positions)
        return k, n

    header = ["id", "species", "reads", "ct_5p1", "ct_5p1_n", "ct_5p1_low", "ct_5p1_high",
              f"ct_5p_{first}_{last}", "ga_3p1", "ga_3p1_n", "ga_3p1_low", "ga_3p1_high",
              f"ga_3p_{first}_{last}", "other_per_base"]
    with open(args.damage_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for group in order:
            for species in args.taxa:
                if not reads[group, species]:
                    continue
                ct, ct_n = window(group, species, "5p", [1])
                ga, ga_n = window(group, species, "3p", [1])
                ct_in = window(group, species, "5p", range(first, last + 1))
                ga_in = window(group, species, "3p", range(first, last + 1))
                other = damage[group, species, "other", 0]
                writer.writerow([group, species, reads[group, species],
                                 rate(ct, ct_n), ct_n, *wilson(ct, ct_n), rate(*ct_in),
                                 rate(ga, ga_n), ga_n, *wilson(ga, ga_n), rate(*ga_in),
                                 rate(other[1], other[0])])

    # Left: the position profile, only for the taxon with the most reads; thin taxa are noise.
    # Right: the terminal rate with its 95% interval for every taxon, against the interior.
    # Below about 50 reads the interval spans most of the axis and hides the others.
    shown = [s for s in args.taxa if reads[POOLED, s] >= MIN_PLOTTED]
    if shown:
        figure, (left, right) = plt.subplots(
            1, 2, figsize=(12, 3.8), facecolor=SURFACE, gridspec_kw={"width_ratios": [1, 1.4]})
        top = max(shown, key=lambda s: reads[POOLED, s])
        positions = range(1, 1 + max(p for (g, s, e, p) in damage if g == POOLED and e == "5p"))
        for end, colour, label in (("5p", FIVE, "C→T from 5′ end"), ("3p", THREE, "G→A from 3′ end")):
            values = [rate(*window(POOLED, top, end, [p])) or 0 for p in positions]
            left.plot(list(positions), values, color=colour, linewidth=2, marker="o",
                      markersize=3, label=label)
        style(left)
        left.set_ylim(bottom=0)
        left.set_title(f"{top}: by position\n{reads[POOLED, top]:,} reads", fontsize=9.5,
                       color=INK, loc="left")
        left.set_xlabel("position from read end", fontsize=8.5, color=MUTED)
        left.set_ylabel("substitution frequency", fontsize=8.5, color=MUTED)
        left.legend(frameon=False, fontsize=8)

        for i, species in enumerate(shown):
            for offset, (end, colour) in ((-0.12, ("5p", FIVE)), (0.12, ("3p", THREE))):
                k, n = window(POOLED, species, end, [1])
                if not n:
                    continue
                low, high = wilson(k, n)
                right.errorbar(i + offset, k / n, yerr=[[k / n - low], [high - k / n]],
                               fmt="o", color=colour, markersize=5, capsize=3, linewidth=1.5)
                interior = rate(*window(POOLED, species, end, range(first, last + 1))) or 0
                right.hlines(interior, i + offset - 0.1, i + offset + 0.1, color=OUTLINE,
                             linewidth=1.5)
        style(right)
        right.set_xticks(range(len(shown)))
        right.set_xticklabels([f"{s}\n{reads[POOLED, s]:,} reads" for s in shown], fontsize=8,
                              color=MUTED)
        right.set_ylim(bottom=0)
        right.set_ylabel("rate at the terminal base", fontsize=8.5, color=MUTED)
        right.set_title("Terminal base, 95% interval; grey bar = interior positions "
                        f"{first}–{last}", fontsize=9.5, color=INK, loc="left")
        figure.suptitle("Step 5: deamination-type substitutions (samples pooled)", fontsize=11,
                        color=INK, x=0.01, ha="left")
        figure.tight_layout(rect=(0, 0, 1, 0.93))
        figure.savefig(args.damage_plot, dpi=150, facecolor=SURFACE)

    # Lengths: assigned taxa, and every merged read of the library as the reference.
    lengths = defaultdict(Counter)
    for row in read_rows(args.lengths):
        for group in groups_of(row["library"]):
            lengths[group, row["species"]][int(row["length"])] += int(row["reads"])
    with open(args.library_lengths, newline="") as fh:
        by_id = {row["id"]: row for row in sheet.values()}
        for row in csv.DictReader(fh, delimiter="\t"):
            library = by_id[row["id"]]["library"]
            for group in groups_of(library):
                lengths[group, "all merged reads"][int(row["length"])] += int(row["reads"])

    def describe(counter):
        values = sorted(counter.elements())
        return {
            "reads": len(values),
            "median": statistics.median(values),
            "mean": round(statistics.mean(values), 1),
            "q25": values[len(values) // 4],
            "q75": values[(3 * len(values)) // 4],
            "below_50": round(sum(v < 50 for v in values) / len(values), 3),
        }

    with open(args.lengths_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "species", "reads", "median", "mean", "q25", "q75", "below_50"])
        for group in order:
            for species in [*args.taxa, "all merged reads"]:
                if lengths[group, species]:
                    writer.writerow([group, species, *describe(lengths[group, species]).values()])

    panels = [(POOLED, s) for s in args.taxa if lengths[POOLED, s]][:4]
    panels.append((args.blank, "all merged reads"))
    figure, axes = plt.subplots(1, len(panels), figsize=(3.3 * len(panels), 3.4),
                                facecolor=SURFACE, squeeze=False)
    bins = list(range(25, 300, 10))
    reference = lengths[POOLED, "all merged reads"]
    for ax, (group, species) in zip(axes[0], panels):
        counter = lengths[group, species]
        ax.hist(list(counter), bins=bins, weights=list(counter.values()), density=True,
                color=FIVE if group == POOLED else THREE, edgecolor=SURFACE, linewidth=0.6)
        ax.hist(list(reference), bins=bins, weights=list(reference.values()), density=True,
                histtype="step", color=OUTLINE, linewidth=1.2, label="all merged reads, samples pooled")
        style(ax)
        stats = describe(counter)
        label = species if group == POOLED else f"{group}: all merged reads"
        ax.set_title(f"{label}\n{stats['reads']:,} reads · median {stats['median']:g} bp",
                     fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("fragment length (bp)", fontsize=8.5, color=MUTED)
        ax.set_yticks([])
    axes[0][0].legend(frameon=False, fontsize=8)
    figure.suptitle("Step 5: fragment length by taxon (samples pooled, deduplicated), against all reads",
                    fontsize=11, color=INK, x=0.01, ha="left")
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(args.lengths_plot, dpi=150, facecolor=SURFACE)

    # Complexity: preseq's extrapolated distinct molecules at more sequencing.
    with open(args.complexity_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "role", "reads", "distinct",
                         *[f"distinct_at_{m}x" for m in MULTIPLIERS], "note"])
        for histogram, curve in zip(sorted(args.histograms), sorted(args.curves)):
            library = histogram.split("/")[-1].removesuffix(".duplicate_histogram.tsv")
            info = sheet[library]
            with open(histogram) as hf:
                pairs = [tuple(map(int, line.split())) for line in hf if line.strip()]
            total = sum(c * n for c, n in pairs)
            distinct = sum(n for _, n in pairs)
            with open(curve, newline="") as cf:
                points = [(float(r["TOTAL_READS"]), float(r["EXPECTED_DISTINCT"]))
                          for r in csv.DictReader(cf, delimiter="\t")]
            projected = []
            for m in MULTIPLIERS:
                projected.append(interpolate(points, m * total))
            note = "" if points else "preseq could not extrapolate"
            writer.writerow([info["id"], info["role"], total, distinct, *projected, note])


if __name__ == "__main__":
    main()
