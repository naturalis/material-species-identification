#!/usr/bin/env python
"""Is a taxon's DNA above lab background, given how deeply each control was sequenced?

Gate 1 asks whether a taxon is >=5x the blank. That question has no answer when the blank
holds zero reads of it: the fold change is undefined, and filling the zero with a Poisson
upper bound makes the bar depend on an assumption rather than on data. This asks the question
the other way round, which the depths can answer:

    if this taxon were present in the control at the same rate per molecule as in the
    samples, how many reads would the control hold, and how unlikely is what we saw?

Counts are step 4's position-deduplicated assigned reads. The denominator is the library's
distinct molecule count (step 5's duplicate histogram), because a deduplicated count scales
with molecules sampled, not with reads. The two deduplications are not identical -- exact
sequence for the denominator, mapped position for the numerator -- so the rates are
comparable between libraries, not absolute.

The test is the exact conditional (binomial) test for two Poisson rates: given K reads in
N sample molecules and k in n control molecules, k | K+k ~ Binomial(K+k, n/(N+n)) under
equal rates. One-sided, testing whether the control is *depleted* relative to the samples.

Where several sample libraries come from ONE extract, pooling them measures that extract more
deeply. It is not several independent observations, and the per-library rows are there to show
the spread, not to be combined as replicates. Check how your libraries relate before reading
the pooled row as replication.
"""

import argparse
import csv
import math
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, SURFACE, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e4e3df"
SAMPLE_COLOUR, CONTROL_COLOUR = "#2a78d6", "#eb6834"
# 95% upper bound on a Poisson mean that produced zero: -ln(0.05).
ZERO_UPPER = 3.0
PER = 1e6


def log_binomial_cdf(k, n, p):
    """log P(X <= k) for X ~ Binomial(n, p), summed in log space."""
    if k >= n:
        return 0.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return -math.inf
    terms = [
        math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
        + i * math.log(p) + (n - i) * math.log1p(-p)
        for i in range(k + 1)
    ]
    top = max(terms)
    return top + math.log(sum(math.exp(t - top) for t in terms))


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metrics", required=True, help="per_reference_metrics.tsv")
    parser.add_argument("--complexity", required=True)
    parser.add_argument("--blank", required=True, help="the extraction blank's id")
    parser.add_argument("--min-sample-reads", type=int, default=1)
    parser.add_argument("--out", required=True)
    parser.add_argument("--figure", required=True)
    args = parser.parse_args()

    molecules, role = {}, {}
    for row in read_tsv(args.complexity):
        molecules[row["id"]] = int(row["distinct"])
        role[row["id"]] = row["role"]
    samples = [i for i in molecules if role[i] == "sample"]
    controls = [i for i in molecules if role[i] == "control"]
    controls.sort(key=lambda i: (i != args.blank, i))

    counts = defaultdict(dict)
    for row in read_tsv(args.metrics):
        counts[row["species"]][row["id"]] = int(row["dedup_reads"])
    taxa = {s: sum(counts[s].get(i, 0) for i in samples) for s in counts}
    taxa = {s: n for s, n in taxa.items() if n >= args.min_sample_reads}
    order = sorted(taxa, key=lambda s: (-taxa[s], s))

    sample_molecules = sum(molecules[i] for i in samples)
    rows = []
    for species in order:
        pooled = taxa[species]
        for group, ids in [("samples pooled", samples), *((i, [i]) for i in samples)]:
            k = sum(counts[species].get(i, 0) for i in ids)
            n = sum(molecules[i] for i in ids)
            rows.append({
                "species": species, "id": group, "role": "sample",
                "dedup_reads": k, "distinct_molecules": n,
                "per_million": round(k / n * PER, 3),
                "compared_with": "", "expected_if_background": "", "p_depleted": "",
                "fold_over_control": "", "control_detection_limit_per_million": "",
            })
        for control in controls:
            k, n = counts[species].get(control, 0), molecules[control]
            expected = pooled / sample_molecules * n
            p = math.exp(log_binomial_cdf(k, pooled + k, n / (sample_molecules + n)))
            # With k = 0 the fold change is undefined; bound the control rate instead.
            bound = max(k, ZERO_UPPER)
            rows.append({
                "species": species, "id": control, "role": "control",
                "dedup_reads": k, "distinct_molecules": n,
                "per_million": round(k / n * PER, 3),
                "compared_with": "samples pooled",
                "expected_if_background": round(expected, 2),
                "p_depleted": f"{p:.3g}",
                "fold_over_control": round((pooled / sample_molecules) / (bound / n), 1),
                "control_detection_limit_per_million": round(ZERO_UPPER / n * PER, 3),
            })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    plotted = order[:8]
    fig, ax = plt.subplots(figsize=(10, 0.62 * len(plotted) + 2.4))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for offset, species in enumerate(plotted):
        y = len(plotted) - offset
        rate = taxa[species] / sample_molecules * PER
        ax.plot([rate], [y + 0.16], "o", color=SAMPLE_COLOUR, markersize=9, zorder=3)
        blank_count = counts[species].get(args.blank, 0)
        blank_rate = blank_count / molecules[args.blank] * PER
        limit = ZERO_UPPER / molecules[args.blank] * PER
        marker = "o" if blank_count else "<"
        ax.plot([blank_rate if blank_count else limit], [y - 0.16], marker,
                color=CONTROL_COLOUR, markersize=9, zorder=3)
        ax.plot([min(rate, blank_rate or limit), max(rate, blank_rate or limit)],
                [y, y], color=GRID, linewidth=8, zorder=1, solid_capstyle="round")
    ax.set_yticks(range(len(plotted), 0, -1),
                  [f"{s}\n{taxa[s]} reads" for s in plotted], fontsize=8, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel("deduplicated reads per million distinct molecules (log scale)",
                  color=MUTED, fontsize=9)
    ax.set_title("Every taxon against the blank, corrected for sequencing depth\n"
                 "A left arrow means the blank held none; the point is all the depth could "
                 "have excluded", color=INK, fontsize=10, loc="left")
    handles = [
        plt.Line2D([], [], marker="o", color=SAMPLE_COLOUR, linestyle="none", markersize=9,
                   label="samples pooled"),
        plt.Line2D([], [], marker="o", color=CONTROL_COLOUR, linestyle="none", markersize=9,
                   label=f"{args.blank}, reads found"),
        plt.Line2D([], [], marker="<", color=CONTROL_COLOUR, linestyle="none", markersize=9,
                   label=f"{args.blank}, none found: 95% upper limit, so the true rate is "
                         "somewhere to the left"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower right")
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    fig.tight_layout()
    fig.savefig(args.figure, dpi=150, facecolor=SURFACE)

    for species in order[:6]:
        blank = next(r for r in rows if r["species"] == species and r["id"] == args.blank)
        print(f"{species}: {taxa[species]} in samples, {blank['dedup_reads']} in {args.blank}, "
              f"expected {blank['expected_if_background']}, p={blank['p_depleted']}")


if __name__ == "__main__":
    main()
