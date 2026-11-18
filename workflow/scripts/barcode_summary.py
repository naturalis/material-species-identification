#!/usr/bin/env python
"""Step 6 summary: BOLD hits by species, and one figure for lineage placement and recovery.

The query is a consensus with N where it is masked, so BLAST's own identity counts every N as
a mismatch, and megablast breaks the query into several segments at runs of N. So each BOLD
record is scored over all its segments, on called query bases only: matches and called bases
aligned, summed. Species are ranked by their best record, on matching bases first. A top
species from the uncurated release still needs the §3b checks by hand (length, real name,
stop codons, contaminants).
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1f9e89", "#8a5cc2"]
TOP = 5
COLUMNS = ["db", "qseqid", "sseqid", "stitle", "pident", "length", "mismatch", "gapopen", "qstart",
           "qend", "sstart", "send", "evalue", "bitscore", "slen", "qseq", "sseq"]
BASES = set("ACGT")


def species_of(title):
    match = re.search(r"s__([^;]*)", title)
    name = match.group(1) if match else ""
    return name if name and name != "None" else "(no species name)"


def summarise_bold(paths, taxa, out):
    records = defaultdict(lambda: {"called": 0, "matches": 0, "segments": 0})
    for taxon, path in zip(taxa, paths):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                row = dict(zip(COLUMNS, line.rstrip("\n").split("\t")))
                key = (taxon, row["qseqid"], row["db"], row["sseqid"])
                record = records[key]
                pairs = [(q, t) for q, t in zip(row["qseq"].upper(), row["sseq"].upper())
                         if q in BASES]
                record["called"] += len(pairs)
                record["matches"] += sum(q == t for q, t in pairs)
                record["segments"] += 1
                record["title"], record["slen"] = row["stitle"], row["slen"]
    best = defaultdict(dict)
    counts = defaultdict(lambda: defaultdict(int))
    for (taxon, query, db, subject), record in records.items():
        species = species_of(record["title"])
        counts[taxon, query, db][species] += 1
        held = best[taxon, query, db].get(species)
        if held is None or (record["matches"], -record["called"]) > (held["matches"], -held["called"]):
            best[taxon, query, db][species] = {**record, "subject": subject}
    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon", "query", "library", "rank", "species", "records_hit",
                         "matching_bases", "called_bases_aligned", "identity_on_called",
                         "segments", "subject_length", "marker", "record"])
        for key in sorted(best):
            ordered = sorted(best[key].values(),
                             key=lambda r: (-r["matches"], r["called"] - r["matches"]))
            for position, record in enumerate(ordered[:TOP], 1):
                parts = record["subject"].split("|")
                species = species_of(record["title"])
                writer.writerow([key[0], key[1].split("|", 1)[-1], key[2], position, species,
                                 counts[key][species], record["matches"], record["called"],
                                 round(record["matches"] / record["called"], 4),
                                 record["segments"], record["slen"],
                                 parts[2] if len(parts) > 2 else "",
                                 parts[1] if len(parts) > 1 else record["subject"]])


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_facecolor(SURFACE)


def figure(lineage_paths, recovery_paths, predicted_paths, out):
    panels = max(1, len(lineage_paths))
    fig, axes = plt.subplots(panels, 2, figsize=(11, 4.4 * panels), squeeze=False)
    fig.patch.set_facecolor(SURFACE)
    for row, (lineage, recovery, predicted) in enumerate(
            zip(lineage_paths, recovery_paths, predicted_paths)):
        left, right = axes[row]
        all_clusters = read_tsv(lineage)
        clusters = [r for r in all_clusters
                    if int(r["reads_favouring"]) > 0 or int(r["genomes"]) >= 5]
        labels = [f"{r['cluster']} ({r['genomes']} genomes)\n{r['organisms'].split(';')[0]}"
                  + ("\ncontains the panel genome" if r["contains_reference"] == "True" else "")
                  for r in clusters]
        y = range(len(clusters))
        left.barh(y, [int(r["reads_favouring"]) for r in clusters], color=SERIES[0], height=0.6,
                  label="fewer edits than any other cluster")
        left.barh(y, [int(r["reads_favouring_margin2"]) for r in clusters], color=SERIES[1],
                  height=0.6, label="by 2 or more edits")
        left.set_yticks(list(y), labels, fontsize=7, color=INK)
        left.invert_yaxis()
        decisive = sum(int(r["reads_favouring"]) for r in all_clusters)
        left.set_xlabel("deduplicated reads favouring the cluster", color=MUTED, fontsize=8)
        left.set_title(f"Which lineage do the reads fit best?\n{decisive} of the reads tell the "
                       f"lineages apart at all", color=INK, fontsize=9, loc="left")
        left.legend(frameon=False, fontsize=7, loc="lower right")
        style(left)

        curves = defaultdict(list)
        for r in read_tsv(recovery):
            curves[r["truth"]].append((int(r["length"]), float(r["recovery"])))
        notes = {r["truth"]: r for r in read_tsv(predicted)}
        for colour, (truth, points) in zip(SERIES, curves.items()):
            xs, ys = zip(*sorted(points))
            note = notes[truth]
            right.plot(xs, ys, color=colour, linewidth=2,
                       label=f"{truth} ({note['accession']}): predicted median "
                             f"{note['predicted_median_if_human_lengths']} bp")
        any_note = next(iter(notes.values()))
        right.set_ylim(0, 1.05)
        right.set_xlabel("simulated read length (bp)", color=MUTED, fontsize=8)
        right.set_ylabel("share assigned to the panel genome", color=MUTED, fontsize=8)
        right.set_title("Does mapping lose this taxon's short reads?\nWith the human length "
                        f"distribution (median {any_note['human_median']} bp) the recovered\n"
                        f"median would be as below. Observed: {any_note['observed_median']} bp",
                        color=INK, fontsize=9, loc="left")
        right.legend(frameon=False, fontsize=7, loc="lower right")
        right.grid(axis="y", color="#e4e3df", linewidth=0.6)
        style(right)
    fig.tight_layout()
    fig.savefig(out, dpi=150, facecolor=SURFACE)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--taxa", nargs="+", required=True)
    parser.add_argument("--hits", nargs="+", required=True)
    parser.add_argument("--summaries", nargs="+", required=True)
    parser.add_argument("--lineage", nargs="*", default=[])
    parser.add_argument("--genome-scores", nargs="*", default=[])
    parser.add_argument("--recovery", nargs="*", default=[])
    parser.add_argument("--predicted", nargs="*", default=[])
    parser.add_argument("--bold-out", required=True)
    parser.add_argument("--figure", required=True)
    args = parser.parse_args()
    summarise_bold(args.hits, args.taxa, args.bold_out)
    figure(args.lineage, args.recovery, args.predicted, args.figure)


if __name__ == "__main__":
    main()
