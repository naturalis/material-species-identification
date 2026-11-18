#!/usr/bin/env python
"""Step 4: assign each read to one panel genome by edit distance, or to the ambiguous bin.

Every read was mapped to every panel genome on its own. A read goes to the genome where its
edit distance (NM) is lowest, when every other genome is at least --margin edits worse and
the read has at least --min-mapq on the winner. MAPQ is judged on that genome alone (§6),
never on the combined panel, where it would discard every read two congeners fit equally.
Reads within the margin go to the ambiguous bin, under the species they tie between.
Assigned reads are then deduplicated by position: same genome, start, end and strand.

Query names u<n>_<count> carry how many merged reads share that exact sequence.
"""

import argparse
import csv
import gzip
import re
import statistics
from collections import Counter, defaultdict

import pysam

COUNT = re.compile(r"^u\d+_(\d+)$")


def reads_in(name):
    match = COUNT.match(name)
    return int(match.group(1)) if match else 1


def species_names(panel):
    names = {}
    with open(panel) as fh:
        for line in fh:
            if line.startswith(">"):
                accession, *words = line[1:].split()
                names[accession] = " ".join(words[:2])
    return names


def load_hits(bams):
    """read -> [(accession, nm, mapq, start, end, strand, length)], plus genome lengths."""
    hits, lengths = defaultdict(list), {}
    for path in bams:
        with pysam.AlignmentFile(path) as bam:
            lengths.update(zip(bam.references, bam.lengths))
            for read in bam:
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                hits[read.query_name].append((
                    read.reference_name, read.get_tag("NM"), read.mapping_quality,
                    read.reference_start + 1, read.reference_end,
                    "-" if read.is_reverse else "+", read.query_length,
                ))
    return hits, lengths


def depth_profile(intervals, length):
    change = [0] * (length + 1)
    for start, end in intervals:
        change[start - 1] += 1
        change[end] -= 1
    depth, running = [], 0
    for step in change[:-1]:
        running += step
        depth.append(running)
    return depth


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--margin", type=int, required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--assigned", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--ambiguous", required=True)
    args = parser.parse_args()

    species = species_names(args.panel)
    hits, lengths = load_hits(args.bams)
    tally = defaultdict(Counter)
    ambiguous = defaultdict(Counter)
    assigned = []
    for name, rows in hits.items():
        reads = reads_in(name)
        for accession, _, mapq, *_ in rows:
            if mapq >= args.min_mapq:
                tally[accession]["single_reads"] += reads
                tally[accession]["single_distinct"] += 1
        best = min(row[1] for row in rows)
        tied = [row for row in rows if row[1] < best + args.margin]
        if len(tied) > 1:
            key = " + ".join(sorted({species[row[0]] for row in tied}))
            ambiguous[key]["reads"] += reads
            ambiguous[key]["distinct"] += 1
            for row in tied:
                tally[row[0]]["ambiguous_distinct"] += 1
        elif tied[0][2] < args.min_mapq:
            tally[tied[0][0]]["low_mapq_distinct"] += 1
        else:
            assigned.append((name, reads, *tied[0]))

    # Position duplicates: keep the copy with the fewest edits, then the lowest name.
    assigned.sort(key=lambda row: (row[2], row[5], row[6], row[7], row[3], row[0]))
    seen, table, kept = set(), [], defaultdict(list)
    for name, reads, accession, nm, mapq, start, end, strand, length in assigned:
        duplicate = (accession, start, end, strand) in seen
        seen.add((accession, start, end, strand))
        table.append((name, accession, species[accession], start, end, strand, nm, mapq,
                      length, reads, int(duplicate)))
        tally[accession]["assigned_reads"] += reads
        tally[accession]["assigned_distinct"] += 1
        if not duplicate:
            kept[accession].append((start, end, nm, length))

    with gzip.open(args.assigned, "wt", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["read", "accession", "species", "start", "end", "strand", "nm", "mapq",
                         "length", "reads", "duplicate"])
        writer.writerows(table)

    header = [
        "library", "accession", "species", "genome_length", "single_reads", "single_distinct",
        "ambiguous_distinct", "low_mapq_distinct", "assigned_reads", "assigned_distinct",
        "dedup_reads", "breadth", "mean_depth", "top5pct_share", "mean_nm", "nm0_fraction",
        "median_length",
    ]
    with open(args.metrics, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for accession in sorted(lengths, key=lambda a: species[a]):
            length, unique, counts = lengths[accession], kept[accession], tally[accession]
            row = [args.library, accession, species[accession], length,
                   *(counts[k] for k in ("single_reads", "single_distinct", "ambiguous_distinct",
                                         "low_mapq_distinct", "assigned_reads",
                                         "assigned_distinct")),
                   len(unique)]
            if unique:
                depth = depth_profile([(s, e) for s, e, _, _ in unique], length)
                busiest = sorted(depth, reverse=True)[: max(1, length // 20)]
                row += [
                    round(sum(d > 0 for d in depth) / length, 4),
                    round(sum(depth) / length, 4),
                    round(sum(busiest) / sum(depth), 3),
                    round(statistics.mean(nm for _, _, nm, _ in unique), 3),
                    round(sum(nm == 0 for _, _, nm, _ in unique) / len(unique), 3),
                    statistics.median(n for _, _, _, n in unique),
                ]
            else:
                row += [0, 0, "", "", "", ""]
            writer.writerow(row)

    with open(args.ambiguous, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["library", "tied_species", "reads", "distinct"])
        for key, counts in sorted(ambiguous.items(), key=lambda kv: -kv[1]["distinct"]):
            writer.writerow([args.library, key, counts["reads"], counts["distinct"]])

    print(f"{args.library}: {len(hits)} distinct reads mapped to at least one genome, "
          f"{len(table)} assigned, {sum(c['distinct'] for c in ambiguous.values())} ambiguous, "
          f"{sum(c['low_mapq_distinct'] for c in tally.values())} below MAPQ {args.min_mapq}")


if __name__ == "__main__":
    main()
