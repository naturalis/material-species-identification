#!/usr/bin/env python
"""Step 5: misincorporation by read position, and fragment lengths, per assigned taxon.

Only reads step 4 assigned to one genome and kept after position deduplication are used,
each once, aligned to the genome they were assigned to. Positions are counted in the read's
own orientation: C->T from the 5' end and G->A from the 3' end, the deamination signal of a
double-stranded library. `other` counts every other mismatch, as the noise floor.
"""

import argparse
import csv
import gzip
from collections import Counter, defaultdict

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bam", required=True)
    parser.add_argument("--assigned", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--positions", type=int, required=True)
    parser.add_argument("--damage", required=True)
    parser.add_argument("--lengths", required=True)
    args = parser.parse_args()

    kept = {}
    lengths = defaultdict(Counter)
    with gzip.open(args.assigned, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["duplicate"] == "0":
                kept[row["read"]] = (row["accession"], row["species"])
                lengths[row["species"]][int(row["length"])] += 1

    genomes = read_fasta(args.panel)
    counts = defaultdict(Counter)
    used, skipped = Counter(), 0
    with pysam.AlignmentFile(args.bam) as bam:
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            target = kept.get(read.query_name)
            if target is None:
                continue
            accession, species = target
            # The combined-panel alignment is used only where it is to the assigned genome.
            if read.reference_name != accession:
                skipped += 1
                continue
            used[species] += 1
            genome, sequence = genomes[accession], read.query_sequence.upper()
            length, reverse = len(sequence), read.is_reverse
            tally = counts[species]
            for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                ref_base, read_base = genome[ref_pos], sequence[query_pos]
                if reverse:
                    ref_base, read_base = ref_base.translate(COMPLEMENT), read_base.translate(COMPLEMENT)
                    from_5 = length - 1 - query_pos
                else:
                    from_5 = query_pos
                from_3 = length - 1 - from_5
                if from_5 < args.positions and ref_base == "C":
                    tally["5p", from_5 + 1, "C"] += 1
                    tally["5p", from_5 + 1, "C>T"] += read_base == "T"
                if from_3 < args.positions and ref_base == "G":
                    tally["3p", from_3 + 1, "G"] += 1
                    tally["3p", from_3 + 1, "G>A"] += read_base == "A"
                if ref_base in "ACGT" and read_base in "ACGT":
                    tally["all", 0, "bases"] += 1
                    if read_base != ref_base and (ref_base, read_base) not in (("C", "T"), ("G", "A")):
                        tally["all", 0, "other"] += 1

    with open(args.damage, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["library", "species", "reads", "end", "position", "reference_bases",
                         "substitutions"])
        for species in sorted(counts):
            tally = counts[species]
            for end, base, change in (("5p", "C", "C>T"), ("3p", "G", "G>A")):
                for position in range(1, args.positions + 1):
                    writer.writerow([args.library, species, used[species], end, position,
                                     tally[end, position, base], tally[end, position, change]])
            writer.writerow([args.library, species, used[species], "other", 0,
                             tally["all", 0, "bases"], tally["all", 0, "other"]])

    with open(args.lengths, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["library", "species", "length", "reads"])
        for species in sorted(lengths):
            for length in sorted(lengths[species]):
                writer.writerow([args.library, species, length, lengths[species][length]])

    print(f"{args.library}: {sum(used.values())} deduplicated reads profiled, {skipped} skipped "
          "because their combined-panel alignment was to another genome")


if __name__ == "__main__":
    main()
