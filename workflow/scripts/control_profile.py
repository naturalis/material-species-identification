#!/usr/bin/env python
"""Read-level profile of one control's merged fragments, from sequences on stdin.

Needs no reference. Taxonomic composition is the rule control_composition, which shares
the step 3 classifier.
"""
import argparse
import csv
import sys
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--id", required=True)
    parser.add_argument("--top", type=int, required=True)
    parser.add_argument("--low-complexity-fraction", type=float, required=True)
    parser.add_argument("--adapter-motifs", nargs="+", required=True)
    parser.add_argument("--out-profile", required=True)
    parser.add_argument("--out-top", required=True)
    args = parser.parse_args()

    counts = Counter(line.strip() for line in sys.stdin)
    counts.pop("", None)
    reads = sum(counts.values())
    gc_bases = bases = low_complexity = adapter = 0
    for seq, n in counts.items():
        length = len(seq)
        bases += n * length
        gc_bases += n * (seq.count("G") + seq.count("C"))
        if max(seq.count(b) for b in "ACGT") >= args.low_complexity_fraction * length:
            low_complexity += n
        if any(motif in seq for motif in args.adapter_motifs):
            adapter += n

    profile = {
        "id": args.id,
        "merged_reads": reads,
        "distinct_fragments": len(counts),
        "duplicate_fraction": round(1 - len(counts) / reads, 4) if reads else "",
        "gc_fraction": round(gc_bases / bases, 4) if bases else "",
        "low_complexity_fraction": round(low_complexity / reads, 4) if reads else "",
        "adapter_insert_fraction": round(adapter / reads, 4) if reads else "",
    }
    with open(args.out_profile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(profile), delimiter="\t")
        writer.writeheader()
        writer.writerow(profile)

    with open(args.out_top, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "rank", "reads", "fraction", "length", "sequence"])
        for rank, (seq, n) in enumerate(counts.most_common(args.top), 1):
            writer.writerow([args.id, rank, n, round(n / reads, 6), len(seq), seq])


if __name__ == "__main__":
    main()
