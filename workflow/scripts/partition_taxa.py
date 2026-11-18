#!/usr/bin/env python
"""Step 6: write the deduplicated reads step 4 assigned to each target taxon, one FASTA each.

Only reads assigned to one genome (margin and MAPQ of step 4) and kept after position
deduplication are written, so every molecule enters the consensus once. Names keep step 3's
u<n>_<count> form. Sequences come from step 4's mapped FASTA, in the read's own orientation.
"""

import argparse
import csv
import gzip


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--assigned", required=True)
    parser.add_argument("--mapped", required=True)
    parser.add_argument("--taxa", nargs="+", required=True)
    parser.add_argument("--out", nargs="+", required=True, help="one FASTA per taxon, same order")
    args = parser.parse_args()
    assert len(args.taxa) == len(args.out), "one output per taxon"

    taxon_of = {}
    with gzip.open(args.assigned, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["duplicate"] == "0" and row["species"] in args.taxa:
                taxon_of[row["read"]] = row["species"]

    handles = {taxon: gzip.open(path, "wt") for taxon, path in zip(args.taxa, args.out)}
    written = dict.fromkeys(args.taxa, 0)
    with gzip.open(args.mapped, "rt") as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                continue
            taxon = taxon_of.get(name)
            if taxon is not None:
                handles[taxon].write(f">{name}\n{line}")
                written[taxon] += 1
    for handle in handles.values():
        handle.close()
    for taxon in args.taxa:
        print(f"{taxon}: {written[taxon]} reads")
    missing = len(taxon_of) - sum(written.values())
    assert missing == 0, f"{missing} assigned reads not found in {args.mapped}"


if __name__ == "__main__":
    main()
