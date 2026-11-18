#!/usr/bin/env python
"""Step 6: the structural quality of a recovered barcode -- length, ambiguities, stop codons.

These are the criteria to apply by hand to any hit taken from an uncurated reference library,
an uncurated BOLD release, and the same three that `barcode_validator` applies in its
structural mode.

`barcode_validator` 2.0.8 itself is NOT what runs here. Its structural mode requires a
taxonomy dump and a resolver for every record, used only to look up which translation table
the record's taxon needs. Our sequences are cut at the annotated CDS of a named mitogenome,
so the reading frame starts at base 1 and the table is known (vertebrate mitochondrial, 2).
Supplying a taxonomy to rediscover that is plumbing, not evidence. The criteria and their
thresholds are theirs; the implementation is ours, and it is stated as ours wherever reported.

Reported against two bars, because they differ:
  - BOLD's formal barcode standard, 487 bp (section 3b).
  - `barcode_validator`'s COI-5P defaults: 500 bp, 0 ambiguous bases, 0 stop codons.
"""

import argparse
import csv

from Bio.Seq import Seq

BASES = set("ACGT")


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = []
        else:
            records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fasta", nargs="+", required=True)
    parser.add_argument("--table", type=int, required=True, help="NCBI translation table")
    parser.add_argument("--bold-min-length", type=int, required=True)
    parser.add_argument("--validator-min-length", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    for path in args.fasta:
        for header, sequence in read_fasta(path).items():
            # Deletions relative to the reference are '*'; they are not bases and not gaps.
            sequence = sequence.replace("*", "")
            called = [b for b in sequence if b in BASES]
            # The frame starts at base 1: the sequence is cut at the annotated CDS start.
            codons = [sequence[i:i + 3] for i in range(0, len(sequence) - len(sequence) % 3, 3)]
            translated = "".join(
                str(Seq(c).translate(table=args.table)) if set(c) <= BASES else "X"
                for c in codons
            )
            internal = translated[:-1] if translated.endswith("*") else translated
            stops = internal.count("*")
            taxon, route = header.split()[0].split("|")
            rows.append({
                "taxon": taxon,
                "route": route,
                "length": len(sequence),
                "called_bases": len(called),
                "ambiguous_bases": len(sequence) - len(called),
                "codons_translated": sum(1 for c in codons if set(c) <= BASES),
                "internal_stop_codons": stops,
                "terminal_stop": translated.endswith("*"),
                "meets_bold_487": len(called) >= args.bold_min_length,
                "meets_validator_defaults": (len(called) >= args.validator_min_length
                                             and len(sequence) == len(called) and stops == 0),
            })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['taxon']} {row['route']}: {row['called_bases']} called, "
              f"{row['ambiguous_bases']} ambiguous, {row['internal_stop_codons']} internal stops, "
              f"BOLD 487 {'pass' if row['meets_bold_487'] else 'FAIL'}")


if __name__ == "__main__":
    main()
