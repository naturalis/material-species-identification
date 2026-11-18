#!/usr/bin/env python
"""Step 6: the barcode gene from route 3 (mapping consensus) and route 1 (MitoGeneExtractor).

Route 3 is cut from the reference-coordinate consensus at the gene's annotated location, at
the masking depth and, labelled as such, at depth 1. Route 1 is MitoGeneExtractor's
consensus over the protein, which starts at the start codon and stops before the stop codon,
so both are compared position by position over the coding length.
"""

import argparse
import csv
import itertools

BASES = set("ACGT")


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip()
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def only(path):
    records = read_fasta(path)
    assert len(records) == 1, f"{path}: {len(records)} records"
    return next(iter(records.values()))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--taxon", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--cds", required=True)
    parser.add_argument("--refcoords", required=True)
    parser.add_argument("--low-depth", required=True)
    parser.add_argument("--mge", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args()

    with open(args.cds, newline="") as fh:
        gene = next(r for r in csv.DictReader(fh, delimiter="\t") if r["accession"] == args.accession)
    assert gene["strand"] == "+", "minus-strand barcode genes are not handled"
    start, coding = int(gene["start"]) - 1, 3 * int(gene["protein_length"])

    routes = {
        "route3_masked": only(args.refcoords)[start:start + coding],
        "route3_depth1": only(args.low_depth)[start:start + coding],
        "route1_mge": only(args.mge).replace("-", "N"),
    }
    for name, sequence in routes.items():
        assert len(sequence) == coding, f"{name}: {len(sequence)} bp, expected {coding}"

    with open(args.fasta, "w") as fh:
        for name, sequence in routes.items():
            fh.write(f">{args.taxon}|{name} {gene['gene']} {args.accession}:{gene['start']}+{coding}\n"
                     f"{sequence}\n")

    rows = []
    for name, sequence in routes.items():
        rows.append({"comparison": name, "length": coding,
                     "called": sum(b in BASES for b in sequence),
                     "deletions": sequence.count("*"), "both_called": "", "differ": "",
                     "differ_at": ""})
    for (a, sa), (b, sb) in itertools.combinations(routes.items(), 2):
        both = [i for i in range(coding) if sa[i] in BASES and sb[i] in BASES]
        differ = [i for i in both if sa[i] != sb[i]]
        rows.append({"comparison": f"{a} vs {b}", "length": coding, "called": "",
                     "deletions": "", "both_called": len(both), "differ": len(differ),
                     "differ_at": ",".join(f"{i + 1}{sa[i]}>{sb[i]}" for i in differ[:50])})
    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
