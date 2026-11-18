#!/usr/bin/env python
"""Step 3: join the per-library tables into one long table, with library totals."""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--results", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.sheet, newline="") as fh:
        sheet = list(csv.DictReader(fh, delimiter="\t"))
    columns = ["id", "role", "level", "name", "reads", "distinct", "strict_reads", "strict_distinct"]
    with open(args.out, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in sheet:
            prefix = f"{args.results}/{row['library']}/{row['library']}"
            with open(f"{prefix}.unique.tsv", newline="") as fh:
                totals = next(csv.DictReader(fh, delimiter="\t"))
            writer.writerow([row["id"], row["role"], "library", row["library"],
                             totals["reads"], totals["distinct"], "", ""])
            with open(f"{prefix}.lca.tsv", newline="") as fh:
                for entry in csv.DictReader(fh, delimiter="\t"):
                    writer.writerow([row["id"], row["role"], *(entry[c] for c in columns[2:])])


if __name__ == "__main__":
    main()
