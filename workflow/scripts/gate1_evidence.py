#!/usr/bin/env python
"""Step 3: put the Gate 1 evidence in one table. It reports; the gate is a manual decision.

Both rules are shown per sample. `plan_rule` is §6's starting rule on all non-human vertebrate
reads. `proposed_rule` works on strict distinct reads within one family, as a rate per
million distinct reads, because the libraries differ in duplication.
"""

import argparse
import csv


def read_table(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--profile", required=True)
    parser.add_argument("--blank", required=True)
    parser.add_argument("--min-reads", type=int, required=True)
    parser.add_argument("--fold", type=float, required=True)
    parser.add_argument("--zero-limit", type=float, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = read_table(args.profile)
    libraries = {row["id"]: row for row in rows if row["level"] == "library"}
    blank = libraries[args.blank]
    blank_million = int(blank["distinct"]) / 1e6
    blank_reads_million = int(blank["reads"]) / 1e6

    def find(library, level, name):
        for row in rows:
            if (row["id"], row["level"], row["name"]) == (library, level, name):
                return row
        return {"reads": 0, "distinct": 0, "strict_reads": 0, "strict_distinct": 0}

    header = [
        "id", "role", "rule", "name", "reads", "distinct", "strict_reads", "strict_distinct",
        "rate_per_M", "blank_rate_per_M", "required_rate_per_M", "passes",
    ]
    table = []
    for library, info in libraries.items():
        if info["role"] != "sample":
            continue
        million = int(info["distinct"]) / 1e6
        reads_million = int(info["reads"]) / 1e6

        plan = find(library, "category", "nonhuman_vertebrate")
        rate = int(plan["reads"]) / reads_million
        blank_rate = int(find(args.blank, "category", "nonhuman_vertebrate")["reads"]) / blank_reads_million
        required = args.fold * blank_rate
        table.append([
            library, info["role"], "plan_rule", "all non-human vertebrate", plan["reads"],
            plan["distinct"], plan["strict_reads"], plan["strict_distinct"], f"{rate:.2f}",
            f"{blank_rate:.2f}", f"{required:.2f}",
            int(plan["reads"]) > args.min_reads and rate > required,
        ])

        reference = find(library, "category", "human")
        rate = int(reference["strict_distinct"]) / million
        blank_rate = int(find(args.blank, "category", "human")["strict_distinct"]) / blank_million
        table.append([
            library, info["role"], "reference_only", "human", reference["reads"],
            reference["distinct"], reference["strict_reads"], reference["strict_distinct"],
            f"{rate:.2f}", f"{blank_rate:.2f}", "", "",
        ])

        families = [row for row in rows if row["id"] == library and row["level"] == "family"]
        families.sort(key=lambda row: -int(row["strict_distinct"]))
        for row in families:
            found = int(find(args.blank, "family", row["name"])["strict_distinct"])
            blank_rate = found / blank_million
            # An empty blank is not a rate of zero: use the upper limit of a zero count.
            required = args.fold * max(found, args.zero_limit) / blank_million
            rate = int(row["strict_distinct"]) / million
            table.append([
                library, info["role"], "proposed_rule", row["name"], row["reads"],
                row["distinct"], row["strict_reads"], row["strict_distinct"], f"{rate:.2f}",
                f"{blank_rate:.2f}", f"{required:.2f}",
                int(row["strict_distinct"]) >= args.min_reads and rate >= required,
            ])

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(table)


if __name__ == "__main__":
    main()
