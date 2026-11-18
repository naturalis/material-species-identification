#!/usr/bin/env python
"""Step 4: join the per-library panel metrics and ambiguity bins, with each library's role.

`breadth_meets_gate2` is evidence for Gate 2 only. Gate 2 also asks that depth is not
confined to conserved blocks (`top5pct_share`), and it stays a manual decision.
"""

import argparse
import csv


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--metrics", nargs="+", required=True)
    parser.add_argument("--ambiguous", nargs="+", required=True)
    parser.add_argument("--min-breadth", type=float, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ambiguous-out", required=True)
    args = parser.parse_args()

    with open(args.sheet, newline="") as fh:
        sheet = {row["library"]: row for row in csv.DictReader(fh, delimiter="\t")}

    metrics = read_rows(args.metrics)
    for row in metrics:
        info = sheet[row["library"]]
        row["id"], row["role"] = info["id"], info["role"]
        row["breadth_meets_gate2"] = float(row["breadth"]) >= args.min_breadth
    metrics.sort(key=lambda r: (r["role"] != "sample", r["id"], -int(r["dedup_reads"]),
                                r["species"]))
    columns = ["id", "role", *[c for c in metrics[0] if c not in ("id", "role")]]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(metrics)

    ambiguous = read_rows(args.ambiguous)
    with open(args.ambiguous_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "role", "tied_species", "reads", "distinct"])
        for row in ambiguous:
            info = sheet[row["library"]]
            writer.writerow([info["id"], info["role"], row["tied_species"], row["reads"],
                             row["distinct"]])


if __name__ == "__main__":
    main()
