#!/usr/bin/env python
"""Combine the per-resource manifests, adding a row for each BLAST database file."""
import argparse
import csv
from pathlib import Path

from resource_manifest import COLUMNS, row, utc_now, write_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", nargs="+", required=True)
    parser.add_argument("--blastdb", required=True, help="database prefix")
    parser.add_argument("--blast-version", required=True, help="file with makeblastdb -version")
    parser.add_argument("--records", required=True, help="the database FASTA's MANIFEST row key")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    for part in args.parts:
        with open(part, newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    records = next(r["records"] for r in rows if r["file"] == args.records)
    tool = Path(args.blast_version).read_text().splitlines()[0]
    prefix = Path(args.blastdb)
    for path in sorted(prefix.parent.glob(f"{prefix.name}.*")):
        rows.append(row(path, f"derived: built from {args.records}", utc_now(), tool, records))
    write_rows(args.out, [{c: r[c] for c in COLUMNS} for r in rows])


if __name__ == "__main__":
    main()
