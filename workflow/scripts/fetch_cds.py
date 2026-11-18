#!/usr/bin/env python
"""Step 6 references: the barcode gene of each target reference, as protein and coordinates.

The protein is MitoGeneExtractor's reference (route 1). The coordinates cut the same gene out
of the mapping consensus (route 3), so the two routes can be compared base by base.
"""

import argparse
import csv
import re
import time
from pathlib import Path

from resource_manifest import TOOL, fetch, row, utc_now, write_rows

LOCATION = re.compile(r"\[location=(complement\()?(\d+)\.\.(\d+)\)?\]")


def efetch_url(base, accession, rettype):
    return (f"{base}?db=nuccore&id={accession}&rettype={rettype}&retmode=text"
            "&tool=material-species-identification")


def records(text):
    for block in text.split(">")[1:]:
        header, *lines = block.splitlines()
        yield header, "".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--accessions", nargs="+", required=True)
    parser.add_argument("--gene-names", nargs="+", required=True, help="first is the label")
    parser.add_argument("--efetch-url", required=True)
    parser.add_argument("--delay", type=float, required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    aliases = {name.upper() for name in args.gene_names}
    label = args.gene_names[0]
    table, manifest = [], []
    for accession in args.accessions:
        url = efetch_url(args.efetch_url, accession, "fasta_cds_aa")
        text = fetch(url).decode()
        fetched = utc_now()
        time.sleep(args.delay)
        hits = [(h, s) for h, s in records(text)
                if (m := re.search(r"\[gene=([^\]]+)\]", h)) and m.group(1).upper() in aliases]
        assert len(hits) == 1, f"{accession}: {len(hits)} CDS named {sorted(aliases)}"
        header, protein = hits[0]
        location = LOCATION.search(header)
        assert location, f"{accession}: no simple location in {header!r}"
        path = outdir / f"{accession}.{label}.faa"
        path.write_text(f">{accession}_{label}\n{protein.rstrip('*')}\n")
        manifest.append(row(path, url, fetched, TOOL, 1))
        table.append({
            "accession": accession, "gene": label,
            "start": int(location.group(2)), "end": int(location.group(3)),
            "strand": "-" if location.group(1) else "+",
            "protein_length": len(protein.rstrip("*")),
            "protein_id": (re.search(r"\[protein_id=([^\]]+)\]", header) or [None, ""])[1],
        })

    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(table[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    write_rows(args.manifest, manifest)


if __name__ == "__main__":
    main()
