#!/usr/bin/env python
"""Download one RefSeq mitochondrion release, reading the file list from the directory.

Writes the concatenated, decompressed FASTA and an accession-to-taxonomy table from the
GenBank flatfiles, which carry the organism field that FASTA titles do not reliably give.
"""
import argparse
import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

from resource_manifest import (
    TOOL,
    count_fasta_records,
    count_genbank_records,
    download,
    fetch,
    row,
    utc_now,
    write_rows,
)

FNA = re.compile(r'href="(mitochondrion\.\d+\.\d+\.genomic\.fna\.gz)"')
GBFF = re.compile(r'href="(mitochondrion\.\d+\.genomic\.gbff\.gz)"')


def parse_gbff(paths):
    """Yield (accession, organism, taxid, lineage) for each flatfile record."""
    for path in paths:
        with gzip.open(path, "rt") as fh:
            accession = organism = taxid = None
            lineage, in_lineage = [], False
            for line in fh:
                if line.startswith("VERSION"):
                    accession = line.split()[1]
                elif line.startswith("  ORGANISM"):
                    organism = line[10:].strip()
                    in_lineage = True
                    lineage = []
                elif in_lineage:
                    if line.startswith("            "):
                        lineage.append(line.strip())
                    else:
                        in_lineage = False
                elif taxid is None and '/db_xref="taxon:' in line:
                    taxid = line.split("taxon:")[1].split('"')[0]
                elif line.startswith("//"):
                    yield accession, organism, taxid, " ".join(lineage).rstrip(".")
                    accession = organism = taxid = None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--url", required=True)
    parser.add_argument("--release-url", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--release-file", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    listing = fetch(args.url).decode()
    fna_names, gbff_names = sorted(set(FNA.findall(listing))), sorted(set(GBFF.findall(listing)))
    if not fna_names or not gbff_names:
        sys.exit(f"{args.url}: no FASTA or GenBank files in the listing")
    release = fetch(args.release_url).decode().strip()
    print(f"release {release}: {fna_names + gbff_names}", file=sys.stderr)

    rows, fna_paths, gbff_paths = [], [], []
    for names, paths in ((fna_names, fna_paths), (gbff_names, gbff_paths)):
        for name in names:
            path = outdir / name
            downloaded = download(args.url + name, path)
            with gzip.open(path, "rb") as fh:  # the gzip CRC catches a truncated download
                while fh.read(1 << 20):
                    pass
            counter = count_fasta_records if ".fna." in name else count_genbank_records
            records = counter(path)
            rows.append(row(path, args.url + name, downloaded, TOOL, records))
            paths.append(path)

    if fetch(args.release_url).decode().strip() != release:
        sys.exit("RefSeq release changed during the download; rerun.")
    Path(args.release_file).write_text(f"{release}\n")
    rows.append(row(args.release_file, args.release_url, utc_now(), TOOL, 1))

    with open(args.fasta, "wb") as out:
        for path in fna_paths:
            with gzip.open(path, "rb") as fh:
                shutil.copyfileobj(fh, out)
    records = count_fasta_records(args.fasta)
    rows.append(
        row(args.fasta, "derived: concatenated *.genomic.fna.gz", utc_now(), TOOL, records)
    )

    taxonomy = list(parse_gbff(gbff_paths))
    with open(args.taxonomy, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "organism", "taxid", "lineage"])
        writer.writerows(taxonomy)
    rows.append(
        row(args.taxonomy, "derived: parsed *.genomic.gbff.gz", utc_now(), TOOL, len(taxonomy))
    )

    with open(args.fasta) as fh:
        fasta_ids = {line[1:].split()[0] for line in fh if line.startswith(">")}
    gbff_ids = {t[0] for t in taxonomy}
    if fasta_ids != gbff_ids:
        sys.exit(
            f"FASTA and GenBank disagree: {len(fasta_ids - gbff_ids)} accessions only in FASTA, "
            f"{len(gbff_ids - fasta_ids)} only in GenBank"
        )
    write_rows(args.manifest, rows)
    print(f"release {release}: {records} records", file=sys.stderr)


if __name__ == "__main__":
    main()
