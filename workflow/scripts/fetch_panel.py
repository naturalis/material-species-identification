#!/usr/bin/env python
"""Fetch each panel accession individually from NCBI and verify it before accepting it.

The check is against the record's /organism field, not its title: a title can read right
while the organism is something else. Anything other than an exact match is reported and
kept out of panel.fasta. Nothing is substituted or fuzzily accepted. Where RefSeq holds
only subspecies records, the request names the subspecies in `expected_organism`.
"""
import argparse
import csv
import re
import sys
import time
from pathlib import Path

from resource_manifest import TOOL, count_fasta_records, download, row, utc_now, write_rows

COMPLETE = re.compile(r"complete (genome|sequence|mitochondrial genome)", re.IGNORECASE)
COLUMNS = [
    "species", "accession", "use", "role", "status", "problems", "organism", "taxid", "definition",
    "length", "fasta_length", "topology", "ambiguous_bases", "in_refseq_release",
]


def efetch_url(base, accession, rettype):
    return f"{base}?db=nuccore&id={accession}&rettype={rettype}&retmode=text&tool=material-species-identification"


def parse_genbank(text):
    """Return the fields needed for verification from one GenBank flatfile record."""
    locus = re.search(r"^LOCUS\s+\S+\s+(\d+) bp\s+\S+\s+(\w+)", text, re.MULTILINE)
    definition = re.search(r"^DEFINITION\s+(.+?)\n(?=\S)", text, re.MULTILINE | re.DOTALL)
    version = re.search(r"^VERSION\s+(\S+)", text, re.MULTILINE)
    organism = re.search(r'/organism="([^"]+)"', text)
    taxid = re.search(r'/db_xref="taxon:(\d+)"', text)
    organelle = re.search(r'/organelle="([^"]+)"', text)
    return {
        "length": int(locus.group(1)) if locus else None,
        "topology": locus.group(2) if locus else "",
        "definition": " ".join(definition.group(1).split()) if definition else "",
        "version": version.group(1) if version else "",
        "organism": organism.group(1) if organism else "",
        "taxid": taxid.group(1) if taxid else "",
        "organelle": organelle.group(1) if organelle else "",
    }


def fasta_sequence(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    return lines[0], "".join(lines[1:]).upper()


def verify(request, record, header, sequence):
    problems = []
    expected = request.get("expected_organism") or request["species"]
    if record["organism"] != expected:
        problems.append(f"organism field is {record['organism']!r}")
    if record["version"] != request["accession"]:
        problems.append(f"returned version {record['version']!r}")
    if record["organelle"] != "mitochondrion":
        problems.append(f"organelle is {record['organelle']!r}")
    if not COMPLETE.search(record["definition"]):
        problems.append("definition does not say the record is complete")
    if not header.startswith(f">{request['accession']} "):
        problems.append(f"FASTA header is {header[:40]!r}")
    if record["length"] != len(sequence):
        problems.append(f"GenBank length {record['length']}, FASTA length {len(sequence)}")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--request", required=True)
    parser.add_argument("--efetch-url", required=True)
    parser.add_argument("--delay", type=float, required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--refseq-taxonomy", required=True)
    parser.add_argument("--panel-fasta", required=True)
    parser.add_argument("--panel-tsv", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    (outdir / "fasta").mkdir(parents=True, exist_ok=True)
    (outdir / "genbank").mkdir(parents=True, exist_ok=True)
    with open(args.refseq_taxonomy) as fh:
        refseq_ids = {r["accession"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(args.request, newline="") as fh:
        requests = list(csv.DictReader(fh))

    table, manifest, accepted = [], [], []
    for request in requests:
        accession = request["accession"]
        stem = f"{request['species'].replace(' ', '_')}.{accession}"
        gb_path = outdir / "genbank" / f"{stem}.gb"
        fasta_path = outdir / "fasta" / f"{stem}.fasta"
        for path, rettype, records in ((gb_path, "gb", 1), (fasta_path, "fasta", 1)):
            url = efetch_url(args.efetch_url, accession, rettype)
            downloaded = download(url, path)
            time.sleep(args.delay)
            manifest.append(row(path, url, downloaded, TOOL, records))

        record = parse_genbank(gb_path.read_text())
        header, sequence = fasta_sequence(fasta_path)
        problems = verify(request, record, header, sequence)
        table.append(
            {
                "species": request["species"],
                "accession": accession,
                "use": request["use"],
                "role": request["role"],
                "status": "MISMATCH" if problems else "verified",
                "problems": "; ".join(problems),
                "organism": record["organism"],
                "taxid": record["taxid"],
                "definition": record["definition"],
                "length": record["length"],
                "fasta_length": len(sequence),
                "topology": record["topology"],
                "ambiguous_bases": sum(1 for base in sequence if base not in "ACGT"),
                "in_refseq_release": accession in refseq_ids,
            }
        )
        if problems:
            print(f"MISMATCH {request['species']} {accession}: {problems}", file=sys.stderr)
        elif request["use"] == "panel":
            accepted.append(fasta_path)

    with open(args.panel_fasta, "w") as out:
        for path in accepted:
            out.write(Path(path).read_text())
    manifest.append(
        row(
            args.panel_fasta, "derived: concatenated verified records", utc_now(), TOOL,
            count_fasta_records(args.panel_fasta),
        )
    )
    with open(args.panel_tsv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    manifest.append(row(args.panel_tsv, "derived: verification table", utc_now(), TOOL, len(table)))
    write_rows(args.manifest, manifest)
    print(f"{len(accepted)} of {len(requests)} requests went into panel.fasta", file=sys.stderr)


if __name__ == "__main__":
    main()
