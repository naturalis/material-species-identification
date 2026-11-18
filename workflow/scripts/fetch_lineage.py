#!/usr/bin/env python
"""Step 6 references: complete mitogenomes of one species, to place its reads by lineage.

Fetches every nuccore record matching an Entrez query as GenBank, keeps the fields that say
where the animal came from (country, breed, strain, isolate, subspecies, note), and adds the
records of the species' family already in the local RefSeq mitochondrion release, so a read
set closer to a related species than to any lineage of this one shows it. Only the query is
sent to NCBI; nothing from the project's data leaves the machine.
"""

import argparse
import csv
import gzip
import re
import time
import urllib.parse
from pathlib import Path

from resource_manifest import TOOL, fetch, row, utc_now, write_rows

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH = 100
QUALIFIERS = ["country", "geo_loc_name", "breed", "strain", "isolate", "sub_species", "note"]
COLUMNS = ["accession", "source", "organism", "definition", "length", "country", "breed",
           "strain", "isolate", "sub_species", "note"]


def eutils(endpoint, **params):
    params["tool"] = "material-species-identification"
    return f"{EUTILS}/{endpoint}.fcgi?{urllib.parse.urlencode(params)}"


def qualifier(source, name):
    match = re.search(rf'/{name}="([^"]*)"', source, re.DOTALL)
    return " ".join(match.group(1).split()) if match else ""


def parse_records(text):
    for record in text.split("\n//"):
        version = re.search(r"^VERSION\s+(\S+)", record, re.MULTILINE)
        if not version:
            continue
        definition = re.search(r"^DEFINITION\s+(.+?)\n(?=\S)", record, re.MULTILINE | re.DOTALL)
        source = re.search(r"^     source .*?(?=^     \S)", record, re.MULTILINE | re.DOTALL)
        source = source.group(0) if source else ""
        origin = record.split("\nORIGIN", 1)[1] if "\nORIGIN" in record else ""
        sequence = re.sub(r"[^A-Za-z]", "", origin).upper()
        fields = {name: qualifier(source, name) for name in QUALIFIERS}
        fields["country"] = fields["country"] or fields.pop("geo_loc_name")
        fields.pop("geo_loc_name", None)
        yield {
            "accession": version.group(1),
            "source": "genbank_query",
            "organism": qualifier(source, "organism"),
            "definition": " ".join(definition.group(1).split()) if definition else "",
            "length": len(sequence),
            **fields,
        }, sequence


def read_fasta(path, wanted):
    records, name, keep = {}, None, False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                keep = name in wanted
                if keep:
                    records[name] = []
            elif keep:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--query", required=True)
    parser.add_argument("--family", required=True, help="family name in the RefSeq lineage")
    parser.add_argument("--refseq-fasta", required=True)
    parser.add_argument("--refseq-taxonomy", required=True)
    parser.add_argument("--delay", type=float, required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--genbank", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    search_url = eutils("esearch", db="nuccore", term=args.query, usehistory="y", retmax=0)
    search = fetch(search_url).decode()
    count = int(re.search(r"<Count>(\d+)</Count>", search).group(1))
    webenv = re.search(r"<WebEnv>(\S+?)</WebEnv>", search).group(1)
    query_key = re.search(r"<QueryKey>(\d+)</QueryKey>", search).group(1)
    print(f"{count} records match {args.query!r}")

    rows, sequences = [], {}
    Path(args.genbank).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.genbank, "wt") as gb:
        for start in range(0, count, BATCH):
            url = eutils("efetch", db="nuccore", WebEnv=webenv, query_key=query_key,
                         retstart=start, retmax=BATCH, rettype="gb", retmode="text")
            text = fetch(url).decode()
            gb.write(text)
            for record, sequence in parse_records(text):
                rows.append(record)
                sequences[record["accession"]] = sequence
            time.sleep(args.delay)
    assert len(rows) == count, f"expected {count} records, parsed {len(rows)}"
    fetched = utc_now()

    with open(args.refseq_taxonomy) as fh:
        family = {r["accession"]: r["organism"] for r in csv.DictReader(fh, delimiter="\t")
                  if args.family in [t.strip() for t in r["lineage"].split(";")]}
    for accession, sequence in read_fasta(args.refseq_fasta, family).items():
        if accession in sequences:
            continue
        rows.append({"accession": accession, "source": f"refseq_{args.family}",
                     "organism": family[accession], "length": len(sequence)})
        sequences[accession] = sequence

    with open(args.fasta, "w") as fh:
        for record in rows:
            fh.write(f">{record['accession']} {record['organism']}\n{sequences[record['accession']]}\n")
    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n", restval="")
        writer.writeheader()
        writer.writerows(rows)
    write_rows(args.manifest, [
        row(args.genbank, f"{EUTILS}/esearch+efetch term={args.query}", fetched, TOOL, count),
        row(args.fasta, f"derived: {count} query records + RefSeq {args.family}", utc_now(),
            TOOL, len(rows)),
    ])
    print(f"wrote {len(rows)} genomes ({len(rows) - count} from RefSeq {args.family})")


if __name__ == "__main__":
    main()
