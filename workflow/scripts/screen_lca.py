#!/usr/bin/env python
"""Step 3: summarise one library's blastn hits by category, taxon and vertebrate family.

Query names u<n>_<count> carry how many merged reads share that exact sequence. Other names
count as one read, so the Day 1 output can be used as a regression test. Same best-hit, tie
and strict rules as analysis/day1_provisional.py, which is the reference for this logic.
"""

import argparse
import csv
import gzip
import re
from collections import Counter

COUNT = re.compile(r"^u\d+_(\d+)$")
FIELDS = ("reads", "distinct", "strict_reads", "strict_distinct")


def species_name(organism):
    """Binomial, so subspecies and strains collapse into one species count."""
    parts = organism.split()
    return " ".join(parts[:2]) if len(parts) > 1 else organism


def family_of(lineage):
    for taxon in reversed(lineage):
        if taxon.endswith("idae"):
            return taxon
    return None


def load_taxonomy(path):
    taxonomy = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            lineage = [t.strip() for t in row["lineage"].split(";") if t.strip()]
            taxonomy[row["accession"]] = (species_name(row["organism"]), lineage)
    return taxonomy


def best_hits(path):
    """Yield (query, rows at its top bitscore)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    current, top, rows = None, 0.0, []
    with opener(path, "rt") as fh:
        for line in fh:
            qseqid, accession, pident, length, bitscore, qlen = line.rstrip("\n").split("\t")
            bitscore = float(bitscore)
            if qseqid != current:
                if current is not None:
                    yield current, rows
                current, top, rows = qseqid, bitscore, []
            if bitscore >= top:
                top = max(top, bitscore)
                rows.append((accession, float(pident), int(length), int(qlen)))
    if current is not None:
        yield current, rows


def lowest_common_taxon(lineages):
    shared = None
    for lineage in lineages:
        shared = lineage if shared is None else [a for a, b in zip(shared, lineage) if a == b]
    return shared[-1] if shared else "root"


def categorise(species, lineages):
    homo = {s for s in species if s.split()[0] == "Homo"}
    vertebrate = all("Vertebrata" in lineage for lineage in lineages)
    if homo and homo == species:
        return "human"
    if homo:
        return "human_tied_with_other"
    if vertebrate:
        return "nonhuman_vertebrate"
    if any("Vertebrata" in lineage for lineage in lineages):
        return "vertebrate_tied_with_nonvertebrate"
    return "non_vertebrate"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--hits", required=True)
    parser.add_argument("--taxonomy", required=True)
    parser.add_argument("--strict-identity", type=float, required=True)
    parser.add_argument("--strict-coverage", type=float, required=True)
    parser.add_argument("--artefact-records", nargs="*", default=[])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    taxonomy = load_taxonomy(args.taxonomy)
    artefacts = set(args.artefact_records)
    table = {}

    def add(level, name, reads, strict):
        counts = table.setdefault((level, name), Counter())
        counts["reads"] += reads
        counts["distinct"] += 1
        if strict:
            counts["strict_reads"] += reads
            counts["strict_distinct"] += 1

    for query, rows in best_hits(args.hits):
        match = COUNT.match(query)
        reads = int(match.group(1)) if match else 1
        hits = [taxonomy[accession] for accession, *_ in rows]
        species = {name for name, _ in hits}
        lineages = [lineage for _, lineage in hits]
        # Species plus its own name as a leaf, so reads on one species keep that species.
        taxon = lowest_common_taxon([lineage + [name] for name, lineage in hits])
        category = categorise(species, lineages)
        strict = max(
            pident >= args.strict_identity and length >= args.strict_coverage * qlen
            for _, pident, length, qlen in rows
        )
        add("total", "reads_with_hit", reads, strict)
        add("category", category, reads, strict)
        add("lca", taxon, reads, strict)
        for name in species:
            add("species_alone" if len(species) == 1 else "species_tied", name, reads, strict)

        # A flag for the family table only: nothing is removed from the counts above.
        flagged = artefacts & {accession for accession, *_ in rows}
        for accession in sorted(flagged):
            add("artefact_record", accession, reads, strict)
        families = {family_of(lineage) for lineage in lineages}
        if category == "nonhuman_vertebrate" and not flagged and len(families) == 1:
            add("family", families.pop() or "no_family", reads, strict)

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["level", "name", *FIELDS])
        for (level, name), counts in sorted(table.items()):
            writer.writerow([level, name, *(counts[field] for field in FIELDS)])


if __name__ == "__main__":
    main()
