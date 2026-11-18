#!/usr/bin/env python
"""Step 6: place a taxon's reads among complete mitogenomes of its species and family.

1. Lineages are found from the sequences, not the labels: pairwise p-distances over the
   MAFFT alignment (columns where at least 90% of genomes have a base, so length-variable
   repeats do not count), average-linkage clustering, cut at --cut.
2. Every read is mapped to every genome on its own, with the step 4 bwa settings.
3. A genome's score is its total edit distance over the reads that map to nearly all genomes;
   a read that fails to map to one genome counts one edit more than its worst mapping.
4. A read supports a cluster when that cluster's best genome fits it with fewer edits than
   any genome of any other cluster. Most reads fit every lineage equally; only the
   difference counts.
Origin metadata is shown as given in GenBank. It is often missing or a breed name.
"""

import argparse
import csv
import gzip
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ORIGIN_FIELDS = ["country", "breed", "strain", "isolate", "sub_species", "note"]


def read_fasta(path, opener=open):
    records, name = {}, None
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def distances(alignment, names):
    matrix = np.frombuffer("".join(alignment[n] for n in names).encode(), dtype=np.uint8)
    matrix = matrix.reshape(len(names), -1)
    base = np.isin(matrix, np.frombuffer(b"ACGT", dtype=np.uint8))
    keep = base.mean(axis=0) >= 0.9
    matrix, base = matrix[:, keep], base[:, keep]
    out = np.zeros((len(names), len(names)))
    for i in range(len(names)):
        valid = base & base[i]
        out[i] = ((matrix != matrix[i]) & valid).sum(axis=1) / np.maximum(valid.sum(axis=1), 1)
    return (out + out.T) / 2, int(keep.sum())


def map_all(reads, genomes, work, aln_params, min_mapq, threads):
    """read -> {genome: nm}, for mappings at min_mapq or above."""
    work.mkdir(parents=True, exist_ok=True)
    hits = defaultdict(dict)
    for accession, sequence in genomes.items():
        reference = work / "genome.fasta"
        reference.write_text(f">{accession}\n{sequence}\n")
        subprocess.run(["bwa", "index", str(reference)], check=True, capture_output=True)
        with open(work / "reads.sai", "wb") as sai:
            subprocess.run(["bwa", "aln", "-t", str(threads), *aln_params.split(), str(reference),
                            str(reads)], check=True, stdout=sai, stderr=subprocess.DEVNULL)
        with open(work / "reads.sam", "wb") as sam:
            subprocess.run(["bwa", "samse", str(reference), str(work / "reads.sai"), str(reads)],
                           check=True, stdout=sam, stderr=subprocess.DEVNULL)
        with pysam.AlignmentFile(str(work / "reads.sam")) as bam:
            for read in bam:
                if not read.is_unmapped and read.mapping_quality >= min_mapq:
                    hits[read.query_name][accession] = read.get_tag("NM")
    return hits


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reads", nargs="+", required=True)
    parser.add_argument("--genomes", required=True)
    parser.add_argument("--table", required=True)
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--reference", required=True, help="the panel genome of this species")
    parser.add_argument("--cut", type=float, required=True)
    parser.add_argument("--aln-params", required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--work", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--genome-scores", required=True)
    parser.add_argument("--read-support", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    genomes = read_fasta(args.genomes)
    with open(args.table, newline="") as fh:
        meta = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for record in meta.values():
        parts = [record.get(f, "") for f in ORIGIN_FIELDS]
        record["origin"] = "; ".join(p for p in parts if p) or record.get("definition", "")
    names = list(genomes)
    assert args.reference in genomes, f"{args.reference} is not in {args.genomes}"

    # 1. Lineages from the sequences.
    alignment = read_fasta(args.alignment)
    matrix, columns = distances(alignment, names)
    tree = linkage(squareform(matrix, checks=False), method="average")
    raw = fcluster(tree, t=args.cut, criterion="distance")
    order = [c for c, _ in Counter(raw).most_common()]
    cluster = {name: f"L{order.index(c) + 1}" for name, c in zip(names, raw)}
    ref_index = names.index(args.reference)
    print(f"{len(names)} genomes, {columns} alignment columns, {len(order)} clusters at {args.cut}")

    # 2. Reads against every genome.
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    reads_fasta = work / "reads.fasta"
    lengths = {}
    with open(reads_fasta, "w") as out:
        for path in args.reads:
            for name, sequence in read_fasta(path, gzip.open).items():
                out.write(f">{name}\n{sequence}\n")
                lengths[name] = len(sequence)
    hits = map_all(reads_fasta, genomes, work, args.aln_params, args.min_mapq, args.threads)
    print(f"{len(lengths)} reads, {len(hits)} map to at least one genome")

    # 3. Genome scores over reads that map nearly everywhere.
    common = [r for r in hits if len(hits[r]) >= 0.95 * len(names)]
    score = {g: sum(hits[r].get(g, max(hits[r].values()) + 1) for r in common) for g in names}
    nm0 = {g: sum(hits[r].get(g) == 0 for r in hits) for g in names}
    mapped = {g: sum(g in hits[r] for r in hits) for g in names}
    ranked = sorted(names, key=lambda g: (score[g], -mapped[g], g))
    with open(args.genome_scores, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["rank", "accession", "cluster", "organism", "origin", "reads_mapped",
                         "common_reads", "edits_on_common", "nm0_reads", "distance_to_reference"])
        for rank, g in enumerate(ranked, 1):
            writer.writerow([rank, g, cluster[g], meta[g]["organism"], meta[g]["origin"],
                             mapped[g], len(common), score[g], nm0[g],
                             round(matrix[ref_index, names.index(g)], 4)])

    # 4. Per-read support for clusters.
    clusters = sorted(set(cluster.values()), key=lambda c: int(c[1:]))
    support, support2 = Counter(), Counter()
    with open(args.read_support, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["read", "length", "genomes_mapped", "best_nm", "best_clusters", "margin",
                         *clusters])
        for read in sorted(hits):
            best = {c: min((nm for g, nm in hits[read].items() if cluster[g] == c), default=None)
                    for c in clusters}
            values = sorted(v for v in best.values() if v is not None)
            top = values[0]
            winners = [c for c in clusters if best[c] == top]
            margin = (values[1] - top) if len(values) > 1 else ""
            if len(winners) == 1 and margin != "":
                support[winners[0]] += 1
                if margin >= 2:
                    support2[winners[0]] += 1
            writer.writerow([read, lengths[read], len(hits[read]), top, "+".join(winners), margin,
                             *("" if best[c] is None else best[c] for c in clusters)])

    with open(args.clusters, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "cluster", "organism", "origin", "distance_to_reference"])
        for g in sorted(names, key=lambda g: (int(cluster[g][1:]), g)):
            writer.writerow([g, cluster[g], meta[g]["organism"], meta[g]["origin"],
                             round(matrix[ref_index, names.index(g)], 4)])

    with open(args.summary, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["cluster", "genomes", "organisms", "countries", "contains_reference",
                         "distance_to_reference", "reads_favouring", "reads_favouring_margin2",
                         "best_genome", "best_genome_rank", "best_genome_edits"])
        for c in clusters:
            members = [g for g in names if cluster[g] == c]
            organisms = Counter(meta[g]["organism"] for g in members)
            countries = Counter((meta[g].get("country") or "not given").split(":")[0]
                                for g in members)
            best = min(members, key=ranked.index)
            writer.writerow([
                c, len(members),
                "; ".join(f"{o} {n}" for o, n in organisms.most_common()),
                "; ".join(f"{o} {n}" for o, n in countries.most_common(6)),
                args.reference in members,
                round(float(np.median([matrix[ref_index, names.index(g)] for g in members])), 4),
                support[c], support2[c], best, ranked.index(best) + 1, score[best],
            ])
    for c in clusters:
        print(c, support[c], support2[c])


if __name__ == "__main__":
    main()
