#!/usr/bin/env python
"""Step 6: how complete the mapping consensus is, whether libraries agree, and mixture.

Everything is in reference coordinates. Base counts come from the per-library BAMs (MAPQ
filtered at mapping). A consensus called from one individual leaves only sequencing error at
each site; DNA from more than one individual leaves a minority allele at the sites where they
differ. Two measures of that:

- minor-allele sites: sites at --mixture-depth or more where the second allele has at least
  --minor-fraction of reads. Its count and the median minority share are reported.
- consensus-defining sites: sites where the pooled consensus differs from the reference. The
  share of reads there that do not carry the consensus allele estimates how much DNA comes
  from other individuals. It is a lower bound: another individual sharing the allele is
  invisible.
"""

import argparse
import csv
import gzip
import statistics

import pysam

BASES = "ACGT"


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def base_counts(path, contig, length):
    with pysam.AlignmentFile(path) as bam:
        reads = bam.count(contig)
        a, c, g, t = bam.count_coverage(contig, 0, length, quality_threshold=0,
                                       read_callback="all")
    return reads, list(zip(a, c, g, t))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--bams", nargs="+", required=True)
    parser.add_argument("--libraries", nargs="+", required=True, help="ids, same order as --bams")
    parser.add_argument("--samples", nargs="+", required=True, help="ids pooled as samples")
    parser.add_argument("--pooled", required=True, help="pooled samples consensus, ref coords")
    parser.add_argument("--per-library", required=True, help="per-library consensus, ref coords")
    parser.add_argument("--min-depth", type=int, required=True)
    parser.add_argument("--mixture-depth", type=int, required=True)
    parser.add_argument("--minor-fraction", type=float, required=True)
    parser.add_argument("--site-counts", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--mixture", required=True)
    args = parser.parse_args()

    (contig, reference), = read_fasta(args.reference).items()
    length = len(reference)
    pooled = next(iter(read_fasta(args.pooled).values()))
    # samtools consensus writes no record for a library without reads: fully masked.
    per_library = read_fasta(args.per_library)
    for library in args.libraries:
        per_library.setdefault(library, "N" * length)
    assert len(pooled) == length, "pooled consensus is not in reference coordinates"

    counts, reads = {}, {}
    for library, path in zip(args.libraries, args.bams):
        reads[library], counts[library] = base_counts(path, contig, length)
    pooled_id = "samples pooled"
    reads[pooled_id] = sum(reads[s] for s in args.samples)
    counts[pooled_id] = [tuple(map(sum, zip(*(counts[s][i] for s in args.samples))))
                         for i in range(length)]
    per_library[pooled_id] = pooled

    with gzip.open(args.site_counts, "wt") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "position", "reference", "consensus", *BASES])
        for library in [*args.libraries, pooled_id]:
            for i, row in enumerate(counts[library]):
                if sum(row):
                    writer.writerow([library, i + 1, reference[i], pooled[i], *row])

    # Sites where the pooled consensus is called and differs from the reference.
    defining = [i for i in range(length)
                if pooled[i] in BASES and reference[i] in BASES and pooled[i] != reference[i]]

    summary, mixture = [], []
    for library in [*args.libraries, pooled_id]:
        depth = [sum(row) for row in counts[library]]
        consensus = per_library[library]
        called = [i for i in range(length) if consensus[i] in BASES]
        both = [i for i in called if pooled[i] in BASES]
        summary.append({
            "id": library,
            "reads": reads[library],
            "mean_depth": round(sum(depth) / length, 2),
            "breadth_1x": round(sum(d >= 1 for d in depth) / length, 4),
            f"breadth_{args.min_depth}x": round(sum(d >= args.min_depth for d in depth) / length, 4),
            "called": len(called),
            "masked_fraction": round(1 - len(called) / length, 4),
            "differs_from_reference": sum(consensus[i] != reference[i] for i in called
                                          if reference[i] in BASES),
            "called_with_pooled": len(both),
            "differs_from_pooled": sum(consensus[i] != pooled[i] for i in both),
        })

        deep = [row for row in counts[library] if sum(row) >= args.mixture_depth]
        minor_shares = [sorted(row)[-2] / sum(row) for row in deep]
        mixed = [s for s in minor_shares if s >= args.minor_fraction]
        at_defining = [counts[library][i] for i in defining]
        consensus_reads = sum(row[BASES.index(pooled[i])] for i, row in zip(defining, at_defining))
        reference_reads = sum(row[BASES.index(reference[i])] for i, row in zip(defining, at_defining))
        total = sum(map(sum, at_defining))
        mixture.append({
            "id": library,
            f"sites_depth_{args.mixture_depth}": len(deep),
            "minor_sites": len(mixed),
            "median_minor_share_all_sites": round(statistics.median(minor_shares), 4) if deep else "",
            "median_minor_share_minor_sites": round(statistics.median(mixed), 3) if mixed else "",
            "defining_sites": len(defining),
            "defining_sites_covered": sum(sum(row) > 0 for row in at_defining),
            "reads_at_defining_sites": total,
            "consensus_allele": consensus_reads,
            "reference_allele": reference_reads,
            "other_allele": total - consensus_reads - reference_reads,
            "non_consensus_share": round(1 - consensus_reads / total, 4) if total else "",
        })

    for path, rows in ((args.summary, summary), (args.mixture, mixture)):
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
