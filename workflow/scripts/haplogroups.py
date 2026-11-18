#!/usr/bin/env python
"""Step 6: human haplogroups of the consensus, of each library, and of each contributor.

A consensus from DNA of two people is masked wherever they differ, and haplogrep then fits a
haplogroup to whatever is left: libraries of the same extract can come out as different
haplogroups. So the mixture is separated before it is classified:

1. Mixed sites: pooled-sample depth of --depth or more, second allele at --fraction or more.
2. Haplotype A takes the non-reference allele at every mixed site, haplotype B the reference
   allele (rCRS). Elsewhere both are the pooled consensus.
3. Both are classified. Alleles at mixed sites that haplogrep lists as private mutations of
   A, not part of its haplogroup, are moved to B, and both are classified again.
4. Each library's share of A is the share of reads carrying A's allele at the mixed sites.

This assumes two contributors. Sites with a second allele below --fraction are counted, not
phased: a third, smaller contributor would show there. Everything runs locally.
"""

import argparse
import csv
import gzip
import math
import subprocess
from collections import defaultdict
from pathlib import Path

BASES = "ACGT"
FIELDS = ["SampleID", "Haplogroup", "Rank", "Quality", "Range", "Not_Found_Polys",
          "Found_Polys", "Remaining_Polys"]


def read_fasta(path):
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().replace(" ", "_")
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    return {name: "".join(parts) for name, parts in records.items()}


def classify(sequences, tree, work, tag):
    fasta, out = work / f"{tag}.fasta", work / f"{tag}.txt"
    with open(fasta, "w") as fh:
        for name, sequence in sequences.items():
            fh.write(f">{name}\n{sequence.replace('*', '')}\n")
    subprocess.run(["haplogrep3", "classify", "--in", str(fasta), "--tree", tree, "--out", str(out),
                    "--extend-report", "--hits", "3"], check=True, capture_output=True)
    with open(out, newline="") as fh:
        rows = [{k.strip('"'): v.strip('"') for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]
    return rows


def private(row):
    """Positions of private mutations in haplogrep's Remaining_Polys."""
    out = {}
    for token, kind in zip(row["Remaining_Polys"].split()[::2], row["Remaining_Polys"].split()[1::2]):
        if "PrivateMutation" in kind and token[:-1].isdigit():
            out[int(token[:-1])] = token[-1]
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return "", ""
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return round(centre - half, 3), round(centre + half, 3)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference", required=True)
    parser.add_argument("--pooled", required=True, help="pooled consensus, reference coordinates")
    parser.add_argument("--consensus", required=True, help="pooled consensus with insertions")
    parser.add_argument("--libraries", required=True, help="per-library consensus with insertions")
    parser.add_argument("--summary", required=True, help="consensus_summary.tsv, for coverage")
    parser.add_argument("--site-counts", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--depth", type=int, required=True)
    parser.add_argument("--fraction", type=float, required=True)
    parser.add_argument("--noise", type=float, required=True, help="lower bound for 'minor' sites")
    parser.add_argument("--min-called", type=float, default=0.5,
                        help="classify a library only if this share of its consensus is called")
    parser.add_argument("--work", required=True)
    parser.add_argument("--haplogroups", required=True)
    parser.add_argument("--contributors", required=True)
    parser.add_argument("--sites", required=True)
    parser.add_argument("--haplotypes", required=True)
    args = parser.parse_args()

    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    reference = next(iter(read_fasta(args.reference).values()))
    pooled = next(iter(read_fasta(args.pooled).values()))
    assert len(pooled) == len(reference)

    # Consensus and per-library classification, libraries only where mostly called.
    with open(args.summary, newline="") as fh:
        called = {r["id"]: 1 - float(r["masked_fraction"]) for r in csv.DictReader(fh, delimiter="\t")}
    sequences = read_fasta(args.consensus)
    for name, sequence in read_fasta(args.libraries).items():
        if called.get(name, 0) >= args.min_called:
            sequences[f"library_{name}"] = sequence
    rows = classify(sequences, args.tree, work, "consensus")

    # Mixed sites from the pooled samples.
    counts = defaultdict(dict)
    with gzip.open(args.site_counts, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            counts[r["id"]][int(r["position"])] = [int(r[b]) for b in BASES]
    mixed, minor_only = {}, []
    for position, row in counts["samples pooled"].items():
        depth = sum(row)
        if depth < args.depth:
            continue
        order = sorted(range(4), key=lambda i: -row[i])
        share = row[order[1]] / depth
        if share >= args.fraction:
            mixed[position] = (BASES[order[0]], BASES[order[1]])
        elif share >= args.noise:
            minor_only.append(position)

    a, b = list(pooled), list(pooled)
    for position, (first, second) in mixed.items():
        ref = reference[position - 1]
        other = second if first == ref else first
        a[position - 1] = other if ref in (first, second) else first
        b[position - 1] = ref if ref in (first, second) else second
    first_pass = classify({"contributor_A": "".join(a), "contributor_B": "".join(b)},
                          args.tree, work, "phase1")
    best_a = next(r for r in first_pass if r["SampleID"] == "contributor_A" and r["Rank"] == "1")
    moved = {p: base for p, base in private(best_a).items()
             if p in mixed and a[p - 1] == base}
    for position, base in moved.items():
        b[position - 1], a[position - 1] = base, next(x for x in mixed[position] if x != base)
    final = classify({"contributor_A": "".join(a), "contributor_B": "".join(b)},
                     args.tree, work, "phase2")
    rows.extend(final)

    with open(args.haplogroups, "w", newline="") as fh:
        writer = csv.DictWriter(fh, FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with open(args.haplotypes, "w") as fh:
        for name, sequence in (("contributor_A", a), ("contributor_B", b)):
            fh.write(f">{name}\n{''.join(sequence).replace('*', '')}\n")

    libraries = [k for k in counts if k != "samples pooled"] + ["samples pooled"]
    with open(args.sites, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["position", "reference", "allele_A", "allele_B", "moved_after_phase1",
                         *(f"{lib}_{x}" for lib in libraries for x in ("A", "B", "other"))])
        for position in sorted(mixed):
            ia, ib = BASES.index(a[position - 1]), BASES.index(b[position - 1])
            values = []
            for lib in libraries:
                row = counts[lib].get(position, [0, 0, 0, 0])
                values += [row[ia], row[ib], sum(row) - row[ia] - row[ib]]
            writer.writerow([position, reference[position - 1], a[position - 1], b[position - 1],
                             position in moved, *values])

    best = {r["SampleID"]: r for r in final if r["Rank"] == "1"}
    with open(args.contributors, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "mixed_sites", "minor_only_sites", "reads_A", "reads_B",
                         "share_A", "share_A_low", "share_A_high", "haplogroup_A", "quality_A",
                         "haplogroup_B", "quality_B"])
        for lib in libraries:
            ka = sum(counts[lib].get(p, [0] * 4)[BASES.index(a[p - 1])] for p in mixed)
            kb = sum(counts[lib].get(p, [0] * 4)[BASES.index(b[p - 1])] for p in mixed)
            low, high = wilson(ka, ka + kb)
            writer.writerow([lib, len(mixed), len(minor_only), ka, kb,
                             round(ka / (ka + kb), 3) if ka + kb else "", low, high,
                             best["contributor_A"]["Haplogroup"], best["contributor_A"]["Quality"],
                             best["contributor_B"]["Haplogroup"], best["contributor_B"]["Quality"]])
    print(f"{len(mixed)} mixed sites, {len(moved)} moved from A to B: {sorted(moved)}")
    print(f"{len(minor_only)} sites with a second allele between {args.noise} and {args.fraction}")


if __name__ == "__main__":
    main()
