#!/usr/bin/env python
"""Step 6: how step 4's mapping and assignment recover reads of a taxon, by read length.

Fragments of fixed lengths are cut from real genomes of the taxon (no damage, no sequencing
error), then put through step 4 exactly: bwa aln against the combined panel, every mapped
read against every panel genome alone, assign_panel.py with the same margin and MAPQ.
Genomes used as the truth:
  - the panel genome itself (no divergence),
  - the lineage genome that fits the real reads best (step 6 placement),
  - the genome of the species most distant from the panel genome.
If recovery falls for short reads of a divergent lineage, the observed length distribution of
that taxon is shifted long by the method. The shift is predicted by weighting the observed
human length distribution of the sample libraries with each recovery curve.
"""

import argparse
import csv
import gzip
import random
import statistics
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


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


def run(command, **kwargs):
    subprocess.run(command, check=True, **kwargs)


def bwa_bam(reference, reads, bam, aln_params, threads, work):
    sai = work / "reads.sai"
    with open(sai, "wb") as fh:
        run(["bwa", "aln", "-t", str(threads), *aln_params.split(), reference, str(reads)],
            stdout=fh, stderr=subprocess.DEVNULL)
    samse = subprocess.Popen(["bwa", "samse", reference, str(sai), str(reads)],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    run(["samtools", "view", "-b", "-F", "4", "-o", str(bam), "-"], stdin=samse.stdout)
    assert samse.wait() == 0
    sai.unlink()


def weighted_median(weights):
    total, running = sum(weights.values()), 0
    for length in sorted(weights):
        running += weights[length]
        if running >= total / 2:
            return length


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--genomes", required=True)
    parser.add_argument("--genome-scores", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--panel", required=True)
    parser.add_argument("--single", required=True)
    parser.add_argument("--aln-params", required=True)
    parser.add_argument("--margin", type=int, required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--max-length", type=int, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--per-length", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--observed-lengths", nargs="+", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--predicted", required=True)
    args = parser.parse_args()

    genomes = read_fasta(args.genomes)
    with open(args.genome_scores, newline="") as fh:
        scores = list(csv.DictReader(fh, delimiter="\t"))
    own = [r for r in scores if r["organism"].startswith(args.species)]
    truths = {"panel genome": args.reference, "best-fitting lineage genome": own[0]["accession"]}
    farthest = max(own, key=lambda r: float(r["distance_to_reference"]))
    truths["most distant genome of the species"] = farthest["accession"]
    distance = {r["accession"]: r["distance_to_reference"] for r in scores}
    cluster = {r["accession"]: r["cluster"] for r in scores}

    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    lengths = range(args.min_length, args.max_length + 1, args.step)
    origin = {}
    reads = work / "simulated.fasta"
    with open(reads, "w") as fh:
        index = 0
        for label, accession in truths.items():
            genome = genomes[accession].replace("N", "")
            for length in lengths:
                for _ in range(args.per_length):
                    start = rng.randrange(len(genome) - length)
                    fragment = genome[start:start + length]
                    if rng.random() < 0.5:
                        fragment = fragment.translate(COMPLEMENT)[::-1]
                    index += 1
                    name = f"u{index}_1"
                    origin[name] = (label, length)
                    fh.write(f">{name}\n{fragment}\n")
    print(f"{index} fragments from {len(truths)} genomes: {truths}")

    # Step 4, as in rules map_panel and assign_panel.
    bwa_bam(args.panel, reads, work / "panel.bam", args.aln_params, args.threads, work)
    mapped = work / "mapped.fasta"
    with open(mapped, "w") as fh:
        run(["samtools", "fasta", str(work / "panel.bam")], stdout=fh, stderr=subprocess.DEVNULL)
    prefilter = {line[1:].strip() for line in open(mapped) if line.startswith(">")}
    single = sorted(Path(args.single).glob("*.fasta"))
    bams = []
    for reference in single:
        bam = work / f"{reference.stem}.bam"
        bwa_bam(str(reference), mapped, bam, args.aln_params, args.threads, work)
        bams.append(str(bam))
    assigned_path = work / "assigned.tsv.gz"
    run(["python", str(Path(__file__).with_name("assign_panel.py")), "--bams", *bams,
         "--panel", args.panel, "--library", "simulated", "--margin", str(args.margin),
         "--min-mapq", str(args.min_mapq), "--assigned", str(assigned_path),
         "--metrics", str(work / "metrics.tsv"), "--ambiguous", str(work / "ambiguous.tsv")],
        stdout=subprocess.DEVNULL)
    assigned = {}
    with gzip.open(assigned_path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            assigned[row["read"]] = row["accession"]

    tally = defaultdict(Counter)
    for name, key in origin.items():
        tally[key]["simulated"] += 1
        tally[key]["prefilter_mapped"] += name in prefilter
        hit = assigned.get(name)
        tally[key]["assigned_target"] += hit == args.reference
        tally[key]["assigned_other"] += hit is not None and hit != args.reference
    recovery = defaultdict(dict)
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["truth", "accession", "cluster", "distance_to_reference", "length",
                         "simulated", "prefilter_mapped", "assigned_target", "assigned_other",
                         "recovery"])
        for label, accession in truths.items():
            for length in lengths:
                t = tally[label, length]
                recovery[label][length] = t["assigned_target"] / t["simulated"]
                writer.writerow([label, accession, cluster.get(accession, ""),
                                 distance.get(accession, ""), length, t["simulated"],
                                 t["prefilter_mapped"], t["assigned_target"], t["assigned_other"],
                                 round(recovery[label][length], 4)])

    observed = defaultdict(Counter)
    for path in args.observed_lengths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                observed[row["species"]][int(row["length"])] += int(row["reads"])
    human, own_lengths = observed["Homo sapiens"], observed[args.species]
    grid = sorted(lengths)

    def rate(curve, length):
        if length <= grid[0]:
            return curve[grid[0]]
        if length >= grid[-1]:
            return curve[grid[-1]]
        low = max(g for g in grid if g <= length)
        high = min(g for g in grid if g >= length)
        if low == high:
            return curve[low]
        return curve[low] + (curve[high] - curve[low]) * (length - low) / (high - low)

    with open(args.predicted, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["truth", "accession", "human_reads", "human_median",
                         "predicted_median_if_human_lengths", "observed_median",
                         "observed_reads", "recovery_below_60", "recovery_60_to_120",
                         "recovery_above_120"])
        for label, accession in truths.items():
            curve = recovery[label]
            weights = {length: n * rate(curve, length) for length, n in human.items()}
            bands = [[g for g in grid if g < 60], [g for g in grid if 60 <= g <= 120],
                     [g for g in grid if g > 120]]
            writer.writerow([
                label, accession, sum(human.values()), weighted_median(human),
                weighted_median(weights), weighted_median(own_lengths),
                sum(own_lengths.values()),
                *(round(statistics.mean(curve[g] for g in band), 3) for band in bands),
            ])


if __name__ == "__main__":
    main()
