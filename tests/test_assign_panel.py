#!/usr/bin/env python
"""Step 4 assignment against simulated reads from the fixture panel, where truth is known.

Run from the repo root in the step 4 environment (bwa and samtools on PATH):
    PATH=.snakemake/conda/<hash>/bin:$PATH python tests/test_assign_panel.py

Reads are cut from the public reference genomes only, never from sequencing data. Pig and
human are far from everything else, so their reads must be assigned to them. Cattle and zebu
differ by about 1%, so cattle reads must be assigned to cattle or tie with zebu, and no read
may be assigned to a genome it did not come from. Every read is sent twice under two names, so half are position duplicates.
"""

import csv
import gzip
import random
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

PANEL = Path("tests/fixture/references/panel.fasta")
ALN = "-l 1024 -n 0.01 -o 2"
SOURCES = {"NC_000845.1": 150, "NC_012920.1": 150, "NC_006853.1": 150}


def read_fasta(path):
    records, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        else:
            records[name].append(line.strip())
    return {name: "".join(parts) for name, parts in records.items()}


def simulate(genomes, out, seed=20261109):
    rng, truth, n = random.Random(seed), {}, 0
    with open(out, "w") as fh:
        for accession, count in SOURCES.items():
            genome = genomes[accession]
            for _ in range(count):
                length = rng.randint(35, 120)
                start = rng.randrange(0, len(genome) - length)
                read = list(genome[start:start + length])
                position = rng.randrange(length)
                read[position] = rng.choice([b for b in "ACGT" if b != read[position]])
                for _ in range(2):
                    n += 1
                    truth[f"u{n}_1"] = accession
                    fh.write(f">u{n}_1\n{''.join(read)}\n")
    return truth


def run(command):
    subprocess.run(command, shell=True, check=True, stderr=subprocess.DEVNULL)


def main():
    genomes = read_fasta(PANEL)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        reads = tmp / "reads.fasta"
        truth = simulate(genomes, reads)
        for accession, sequence in genomes.items():
            reference = tmp / f"{accession}.fasta"
            reference.write_text(f">{accession}\n{sequence}\n")
            run(f"bwa index {reference}")
            run(f"bwa aln {ALN} {reference} {reads} > {tmp}/{accession}.sai")
            run(f"bwa samse {reference} {tmp}/{accession}.sai {reads} "
                f"| samtools view -b -F 4 -o {tmp}/{accession}.bam -")
        run(f"{sys.executable} workflow/scripts/assign_panel.py --bams {tmp}/NC_*.bam "
            f"--panel {PANEL} --library test --margin 2 --min-mapq 25 "
            f"--assigned {tmp}/assigned.tsv.gz --metrics {tmp}/metrics.tsv "
            f"--ambiguous {tmp}/ambiguous.tsv")
        with gzip.open(tmp / "assigned.tsv.gz", "rt") as fh:
            assigned = list(csv.DictReader(fh, delimiter="\t"))
        ambiguous = list(csv.DictReader(open(tmp / "ambiguous.tsv"), delimiter="\t"))

    outcome = Counter((truth[row["read"]], row["accession"]) for row in assigned)
    print("assigned (true, called):", dict(outcome))
    print("ambiguous:", [(row["tied_species"], row["distinct"]) for row in ambiguous])
    sent = {accession: 2 * count for accession, count in SOURCES.items()}

    wrong = {pair: n for pair, n in outcome.items() if pair[0] != pair[1]}
    assert not wrong, f"reads assigned to the wrong genome: {wrong}"
    for accession in ("NC_000845.1", "NC_012920.1"):
        assert outcome[accession, accession] >= 0.95 * sent[accession], accession
    # Only cattle has a close relative, so every tie must include cattle. A short read in a
    # conserved gene can also tie with pig or sheep; that is the bin doing its job.
    ties = {row["tied_species"] for row in ambiguous}
    assert all("Bos taurus" in tie for tie in ties), f"ties without the source: {ties}"
    tied_to_zebu = sum(int(r["distinct"]) for r in ambiguous
                       if r["tied_species"] == "Bos indicus + Bos taurus")
    cattle = outcome["NC_006853.1", "NC_006853.1"] + tied_to_zebu
    assert cattle >= 0.95 * sent["NC_006853.1"], f"{cattle} cattle reads kept to Bos"
    duplicates = sum(row["duplicate"] == "1" for row in assigned)
    assert duplicates == len(assigned) // 2, f"{duplicates} duplicates of {len(assigned)}"
    print("ok")


if __name__ == "__main__":
    main()
