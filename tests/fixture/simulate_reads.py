#!/usr/bin/env python
"""Build a synthetic test run from the public reference genomes in tests/fixture/references.

The pipeline's own fixture must not carry anyone's sequencing data, so the reads here are
simulated: fragments are cut from RefSeq mitogenomes, given a length distribution like a real
degraded library, sequenced from both ends at 150 bp with adapter read-through where the
fragment is shorter than the read, and written as gzipped paired FASTQ next to a manifest in
the layout step 0 expects.

It is a fixture, not a simulator anyone should trust for benchmarking: error and damage models
are crude and there is no PCR duplication beyond what resampling produces.

    tests/fixture/simulate_reads.py            # from the repository root
"""

import argparse
import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
READ_LENGTH = 150
# TruSeq-style, with the 10 nt index masked out; matches config adapter1/adapter2.
ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACGTACGTACATCTCGTATGCCGTCTTCTGCTTG"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTACGTACGTACGTGTAGATCTCGGTGGTCGCCGTATC"

# id -> (role, number of fragments, share drawn from the panel rather than the background)
LIBRARIES = {
    "sample-1": ("sample", 20000, 0.55),
    "sample-2": ("sample", 20000, 0.55),
    "blank-1": ("control", 4000, 0.02),
}
# Deliberately artificial index sequences: these are a fixture, not a real kit.
INDEX1 = "ACGTACGTAC"
INDEX2 = {"sample-1": "TGCATGCATG", "sample-2": "GGCCTTAAGG", "blank-1": "CCGGAATTCC"}
WELL = {"sample-1": "A01", "sample-2": "B01", "blank-1": "C01"}


def read_fasta(path, rng=None):
    """Read a FASTA, replacing IUPAC ambiguity codes with a concrete base.

    Reference mitogenomes carry R, Y, S and friends. Real reads never do, and AdapterRemoval
    refuses them, so they are resolved here rather than propagated into the fixture.
    """
    records, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            else:
                records[name].append(line.strip().upper())
    resolved = {}
    for name, parts in records.items():
        sequence = "".join(parts)
        if set(sequence) - set("ACGT"):
            picker = rng or random
            sequence = "".join(b if b in "ACGT" else picker.choice("ACGT") for b in sequence)
        resolved[name] = sequence
    return resolved


def quality(length, rng):
    return "".join(chr(33 + min(41, max(20, int(rng.gauss(37, 3))))) for _ in range(length))


def sequence_errors(fragment, rng, rate=0.001):
    if rng.random() > rate * len(fragment):
        return fragment
    bases = list(fragment)
    i = rng.randrange(len(bases))
    bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--outdir", default=str(HERE))
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--median-length", type=int, default=95)
    args = parser.parse_args()

    out = Path(args.outdir)
    reads = out / "reads"
    reads.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    panel = read_fasta(out / "references" / "panel.fasta", rng)
    everything = read_fasta(out / "references" / "mito_subset.fasta", rng)
    background = {k: v for k, v in everything.items() if k not in panel}
    assert panel and background, "run tests/fixture/build_references.sh first"

    manifest_rows = []
    for library, (role, count, panel_share) in LIBRARIES.items():
        name = f"DEMO_{library}_IDX_i7_1_IDX_i5_1_{WELL[library]}"
        r1_path = reads / f"{name}_R1.fastq.gz"
        r2_path = reads / f"{name}_R2.fastq.gz"
        with gzip.open(r1_path, "wt") as r1, gzip.open(r2_path, "wt") as r2:
            for i in range(count):
                source = panel if rng.random() < panel_share else background
                accession = rng.choice(sorted(source))
                genome = source[accession]
                # Lognormal fragment lengths: short, long-tailed, like a degraded library.
                length = int(rng.lognormvariate(0, 0.45) * args.median_length)
                length = max(30, min(length, len(genome) - 1, 2 * READ_LENGTH + 40))
                start = rng.randrange(len(genome) - length)
                fragment = genome[start:start + length]
                if rng.random() < 0.5:
                    fragment = fragment.translate(COMPLEMENT)[::-1]
                fragment = sequence_errors(fragment, rng)
                # Read through into the adapter when the fragment is shorter than the read.
                forward = (fragment + ADAPTER1)[:READ_LENGTH]
                reverse = (fragment.translate(COMPLEMENT)[::-1] + ADAPTER2)[:READ_LENGTH]
                read_name = f"SIM:1:FIXTURE:1:1:{i}:{i} 1:N:0:{INDEX1}+{INDEX2[library]}"
                r1.write(f"@{read_name}\n{forward}\n+\n{quality(len(forward), rng)}\n")
                r2.write(f"@{read_name.replace(' 1:N', ' 2:N')}\n{reverse}\n+\n"
                         f"{quality(len(reverse), rng)}\n")
        manifest_rows.append((name, INDEX1, INDEX2[library]))
        print(f"{name}: {count} pairs, role {role}")

    manifest = out / "RunManifest.csv"
    with open(manifest, "w", newline="") as fh:
        fh.write("[SETTINGS],,,,,,\nSettingName,Value,Lane,,,,\n")
        fh.write("R1Adapter,AAAAAAAAAAAAAAAAAAA,1+2,,,,\nR1AdapterTrim,FALSE,1+2,,,,\n")
        fh.write("R2Adapter,TTTTTTTTTTTTTTTTTTT,1+2,,,,\nR2AdapterTrim,FALSE,1+2,,,,\n")
        fh.write("[RunValues],,,,,,\nkeyname,Value,,,,,\n")
        fh.write("RUN_NAME,FIXTURE-SIMULATED,,,,,\nSEQUENCING_TYPE,Skimming,,,,,\n")
        fh.write("[SAMPLES],,,,,,\n")
        fh.write("SampleName,Index1,Index2,Lane,PROJECT_NAME,COLLECTION_ID,TAXA\n")
        fh.write("PhiX,AAACCCGGGT,TTTGGGCCCA,1+2,,,\n")
        for name, index1, index2 in manifest_rows:
            fh.write(f"{name},{index1},{index2},1+2,FIXTURE,,\n")
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
