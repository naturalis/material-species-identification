#!/usr/bin/env python
"""Does bleed_scan.py recover a known spillover fraction?

The fixture cannot show this: a 50k subsample has too little self-sharing to call. So build
a source library from a molecule pool with PCR-like duplication, and controls that draw a
known fraction of their reads from that same pool, on either strand. Run from the repo root
inside the env: python tests/test_bleed_scan.py
"""
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = "workflow/scripts/bleed_scan.py"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MIN_LENGTH = 35
ADAPTER_MOTIFS = ["GATCGGAAGAGC", "GCTCTTCCGATC"]
ADAPTER_CHIMERA = "GATCGGAAGAGCACACGTCTGAACTCCAGTCTACACTCTTTCCCTACACGACGCTCTTCCGATC"
SOURCE_READS = 120_000
CONTROL_READS = 8_000
TOLERANCE = 0.08


def molecules(rng, n):
    return [
        "".join(rng.choices("ACGT", k=rng.randint(MIN_LENGTH, 120))) for _ in range(n)
    ]


def draw(rng, pool, weights, n):
    reads = rng.choices(pool, weights=weights, k=n)
    return [r if rng.random() < 0.5 else r.translate(COMPLEMENT)[::-1] for r in reads]


def write_fastq(path, reads):
    with gzip.open(path, "wt", compresslevel=1) as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")


def run_scan(directory, source, controls):
    out = directory / "bleed.tsv"
    command = [
        sys.executable, SCRIPT, "--source-id", "S", "--source", str(source),
        "--min-length", str(MIN_LENGTH), "--probe-reads", "20000", "--min-self-hits", "50",
        "--min-expected-shared", "10", "--exclude-motifs", *ADAPTER_MOTIFS,
        "--flag-fraction", "0.10", "--top-shared", "5", "--seed", "1",
        "--out", str(out), "--out-top", str(directory / "top.tsv"),
    ]
    for name, path in controls.items():
        command += ["--control", f"{name}={path}"]
    subprocess.run(command, check=True)
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    return {r["control"]: r for r in (dict(zip(header, l.split("\t"))) for l in lines[1:])}


def main():
    rng = random.Random(7)
    pool = molecules(rng, 60_000)
    weights = [rng.lognormvariate(0, 1.5) for _ in pool]
    unrelated = molecules(rng, 60_000)

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        source = directory / "source.fastq.gz"
        write_fastq(source, draw(rng, pool, weights, SOURCE_READS) + [ADAPTER_CHIMERA] * 500)

        truth = {"none": 0.0, "low": 0.05, "mid": 0.30, "all": 1.0}
        controls = {}
        for name, fraction in truth.items():
            spilled = int(CONTROL_READS * fraction)
            reads = draw(rng, pool, weights, spilled) + draw(
                rng, unrelated, [1] * len(unrelated), CONTROL_READS - spilled
            )
            controls[name] = directory / f"{name}.fastq.gz"
            write_fastq(controls[name], reads)

        # Shares only adapter chimeras with the source: must not read as spillover.
        controls["adapter"] = directory / "adapter.fastq.gz"
        write_fastq(
            controls["adapter"],
            [ADAPTER_CHIMERA] * 3000 + draw(rng, unrelated, [1] * len(unrelated), 5000),
        )
        # Too few reads to show anything: must not read as cleared.
        controls["tiny"] = directory / "tiny.fastq.gz"
        write_fastq(controls["tiny"], draw(rng, unrelated, [1] * len(unrelated), 50))

        results = run_scan(directory, source, controls)

    failed = False
    print("control  truth  estimate  verdict")
    for name, fraction in truth.items():
        estimate = float(results[name]["est_fraction"])
        print(f"{name:7}  {fraction:5.2f}  {estimate:8.3f}  {results[name]['verdict']}")
        failed |= abs(estimate - fraction) > TOLERANCE
    print(f"adapter  shares {results['adapter']['shared_reads']} reads, "
          f"{results['adapter']['verdict']}; tiny {results['tiny']['verdict']}")
    expected = {
        "none": "not_flagged", "low": "not_flagged", "mid": "flagged", "all": "flagged",
        "adapter": "not_flagged", "tiny": "underpowered",
    }
    failed |= any(results[k]["verdict"] != v for k, v in expected.items())
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
