# Test fixture

A small, fully synthetic run for developing and smoke-testing the pipeline. Steps 0–2 complete
in about 7 seconds.

**No real sequencing data is used or shipped here.** Reads are simulated from public RefSeq
mitochondrial genomes, so the fixture works in a fresh clone with nothing of your own.

## Building it

```bash
tests/fixture/simulate_reads.py     # writes reads/ and RunManifest.csv
tests/run_smoke.sh                  # does the above if needed, then runs steps 0-2
```

`reads/` and the generated `RunManifest.csv` and `samples.tsv` are not tracked in Git — they are
regenerated deterministically from the seed in `simulate_reads.py`.

## What it contains

| Library | Role | Pairs | Panel share |
|---|---|---|---|
| `sample-1` | sample | 20,000 | 55% |
| `sample-2` | sample | 20,000 | 55% |
| `blank-1` | control | 4,000 | 2% |

"Panel share" is the proportion of fragments drawn from the five panel genomes rather than from
the background mitogenomes — so the samples carry a signal and the blank is nearly clean, which
is the situation the pipeline is built to tell apart.

Fragments have lognormal lengths around 95 bp, are sequenced 150 bp from each end with adapter
read-through, and carry a low rate of substitution errors. IUPAC ambiguity codes in the
reference genomes are resolved to concrete bases, because real reads do not contain them and
AdapterRemoval rejects them.

## References

Built by `build_references.sh` from one RefSeq `mitochondrion` release:

- `panel.fasta` / `panel.csv` — 5 genomes: human, pig, cattle, zebu, sheep. The cattle/zebu pair
  is deliberate: they are close enough that many reads cannot be told apart, so the ambiguous
  bin has something real to do and `test_assign_panel.py` can check it behaves.
- `mito_subset.fasta` — the panel plus 45 random mitogenomes under 100 kb, as background.

To rebuild them you need a RefSeq release file:

```bash
tests/fixture/build_references.sh path/to/mitochondrion.1.1.genomic.fna.gz
```

## Caveats

This is a fixture, not a benchmark. The error model is crude, there is no damage model, and PCR
duplication only arises from resampling the same positions. Numbers produced from it mean
nothing biologically — they exist to prove the code runs and the plumbing is connected.
