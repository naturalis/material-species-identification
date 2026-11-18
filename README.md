# material-species-identification

A Snakemake pipeline for identifying **which animal a material is made of** — leather, parchment,
horn, bone, fur — from shallow shotgun sequencing (genome skimming) of a sample taken from the
object.

It is built for the hard case: a handful of degraded mitochondrial fragments buried in modern
handling contamination, where the honest answer is often "not enough signal". The pipeline is
therefore organised around **telling a real signal from background**, not around producing an
identification regardless.

## What it does

| Step | Question it answers |
|---|---|
| 0–1 | Did the run demultiplex cleanly, and how long are the fragments? |
| 2 | Are the controls clean, and is there index bleed between libraries? |
| 3 | What is in here at all? (`blastn` against RefSeq mitochondrion + LCA) |
| 4 | Which species do reads really belong to, when close relatives compete? |
| 5 | Is the DNA old, how fragmented is it, and would more sequencing help? |
| 6 | Reconstruct consensus genomes, barcodes, lineage placement, contributor mixture |
| 7–8 | Phylogenetic placement and report (stubs — not implemented) |

Two **manual decision points** ("gates") sit between the steps. The pipeline computes the
evidence and stops; a human makes the call and records it. This is deliberate — these are
judgement calls that should be attributable, not thresholds buried in a config file.

## Design choices worth knowing

- **Competitive mapping, not best-hit.** A short fragment from a conserved region matches cattle,
  pig and sheep equally well. Reads are mapped against every panel genome separately and assigned
  only when one genome beats the rest by a margin. Reads that can't be separated go to a reported
  **ambiguous bin**, never silently to the best scorer.
- **MAPQ is applied per-reference, never on the combined panel.** On a combined panel MAPQ
  deletes exactly the reads that distinguish congeners.
- **Nothing is filtered by taxon.** No contaminant blacklists. Filtering is on alignment quality
  and comparison against the blank.
- **Depth-aware background testing.** "Zero reads in the blank" means nothing if the blank was
  sequenced ten times less deeply. `background_test.py` asks how many reads the blank *should*
  hold if a taxon were background, and tests the shortfall exactly.
- **Exact duplicates collapsed before searching**, counts carried in the read name; position
  deduplication after mapping. Counts are distinct molecules, not sequencer output.

## Requirements

- Linux, conda/mamba, and [Snakemake](https://snakemake.readthedocs.io) ≥ 8
- Per-rule tool environments are created automatically with `--sdm conda`
- Internet access for the reference downloads (NCBI), once
- Roughly 50 GB for the RefSeq mitochondrion BLAST database

## Quick start

```bash
# 1. Base environment
conda create -n material-species-identification --override-channels \
    -c conda-forge -c bioconda --strict-channel-priority \
    python=3.12 snakemake pandas matplotlib seqkit 'adapterremoval>=2.3,<3'
conda activate material-species-identification
```

Keep the `<3` pin: bioconda's current `adapterremoval` is 3.x with a different CLI.

```bash
# 2. Configure
cp config/config.example.yaml config/config.yaml
$EDITOR config/config.yaml          # run_dir, library_id_regex, roles

# 3. References (downloads RefSeq mitochondrion + your panel, verifies every record)
snakemake make_references --sdm conda --cores 8

# 4. Sample sheet, generated from the run manifest — do not hand-edit
python workflow/scripts/make_sample_sheet.py --config config/config.yaml

# 5. Run
snakemake --cores 24 --resources mem_mb=64000                          # steps 0-2
snakemake screen       --sdm conda --rerun-triggers mtime --cores 48 --resources mem_mb=120000  # step 3, hours
snakemake mapping      --sdm conda --rerun-triggers mtime --cores 48 --resources mem_mb=120000  # step 4
snakemake authenticate --sdm conda --rerun-triggers mtime --cores 24 --resources mem_mb=64000   # step 5
snakemake barcode      --sdm conda --rerun-triggers mtime --cores 48 --resources mem_mb=120000  # step 6
```

Always pass `--rerun-triggers mtime`. Snakemake's default includes a software-environment
trigger whose hash contains the absolute environment path, so **moving the project directory
would otherwise rerun everything**.

To see the intended shape without running anything:

```bash
snakemake -n full_plan
snakemake full_plan --rulegraph | dot -Tpng > dag.png
```

## Configuring the reference panel

`config/reference_panel_mitogenomes.csv` lists the complete mitochondrial genomes that compete
in step 4. The shipped panel has 64 records aimed at materials — the usual leather and parchment
animals, their close relatives (so the pipeline can *fail visibly* rather than pick one of a
congener pair), and likely contaminants.

Add a row per genome:

```csv
species,accession,role,use,expected_organism
Sus scrofa,NC_000845.1,panel,domestic pig,
```

Every accession is fetched individually and checked against the record's `/organism` field, not
its title. Mismatches are reported and excluded rather than silently accepted. Where RefSeq holds
only a subspecies record, name it in `expected_organism`.

## Tests

```bash
tests/fixture/simulate_reads.py          # build a synthetic test run from public genomes
tests/run_smoke.sh                       # steps 0-2 end to end, ~10 s
python tests/test_assign_panel.py        # read assignment, with a known truth
python tests/test_bleed_scan.py          # index-bleed estimator
```

The fixture is simulated from public RefSeq mitogenomes, so it ships with the repository and
needs no data of your own.

## Interpreting the output — three cautions

1. **Gate outcomes are not conclusions.** Passing gate 2 means a taxon has plausible coverage
   geometry. It does not mean the object is made of it. Check `background_test.tsv` as well:
   taxa can pass the geometry test and still be indistinguishable from lab background.
2. **An extraction blank controls for extraction onward.** It does not control for anything the
   object picked up before sampling — surface contamination from storage, handling or
   environment. "Above the blank" is not "from the material".
3. **If libraries were built with the USER enzyme**, deamination damage is removed by design and
   damage-based dating is uninformative. The pipeline still reports it; read the result as a
   confirmation, not a dating.

## Status

Steps 0–6 are implemented. Steps 7 (phylogenetic placement) and 8 (report) are stubs that fail
if executed. `contamMix` and the Krona plots are also stubs.

## Licence

MIT — see [LICENSE](LICENSE).

Developed at [Naturalis Biodiversity Center](https://www.naturalis.nl/). The pipeline was built
for a specific forensic material-identification case; that analysis and its data are not part of
this repository.
