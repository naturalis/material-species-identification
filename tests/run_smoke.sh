#!/usr/bin/env bash
# Steps 0-2 on the simulated fixture, from clean. Run from the repo root inside the env.
set -euo pipefail

# The fixture config is an overlay on config/config.yaml, so that has to exist.
[ -f config/config.yaml ] || cp config/config.example.yaml config/config.yaml
[ -f tests/fixture/RunManifest.csv ] || tests/fixture/simulate_reads.py

rm -rf tests/fixture/results tests/fixture/logs
python workflow/scripts/make_sample_sheet.py --config config/config.yaml --config tests/fixture/config.yaml
time snakemake --configfile tests/fixture/config.yaml --rerun-triggers mtime \
    --cores 8 --resources mem_mb=16000 --quiet progress
