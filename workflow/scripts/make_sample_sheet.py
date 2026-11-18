#!/usr/bin/env python
"""Generate the sample sheet from RunManifest.csv. The sheet is never hand-edited.

Snakemake reads the sheet at parse time, so this runs before the first snakemake call:

    python workflow/scripts/make_sample_sheet.py --config config/config.yaml

Later --config files overlay earlier ones, like snakemake --configfile. The ingest rule
calls build_rows() again and fails if the sheet on disk no longer matches the manifest.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

COLUMNS = ["library", "id", "role", "index1", "index2", "well", "R1", "R2"]
ROLES = {"sample", "control"}


def read_manifest_samples(path):
    """Return the rows of the [SAMPLES] section of an AVITI RunManifest.csv."""
    rows, header, in_samples = [], None, False
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            first = row[0].strip() if row else ""
            if first.startswith("["):
                in_samples = first == "[SAMPLES]"
                continue
            if not in_samples or not first or first.startswith("#"):
                continue
            if header is None:
                header = [c.strip() for c in row]
                continue
            rows.append(dict(zip(header, (c.strip() for c in row))))
    if not rows:
        sys.exit(f"{path}: no [SAMPLES] rows found")
    return rows


def build_rows(cfg):
    run_dir = Path(cfg["run_dir"])
    fastq_dir = run_dir / cfg["fastq_dir"]
    id_re, well_re = re.compile(cfg["library_id_regex"]), re.compile(cfg["well_regex"])
    errors, rows = [], []

    for entry in read_manifest_samples(run_dir / cfg["manifest"]):
        library = entry["SampleName"]
        if library in cfg["non_libraries"]:
            continue
        id_match, well_match = id_re.search(library), well_re.search(library)
        if not id_match or not well_match:
            errors.append(f"{library}: id or well not parseable from the name")
            continue
        library_id = id_match.group(1)
        role = str(cfg["roles"].get(library_id))
        if role not in ROLES:
            errors.append(f"{library}: id {library_id!r} has no valid role in config 'roles'")
            continue
        rows.append(
            {
                "library": library,
                "id": library_id,
                "role": role,
                "index1": entry["Index1"],
                "index2": entry["Index2"],
                "well": well_match.group(1),
                "R1": str(fastq_dir / f"{library}_R1.fastq.gz"),
                "R2": str(fastq_dir / f"{library}_R2.fastq.gz"),
            }
        )

    libraries = [r["library"] for r in rows]
    for library in sorted({x for x in libraries if libraries.count(x) > 1}):
        errors.append(f"{library}: listed more than once in the manifest")
    if errors:
        sys.exit("Sample sheet not generated:\n  " + "\n  ".join(errors))
    return rows


def format_sheet(rows):
    lines = ["\t".join(COLUMNS)]
    lines += ["\t".join(r[c] for c in COLUMNS) for r in rows]
    return "\n".join(lines) + "\n"


def load_config(paths):
    cfg = {}

    def overlay(base, new):
        for key, value in new.items():
            if isinstance(value, dict) and isinstance(base.get(key), dict):
                overlay(base[key], value)
            else:
                base[key] = value

    for path in paths:
        overlay(cfg, yaml.safe_load(Path(path).read_text()))
    return cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", action="append", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    Path(cfg["sample_sheet"]).write_text(format_sheet(build_rows(cfg)))
    print(f"wrote {cfg['sample_sheet']}", file=sys.stderr)


if __name__ == "__main__":
    main()
