"""Step 0 control point: the sheet matches the manifest, and every library resolves to two
readable, paired files. Runs from a snakemake script: directive."""
import csv
import gzip
import os
import sys
from pathlib import Path

sys.path.insert(0, snakemake.scriptdir)
from make_sample_sheet import build_rows, format_sheet  # noqa: E402

cfg = snakemake.config
sheet_path = Path(cfg["sample_sheet"])
if sheet_path.read_text() != format_sheet(build_rows(cfg)):
    sys.exit(f"{sheet_path} is stale: it no longer matches {cfg['manifest']}. Regenerate it.")


def first_record(path):
    with gzip.open(path, "rt") as fh:
        return [fh.readline().rstrip("\n") for _ in range(4)]


def polonies_by_library():
    if not cfg.get("index_assignment"):
        return {}
    with open(Path(cfg["run_dir"]) / cfg["index_assignment"], newline="") as fh:
        return {
            row["SampleName"]: row["NumPoloniesAssigned"]
            for row in csv.DictReader(fh)
            if row["Lane"] == "1+2"
        }


polonies = polonies_by_library()
failures, report = [], []
with open(sheet_path, newline="") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        library = row["library"]
        records = {}
        for mate in ("R1", "R2"):
            path = Path(row[mate])
            if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size == 0:
                failures.append(f"{library} {mate}: missing, unreadable or empty: {path}")
                continue
            record = first_record(path)
            if not record[0].startswith("@") or len(record[1]) != len(record[3]):
                failures.append(f"{library} {mate}: first record is not FASTQ: {path}")
            records[mate] = record
        if len(records) == 2 and records["R1"][0].split()[0] != records["R2"][0].split()[0]:
            failures.append(f"{library}: first R1 and R2 read names differ")
        if polonies and library not in polonies:
            failures.append(f"{library}: not in {cfg['index_assignment']}")
        report.append(
            {
                "library": library,
                "id": row["id"],
                "role": row["role"],
                "R1_bytes": Path(row["R1"]).stat().st_size if "R1" in records else "",
                "R2_bytes": Path(row["R2"]).stat().st_size if "R2" in records else "",
                "polonies": polonies.get(library, ""),
            }
        )

if failures:
    sys.exit("Ingest failed:\n  " + "\n  ".join(failures))

with open(snakemake.output.report, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(report[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(report)
