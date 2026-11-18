#!/usr/bin/env python
"""Read-survival table and fragment-length table. Each library stays its own row.

Counts come from the AdapterRemoval .settings file, fragment lengths from seqkit. Every
input read must be accounted for, seqkit and AdapterRemoval must agree on the merged count,
and raw pairs must equal the demultiplexer's polonies when known. Any disagreement is a hard
error. An implausible merge rate is only flagged: the gates are manual, and a control is
allowed to look strange.
"""
import argparse
import re
import sys
from pathlib import Path

import pandas as pd

SETTINGS_FIELDS = {
    "raw_pairs": "Total number of read pairs",
    "aligned_pairs": "Number of well aligned read pairs",
    "reads_with_adapter": "Number of reads with adapters[1]",
    "discarded_mate1": "Number of discarded mate 1 reads",
    "discarded_mate2": "Number of discarded mate 2 reads",
    "singleton_mate1": "Number of singleton mate 1 reads",
    "singleton_mate2": "Number of singleton mate 2 reads",
    "collapsed_full": "Number of full-length collapsed pairs",
    "collapsed_truncated": "Number of truncated collapsed pairs",
    "retained_reads": "Number of retained reads",
}


def parse_settings(path):
    text = Path(path).read_text()
    values = {}
    for field, label in SETTINGS_FIELDS.items():
        match = re.search(rf"^{re.escape(label)}: (\d+)$", text, re.MULTILINE)
        if not match:
            sys.exit(f"{path}: no line for {label!r}; AdapterRemoval format changed?")
        values[field] = int(match.group(1))
    return values


def weighted_median(lengths):
    cumulative = lengths["reads"].cumsum()
    return int(lengths.loc[cumulative >= lengths["reads"].sum() / 2, "length"].iloc[0])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sheet", required=True)
    parser.add_argument("--ingest", required=True)
    parser.add_argument("--directory", required=True)
    parser.add_argument("--merge-rate-min", type=float, required=True)
    parser.add_argument("--merge-rate-max", type=float, required=True)
    parser.add_argument("--out-stats", required=True)
    parser.add_argument("--out-lengths", required=True)
    args = parser.parse_args()

    sheet = pd.read_csv(args.sheet, sep="\t")
    ingest = pd.read_csv(args.ingest, sep="\t").set_index("library")

    rows, length_tables, errors = [], [], []
    for lib in sheet.itertuples():
        prefix = Path(args.directory) / lib.library / lib.library
        row = {"library": lib.library, "id": lib.id, "role": lib.role}
        row.update(parse_settings(f"{prefix}.settings"))

        lengths = pd.read_csv(f"{prefix}.lengths.tsv", sep="\t")
        merged = row["collapsed_full"] + row["collapsed_truncated"]
        if lengths["reads"].sum() != merged:
            errors.append(
                f"{lib.id}: seqkit counts {lengths['reads'].sum()} merged, settings {merged}"
            )
        polonies = ingest.loc[lib.library, "polonies"]
        if pd.notna(polonies) and int(polonies) != row["raw_pairs"]:
            errors.append(f"{lib.id}: {row['raw_pairs']} raw pairs, {int(polonies)} polonies")

        raw = row["raw_pairs"]
        singletons = row["singleton_mate1"] + row["singleton_mate2"]
        discarded = row["discarded_mate1"] + row["discarded_mate2"]
        unmerged = (row["retained_reads"] - merged - singletons) // 2
        if 2 * raw != 2 * merged + 2 * unmerged + singletons + discarded:
            errors.append(f"{lib.id}: reads in and reads out do not balance")

        row["polonies"] = polonies
        row["merged"] = merged
        row["unmerged_pairs"] = unmerged
        row["merge_rate"] = round(merged / raw, 4) if raw else float("nan")
        row["discarded_fraction"] = round(discarded / (2 * raw), 4) if raw else float("nan")
        row["merge_rate_ok"] = args.merge_rate_min <= row["merge_rate"] <= args.merge_rate_max
        if len(lengths):
            row["length_mode"] = int(lengths.loc[lengths["reads"].idxmax(), "length"])
            row["length_median"] = weighted_median(lengths)
            total_bases = (lengths["length"] * lengths["reads"]).sum()
            row["length_mean"] = round(total_bases / lengths["reads"].sum(), 1)
        rows.append(row)
        length_tables.append(lengths.assign(id=lib.id, role=lib.role))

    if errors:
        sys.exit("Read accounting failed:\n  " + "\n  ".join(errors))

    stats = pd.DataFrame(rows)
    stats.to_csv(args.out_stats, sep="\t", index=False)
    pd.concat(length_tables)[["id", "role", "length", "reads"]].to_csv(
        args.out_lengths, sep="\t", index=False
    )
    for row in stats.loc[~stats["merge_rate_ok"]].itertuples():
        print(f"WARNING {row.id}: merge rate {row.merge_rate} outside the plausible range")


if __name__ == "__main__":
    main()
