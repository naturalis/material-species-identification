#!/usr/bin/env python
"""Join per-control read-level profiles with their step 1 numbers. One row per control."""
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readstats", required=True)
    parser.add_argument("--profiles", nargs="+", required=True)
    parser.add_argument("--top", nargs="+", required=True)
    parser.add_argument("--out-profile", required=True)
    parser.add_argument("--out-top", required=True)
    args = parser.parse_args()

    stats = pd.read_csv(args.readstats, sep="\t")
    keep = ["id", "role", "raw_pairs", "merge_rate", "discarded_fraction",
            "length_mode", "length_median"]
    profiles = pd.concat(pd.read_csv(p, sep="\t") for p in args.profiles)
    joined = stats.loc[stats["role"] == "control", keep].merge(profiles, on="id", validate="1:1")
    joined.to_csv(args.out_profile, sep="\t", index=False)

    pd.concat(pd.read_csv(p, sep="\t") for p in args.top).to_csv(
        args.out_top, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
