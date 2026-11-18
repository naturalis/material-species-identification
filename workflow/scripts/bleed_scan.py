#!/usr/bin/env python
"""Exact-sequence bleed test: how much of each control looks like spillover from one source.

All libraries share an i7, so a misassigned read carries a molecule from the source
library into the control with its exact fragment boundaries intact. Merged fragments are
therefore compared by identity, on either strand.

Raw sharing is not the answer. A spilled read only matches if another copy of the same
molecule was also sequenced in the source, and that chance depends on how deeply the
source is sequenced. So the source is calibrated against itself: a random probe of its
own reads gives self_share, the fraction having another copy in the source. A control
that is entirely spillover would show the same fraction, hence

    est_fraction = control_share / self_share

Fragments carrying an adapter motif are excluded, since they are shared without any bleed.
Other sequences shared for another reason, such as reagent contaminants, still inflate the
estimate, so est_fraction is an upper bound and can exceed 1. The top shared sequences are
written out so a flag can be read rather than trusted.

No shared reads only clears a control if it could have shown some. expected_at_flag is the
number of shared reads a control would show at exactly flag_fraction; if that is too low, or
the source baseline is too thin, the verdict is underpowered rather than not_flagged.

Memory scales with the reads in the controls plus the probe, not with the source.
"""
import argparse
import csv
import random
import subprocess
from collections import Counter

COMPLEMENT = bytes.maketrans(b"ACGTN", b"TGCAN")


def canonical(seq):
    return min(seq, seq.translate(COMPLEMENT)[::-1])


def read_fragments(paths, min_length, exclude):
    """Yield the canonical form of every sequence of at least min_length with no excluded motif."""
    proc = subprocess.Popen(["seqkit", "seq", "-s", "-w", "0", *paths], stdout=subprocess.PIPE)
    for line in proc.stdout:
        seq = line.rstrip()
        if len(seq) >= min_length and not any(motif in seq for motif in exclude):
            yield canonical(seq)
    if proc.wait() != 0:
        raise SystemExit(f"seqkit failed on {paths}")


def reservoir(items, size, rng):
    sample = []
    for i, item in enumerate(items):
        if i < size:
            sample.append(item)
        else:
            j = rng.randrange(i + 1)
            if j < size:
                sample[j] = item
    return sample


def verdict(est_fraction, flag_fraction, powered):
    if not powered:
        return "underpowered"
    return "flagged" if est_fraction >= flag_fraction else "not_flagged"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--source", nargs="+", required=True)
    parser.add_argument("--control", action="append", required=True, metavar="ID=FILE[,FILE]")
    parser.add_argument("--min-length", type=int, required=True)
    parser.add_argument("--probe-reads", type=int, required=True)
    parser.add_argument("--exclude-motifs", nargs="*", default=[])
    parser.add_argument("--min-self-hits", type=int, required=True)
    parser.add_argument("--min-expected-shared", type=float, required=True)
    parser.add_argument("--flag-fraction", type=float, required=True)
    parser.add_argument("--top-shared", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--out-top", required=True)
    args = parser.parse_args()

    exclude = [motif.encode() for motif in args.exclude_motifs]
    controls = {}
    for spec in args.control:
        control_id, paths = spec.split("=", 1)
        controls[control_id] = Counter(read_fragments(paths.split(","), args.min_length, exclude))

    probe = reservoir(
        read_fragments(args.source, args.min_length, exclude),
        args.probe_reads,
        random.Random(args.seed),
    )

    # Second pass: occurrences in the source of everything we need to look up.
    source_counts = dict.fromkeys(probe, 0)
    for counts in controls.values():
        source_counts.update(dict.fromkeys(counts, 0))
    source_reads = 0
    for fragment in read_fragments(args.source, args.min_length, exclude):
        source_reads += 1
        if fragment in source_counts:
            source_counts[fragment] += 1

    # A probe read always finds itself, so another copy means a count of at least 2.
    probe_hits = sum(1 for fragment in probe if source_counts[fragment] >= 2)
    self_share = probe_hits / len(probe) if probe else 0.0

    rows, top = [], []
    for control_id, counts in controls.items():
        control_reads = sum(counts.values())
        shared = {f: n for f, n in counts.items() if source_counts[f] >= 1}
        shared_reads = sum(shared.values())
        control_share = shared_reads / control_reads if control_reads else 0.0
        est_fraction = control_share / self_share if self_share else float("nan")
        expected_at_flag = control_reads * self_share * args.flag_fraction
        powered = (
            probe_hits >= args.min_self_hits and expected_at_flag >= args.min_expected_shared
        )
        rows.append(
            {
                "source": args.source_id,
                "control": control_id,
                "source_reads": source_reads,
                "control_reads": control_reads,
                "control_distinct": len(counts),
                "shared_reads": shared_reads,
                "shared_distinct": len(shared),
                "control_share": round(control_share, 6),
                "probe_reads": len(probe),
                "probe_hits": probe_hits,
                "self_share": round(self_share, 6),
                "est_fraction": round(est_fraction, 4),
                "expected_at_flag": round(expected_at_flag, 1),
                "verdict": verdict(est_fraction, args.flag_fraction, powered),
            }
        )
        ranked = sorted(shared.items(), key=lambda kv: kv[1], reverse=True)
        for fragment, n in ranked[: args.top_shared]:
            top.append([args.source_id, control_id, n, source_counts[fragment], fragment.decode()])

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with open(args.out_top, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["source", "control", "control_reads", "source_reads", "sequence"])
        writer.writerows(top)


if __name__ == "__main__":
    main()
