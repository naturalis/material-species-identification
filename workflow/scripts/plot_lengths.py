#!/usr/bin/env python
"""Fragment-length distribution of merged pairs: one panel per library, never pooled."""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROLE_COLOURS = {"sample": "#2a78d6", "control": "#eb6834"}
SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e4e3df"
COLUMNS = 4


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lengths", required=True)
    parser.add_argument("--read-length", type=int, required=True)
    parser.add_argument("--min-overlap", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.lengths, sep="\t")
    # Samples first, then controls, each in sheet order.
    libraries = table[["id", "role"]].drop_duplicates()
    libraries = pd.concat(libraries[libraries["role"] == r] for r in ROLE_COLOURS)
    ceiling = 2 * args.read_length - args.min_overlap

    columns = min(COLUMNS, len(libraries))
    rows = -(-len(libraries) // columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(3.2 * columns + 0.6, 2.5 * rows + 0.9), sharex=True, sharey=True,
        squeeze=False, facecolor=SURFACE,
    )
    axes = axes.flatten()

    for ax, library in zip(axes, libraries.itertuples()):
        ax.set_facecolor(SURFACE)
        data = table[table["id"] == library.id].sort_values("length")
        total = data["reads"].sum()
        if total:
            share = data["reads"] / total
            ax.plot(data["length"], share, color=ROLE_COLOURS[library.role], linewidth=1.5,
                    solid_capstyle="round", solid_joinstyle="round")
            mode = data.loc[data["reads"].idxmax(), "length"]
            ax.text(0.88, 0.92, f"mode {mode} bp\n{total:,} merged", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8, color=MUTED)
        else:
            ax.text(0.5, 0.5, "no merged reads", transform=ax.transAxes, ha="center",
                    fontsize=8, color=MUTED)
        ax.set_title(f"{library.id} ({library.role})", loc="left", fontsize=10, color=TEXT)
        ax.axvline(ceiling, color=MUTED, linewidth=0.8, linestyle=(0, (2, 3)))
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=8, length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(GRID)

    for ax in axes[len(libraries):]:
        ax.set_visible(False)

    fig.supxlabel("Merged fragment length (bp)", fontsize=9, color=MUTED)
    fig.supylabel("Share of the library's merged reads", fontsize=9, color=MUTED)
    fig.suptitle("Fragment lengths by library", x=0.01, y=0.995, ha="left", fontsize=12, color=TEXT)
    handles = [plt.Line2D([], [], color=c, linewidth=2, label=r) for r, c in ROLE_COLOURS.items()]
    handles.append(plt.Line2D([], [], color=MUTED, linewidth=0.8, linestyle=(0, (2, 3)),
                              label=f"merge ceiling, {ceiling} bp"))
    fig.legend(handles=handles, loc="upper right", ncol=3, frameon=False, fontsize=8,
               labelcolor=MUTED)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)


if __name__ == "__main__":
    main()
