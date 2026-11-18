#!/usr/bin/env bash
# Cuts the fixture references from one RefSeq mitochondrion release.
# Usage, from the repo root: tests/fixture/build_references.sh path/to/mitochondrion.1.1.genomic.fna.gz
#
# Panel: 5 taxa including the Bos congener pair, so the ambiguous bin has something to do.
# Subset: the panel plus 45 random mitogenomes under 100 kb (plants would dominate the bases).
set -euo pipefail

RELEASE_FASTA=$1
SEED=20261030
OUT=tests/fixture/references
PANEL="NC_012920.1 NC_006853.1 NC_005971.1 NC_000845.1 NC_001941.1"

mkdir -p "$OUT"
printf '%s\n' $PANEL > "$OUT/panel_accessions.txt"
seqkit grep -n -r -f <(sed 's/^/^/; s/$/ /' "$OUT/panel_accessions.txt") "$RELEASE_FASTA" \
    -w 0 -o "$OUT/panel.fasta"

seqkit grep -n -v -r -f <(sed 's/^/^/; s/$/ /' "$OUT/panel_accessions.txt") "$RELEASE_FASTA" -w 0 \
    | seqkit seq -M 100000 -w 0 \
    | seqkit shuffle -s "$SEED" -w 0 -o "$OUT/background_shuffled.fasta" 2> /dev/null
# Not piped into head: head closing the pipe would kill the shuffle under pipefail.
seqkit head -n 45 -w 0 "$OUT/background_shuffled.fasta" -o "$OUT/background.fasta"
cat "$OUT/panel.fasta" "$OUT/background.fasta" > "$OUT/mito_subset.fasta"
rm "$OUT/background.fasta" "$OUT/background_shuffled.fasta"

{
    printf 'accession\ttaxon\tin_panel\n'
    seqkit seq -n "$OUT/mito_subset.fasta" | awk -v p="$PANEL" '
        BEGIN { n = split(p, a, " "); for (i = 1; i <= n; i++) panel[a[i]] = 1 }
        { acc = $1; sub(/^[^ ]+ /, ""); sub(/ mitochondrion.*/, ""); print acc "\t" $0 "\t" (acc in panel ? "yes" : "no") }'
} > "$OUT/mito_subset.tsv"
