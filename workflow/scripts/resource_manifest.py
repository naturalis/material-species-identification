"""Shared helpers for reference acquisition: downloads and resources/MANIFEST.tsv rows."""
import csv
import gzip
import hashlib
import http.client
import platform
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

COLUMNS = ["file", "source", "downloaded_utc", "tool", "bytes", "records", "sha256"]
TOOL = f"python {platform.python_version()} urllib"
RETRIES = 8
# NCBI drops TLS connections now and then, sometimes for a minute or more. OSError covers
# URLError, timeouts and ssl.SSLError; HTTPException covers a body cut short.
TRANSIENT = (OSError, http.client.HTTPException)
USER_AGENT = "material-species-identification reference acquisition"


def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def backoff(url, error, attempt):
    print(f"retrying {url}: {error}", file=sys.stderr)
    time.sleep(15 * (attempt + 1))


def fetch(url):
    """Return the bytes at url, backing off on transient failures."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                return response.read()
        except TRANSIENT as error:
            if attempt == RETRIES - 1:
                raise
            backoff(url, error, attempt)


def download(url, path):
    """Stream url to path, return the UTC time the download finished."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    partial = Path(f"{path}.partial")
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(request, timeout=300) as response, open(partial, "wb") as fh:
                while chunk := response.read(1 << 20):
                    fh.write(chunk)
            partial.rename(path)
            return utc_now()
        except TRANSIENT as error:
            if attempt == RETRIES - 1:
                raise
            backoff(url, error, attempt)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def count_fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith(">"))


def count_genbank_records(path):
    with gzip.open(path, "rt") as fh:
        return sum(1 for line in fh if line.startswith("LOCUS"))


def row(path, source, downloaded, tool, records):
    path = Path(path)
    return {
        "file": str(path),
        "source": source,
        "downloaded_utc": downloaded,
        "tool": tool,
        "bytes": path.stat().st_size,
        "records": records,
        "sha256": sha256(path),
    }


def write_rows(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
