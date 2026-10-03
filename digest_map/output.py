"""CSV and JSON output writers."""

from __future__ import annotations

import csv
import json
import os
from collections import Counter

from .core import DigestResult


def fragments_csv_rows(result: DigestResult) -> list[list]:
    rows = []
    for index, frag in enumerate(result.fragments):
        rows.append([index, frag.start, frag.length, frag.sequence,
                     "yes" if frag.wraps_origin else "no"])
    return rows


CSV_HEADER = ["fragment_index", "start", "length", "sequence", "wraps_origin"]


def result_to_summary(result: DigestResult) -> dict:
    histogram = Counter(f.length for f in result.fragments)
    return {
        "digest_id": result.digest_id,
        "dna_id": result.dna_id,
        "topology": result.topology,
        "dna_length": result.dna_length,
        "cut_sites": [
            {"position": site.position,
             "enzymes": list(site.enzymes),
             "match_starts": list(site.match_starts)}
            for site in result.cut_sites
        ],
        "fragment_count": len(result.fragments),
        "total_fragment_length": sum(f.length for f in result.fragments),
        "length_histogram": {str(k): histogram[k] for k in sorted(histogram)},
        "fragments_file": f"fragments_{result.digest_id}.csv",
    }


def write_outputs(results: list[DigestResult], out_dir: str) -> None:
    """Write all outputs atomically-ish: build in a temp dir, then move into place."""
    tmp_dir = out_dir + ".tmp"
    if os.path.exists(tmp_dir):
        _rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    for result in results:
        csv_path = os.path.join(tmp_dir, f"fragments_{result.digest_id}.csv")
        with open(csv_path, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(CSV_HEADER)
            writer.writerows(fragments_csv_rows(result))
    summary = {"digests": [result_to_summary(r) for r in results]}
    with open(os.path.join(tmp_dir, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    if os.path.exists(out_dir):
        _rmtree(out_dir)
    os.replace(tmp_dir, out_dir)


def _rmtree(path: str) -> None:
    for root, dirs, files in os.walk(path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(path)
