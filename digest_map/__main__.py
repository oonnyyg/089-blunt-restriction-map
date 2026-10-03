"""Command line interface: python -m digest_map input.json --out output_dir"""

from __future__ import annotations

import argparse
import sys

from .core import run_digest
from .errors import InputError
from .output import write_outputs
from .validate import load_task


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m digest_map",
        description="Predict blunt-end restriction fragments for linear and "
                    "circular DNA.")
    parser.add_argument("input", help="task JSON file")
    parser.add_argument("--out", required=True, help="output directory")
    args = parser.parse_args(argv)

    try:
        dnas, enzymes, digests = load_task(args.input)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"error: cannot read input: {exc}", file=sys.stderr)
        return 1

    dna_by_id = {d.id: d for d in dnas}
    enzyme_by_id = {e.id: e for e in enzymes}
    results = [run_digest(digest.id, dna_by_id[digest.dna_id],
                          [enzyme_by_id[e] for e in digest.enzyme_ids])
               for digest in digests]
    write_outputs(results, args.out)
    print(f"wrote {len(results)} digest(s) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
