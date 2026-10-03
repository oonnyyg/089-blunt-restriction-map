"""CLI: python -m digest_map 输入.json --out 输出目录"""

from __future__ import annotations

import argparse
import sys

from .core import DigestError, load_task, run_task


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m digest_map",
        description="Predict blunt-end restriction fragments for linear and circular DNA.")
    parser.add_argument("task", help="task JSON file (dnas, enzymes, digests)")
    parser.add_argument("--out", required=True, help="output directory")
    args = parser.parse_args(argv)

    try:
        task = load_task(args.task)
        summary = run_task(task, args.out)
    except DigestError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"error: cannot write output: {exc}", file=sys.stderr)
        return 1

    count = len(summary["digests"])
    print(f"wrote {count} digest(s) to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
