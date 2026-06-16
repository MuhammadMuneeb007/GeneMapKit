#!/usr/bin/env python3
"""Validate GeneMapKit conversion relationships against raw truth sources."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from genemapkit.core.converter import MAPPING_POLICIES
from genemapkit.core.validation import VALIDATION_CONVERSIONS, validate_database


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description="Benchmark a built GeneMapKit SQLite database against validation sources."
    )
    parser.add_argument("--db", default="data/build/genemapkit.db", help="SQLite database")
    parser.add_argument("--raw-dir", default="data/raw", help="Raw snapshot directory")
    parser.add_argument(
        "--out", default="results/validation", help="Validation report directory"
    )
    parser.add_argument(
        "--source",
        action="append",
        choices=sorted(VALIDATION_CONVERSIONS),
        help="Truth source to evaluate; repeat to select multiple",
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=1000,
        help="Maximum unique truth inputs per conversion; use 0 for complete sources",
    )
    parser.add_argument(
        "--max-examples", type=int, default=25, help="Disagreement examples per conversion"
    )
    parser.add_argument(
        "--require-held-out",
        action="store_true",
        help="Skip truth sources that were included in the database build",
    )
    parser.add_argument(
        "--policy",
        default="all-supported",
        choices=MAPPING_POLICIES,
        help="Mapping-selection policy to benchmark",
    )
    parser.add_argument(
        "--no-progress", action="store_false", dest="progress", default=True,
        help="disable progress bars",
    )
    args = parser.parse_args(argv)
    if args.sample < 0:
        parser.error("--sample must be zero or a positive integer")
    report = validate_database(
        Path(args.db),
        Path(args.raw_dir),
        sources=args.source,
        sample=None if args.sample == 0 else args.sample,
        max_examples=args.max_examples,
        require_held_out=args.require_held_out,
        mapping_policy=args.policy,
        output_dir=Path(args.out),
        progress=args.progress,
    )
    print(json.dumps({"diagnostics": report["diagnostics"], "benchmarks": report["benchmarks"]}, indent=2))
    print(f"\nReports: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
