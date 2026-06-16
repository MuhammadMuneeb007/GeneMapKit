#!/usr/bin/env python3
"""Run GeneMapKit all-to-all and conversion-path consistency benchmarks."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from genemapkit.core.path_benchmark import (
    DEFAULT_PATH_ID_TYPES,
    benchmark_conversion_paths,
)


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description=(
            "Measure direct all-to-all conversion coverage, round-trip recovery, "
            "and direct-versus-pivot path consistency."
        )
    )
    parser.add_argument("--db", default="data/build/genemapkit_final_v21.db")
    parser.add_argument(
        "--id-type",
        action="append",
        choices=DEFAULT_PATH_ID_TYPES,
        help="Namespace to include; repeatable",
    )
    parser.add_argument(
        "--pivot-type", default="entrez_id", choices=DEFAULT_PATH_ID_TYPES
    )
    parser.add_argument(
        "--sample",
        type=int,
        default=1000,
        help="Identifiers per source namespace for path tests; 0 means complete",
    )
    parser.add_argument("--max-examples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260616)
    parser.add_argument("--out", default="results/path_benchmark")
    parser.add_argument(
        "--no-progress", action="store_false", dest="progress", default=True,
        help="disable progress bars",
    )
    args = parser.parse_args(argv)
    if args.sample < 0:
        parser.error("--sample must be zero or a positive integer")
    if args.max_examples < 0:
        parser.error("--max-examples must be zero or a positive integer")
    report = benchmark_conversion_paths(
        Path(args.db),
        id_types=args.id_type or DEFAULT_PATH_ID_TYPES,
        pivot_type=args.pivot_type,
        path_sample=None if args.sample == 0 else args.sample,
        max_examples=args.max_examples,
        random_seed=args.seed,
        output_dir=Path(args.out),
        progress=args.progress,
    )
    print(f"Path benchmark runtime: {report['elapsed_seconds']} seconds")
    print(f"Pivot namespace: {report['pivot_type']}")
    print(f"Coverage rows: {len(report['conversion_coverage'])}")
    print(f"Round-trip rows: {len(report['round_trip_consistency'])}")
    print(f"Pivot-path rows: {len(report['pivot_path_consistency'])}")
    print(f"Reports: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
