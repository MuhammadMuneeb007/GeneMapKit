#!/usr/bin/env python3
"""Run the complete GeneMapKit publication benchmark suite."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from genemapkit.core.converter import MAPPING_POLICIES
from genemapkit.core.external_benchmark import DEFAULT_OUTPUTS, DEFAULT_PROVIDERS
from genemapkit.core.path_benchmark import DEFAULT_PATH_ID_TYPES
from genemapkit.core.publication_benchmark import run_publication_benchmark


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description=(
            "Run local performance, coordinate concordance, multi-API concordance, "
            "and held-out validation aggregation."
        )
    )
    parser.add_argument("input_file", help="CSV/TSV file containing symbols")
    parser.add_argument("--column", default="symbol")
    parser.add_argument("--delimiter", default=",", help=r"Use \t for TSV")
    parser.add_argument("--db", default="data/build/genemapkit_final_v21.db")
    parser.add_argument("--gtf", default="data/raw/Homo_sapiens.GRCh38.116.gtf.gz")
    parser.add_argument(
        "--held-out-report",
        action="append",
        default=[],
        help="Directory containing validation_metrics.json; repeatable",
    )
    parser.add_argument(
        "--batch-size", action="append", type=int, help="Local scaling size; repeatable"
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--policy", default="all-supported", choices=MAPPING_POLICIES)
    parser.add_argument("--out", default="results/publication_benchmark")
    parser.add_argument("--refresh-external-cache", action="store_true")
    parser.add_argument(
        "--coordinate-batch-size",
        type=int,
        default=1000,
        help="Number of symbols per Ensembl coordinate API request",
    )
    parser.add_argument(
        "--path-pivot-type", default="entrez_id", choices=DEFAULT_PATH_ID_TYPES
    )
    parser.add_argument(
        "--path-sample",
        type=int,
        default=1000,
        help="Identifiers per namespace for round-trip/path analysis; 0 means complete",
    )
    parser.add_argument("--path-seed", type=int, default=20260616)
    parser.add_argument(
        "--no-progress", action="store_false", dest="progress", default=True,
        help="disable progress bars",
    )
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error("--repeats must be a positive integer")
    if args.path_sample < 0:
        parser.error("--path-sample must be zero or a positive integer")
    if args.batch_size and any(value < 1 for value in args.batch_size):
        parser.error("--batch-size values must be positive integers")
    if args.coordinate_batch_size < 1:
        parser.error("--coordinate-batch-size must be a positive integer")
    delimiter = "\t" if args.delimiter == r"\t" else args.delimiter
    frame = pd.read_csv(args.input_file, delimiter=delimiter, dtype=str)
    if args.column not in frame:
        parser.error(
            f"Column '{args.column}' not found. Available columns: {', '.join(frame.columns)}"
        )
    held_out = args.held_out_report or [
        "results/mane_heldout_all_supported",
        "results/ncbi_heldout_all_supported",
    ]
    report = run_publication_benchmark(
        Path(args.db),
        Path(args.gtf),
        frame[args.column].dropna(),
        output_types=DEFAULT_OUTPUTS,
        providers=DEFAULT_PROVIDERS,
        held_out_report_dirs=[Path(value) for value in held_out],
        batch_sizes=args.batch_size or [1, 10, 100, 1000],
        repeats=args.repeats,
        mapping_policy=args.policy,
        output_dir=Path(args.out),
        refresh_external_cache=args.refresh_external_cache,
        coordinate_batch_size=args.coordinate_batch_size,
        path_pivot_type=args.path_pivot_type,
        path_sample=None if args.path_sample == 0 else args.path_sample,
        path_random_seed=args.path_seed,
        progress=args.progress,
    )
    print(f"Publication benchmark: {args.out}/publication_benchmark_summary.md")
    print(
        "Offline local conversion:",
        report["local_performance"]["offline_operation_passed"],
    )
    print(
        "Coordinate exact matches:",
        report["coordinate_concordance"]["summary"]["exact"],
        "/",
        report["coordinate_concordance"]["summary"]["tested_symbols"],
    )
    print(
        "Independent held-out rows:",
        len(report["held_out_validation"]["independent_benchmarks"]),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
