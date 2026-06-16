#!/usr/bin/env python3
"""Benchmark final GeneMapKit symbol mappings against multiple live APIs."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

from genemapkit.core.converter import MAPPING_POLICIES
from genemapkit.core.external_benchmark import (
    DEFAULT_OUTPUTS,
    DEFAULT_PROVIDERS,
    PROVIDER_OUTPUTS,
    benchmark_external_services,
)


def _print_table(headers, rows) -> None:
    values = [[str(value) for value in row] for row in rows]
    widths = [
        max(len(str(header)), *(len(row[index]) for row in values))
        for index, header in enumerate(headers)
    ]
    print("  ".join(str(header).ljust(widths[index]) for index, header in enumerate(headers)))
    print("  ".join("-" * width for width in widths))
    for row in values:
        print("  ".join(value.ljust(widths[index]) for index, value in enumerate(row)))


def _pct(value) -> str:
    return "NA" if value is None else f"{value}%"


def _pct_ci(row, prefix) -> str:
    value = row[f"{prefix}_pct"]
    if value is None:
        return "NA"
    return (
        f"{value}% [{row[f'{prefix}_ci_low_pct']}, "
        f"{row[f'{prefix}_ci_high_pct']}]"
    )


def main(argv=None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(
        description=(
            "Compare GeneMapKit symbol mappings with live MyGene, Ensembl, NCBI, "
            "UniProt, g:Profiler, and BridgeDb services. Results measure "
            "concordance, not accuracy."
        )
    )
    parser.add_argument("input_file", help="CSV/TSV file containing symbols")
    parser.add_argument("--column", default="symbol", help="Symbol column name")
    parser.add_argument("--delimiter", default=",", help=r"Input delimiter; use \t for TSV")
    parser.add_argument(
        "--db", default="data/build/genemapkit_final_v21.db", help="Final SQLite database"
    )
    parser.add_argument(
        "--provider",
        action="append",
        choices=sorted(PROVIDER_OUTPUTS),
        help="External service; repeatable (default: all)",
    )
    parser.add_argument(
        "--output-type",
        action="append",
        choices=sorted(set().union(*PROVIDER_OUTPUTS.values())),
        help="Output namespace; repeatable (default: common symbol outputs)",
    )
    parser.add_argument(
        "--policy", default="all-supported", choices=MAPPING_POLICIES
    )
    parser.add_argument(
        "--out", default="results/external_benchmark", help="Output report directory"
    )
    parser.add_argument(
        "--refresh-cache", action="store_true", help="Ignore cached API results"
    )
    parser.add_argument(
        "--no-progress", action="store_false", dest="progress", default=True,
        help="disable progress bars",
    )
    args = parser.parse_args(argv)
    delimiter = "\t" if args.delimiter == r"\t" else args.delimiter
    frame = pd.read_csv(args.input_file, delimiter=delimiter, dtype=str)
    if args.column not in frame:
        parser.error(
            f"Column '{args.column}' not found. Available columns: {', '.join(frame.columns)}"
        )
    report = benchmark_external_services(
        Path(args.db),
        frame[args.column].dropna(),
        args.output_type or DEFAULT_OUTPUTS,
        providers=args.provider or DEFAULT_PROVIDERS,
        mapping_policy=args.policy,
        output_dir=Path(args.out),
        refresh_cache=args.refresh_cache,
        progress=args.progress,
    )
    print("\nCombined external-union concordance")
    _print_table(
        [
            "Output",
            "Exact",
            "Supported+extras",
            "Partial",
            "Conflict",
            "GMK-only",
            "External-only",
            "Precision (95% CI)",
            "Recall (95% CI)",
        ],
        [
            [
                row["output_type"],
                row["exact_external_union"],
                row["supported_with_external_extras"],
                row["partial_support"],
                row["conflict"],
                row["genemapkit_only"],
                row["external_only"],
                _pct_ci(row, "concordance_precision"),
                _pct_ci(row, "concordance_recall"),
            ]
            for row in report["summaries"]
        ],
    )
    print("\nPer-provider concordance")
    _print_table(
        [
            "Provider",
            "Output",
            "GMK",
            "External",
            "Common",
            "Precision (95% CI)",
            "Recall (95% CI)",
        ],
        [
            [
                row["provider"],
                row["output_type"],
                row["genemapkit_mappings"],
                row["provider_mappings"],
                row["common_mappings"],
                _pct_ci(row, "concordance_precision"),
                _pct_ci(row, "concordance_recall"),
            ]
            for row in report["provider_summaries"]
        ],
    )
    if report["provider_errors"]:
        print("\nProvider errors")
        for provider, error in report["provider_errors"].items():
            print(f"- {provider}: {error}")
    print("\nImportant: external API agreement measures concordance, not independent accuracy.")
    print(f"Reports and API cache: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
