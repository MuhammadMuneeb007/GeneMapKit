#!/usr/bin/env python3
"""Build the GeneMapKit SQLite identifier index from a raw snapshot."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from genemapkit.core.builder import PARSERS, build_database, convert_identifier


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "GeneMapKit Stage 1: verify a raw snapshot and build a normalized, "
            "provenance-aware SQLite identifier index."
        )
    )
    parser.add_argument(
        "--raw-dir",
        "--raw",
        dest="raw_dir",
        default="data/raw",
        help="snapshot directory containing raw files and manifest.json",
    )
    parser.add_argument(
        "--out",
        default="data/build/genemapkit.db",
        help="output SQLite database path",
    )
    parser.add_argument(
        "--report",
        default="data/build/build_report.json",
        help="output JSON build report path",
    )
    parser.add_argument(
        "--source",
        action="append",
        choices=sorted(PARSERS),
        help="build from this source; repeat to select multiple sources",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        choices=sorted(PARSERS),
        help="exclude this source; repeat to exclude multiple sources",
    )
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="skip manifest checksum verification; not recommended for publication",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="continue when selected sources are missing or empty",
    )
    parser.add_argument("--force", action="store_true", help="replace an existing output")
    parser.add_argument(
        "--max-records",
        type=int,
        help="development-only record limit per source; creates an incomplete index",
    )
    parser.add_argument(
        "--check",
        nargs=3,
        metavar=("VALUE", "INPUT_TYPE", "OUTPUT_TYPE"),
        help="run one conversion check after building",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="suppress per-source tqdm progress bars",
    )
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    try:
        report = build_database(
            Path(args.raw_dir),
            Path(args.out),
            sources=args.source,
            exclude=args.exclude,
            verify=not args.no_verify,
            strict=not args.allow_missing,
            force=args.force,
            max_records=args.max_records,
            report_path=Path(args.report),
            progress=not args.no_progress,
        )
    except Exception as error:
        logging.error("Build failed: %s", error)
        return 1

    print(json.dumps(report, indent=2))
    if args.check:
        value, input_type, output_type = args.check
        results = convert_identifier(Path(args.out), value, input_type, output_type)
        print(f"\n{value} ({input_type}) -> {output_type}:")
        print("\n".join(results) if results else "(no mappings)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
