#!/usr/bin/env python3
"""Download reproducible raw database snapshots for GeneMapKit."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from genemapkit.core.downloader import (
    DEFAULT_ENSEMBL_RELEASE,
    TIERS,
    VALIDATION_SOURCES,
    DatabaseDownloader,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "GeneMapKit Stage 0: download immutable human gene-identifier "
            "source files and record checksums in a manifest."
        )
    )
    parser.add_argument("--out", default="data/raw", help="raw snapshot directory")
    parser.add_argument(
        "--ensembl-release",
        type=int,
        default=DEFAULT_ENSEMBL_RELEASE,
        help="pinned Ensembl release (default: %(default)s)",
    )
    parser.add_argument(
        "--core-only", action="store_true", help="download only essential core sources"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="download every registered source, including large/niche sources",
    )
    parser.add_argument(
        "--validation",
        action="store_true",
        help=(
            "also download independent benchmark-validation sources "
            "(NCBI gene2ensembl and large gene2refseq)"
        ),
    )
    parser.add_argument(
        "--sources",
        nargs="+",
        metavar="SOURCE",
        help="download only explicitly named sources",
    )
    parser.add_argument(
        "--for-conversion",
        nargs="+",
        metavar=("INPUT_TYPE", "OUTPUT_TYPE"),
        help="download sources required for INPUT_TYPE to OUTPUT_TYPE(s)",
    )
    parser.add_argument(
        "--minimal",
        action="store_true",
        help="with --for-conversion, select the smallest direct source set",
    )
    parser.add_argument("--list", action="store_true", help="list sources and exit")
    parser.add_argument("--coverage", action="store_true", help="print ID coverage and exit")
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    return parser


def print_sources(downloader: DatabaseDownloader) -> None:
    for tier in TIERS:
        print(f"\n{tier.upper()}")
        for source in downloader.sources_for_tiers((tier,)):
            info = downloader.resolve_source(source)
            flags = []
            if source in VALIDATION_SOURCES:
                flags.append("validation")
            if info.get("large"):
                flags.append("large")
            if info.get("taxon_filter"):
                flags.append("filter taxid 9606 during build")
            suffix = f" [{'; '.join(flags)}]" if flags else ""
            print(f"  {source:20} {info['description']}{suffix}")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    downloader = DatabaseDownloader(
        data_dir=args.out,
        ensembl_release=args.ensembl_release,
    )

    if args.list:
        print_sources(downloader)
        return 0
    if args.coverage:
        downloader.print_coverage()
        return 0
    if args.for_conversion and len(args.for_conversion) < 2:
        build_parser().error(
            "--for-conversion requires one input type and at least one output type"
        )
    if args.minimal and not args.for_conversion:
        build_parser().error("--minimal must be used with --for-conversion")
    if args.core_only and (args.all or args.validation or args.sources or args.for_conversion):
        build_parser().error(
            "--core-only cannot be combined with --all, --validation, --sources, "
            "or --for-conversion"
        )
    if args.sources and (args.all or args.validation or args.for_conversion):
        build_parser().error(
            "--sources cannot be combined with --all, --validation, or --for-conversion"
        )
    if args.for_conversion and (args.all or args.validation):
        build_parser().error(
            "--for-conversion cannot be combined with --all or --validation"
        )

    if args.sources:
        results = downloader.download_specific_sources(args.sources, args.force)
    elif args.for_conversion:
        results = downloader.download_for_conversion(
            args.for_conversion[0],
            args.for_conversion[1:],
            force_download=args.force,
            minimal=args.minimal,
        )
    elif args.core_only:
        results = downloader.download_specific_sources(
            downloader.sources_for_tiers(("core",)), args.force
        )
    elif args.all:
        # --validation is intentionally accepted with --all. All registered
        # sources already include validation sources.
        results = downloader.download_all(args.force, include_optional=True)
    elif args.validation:
        results = downloader.download_validation(args.force, include_default=True)
    else:
        # Sensible default: core plus cross-reference validation, excluding
        # all-species multi-gigabyte and niche sources.
        results = downloader.download_default(args.force)

    successful = sum(results.values())
    print(f"\nAvailable sources: {successful}/{len(results)}")
    print(f"Manifest: {Path(args.out) / 'manifest.json'}")
    return 0 if successful == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
