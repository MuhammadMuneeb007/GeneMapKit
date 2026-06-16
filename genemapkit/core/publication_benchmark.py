"""Publication-oriented benchmark suite for GeneMapKit."""

from __future__ import annotations

import csv
import json
import logging
import os
import platform
import statistics
import threading
import time
import tracemalloc
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence
from unittest.mock import patch

import requests
import psutil

from .builder import GTF_ATTRIBUTE, normalize_identifier, open_text
from .converter import GeneConverter
from .external_benchmark import (
    DEFAULT_OUTPUTS,
    DEFAULT_PROVIDERS,
    benchmark_external_services,
)
from .path_benchmark import DEFAULT_PATH_ID_TYPES, benchmark_conversion_paths
from .statistics import interval_fields
from genemapkit.utils.progress import track

LOGGER = logging.getLogger("genemapkit")


def _percent(numerator: int, denominator: int) -> Optional[float]:
    return round(100.0 * numerator / denominator, 4) if denominator else None


def _display_pct(value: object) -> str:
    return "NA" if value is None else f"{value}%"


def _display_pct_ci(row: Mapping[str, object], prefix: str) -> str:
    value = row.get(f"{prefix}_pct")
    low = row.get(f"{prefix}_ci_low_pct")
    high = row.get(f"{prefix}_ci_high_pct")
    return "NA" if value is None else f"{value}% [{low}, {high}]"


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def benchmark_local_performance(
    database: Path,
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_sizes: Sequence[int] = (1, 10, 100, 1000),
    repeats: int = 3,
    mapping_policy: str = "all-supported",
    progress: bool = False,
) -> Dict[str, object]:
    """Measure local converter startup, throughput, Python memory, and offline use."""
    unique = list(dict.fromkeys(symbols))
    if not unique:
        raise ValueError("No symbols provided for local performance benchmark")
    LOGGER.info(
        "Local performance: %d unique symbols, %d output types, batch sizes=%s, repeats=%d",
        len(unique),
        len(output_types),
        ",".join(str(value) for value in batch_sizes),
        repeats,
    )
    startup_times = []
    for _ in track(
        range(repeats),
        enabled=progress,
        description="Local startup timing",
        total=repeats,
        unit="runs",
    ):
        started = time.perf_counter()
        GeneConverter(database, mapping_policy=mapping_policy)
        startup_times.append(time.perf_counter() - started)
    converter = GeneConverter(database, mapping_policy=mapping_policy)
    rows: List[Dict[str, object]] = []
    for size in track(
        batch_sizes,
        enabled=progress,
        description="Local performance batch sizes",
        total=len(batch_sizes),
        unit="sizes",
    ):
        LOGGER.info("Local performance: benchmarking batch size %d", size)
        values = [unique[index % len(unique)] for index in range(size)]
        timings = []
        peaks = []
        rss_peaks = []
        rss_deltas = []
        mapped_rows = 0
        for repeat_index in track(
            range(repeats),
            enabled=progress,
            description=f"Local batch {size}",
            total=repeats,
            unit="runs",
        ):
            LOGGER.info(
                "Local performance: batch size %d repeat %d/%d",
                size,
                repeat_index + 1,
                repeats,
            )
            process = psutil.Process(os.getpid())
            rss_before = process.memory_info().rss
            rss_samples = [rss_before]
            stop_sampling = threading.Event()

            def sample_rss() -> None:
                while not stop_sampling.wait(0.01):
                    rss_samples.append(process.memory_info().rss)

            sampler = threading.Thread(target=sample_rss, daemon=True)
            sampler.start()
            tracemalloc.start()
            started = time.perf_counter()
            frame = converter.convert_batch(
                values,
                "symbol",
                output_types,
                output_format="wide",
                mapping_policy=mapping_policy,
                progress=progress,
            )
            elapsed = time.perf_counter() - started
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            stop_sampling.set()
            sampler.join()
            rss_after = process.memory_info().rss
            rss_samples.append(rss_after)
            timings.append(elapsed)
            peaks.append(peak)
            rss_peaks.append(max(rss_samples))
            rss_deltas.append(max(rss_samples) - rss_before)
            mapped_rows = int(
                frame[[value for value in output_types if value in frame]]
                .ne("")
                .any(axis=1)
                .sum()
            )
            LOGGER.info(
                "Local performance: batch size %d repeat %d/%d finished in %.4f s",
                size,
                repeat_index + 1,
                repeats,
                elapsed,
            )
        median = statistics.median(timings)
        LOGGER.info(
            "Local performance: batch size %d median %.4f s, %.2f identifiers/s",
            size,
            median,
            size / median if median else 0,
        )
        rows.append(
            {
                "batch_size": size,
                "repeats": repeats,
                "median_seconds": round(median, 6),
                "min_seconds": round(min(timings), 6),
                "max_seconds": round(max(timings), 6),
                "identifiers_per_second": round(size / median, 4) if median else None,
                "median_peak_python_memory_mb": round(
                    statistics.median(peaks) / (1024 * 1024), 4
                ),
                "median_peak_process_rss_mb": round(
                    statistics.median(rss_peaks) / (1024 * 1024), 4
                ),
                "median_process_rss_increase_mb": round(
                    statistics.median(rss_deltas) / (1024 * 1024), 4
                ),
                "mapped_rows": mapped_rows,
            }
        )

    def deny_network(*args, **kwargs):
        raise AssertionError("network access attempted during local conversion")

    offline_passed = False
    offline_error = None
    try:
        with patch("requests.sessions.Session.request", side_effect=deny_network):
            converter.convert_single(
                unique[0], "symbol", output_types, mapping_policy=mapping_policy
            )
        offline_passed = True
    except Exception as error:  # benchmark must record failures, not hide them
        offline_error = str(error)

    return {
        "database": str(database),
        "database_size_bytes": database.stat().st_size,
        "database_size_gb": round(database.stat().st_size / (1024**3), 4),
        "mapping_policy": mapping_policy,
        "output_types": list(output_types),
        "startup_repeats": repeats,
        "median_startup_seconds": round(statistics.median(startup_times), 6),
        "offline_operation_passed": offline_passed,
        "offline_operation_error": offline_error,
        "memory_metric": (
            "Process RSS sampled with psutil plus tracemalloc peak Python allocations"
        ),
        "scaling": rows,
    }


def build_system_comparison(
    performance: Mapping[str, object],
    coordinates: Mapping[str, object],
    external: Mapping[str, object],
    tested_symbols: int,
) -> List[Dict[str, object]]:
    """Create a direct local-versus-service runtime and capability table."""
    local = next(
        (
            row
            for row in performance["scaling"]
            if row["batch_size"] == tested_symbols
        ),
        performance["scaling"][-1],
    )
    rows = [
        {
            "system": "GeneMapKit SQLite",
            "operation": "five-namespace identifier conversion",
            "symbols": local["batch_size"],
            "elapsed_seconds": local["median_seconds"],
            "symbols_per_second": local["identifiers_per_second"],
            "network_required": False,
            "offline_capable": performance["offline_operation_passed"],
            "process_peak_rss_mb": local["median_peak_process_rss_mb"],
            "notes": "Local final database; median across repeats",
        },
        {
            "system": "GeneMapKit pinned GTF",
            "operation": "coordinate lookup by full GTF scan",
            "symbols": tested_symbols,
            "elapsed_seconds": coordinates["timings"]["local_gtf_scan_seconds"],
            "symbols_per_second": round(
                tested_symbols / coordinates["timings"]["local_gtf_scan_seconds"], 4
            ),
            "network_required": False,
            "offline_capable": True,
            "process_peak_rss_mb": None,
            "notes": "Coordinates are not indexed in SQLite schema 2.1.0",
        },
    ]
    for timing in external["provider_timings"]:
        elapsed = timing["elapsed_seconds"]
        rows.append(
            {
                "system": timing["provider"],
                "operation": "external identifier API query",
                "symbols": timing["queried_symbols"],
                "elapsed_seconds": elapsed,
                "symbols_per_second": (
                    round(timing["queried_symbols"] / elapsed, 4)
                    if elapsed and timing["queried_symbols"]
                    else None
                ),
                "network_required": True,
                "offline_capable": False,
                "process_peak_rss_mb": None,
                "notes": (
                    "Live API query"
                    if timing["queried_symbols"]
                    else "Cached response; not a live timing"
                ),
            }
        )
    return rows


def load_gtf_gene_coordinates(
    gtf_path: Path, symbols: Sequence[str]
) -> Dict[str, List[Dict[str, object]]]:
    """Load requested gene coordinates from the pinned local Ensembl GTF."""
    LOGGER.info(
        "GTF coordinate scan: looking up %d symbols in %s", len(symbols), gtf_path
    )
    wanted = {symbol.casefold(): symbol for symbol in symbols}
    results: Dict[str, List[Dict[str, object]]] = {symbol: [] for symbol in symbols}
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attributes = dict(GTF_ATTRIBUTE.findall(fields[8]))
            symbol = attributes.get("gene_name", "")
            requested = wanted.get(symbol.casefold())
            if not requested:
                continue
            results[requested].append(
                {
                    "ensembl_gene_id": normalize_identifier(
                        "ensembl_gene_id", attributes.get("gene_id", "")
                    ),
                    "assembly": "GRCh38",
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                }
            )
    found = sum(1 for value in results.values() if value)
    LOGGER.info(
        "GTF coordinate scan: found coordinates for %d/%d requested symbols",
        found, len(symbols),
    )
    return results


def _request_with_retries(
    client: requests.Session,
    method: str,
    url: str,
    *,
    max_retries: int = 5,
    **kwargs,
) -> requests.Response:
    """Issue an HTTP request with bounded retry/backoff for transient failures."""
    response = None
    for attempt in range(max_retries + 1):
        response = client.request(method, url, **kwargs)
        if response.status_code not in {429, 500, 502, 503, 504}:
            return response
        if attempt == max_retries:
            break
        retry_after = response.headers.get("Retry-After")
        wait = float(retry_after) if retry_after else min(2**attempt, 30)
        LOGGER.warning(
            "Coordinate API: %s returned HTTP %s; retrying in %.1f seconds",
            url,
            response.status_code,
            wait,
        )
        time.sleep(wait)
    return response


def query_ensembl_coordinates(
    symbols: Sequence[str],
    *,
    cache_path: Path,
    refresh_cache: bool = False,
    batch_size: int = 1000,
    timeout: float = 120.0,
    delay: float = 0.1,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> tuple[Dict[str, List[Dict[str, object]]], Dict[str, object]]:
    """Query and cache live Ensembl gene coordinates in symbol batches."""
    if batch_size < 1:
        raise ValueError("Coordinate API batch_size must be a positive integer")
    cached = {}
    if cache_path.is_file() and not refresh_cache:
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
    missing = [symbol for symbol in symbols if symbol not in cached]
    total_batches = (len(missing) + batch_size - 1) // batch_size if missing else 0
    LOGGER.info(
        (
            "Coordinate API: %d requested symbols, %d cached, %d live Ensembl "
            "symbols across %d batches of up to %d"
        ),
        len(symbols),
        len(symbols) - len(missing),
        len(missing),
        total_batches,
        batch_size,
    )
    client = session or requests.Session()
    started = time.perf_counter()
    starts = range(0, len(missing), batch_size)
    for batch_number, start in enumerate(
        track(
            starts,
            enabled=progress,
            description="Ensembl coordinate API",
            total=total_batches,
            unit="batches",
        ),
        1,
    ):
        batch = list(missing[start : start + batch_size])
        LOGGER.info(
            "Coordinate API: querying Ensembl batch %d/%d (%d symbols)",
            batch_number,
            total_batches,
            len(batch),
        )
        response = _request_with_retries(
            client,
            "POST",
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            max_retries=max_retries,
            params={"expand": 0},
            json={"symbols": batch},
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": "GeneMapKit publication benchmark",
            },
            timeout=timeout,
        )
        if response.status_code in {400, 404}:
            for symbol in batch:
                cached[symbol] = []
            continue
        response.raise_for_status()
        records = response.json()
        if not isinstance(records, Mapping):
            records = {}
        found = 0
        for symbol in batch:
            record = records.get(symbol)
            if not isinstance(record, Mapping):
                cached[symbol] = []
                continue
            cached[symbol] = [
                {
                    "ensembl_gene_id": normalize_identifier(
                        "ensembl_gene_id", record.get("id", "")
                    ),
                    "assembly": record.get("assembly_name"),
                    "chromosome": str(record.get("seq_region_name", "")),
                    "start": record.get("start"),
                    "end": record.get("end"),
                    "strand": "+" if record.get("strand") == 1 else "-",
                }
            ]
            found += 1
        LOGGER.info(
            "Coordinate API: batch %d/%d returned coordinates for %d/%d symbols",
            batch_number,
            total_batches,
            found,
            len(batch),
        )
        if delay:
            time.sleep(delay)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cached, indent=2) + "\n", encoding="utf-8")
    return (
        {symbol: list(cached.get(symbol, [])) for symbol in symbols},
        {
            "requested_symbols": len(symbols),
            "queried_symbols": len(missing),
            "cached_symbols": len(symbols) - len(missing),
            "elapsed_seconds": round(time.perf_counter() - started, 6),
        },
    )


def benchmark_coordinates(
    database: Path,
    gtf_path: Path,
    symbols: Sequence[str],
    *,
    output_dir: Path,
    refresh_cache: bool = False,
    coordinate_batch_size: int = 1000,
    progress: bool = False,
) -> Dict[str, object]:
    """Compare pinned local-GTF coordinates with live Ensembl coordinates."""
    LOGGER.info(
        "Coordinate concordance: scanning local GTF for %d symbols from %s",
        len(symbols),
        gtf_path,
    )
    started = time.perf_counter()
    local = load_gtf_gene_coordinates(gtf_path, symbols)
    local_seconds = time.perf_counter() - started
    LOGGER.info(
        "Coordinate concordance: local GTF scan finished in %.4f seconds",
        local_seconds,
    )
    external, external_timing = query_ensembl_coordinates(
        symbols,
        cache_path=output_dir / "cache" / "ensembl_coordinates.json",
        refresh_cache=refresh_cache,
        batch_size=coordinate_batch_size,
        progress=progress,
    )
    details = []
    counts: Dict[str, int] = {}
    for symbol in track(
        symbols,
        enabled=progress,
        description="Comparing coordinates",
        total=len(symbols),
        unit="genes",
    ):
        local_rows = local.get(symbol, [])
        external_rows = external.get(symbol, [])
        local_set = {
            (
                row["ensembl_gene_id"],
                row["assembly"],
                row["chromosome"],
                row["start"],
                row["end"],
                row["strand"],
            )
            for row in local_rows
        }
        external_set = {
            (
                row["ensembl_gene_id"],
                row["assembly"],
                row["chromosome"],
                row["start"],
                row["end"],
                row["strand"],
            )
            for row in external_rows
        }
        if local_set == external_set and local_set:
            status = "exact"
        elif local_set and external_set and local_set.intersection(external_set):
            status = "partial"
        elif local_set and not external_set:
            status = "local_only"
        elif external_set and not local_set:
            status = "external_only"
        elif not local_set and not external_set:
            status = "both_unmapped"
        else:
            status = "conflict"
        counts[status] = counts.get(status, 0) + 1
        details.append(
            {
                "symbol": symbol,
                "status": status,
                "local_gtf_coordinates": json.dumps(local_rows, separators=(",", ":")),
                "ensembl_live_coordinates": json.dumps(
                    external_rows, separators=(",", ":")
                ),
            }
        )
    summary = {
        "tested_symbols": len(symbols),
        "exact": counts.get("exact", 0),
        "partial": counts.get("partial", 0),
        "conflict": counts.get("conflict", 0),
        "local_only": counts.get("local_only", 0),
        "external_only": counts.get("external_only", 0),
        "both_unmapped": counts.get("both_unmapped", 0),
        "exact_pct": _percent(counts.get("exact", 0), len(symbols)),
    }
    summary.update(interval_fields("exact", counts.get("exact", 0), len(symbols)))
    return {
        "capability": {
            "sqlite_coordinates_available": False,
            "local_coordinate_source": str(gtf_path),
            "limitation": (
                "Build schema 2.1.0 does not store coordinates; local coordinate "
                "concordance is measured directly from the pinned GTF snapshot."
            ),
        },
        "summary": summary,
        "timings": {
            "local_gtf_scan_seconds": round(local_seconds, 6),
            "ensembl_live": external_timing,
        },
        "details": details,
    }


def aggregate_held_out_reports(
    report_dirs: Sequence[Path],
) -> Dict[str, object]:
    """Collect only independent rows from existing held-out validation reports."""
    rows = []
    sources = []
    for directory in report_dirs:
        path = Path(directory) / "validation_metrics.json"
        if not path.is_file():
            sources.append({"directory": str(directory), "status": "missing"})
            continue
        report = json.loads(path.read_text(encoding="utf-8"))
        independent = [row for row in report.get("benchmarks", []) if row.get("independent")]
        rows.extend(independent)
        sources.append(
            {
                "directory": str(directory),
                "database": report.get("database"),
                "status": "included",
                "independent_rows": len(independent),
            }
        )
    return {"sources": sources, "independent_benchmarks": rows}


def run_publication_benchmark(
    database: Path,
    gtf_path: Path,
    symbols: Iterable[str],
    *,
    output_dir: Path = Path("results/publication_benchmark"),
    output_types: Sequence[str] = DEFAULT_OUTPUTS,
    providers: Sequence[str] = DEFAULT_PROVIDERS,
    held_out_report_dirs: Sequence[Path] = (),
    batch_sizes: Sequence[int] = (1, 10, 100, 1000),
    repeats: int = 3,
    mapping_policy: str = "all-supported",
    refresh_external_cache: bool = False,
    coordinate_batch_size: int = 1000,
    path_pivot_type: str = "entrez_id",
    path_sample: Optional[int] = 1000,
    path_random_seed: int = 20260616,
    progress: bool = False,
) -> Dict[str, object]:
    """Run all publication benchmark sections and write combined tables."""
    database, gtf_path, output_dir = Path(database), Path(gtf_path), Path(output_dir)
    unique_symbols = list(
        dict.fromkeys(str(value).strip() for value in symbols if str(value).strip())
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    scaling_sizes = list(dict.fromkeys([*batch_sizes, len(unique_symbols)]))
    LOGGER.info("Publication benchmark: local performance and memory scaling")
    performance = benchmark_local_performance(
        database,
        unique_symbols,
        output_types,
        batch_sizes=scaling_sizes,
        repeats=repeats,
        mapping_policy=mapping_policy,
        progress=progress,
    )
    LOGGER.info("Publication benchmark: coordinate concordance")
    coordinates = benchmark_coordinates(
        database,
        gtf_path,
        unique_symbols,
        output_dir=output_dir,
        refresh_cache=refresh_external_cache,
        coordinate_batch_size=coordinate_batch_size,
        progress=progress,
    )
    LOGGER.info("Publication benchmark: external API concordance")
    external = benchmark_external_services(
        database,
        unique_symbols,
        output_types,
        providers=providers,
        mapping_policy=mapping_policy,
        output_dir=output_dir / "external",
        refresh_cache=refresh_external_cache,
        progress=progress,
    )
    LOGGER.info("Publication benchmark: aggregating held-out validation reports")
    held_out = aggregate_held_out_reports(held_out_report_dirs)
    LOGGER.info("Publication benchmark: all-to-all paths and round trips")
    paths = benchmark_conversion_paths(
        database,
        id_types=DEFAULT_PATH_ID_TYPES,
        pivot_type=path_pivot_type,
        path_sample=path_sample,
        random_seed=path_random_seed,
        output_dir=output_dir / "paths",
        progress=progress,
    )
    LOGGER.info("Publication benchmark: writing combined reports")
    comparison = build_system_comparison(
        performance, coordinates, external, len(unique_symbols)
    )
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
        },
        "database": str(database),
        "gtf": str(gtf_path),
        "tested_symbols": len(unique_symbols),
        "mapping_policy": mapping_policy,
        "local_performance": performance,
        "coordinate_concordance": coordinates,
        "external_concordance": external,
        "held_out_validation": held_out,
        "conversion_paths": paths,
        "system_comparison": comparison,
    }
    write_publication_benchmark_report(output_dir, report)
    return report


def write_publication_benchmark_report(
    output_dir: Path, report: Mapping[str, object]
) -> None:
    """Write combined publication benchmark JSON, CSV tables, and Markdown."""
    (output_dir / "publication_benchmark.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    performance = report["local_performance"]
    coordinates = report["coordinate_concordance"]
    external = report["external_concordance"]
    held_out = report["held_out_validation"]
    paths = report["conversion_paths"]
    comparison = report["system_comparison"]
    _write_csv(output_dir / "performance_scaling.csv", performance["scaling"])
    _write_csv(output_dir / "coordinate_concordance.csv", coordinates["details"])
    _write_csv(
        output_dir / "external_concordance_summary.csv", external["summaries"]
    )
    _write_csv(output_dir / "external_provider_timings.csv", external["provider_timings"])
    _write_csv(
        output_dir / "held_out_validation.csv", held_out["independent_benchmarks"]
    )
    _write_csv(output_dir / "system_comparison.csv", comparison)
    _write_csv(
        output_dir / "path_conversion_coverage.csv", paths["conversion_coverage"]
    )
    _write_csv(
        output_dir / "path_round_trip_consistency.csv",
        paths["round_trip_consistency"],
    )
    _write_csv(
        output_dir / "path_pivot_consistency.csv", paths["pivot_path_consistency"]
    )

    coordinate_summary = coordinates["summary"]
    lines = [
        "# GeneMapKit Publication Benchmark",
        "",
        "This report separates independent held-out accuracy, external concordance, "
        "coordinate concordance, and local performance.",
        "",
        "## Capability and offline operation",
        "",
        "| Capability | Result |",
        "|---|---|",
        f"| Local SQLite identifier conversion | Available |",
        f"| Offline local conversion test | {performance['offline_operation_passed']} |",
        f"| Direct SQLite coordinate query | {coordinates['capability']['sqlite_coordinates_available']} |",
        f"| Local coordinate source | `{coordinates['capability']['local_coordinate_source']}` |",
        f"| Database size | {performance['database_size_gb']} GiB |",
        f"| Median converter startup | {performance['median_startup_seconds']} seconds |",
        "",
        "## Local performance and batch scalability",
        "",
        "| Batch size | Median seconds | IDs/second | Peak process RSS (MiB) | RSS increase (MiB) | Peak Python memory (MiB) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in performance["scaling"]:
        lines.append(
            f"| {row['batch_size']} | {row['median_seconds']} | "
            f"{row['identifiers_per_second']} | {row['median_peak_process_rss_mb']} | "
            f"{row['median_process_rss_increase_mb']} | "
            f"{row['median_peak_python_memory_mb']} |"
        )
    lines.extend(
        [
            "",
            "## Direct system comparison",
            "",
            "| System | Operation | Symbols | Seconds | Symbols/second | Network required | Offline capable |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in comparison:
        lines.append(
            f"| `{row['system']}` | {row['operation']} | {row['symbols']} | "
            f"{row['elapsed_seconds']} | {row['symbols_per_second']} | "
            f"{row['network_required']} | {row['offline_capable']} |"
        )
    lines.extend(
        [
            "",
            "## Coordinate concordance",
            "",
            coordinates["capability"]["limitation"],
            "",
            "| Symbols | Exact | Partial | Conflict | Local-only | External-only | Exact % (95% CI) | Local scan seconds |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
            f"| {coordinate_summary['tested_symbols']} | {coordinate_summary['exact']} | "
            f"{coordinate_summary['partial']} | {coordinate_summary['conflict']} | "
            f"{coordinate_summary['local_only']} | {coordinate_summary['external_only']} | "
            f"{_display_pct_ci(coordinate_summary, 'exact')} | "
            f"{coordinates['timings']['local_gtf_scan_seconds']} |",
            "",
            "## Multi-API external concordance",
            "",
            "| Output | Precision (95% CI) | Recall (95% CI) | Conflict | GeneMapKit-only |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in external["summaries"]:
        lines.append(
            f"| `{row['output_type']}` | {_display_pct_ci(row, 'concordance_precision')} | "
            f"{_display_pct_ci(row, 'concordance_recall')} | {row['conflict']} | "
            f"{row['genemapkit_only']} |"
        )
    lines.extend(
        [
            "",
            "### External query timings",
            "",
            "| Provider | Queried | Cached | Seconds |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in external["provider_timings"]:
        lines.append(
            f"| `{row['provider']}` | {row['queried_symbols']} | "
            f"{row['cached_symbols']} | {row['elapsed_seconds']} |"
        )
    lines.extend(
        [
            "",
            "## Independent held-out validation",
            "",
            "| Truth source | Conversion | Inputs | Precision (95% CI) | Recall (95% CI) | Exact sets |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in held_out["independent_benchmarks"]:
        lines.append(
            f"| `{row['truth_source']}` | `{row['input_type']} -> {row['output_type']}` | "
            f"{row['tested_inputs']} | {_display_pct_ci(row, 'precision')} | "
            f"{_display_pct_ci(row, 'recall')} | "
            f"{row['exact_set_match_pct']}% |"
        )
    if not held_out["independent_benchmarks"]:
        lines.append("| No independent held-out reports supplied | NA | 0 | NA | NA | NA |")
    lines.extend(
        [
            "",
            "## All-to-all conversion and path consistency",
            "",
            f"Direct all-to-all coverage uses the complete database inventory. "
            f"Round-trip and multi-step paths use `{paths['path_sampling']}` "
            f"deterministically sampled identifiers per namespace with seed "
            f"`{paths['random_seed']}`.",
            "",
            "| Source | Intermediate | Forward mapped | Original recovered | Exact round trip | Ambiguity expanded |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for row in paths["round_trip_consistency"]:
        if row["intermediate_type"] != paths["pivot_type"]:
            continue
        lines.append(
            f"| `{row['source_type']}` | `{row['intermediate_type']}` | "
            f"{_display_pct(row['forward_mapping_pct'])} | "
            f"{_display_pct(row['round_trip_recovery_pct'])} | "
            f"{_display_pct(row['exact_round_trip_pct'])} | "
            f"{_display_pct(row['ambiguity_expansion_pct'])} |"
        )
    lines.extend(
        [
            "",
            f"### Direct versus via `{paths['pivot_type']}`",
            "",
            "| Source | Target | Direct mapped | Path mapped | Exact sets | Overlap | Path precision | Path recall |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in paths["pivot_path_consistency"]:
        lines.append(
            f"| `{row['source_type']}` | `{row['target_type']}` | "
            f"{_display_pct(row['direct_mapping_pct'])} | "
            f"{_display_pct(row['path_mapping_pct'])} | "
            f"{_display_pct(row['exact_nonempty_set_pct'])} | "
            f"{_display_pct(row['overlap_input_pct'])} | "
            f"{_display_pct(row['path_vs_direct_precision_pct'])} | "
            f"{_display_pct(row['path_vs_direct_recall_pct'])} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation boundaries",
            "",
            "- Held-out validation is the independent accuracy analysis.",
            "- External API agreement is concordance because upstream sources overlap.",
            "- Confidence intervals are Wilson 95% intervals. Mapping-level intervals may "
            "be optimistic when mappings from the same input are correlated.",
            "- Coordinate results compare the pinned local GTF with live Ensembl; coordinates "
            "are not currently stored in SQLite schema 2.1.0.",
            "- Round-trip exact equality is not universally expected because one-to-many "
            "gene, transcript, and protein relationships can legitimately expand paths.",
            f"- Memory values are {performance['memory_metric']}.",
        ]
    )
    (output_dir / "publication_benchmark_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
