"""External concordance benchmarks for a built GeneMapKit database."""

from __future__ import annotations

import csv
import json
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set

import requests

from .builder import normalize_identifier
from .converter import GeneConverter, canonical_id_type, canonical_mapping_policy
from .statistics import interval_fields
from genemapkit.utils.progress import track


MYGENE_ENDPOINT = "https://mygene.info/v3/query"
MYGENE_OUTPUT_FIELDS = {
    "symbol": ("symbol",),
    "entrez_id": ("entrezgene",),
    "ensembl_gene_id": ("ensembl.gene",),
    "refseq_rna_id": ("refseq.rna",),
    "refseq_protein_id": ("refseq.protein",),
    "uniprot_id": ("uniprot.Swiss-Prot", "uniprot.TrEMBL"),
}


def _flatten(value: object) -> List[str]:
    if value is None:
        return []
    if isinstance(value, dict):
        values: List[str] = []
        for item in value.values():
            values.extend(_flatten(item))
        return values
    if isinstance(value, (list, tuple, set)):
        values = []
        for item in value:
            values.extend(_flatten(item))
        return values
    return [str(value)]


def _nested(record: Mapping[str, object], field: str) -> object:
    value: object = record
    for part in field.split("."):
        if not isinstance(value, Mapping):
            return None
        value = value.get(part)
    return value


def extract_mygene_values(
    record: Mapping[str, object], output_type: str
) -> Set[str]:
    """Extract and normalize one identifier namespace from a MyGene record."""
    canonical = canonical_id_type(output_type)
    values = set()
    for field in MYGENE_OUTPUT_FIELDS[canonical]:
        for value in _flatten(_nested(record, field)):
            normalized = normalize_identifier(canonical, value)
            if normalized:
                values.add(normalized)
    return values


def query_mygene_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 250,
    timeout: float = 60.0,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> Dict[str, Dict[str, Set[str]]]:
    """Query MyGene.info for exact human symbol matches in deterministic batches."""
    canonical_outputs = [canonical_id_type(value) for value in output_types]
    unsupported = sorted(set(canonical_outputs) - set(MYGENE_OUTPUT_FIELDS))
    if unsupported:
        raise ValueError(
            "MyGene comparison does not support output type(s): "
            + ", ".join(unsupported)
        )
    fields = sorted(
        {field for output in canonical_outputs for field in MYGENE_OUTPUT_FIELDS[output]}
        | {"symbol"}
    )
    client = session or requests.Session()
    results: Dict[str, Dict[str, Set[str]]] = {
        symbol: {output: set() for output in canonical_outputs} for symbol in symbols
    }
    starts = range(0, len(symbols), batch_size)
    for start in track(
        starts,
        enabled=progress,
        description="MyGene API",
        total=len(starts),
        unit="batches",
    ):
        batch = list(symbols[start : start + batch_size])
        response = None
        for attempt in range(max_retries + 1):
            response = client.post(
                MYGENE_ENDPOINT,
                data={
                    "q": ",".join(batch),
                    "scopes": "symbol",
                    "fields": ",".join(fields),
                    "species": "human",
                },
                timeout=timeout,
                headers={"User-Agent": "GeneMapKit external concordance benchmark"},
            )
            if response.status_code not in {429, 500, 502, 503, 504}:
                break
            if attempt == max_retries:
                break
            retry_after = response.headers.get("Retry-After")
            time.sleep(float(retry_after) if retry_after else min(2**attempt, 30))
        response.raise_for_status()
        records = response.json()
        if not isinstance(records, list):
            raise ValueError("Unexpected MyGene response; expected a list of records")
        for record in records:
            if not isinstance(record, Mapping) or record.get("notfound"):
                continue
            query = str(record.get("query", "")).strip()
            approved_symbol = str(record.get("symbol", "")).strip()
            if query not in results or approved_symbol.casefold() != query.casefold():
                continue
            for output in canonical_outputs:
                results[query][output].update(extract_mygene_values(record, output))
    return results


def _comparison_status(gmk: Set[str], external: Set[str]) -> str:
    if not gmk and not external:
        return "both_unmapped"
    if gmk == external:
        return "exact_set"
    if gmk and external and gmk.intersection(external):
        return "overlap"
    if gmk and not external:
        return "genemapkit_only"
    if external and not gmk:
        return "mygene_only"
    return "disjoint"


def compare_mygene(
    database: Path,
    symbols: Iterable[str],
    output_types: Sequence[str],
    *,
    mapping_policy: str = "all-supported",
    output_dir: Optional[Path] = None,
    fetcher: Callable[
        [Sequence[str], Sequence[str]], Dict[str, Dict[str, Set[str]]]
    ] = query_mygene_symbols,
    progress: bool = False,
) -> Dict[str, object]:
    """Compare GeneMapKit conversion sets with live MyGene.info annotations."""
    unique_symbols = list(dict.fromkeys(str(value).strip() for value in symbols if str(value).strip()))
    if not unique_symbols:
        raise ValueError("No non-empty symbols were provided")
    outputs = [canonical_id_type(value) for value in output_types]
    policy = canonical_mapping_policy(mapping_policy)
    converter = GeneConverter(database, mapping_policy=policy)
    mygene = (
        fetcher(unique_symbols, outputs, progress=progress)
        if fetcher is query_mygene_symbols
        else fetcher(unique_symbols, outputs)
    )
    detail_rows: List[Dict[str, object]] = []
    counters: Dict[str, Dict[str, int]] = {
        output: defaultdict(int) for output in outputs
    }

    with converter.open_connection() as connection:
        for symbol in track(
            unique_symbols,
            enabled=progress,
            description="Comparing MyGene mappings",
            total=len(unique_symbols),
            unit="genes",
        ):
            for output in outputs:
                gmk = set(
                    converter.convert_single(
                        symbol,
                        "symbol",
                        [output],
                        mapping_policy=policy,
                        connection=connection,
                    )[output]
                )
                external = set(mygene.get(symbol, {}).get(output, set()))
                common = gmk.intersection(external)
                status = _comparison_status(gmk, external)
                counters[output][status] += 1
                counters[output]["genemapkit_mappings"] += len(gmk)
                counters[output]["mygene_mappings"] += len(external)
                counters[output]["common_mappings"] += len(common)
                detail_rows.append(
                    {
                        "symbol": symbol,
                        "output_type": output,
                        "mapping_policy": policy,
                        "status": status,
                        "genemapkit": ";".join(sorted(gmk)),
                        "mygene": ";".join(sorted(external)),
                        "common": ";".join(sorted(common)),
                        "genemapkit_only": ";".join(sorted(gmk - external)),
                        "mygene_only": ";".join(sorted(external - gmk)),
                    }
                )

    summaries = []
    for output in outputs:
        counts = counters[output]
        common = counts["common_mappings"]
        gmk_total = counts["genemapkit_mappings"]
        mygene_total = counts["mygene_mappings"]
        row = {
                "output_type": output,
                "tested_symbols": len(unique_symbols),
                "exact_set_inputs": counts["exact_set"],
                "overlap_inputs": counts["overlap"],
                "disjoint_inputs": counts["disjoint"],
                "genemapkit_only_inputs": counts["genemapkit_only"],
                "mygene_only_inputs": counts["mygene_only"],
                "both_unmapped_inputs": counts["both_unmapped"],
                "genemapkit_mappings": gmk_total,
                "mygene_mappings": mygene_total,
                "common_mappings": common,
                "concordance_precision_pct": (
                    round(100.0 * common / gmk_total, 4) if gmk_total else None
                ),
                "concordance_recall_pct": (
                    round(100.0 * common / mygene_total, 4) if mygene_total else None
                ),
            }
        row.update(interval_fields("concordance_precision", common, gmk_total))
        row.update(interval_fields("concordance_recall", common, mygene_total))
        summaries.append(row)

    report: Dict[str, object] = {
        "comparison_type": "external_concordance_not_ground_truth",
        "external_service": "MyGene.info",
        "external_endpoint": MYGENE_ENDPOINT,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "database": str(database),
        "database_metadata": converter.metadata(),
        "mapping_policy": policy,
        "input_type": "symbol",
        "output_types": outputs,
        "tested_symbols": len(unique_symbols),
        "summaries": summaries,
        "details": detail_rows,
        "interpretation": (
            "MyGene.info integrates overlapping upstream databases and is not an "
            "independent gold standard. These metrics measure external concordance."
        ),
    }
    if output_dir:
        write_mygene_report(Path(output_dir), report)
    return report


def write_mygene_report(output_dir: Path, report: Mapping[str, object]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "mygene_concordance.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    summaries = list(report["summaries"])
    details = list(report["details"])
    for name, rows in (
        ("mygene_concordance_summary.csv", summaries),
        ("mygene_concordance_details.csv", details),
    ):
        if not rows:
            continue
        with (output_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    lines = [
        "# GeneMapKit / MyGene.info Concordance",
        "",
        "> MyGene.info is an external integrated annotation service, not an independent "
        "gold standard. Report these values as concordance, not accuracy.",
        "",
        f"- Mapping policy: `{report['mapping_policy']}`",
        f"- Tested symbols: `{report['tested_symbols']}`",
        "- Intervals: Wilson 95% confidence intervals.",
        "",
        "| Output type | Exact inputs | Overlap inputs | GeneMapKit-only | MyGene-only | Concordance precision (95% CI) | Concordance recall (95% CI) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            f"| `{row['output_type']}` | {row['exact_set_inputs']} | "
            f"{row['overlap_inputs']} | {row['genemapkit_only_inputs']} | "
            f"{row['mygene_only_inputs']} | {row['concordance_precision_pct']}% "
            f"[{row['concordance_precision_ci_low_pct']}, {row['concordance_precision_ci_high_pct']}] | "
            f"{row['concordance_recall_pct']}% "
            f"[{row['concordance_recall_ci_low_pct']}, {row['concordance_recall_ci_high_pct']}] |"
        )
    (output_dir / "mygene_concordance_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
