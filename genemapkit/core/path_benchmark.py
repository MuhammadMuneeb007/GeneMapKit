"""All-to-all conversion coverage and path-consistency benchmarks."""

from __future__ import annotations

import csv
import json
import logging
import random
import sqlite3
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set

from .converter import canonical_id_type
from .statistics import interval_fields
from genemapkit.utils.progress import track

LOGGER = logging.getLogger("genemapkit")


DEFAULT_PATH_ID_TYPES = (
    "symbol",
    "hgnc_id",
    "entrez_id",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_protein_id",
    "refseq_rna_id",
    "refseq_protein_id",
    "uniprot_id",
)


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


def _connect(database: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(f"file:{database.resolve()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def _sample_inputs(
    connection: sqlite3.Connection,
    id_type: str,
    sample: Optional[int],
    *,
    random_seed: int,
) -> List[str]:
    values = [
        row[0]
        for row in connection.execute(
            "SELECT DISTINCT id_value FROM mappings WHERE id_type = ? ORDER BY id_value",
            (id_type,),
        )
    ]
    if sample is None or sample >= len(values):
        return values
    generator = random.Random(f"{random_seed}:{id_type}")
    return sorted(generator.sample(values, sample))


def _direct_sets(
    connection: sqlite3.Connection,
    source_type: str,
    target_type: str,
    inputs: Sequence[str],
) -> Dict[str, Set[str]]:
    results: Dict[str, Set[str]] = {value: set() for value in inputs}
    for start in range(0, len(inputs), 500):
        batch = list(inputs[start : start + 500])
        placeholders = ",".join("?" for _ in batch)
        sql = f"""
            SELECT DISTINCT source_value, target_value
            FROM relationships
            WHERE source_type = ? AND target_type = ?
              AND source_value IN ({placeholders})
        """
        for row in connection.execute(sql, [source_type, target_type, *batch]):
            results[row["source_value"]].add(row["target_value"])
    return results


def _path_sets(
    connection: sqlite3.Connection,
    source_type: str,
    pivot_type: str,
    target_type: str,
    inputs: Sequence[str],
) -> Dict[str, Set[str]]:
    results: Dict[str, Set[str]] = {value: set() for value in inputs}
    for start in range(0, len(inputs), 250):
        batch = list(inputs[start : start + 250])
        placeholders = ",".join("?" for _ in batch)
        sql = f"""
            SELECT DISTINCT first.source_value, second.target_value
            FROM relationships AS first
            JOIN relationships AS second
              ON second.source_type = ?
             AND second.source_value = first.target_value
             AND second.target_type = ?
            WHERE first.source_type = ?
              AND first.target_type = ?
              AND first.source_value IN ({placeholders})
        """
        parameters = [
            pivot_type,
            target_type,
            source_type,
            pivot_type,
            *batch,
        ]
        for row in connection.execute(sql, parameters):
            results[row["source_value"]].add(row["target_value"])
    return results


def benchmark_conversion_paths(
    database: Path,
    *,
    id_types: Sequence[str] = DEFAULT_PATH_ID_TYPES,
    pivot_type: str = "entrez_id",
    path_sample: Optional[int] = 1000,
    random_seed: int = 20260616,
    max_examples: int = 100,
    output_dir: Optional[Path] = None,
    progress: bool = False,
) -> Dict[str, object]:
    """Benchmark direct all-to-all coverage and multi-step path consistency."""
    database = Path(database)
    types = list(dict.fromkeys(canonical_id_type(value) for value in id_types))
    pivot = canonical_id_type(pivot_type)
    if pivot not in types:
        types.append(pivot)
    started = time.perf_counter()
    LOGGER.info(
        "Path benchmark: %d namespaces, pivot=%s, sample=%s",
        len(types), pivot, "complete" if path_sample is None else path_sample,
    )
    with _connect(database) as connection:
        inventory = {
            row["id_type"]: row["identifier_count"]
            for row in connection.execute(
                """
                SELECT id_type, COUNT(DISTINCT id_value) AS identifier_count
                FROM mappings
                WHERE id_type IN ({})
                GROUP BY id_type
                """.format(",".join("?" for _ in types)),
                types,
            )
        }
        pair_counts = {
            (row["source_type"], row["target_type"]): (
                row["mapped_inputs"],
                row["mapping_edges"],
            )
            for row in connection.execute(
                """
                SELECT source_type, target_type,
                       COUNT(DISTINCT source_value) AS mapped_inputs,
                       COUNT(DISTINCT source_value || char(31) || target_value)
                           AS mapping_edges
                FROM relationships
                WHERE source_type IN ({types}) AND target_type IN ({types})
                GROUP BY source_type, target_type
                """.format(types=",".join("?" for _ in types)),
                [*types, *types],
            )
        }
        coverage_rows: List[Dict[str, object]] = []
        roundtrip_rows: List[Dict[str, object]] = []
        path_rows: List[Dict[str, object]] = []
        examples: List[Dict[str, object]] = []

        for source_type in track(
            types,
            enabled=progress,
            description="Conversion coverage",
            total=len(types),
            unit="namespaces",
        ):
            total = inventory.get(source_type, 0)
            for target_type in types:
                if source_type == target_type:
                    continue
                mapped, edges = pair_counts.get((source_type, target_type), (0, 0))
                row = {
                        "source_type": source_type,
                        "target_type": target_type,
                        "total_source_identifiers": total,
                        "mapped_source_identifiers": mapped,
                        "mapping_edges": edges,
                        "conversion_success_pct": _percent(mapped, total),
                        "mean_outputs_per_mapped_input": (
                            round(edges / mapped, 4) if mapped else None
                        ),
                    }
                row.update(interval_fields("conversion_success", mapped, total))
                coverage_rows.append(row)

        for source_type in track(
            types,
            enabled=progress,
            description="Round trips and paths",
            total=len(types),
            unit="namespaces",
        ):
            inputs = _sample_inputs(
                connection, source_type, path_sample, random_seed=random_seed
            )
            if not inputs:
                continue
            for intermediate_type in types:
                if intermediate_type == source_type:
                    continue
                forward = _direct_sets(
                    connection, source_type, intermediate_type, inputs
                )
                reached = [value for value in inputs if forward[value]]
                if reached:
                    intermediate_values = sorted(
                        set().union(*(forward[value] for value in reached))
                    )
                    backward_by_intermediate = _direct_sets(
                        connection,
                        intermediate_type,
                        source_type,
                        intermediate_values,
                    )
                else:
                    backward_by_intermediate = {}
                recovered = exact = expanded = 0
                for value in inputs:
                    returned = set().union(
                        *(
                            backward_by_intermediate.get(item, set())
                            for item in forward[value]
                        )
                    ) if forward[value] else set()
                    recovered += value in returned
                    exact += returned == {value}
                    expanded += bool(returned - {value})
                    if (
                        returned != {value}
                        and len(examples) < max_examples
                    ):
                        examples.append(
                            {
                                "analysis": "round_trip",
                                "source_type": source_type,
                                "pivot_type": intermediate_type,
                                "target_type": source_type,
                                "input_id": value,
                                "direct": ";".join(sorted(forward[value])),
                                "via_path": ";".join(sorted(returned)),
                                "missing": value if value not in returned else "",
                                "extra": ";".join(sorted(returned - {value})),
                            }
                        )
                row = {
                        "source_type": source_type,
                        "intermediate_type": intermediate_type,
                        "sampled_inputs": len(inputs),
                        "forward_mapped_inputs": len(reached),
                        "forward_mapping_pct": _percent(len(reached), len(inputs)),
                        "original_recovered_inputs": recovered,
                        "round_trip_recovery_pct": _percent(recovered, len(inputs)),
                        "exact_round_trip_inputs": exact,
                        "exact_round_trip_pct": _percent(exact, len(inputs)),
                        "ambiguity_expanded_inputs": expanded,
                        "ambiguity_expansion_pct": _percent(expanded, len(inputs)),
                    }
                row.update(interval_fields("forward_mapping", len(reached), len(inputs)))
                row.update(interval_fields("round_trip_recovery", recovered, len(inputs)))
                row.update(interval_fields("exact_round_trip", exact, len(inputs)))
                row.update(interval_fields("ambiguity_expansion", expanded, len(inputs)))
                roundtrip_rows.append(row)

            if source_type == pivot:
                continue
            for target_type in types:
                if target_type in {source_type, pivot}:
                    continue
                direct = _direct_sets(connection, source_type, target_type, inputs)
                via = _path_sets(
                    connection, source_type, pivot, target_type, inputs
                )
                direct_mapped = via_mapped = exact = overlap = 0
                true_positive = direct_total = path_total = 0
                path_only_inputs = direct_only_inputs = 0
                for value in inputs:
                    direct_set, path_set = direct[value], via[value]
                    common = direct_set.intersection(path_set)
                    direct_mapped += bool(direct_set)
                    via_mapped += bool(path_set)
                    exact += bool(direct_set) and direct_set == path_set
                    overlap += bool(common)
                    true_positive += len(common)
                    direct_total += len(direct_set)
                    path_total += len(path_set)
                    path_only_inputs += bool(path_set and not direct_set)
                    direct_only_inputs += bool(direct_set and not path_set)
                    if (
                        direct_set != path_set
                        and len(examples) < max_examples
                    ):
                        examples.append(
                            {
                                "analysis": "pivot_path",
                                "source_type": source_type,
                                "pivot_type": pivot,
                                "target_type": target_type,
                                "input_id": value,
                                "direct": ";".join(sorted(direct_set)),
                                "via_path": ";".join(sorted(path_set)),
                                "missing": ";".join(sorted(direct_set - path_set)),
                                "extra": ";".join(sorted(path_set - direct_set)),
                            }
                        )
                row = {
                        "source_type": source_type,
                        "pivot_type": pivot,
                        "target_type": target_type,
                        "sampled_inputs": len(inputs),
                        "direct_mapped_inputs": direct_mapped,
                        "path_mapped_inputs": via_mapped,
                        "direct_mapping_pct": _percent(direct_mapped, len(inputs)),
                        "path_mapping_pct": _percent(via_mapped, len(inputs)),
                        "exact_nonempty_set_inputs": exact,
                        "exact_nonempty_set_pct": _percent(exact, len(inputs)),
                        "overlap_inputs": overlap,
                        "overlap_input_pct": _percent(overlap, len(inputs)),
                        "path_vs_direct_precision_pct": _percent(
                            true_positive, path_total
                        ),
                        "path_vs_direct_recall_pct": _percent(
                            true_positive, direct_total
                        ),
                        "path_only_inputs": path_only_inputs,
                        "direct_only_inputs": direct_only_inputs,
                    }
                row.update(interval_fields("direct_mapping", direct_mapped, len(inputs)))
                row.update(interval_fields("path_mapping", via_mapped, len(inputs)))
                row.update(interval_fields("exact_nonempty_set", exact, len(inputs)))
                row.update(interval_fields("overlap_input", overlap, len(inputs)))
                row.update(interval_fields("path_vs_direct_precision", true_positive, path_total))
                row.update(interval_fields("path_vs_direct_recall", true_positive, direct_total))
                path_rows.append(row)

    report: Dict[str, object] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "database": str(database),
        "id_types": types,
        "pivot_type": pivot,
        "path_sampling": "complete" if path_sample is None else path_sample,
        "random_seed": random_seed,
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "interpretation": {
            "conversion_success": "source identifiers with at least one direct target",
            "round_trip_recovery": "original identifier appears after A -> B -> A",
            "exact_round_trip": "A -> B -> A returns only the original identifier",
            "path_consistency": (
                "direct A -> C outputs compared with A -> pivot -> C outputs"
            ),
            "warning": (
                "Exact equality is not universally expected because biological "
                "gene/transcript/protein mappings are often one-to-many."
            ),
        },
        "namespace_inventory": [
            {"id_type": value, "identifier_count": inventory.get(value, 0)}
            for value in types
        ],
        "conversion_coverage": coverage_rows,
        "round_trip_consistency": roundtrip_rows,
        "pivot_path_consistency": path_rows,
        "examples": examples,
    }
    if output_dir:
        write_path_benchmark_report(Path(output_dir), report)
    LOGGER.info("Path benchmark complete in %s seconds", report["elapsed_seconds"])
    return report


def write_path_benchmark_report(
    output_dir: Path, report: Mapping[str, object]
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "path_benchmark.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    _write_csv(output_dir / "namespace_inventory.csv", report["namespace_inventory"])
    _write_csv(output_dir / "conversion_coverage.csv", report["conversion_coverage"])
    _write_csv(
        output_dir / "round_trip_consistency.csv", report["round_trip_consistency"]
    )
    _write_csv(
        output_dir / "pivot_path_consistency.csv", report["pivot_path_consistency"]
    )
    _write_csv(output_dir / "path_disagreement_examples.csv", report["examples"])
    types = list(report["id_types"])
    coverage = {
        (row["source_type"], row["target_type"]): row
        for row in report["conversion_coverage"]
    }
    lines = [
        "# All-to-All Conversion and Path Consistency",
        "",
        f"- Pivot namespace: `{report['pivot_type']}`",
        f"- Path sample per source namespace: `{report['path_sampling']}`",
        f"- Random seed: `{report['random_seed']}`",
        f"- Runtime: `{report['elapsed_seconds']}` seconds",
        "",
        "Exact round-trip equality is not universally expected because gene, transcript, "
        "and protein mappings can be one-to-many.",
        "Reported intervals are Wilson 95% confidence intervals.",
        "",
        "## Direct conversion success matrix (%)",
        "",
        "| Source \\ Target | " + " | ".join(f"`{value}`" for value in types) + " |",
        "|---|" + "---:|" * len(types),
    ]
    for source in types:
        values = []
        for target in types:
            if source == target:
                values.append("100.0")
            else:
                values.append(str(coverage[(source, target)]["conversion_success_pct"]))
        lines.append(f"| `{source}` | " + " | ".join(values) + " |")
    lines.extend(
        [
            "",
            "## Round-trip consistency",
            "",
            "| Source | Intermediate | Inputs | Forward mapped | Original recovered (95% CI) | Exact round trip | Ambiguity expanded |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in report["round_trip_consistency"]:
        lines.append(
            f"| `{row['source_type']}` | `{row['intermediate_type']}` | "
            f"{row['sampled_inputs']} | {_display_pct(row['forward_mapping_pct'])} | "
            f"{_display_pct_ci(row, 'round_trip_recovery')} | "
            f"{_display_pct(row['exact_round_trip_pct'])} | "
            f"{_display_pct(row['ambiguity_expansion_pct'])} |"
        )
    lines.extend(
        [
            "",
            f"## Direct versus via `{report['pivot_type']}`",
            "",
            "| Source | Target | Direct mapped | Path mapped | Exact sets | Overlap | Path precision (95% CI) | Path recall (95% CI) |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in report["pivot_path_consistency"]:
        lines.append(
            f"| `{row['source_type']}` | `{row['target_type']}` | "
            f"{_display_pct(row['direct_mapping_pct'])} | "
            f"{_display_pct(row['path_mapping_pct'])} | "
            f"{_display_pct(row['exact_nonempty_set_pct'])} | "
            f"{_display_pct(row['overlap_input_pct'])} | "
            f"{_display_pct_ci(row, 'path_vs_direct_precision')} | "
            f"{_display_pct_ci(row, 'path_vs_direct_recall')} |"
        )
    (output_dir / "path_benchmark_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
