"""Relationship-level validation for a built GeneMapKit SQLite database."""

from __future__ import annotations

import csv
import json
import logging
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from .builder import PARSERS, load_manifest, normalize_identifier, select_sources
from .converter import GeneConverter
from .statistics import interval_fields
from genemapkit.utils.progress import track

LOGGER = logging.getLogger("genemapkit")


Conversion = Tuple[str, str]
TruthMap = Dict[str, Set[str]]

VALIDATION_CONVERSIONS: Dict[str, Tuple[Conversion, ...]] = {
    "mane_summary": (
        ("ensembl_transcript_id", "ensembl_protein_id"),
        ("ensembl_transcript_id", "refseq_rna_id"),
        ("ensembl_protein_id", "refseq_protein_id"),
        ("refseq_rna_id", "ensembl_transcript_id"),
    ),
    "ncbi_gene2ensembl": (
        ("entrez_id", "ensembl_gene_id"),
        ("ensembl_gene_id", "entrez_id"),
        ("ensembl_transcript_id", "ensembl_protein_id"),
        ("ensembl_transcript_id", "refseq_rna_id"),
        ("ensembl_protein_id", "refseq_protein_id"),
    ),
    "ncbi_gene2refseq": (
        ("entrez_id", "refseq_rna_id"),
        ("entrez_id", "refseq_protein_id"),
        ("refseq_rna_id", "entrez_id"),
        ("refseq_protein_id", "entrez_id"),
    ),
}


def _record_values(record: Iterable[Tuple[str, str]]) -> Dict[str, Set[str]]:
    values: Dict[str, Set[str]] = defaultdict(set)
    for id_type, value in record:
        normalized = normalize_identifier(id_type, value)
        if normalized:
            values[id_type].add(normalized)
    return values


def load_truth(
    source: str,
    path: Path,
    conversions: Sequence[Conversion],
    *,
    sample: Optional[int] = None,
) -> Dict[Conversion, TruthMap]:
    """Load expected mappings directly from one validation source."""
    LOGGER.info(
        "Loading truth source %s for %d conversion(s)%s",
        source,
        len(conversions),
        f" (sample cap: {sample})" if sample is not None else " (complete)",
    )
    truth: Dict[Conversion, TruthMap] = {
        conversion: defaultdict(set) for conversion in conversions
    }
    parser = PARSERS[source]
    records_read = 0
    for record, _ in parser(path):
        records_read += 1
        values = _record_values(record)
        for conversion in conversions:
            input_type, output_type = conversion
            inputs = values.get(input_type, set())
            outputs = values.get(output_type, set())
            if not inputs or not outputs:
                continue
            mapping = truth[conversion]
            for input_id in inputs:
                if sample is None or input_id in mapping or len(mapping) < sample:
                    mapping[input_id].update(outputs)
        if sample is not None and all(len(mapping) >= sample for mapping in truth.values()):
            break
    result = {conversion: dict(mapping) for conversion, mapping in truth.items()}
    for conversion, mapping in result.items():
        LOGGER.info(
            "  %s %s -> %s: %d truth inputs loaded (from %d records parsed)",
            source, conversion[0], conversion[1], len(mapping), records_read,
        )
    return result


def _percent(numerator: int, denominator: int) -> Optional[float]:
    return round(100.0 * numerator / denominator, 4) if denominator else None


def evaluate_truth(
    converter: GeneConverter,
    truth: Mapping[str, Set[str]],
    input_type: str,
    output_type: str,
    *,
    mapping_policy: str = "all-supported",
    max_examples: int = 25,
    progress: bool = False,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    """Compare converter output sets with expected output sets."""
    true_positive = false_positive = false_negative = 0
    exact = overmapped = undermapped = unmapped = 0
    examples: List[Dict[str, object]] = []

    inputs = sorted(truth)
    with converter._connect() as connection:
        for input_id in track(
            inputs,
            enabled=progress,
            description=f"Validating {input_type} -> {output_type}",
            total=len(inputs),
            unit="ids",
        ):
            expected = set(truth[input_id])
            predicted = set(
                converter.convert_single(
                    input_id,
                    input_type,
                    [output_type],
                    mapping_policy=mapping_policy,
                    connection=connection,
                )[output_type]
            )
            extra = predicted - expected
            missing = expected - predicted
            true_positive += len(predicted & expected)
            false_positive += len(extra)
            false_negative += len(missing)
            exact += predicted == expected
            overmapped += bool(extra)
            undermapped += bool(missing)
            unmapped += not predicted
            if (extra or missing) and len(examples) < max_examples:
                examples.append(
                    {
                        "input_type": input_type,
                        "output_type": output_type,
                        "input_id": input_id,
                        "mapping_policy": mapping_policy,
                        "expected": ";".join(sorted(expected)),
                        "predicted": ";".join(sorted(predicted)),
                        "extra": ";".join(sorted(extra)),
                        "missing": ";".join(sorted(missing)),
                    }
                )

    tested = len(truth)
    metrics: Dict[str, object] = {
        "input_type": input_type,
        "output_type": output_type,
        "mapping_policy": mapping_policy,
        "tested_inputs": tested,
        "expected_mappings": true_positive + false_negative,
        "predicted_mappings": true_positive + false_positive,
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "precision_pct": _percent(true_positive, true_positive + false_positive),
        "recall_pct": _percent(true_positive, true_positive + false_negative),
        "exact_set_matches": exact,
        "exact_set_match_pct": _percent(exact, tested),
        "overmapped_inputs": overmapped,
        "overmapped_input_pct": _percent(overmapped, tested),
        "undermapped_inputs": undermapped,
        "undermapped_input_pct": _percent(undermapped, tested),
        "unmapped_inputs": unmapped,
        "unmapped_input_pct": _percent(unmapped, tested),
    }
    metrics.update(
        interval_fields("precision", true_positive, true_positive + false_positive)
    )
    metrics.update(
        interval_fields("recall", true_positive, true_positive + false_negative)
    )
    metrics.update(interval_fields("exact_set_match", exact, tested))
    metrics.update(interval_fields("overmapped_input", overmapped, tested))
    metrics.update(interval_fields("undermapped_input", undermapped, tested))
    metrics.update(interval_fields("unmapped_input", unmapped, tested))
    return metrics, examples


def database_diagnostics(database: Path) -> Dict[str, object]:
    """Measure group shapes that can create broad or conflicting mappings."""
    with sqlite3.connect(database) as connection:
        diagnostics: Dict[str, object] = {}
        diagnostics["direct_relationship_rows"] = connection.execute(
            "SELECT COUNT(*) FROM relationships"
        ).fetchone()[0]
        for input_type, output_type in (
            ("ensembl_transcript_id", "ensembl_protein_id"),
            ("ensembl_transcript_id", "refseq_rna_id"),
            ("ensembl_protein_id", "refseq_protein_id"),
            ("entrez_id", "ensembl_gene_id"),
        ):
            diagnostics[
                f"inputs_with_multiple_direct_{input_type}_to_{output_type}"
            ] = connection.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT source_value
                    FROM relationships
                    WHERE source_type = ? AND target_type = ?
                    GROUP BY source_value
                    HAVING COUNT(DISTINCT target_value) > 1
                )
                """,
                (input_type, output_type),
            ).fetchone()[0]
        for id_type in (
            "hgnc_id",
            "entrez_id",
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "ensembl_protein_id",
            "refseq_rna_id",
            "refseq_protein_id",
            "uniprot_id",
        ):
            diagnostics[f"groups_with_multiple_{id_type}"] = connection.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT group_id
                    FROM mappings
                    WHERE id_type = ?
                    GROUP BY group_id
                    HAVING COUNT(DISTINCT id_value) > 1
                )
                """,
                (id_type,),
            ).fetchone()[0]
        diagnostics["groups_with_gene_transcript_and_protein"] = connection.execute(
            """
            SELECT COUNT(*) FROM (
                SELECT group_id
                FROM mappings
                WHERE id_type IN (
                    'ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_protein_id'
                )
                GROUP BY group_id
                HAVING COUNT(DISTINCT CASE WHEN id_type='ensembl_gene_id'
                                           THEN id_value END) > 0
                   AND COUNT(DISTINCT CASE WHEN id_type='ensembl_transcript_id'
                                           THEN id_value END) > 0
                   AND COUNT(DISTINCT CASE WHEN id_type='ensembl_protein_id'
                                           THEN id_value END) > 0
            )
            """
        ).fetchone()[0]
        diagnostics["groups_with_multiple_genes_and_multiple_proteins"] = (
            connection.execute(
                """
                SELECT COUNT(*) FROM (
                    SELECT group_id
                    FROM mappings
                    WHERE id_type IN ('ensembl_gene_id', 'ensembl_protein_id')
                    GROUP BY group_id
                    HAVING COUNT(DISTINCT CASE WHEN id_type='ensembl_gene_id'
                                               THEN id_value END) > 1
                       AND COUNT(DISTINCT CASE WHEN id_type='ensembl_protein_id'
                                               THEN id_value END) > 1
                )
                """
            ).fetchone()[0]
        )
        return diagnostics


def validate_database(
    database: Path,
    raw_dir: Path,
    *,
    sources: Optional[Sequence[str]] = None,
    sample: Optional[int] = 1000,
    max_examples: int = 25,
    require_held_out: bool = False,
    mapping_policy: str = "all-supported",
    output_dir: Optional[Path] = None,
    progress: bool = False,
) -> Dict[str, object]:
    """Run relationship benchmarks against raw validation sources."""
    database, raw_dir = Path(database), Path(raw_dir)
    converter = GeneConverter(database, mapping_policy=mapping_policy)
    metadata = converter.metadata()
    built_sources = set(metadata.get("selected_sources", []))
    manifest = load_manifest(raw_dir)
    available = select_sources(raw_dir, manifest)
    requested = list(sources or VALIDATION_CONVERSIONS)
    missing = sorted(set(requested) - set(available))
    if missing:
        raise FileNotFoundError(
            f"Validation sources are missing from the raw snapshot: {', '.join(missing)}"
        )

    rows: List[Dict[str, object]] = []
    examples: List[Dict[str, object]] = []
    skipped: List[Dict[str, str]] = []
    for source in track(
        requested,
        enabled=progress,
        description="Validation sources",
        total=len(requested),
        unit="sources",
    ):
        LOGGER.info("Validation: loading truth source %s", source)
        if source not in VALIDATION_CONVERSIONS:
            raise ValueError(f"No validation conversions registered for source: {source}")
        leaked = source in built_sources
        if leaked and require_held_out:
            skipped.append(
                {
                    "source": source,
                    "reason": "source was included in the database build",
                }
            )
            continue
        truth_by_conversion = load_truth(
            source,
            available[source],
            VALIDATION_CONVERSIONS[source],
            sample=sample,
        )
        for conversion, truth in truth_by_conversion.items():
            input_type, output_type = conversion
            LOGGER.info(
                "Validation: %s %s -> %s using %d truth inputs",
                source, input_type, output_type, len(truth),
            )
            metrics, disagreements = evaluate_truth(
                converter,
                truth,
                input_type,
                output_type,
                mapping_policy=mapping_policy,
                max_examples=max_examples,
                progress=progress,
            )
            metrics.update(
                {
                    "truth_source": source,
                    "independent": not leaked,
                    "truth_sampling": (
                        "complete" if sample is None else f"first {sample} unique inputs"
                    ),
                }
            )
            rows.append(metrics)
            for example in disagreements:
                example["truth_source"] = source
                example["independent"] = not leaked
                examples.append(example)

    report: Dict[str, object] = {
        "database": str(database),
        "raw_dir": str(raw_dir),
        "mapping_policy": mapping_policy,
        "database_metadata": metadata,
        "diagnostics": database_diagnostics(database),
        "benchmarks": rows,
        "disagreement_examples": examples,
        "skipped": skipped,
        "interpretation": {
            "independent": "true only when the truth source was excluded from the build",
            "precision": "fraction of returned mappings present in the truth source",
            "recall": "fraction of truth-source mappings returned by GeneMapKit",
            "overmapped": "input returned at least one mapping absent from truth",
            "warning": (
                "A disagreement can reflect source-version or scope differences; "
                "examples require review before publication claims."
            ),
        },
    }
    if output_dir:
        write_validation_report(Path(output_dir), report)
    LOGGER.info("Validation complete: %d benchmark rows", len(rows))
    return report


def write_validation_report(output_dir: Path, report: Mapping[str, object]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "validation_metrics.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    benchmarks = list(report.get("benchmarks", []))
    examples = list(report.get("disagreement_examples", []))
    if benchmarks:
        with (output_dir / "validation_benchmarks.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(benchmarks[0]))
            writer.writeheader()
            writer.writerows(benchmarks)
    if examples:
        with (output_dir / "validation_disagreements.csv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(examples[0]))
            writer.writeheader()
            writer.writerows(examples)
    summary = [
        "# GeneMapKit Validation Summary",
        "",
        f"Mapping policy: `{report.get('mapping_policy', 'all-supported')}`",
        "",
        "Truth sources marked `independent=true` were excluded from the database build.",
        "Other rows diagnose behavior but are not independent accuracy estimates.",
        "Intervals are Wilson 95% confidence intervals.",
        "",
        "| Truth source | Independent | Conversion | Inputs | Precision (95% CI) | Recall (95% CI) | Over-mapped | Exact sets |",
        "|---|---:|---|---:|---:|---:|---:|---:|",
    ]
    for row in benchmarks:
        summary.append(
            f"| {row['truth_source']} | {row['independent']} | "
            f"`{row['input_type']} -> {row['output_type']}` | "
            f"{row['tested_inputs']} | {row['precision_pct']}% "
            f"[{row['precision_ci_low_pct']}, {row['precision_ci_high_pct']}] | "
            f"{row['recall_pct']}% "
            f"[{row['recall_ci_low_pct']}, {row['recall_ci_high_pct']}] | "
            f"{row['overmapped_input_pct']}% | "
            f"{row['exact_set_match_pct']}% |"
        )
    summary.extend(["", "## Structural Diagnostics", ""])
    for key, value in report.get("diagnostics", {}).items():
        summary.append(f"- `{key}`: {value}")
    (output_dir / "validation_summary.md").write_text(
        "\n".join(summary) + "\n", encoding="utf-8"
    )
