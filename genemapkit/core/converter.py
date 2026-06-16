"""SQLite-backed identifier conversion for a built GeneMapKit snapshot."""

from __future__ import annotations

import json
import sqlite3
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import pandas as pd
from genemapkit.utils.progress import track

from .builder import BUILD_SCHEMA_VERSION, normalize_identifier
from .downloader import ID_TYPES, ID_TYPE_ALIASES


SUPPORTED_ID_TYPES = tuple(ID_TYPES)
SYMBOL_POLICIES = {
    "approved": ("symbol",),
    "historical": ("symbol", "previous_symbol"),
    "all": ("symbol", "previous_symbol", "alias_symbol"),
}
MAPPING_POLICIES = (
    "all-supported",
    "direct-only",
    "corroborated",
    "mane-select",
    "authoritative",
)
AUTHORITATIVE_SOURCES = {
    "symbol": ("hgnc", "ncbi_gene_info", "mane_summary", "ccds", "uniprot_human"),
    "alias_symbol": ("hgnc", "ncbi_gene_info"),
    "previous_symbol": ("hgnc",),
    "hgnc_id": ("hgnc", "mane_summary", "ncbi_gene_info"),
    "entrez_id": (
        "ncbi_gene_info",
        "ncbi_gene2ensembl",
        "ncbi_gene2refseq",
        "mane_summary",
        "hgnc",
        "ensembl_entrez",
        "ccds",
        "uniprot_human",
    ),
    "ensembl_gene_id": (
        "ensembl_gtf",
        "ensembl_entrez",
        "ensembl_refseq",
        "ensembl_uniprot",
        "ncbi_gene2ensembl",
        "mane_summary",
        "hgnc",
    ),
    "ensembl_transcript_id": (
        "ensembl_gtf",
        "ensembl_refseq",
        "ensembl_uniprot",
        "ncbi_gene2ensembl",
        "mane_summary",
    ),
    "ensembl_protein_id": (
        "ensembl_refseq",
        "ensembl_uniprot",
        "ncbi_gene2ensembl",
        "mane_summary",
        "uniprot_human",
    ),
    "refseq_rna_id": ("ncbi_gene2refseq", "mane_summary", "ensembl_refseq", "hgnc"),
    "refseq_protein_id": (
        "ncbi_gene2refseq",
        "mane_summary",
        "ensembl_refseq",
        "uniprot_human",
    ),
    "uniprot_id": ("uniprot_human", "ensembl_uniprot", "hgnc", "sifts_pdb"),
    "pdb_id": ("sifts_pdb", "uniprot_human"),
    "go_id": ("ncbi_gene2go",),
    "ccds_id": ("ccds", "hgnc"),
    "omim_id": ("hgnc", "ncbi_gene_info"),
    "ucsc_id": ("hgnc",),
    "vega_id": ("hgnc",),
    "ena_id": ("hgnc",),
}


def canonical_id_type(id_type: str) -> str:
    """Return the canonical Stage 1 namespace name."""
    canonical = ID_TYPE_ALIASES.get(str(id_type).strip(), str(id_type).strip())
    if canonical not in SUPPORTED_ID_TYPES:
        raise ValueError(
            f"Unsupported identifier type: {id_type}. "
            f"Supported types: {', '.join(SUPPORTED_ID_TYPES)}"
        )
    return canonical


def canonical_mapping_policy(policy: str) -> str:
    """Validate and normalize a conversion mapping policy."""
    normalized = str(policy).strip().lower()
    if normalized not in MAPPING_POLICIES:
        raise ValueError(
            f"Unknown mapping policy: {policy}. "
            f"Choose from: {', '.join(MAPPING_POLICIES)}"
        )
    return normalized


@dataclass(frozen=True)
class MappingResult:
    """One supported mapping with its source-level provenance."""

    input_id: str
    normalized_input: str
    requested_input_type: str
    matched_input_type: str
    output_id: Optional[str]
    output_type: str
    mapping_status: str
    mapping_count: int
    mapping_policy: str = "all-supported"
    selection_reason: Optional[str] = None
    source_support_count: int = 0
    supporting_sources: Optional[str] = None
    agreement_status: Optional[str] = None
    input_source: Optional[str] = None
    input_evidence: Optional[str] = None
    output_source: Optional[str] = None
    output_evidence: Optional[str] = None

    def as_dict(self) -> Dict[str, object]:
        return asdict(self)


class GeneConverter:
    """Query a normalized Stage 1 SQLite identifier index."""

    def __init__(
        self,
        database: str | Path = "data/build/genemapkit.db",
        *,
        symbol_policy: str = "all",
        mapping_policy: str = "all-supported",
    ):
        self.db_path = Path(database)
        if not self.db_path.is_file():
            raise FileNotFoundError(
                f"GeneMapKit database not found: {self.db_path}. "
                "Run `python -m genemapkit.cli build` first."
            )
        if symbol_policy not in SYMBOL_POLICIES:
            raise ValueError(
                f"Unknown symbol policy: {symbol_policy}. "
                f"Choose from: {', '.join(SYMBOL_POLICIES)}"
            )
        self.symbol_policy = symbol_policy
        self.mapping_policy = canonical_mapping_policy(mapping_policy)
        self._validate_database()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(
            f"file:{self.db_path.resolve()}?mode=ro&immutable=1",
            uri=True,
            timeout=60.0,
        )
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only = ON")
        return connection

    def open_connection(self) -> sqlite3.Connection:
        """Open a read-only SQLite connection for repeated conversion queries."""
        return self._connect()

    def _validate_database(self) -> None:
        with self._connect() as connection:
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master WHERE type = 'table'"
                )
            }
            required = {"mappings", "relationships", "metadata"}
            if not required.issubset(tables):
                raise ValueError(
                    f"Unsupported GeneMapKit database schema; missing tables: "
                    f"{', '.join(sorted(required - tables))}"
                )
            version_row = connection.execute(
                "SELECT value FROM metadata WHERE key = 'build_schema_version'"
            ).fetchone()
            if not version_row or version_row[0] != BUILD_SCHEMA_VERSION:
                found = version_row[0] if version_row else "(missing)"
                raise ValueError(
                    f"Database schema version {found} is incompatible with converter "
                    f"schema version {BUILD_SCHEMA_VERSION}"
                )

    def metadata(self) -> Dict[str, object]:
        with self._connect() as connection:
            metadata = {}
            for key, value in connection.execute("SELECT key, value FROM metadata"):
                try:
                    metadata[key] = json.loads(value)
                except (json.JSONDecodeError, TypeError):
                    metadata[key] = value
            return metadata

    def _input_types(self, input_type: str, symbol_policy: Optional[str]) -> Sequence[str]:
        canonical = canonical_id_type(input_type)
        if canonical != "symbol":
            return (canonical,)
        policy = symbol_policy or self.symbol_policy
        if policy not in SYMBOL_POLICIES:
            raise ValueError(
                f"Unknown symbol policy: {policy}. Choose from: "
                f"{', '.join(SYMBOL_POLICIES)}"
            )
        return SYMBOL_POLICIES[policy]

    def query(
        self,
        value: object,
        input_type: str,
        output_type: str,
        *,
        symbol_policy: Optional[str] = None,
        mapping_policy: Optional[str] = None,
        include_retired: bool = False,
        connection: Optional[sqlite3.Connection] = None,
    ) -> List[MappingResult]:
        """Return all exact mappings and provenance for one target namespace."""
        requested_type = canonical_id_type(input_type)
        target_type = canonical_id_type(output_type)
        raw_value = str(value).strip()
        if not raw_value:
            return []
        normalized = normalize_identifier(requested_type, raw_value)
        policy = canonical_mapping_policy(mapping_policy or self.mapping_policy)
        from_types = self._input_types(requested_type, symbol_policy)
        placeholders = ",".join("?" for _ in from_types)
        parameters = (*from_types, normalized, target_type)
        retired_filter = "" if include_retired else "AND relationship.evidence != 'retired'"
        sql = f"""
            SELECT DISTINCT
                relationship.source_type AS matched_input_type,
                relationship.target_value AS output_id,
                relationship.source AS input_source,
                relationship.evidence AS input_evidence,
                relationship.source AS output_source,
                relationship.evidence AS output_evidence
            FROM relationships AS relationship
            WHERE relationship.source_type IN ({placeholders})
              AND relationship.source_value = ?
              AND relationship.target_type = ?
              {retired_filter}
            ORDER BY relationship.target_value, relationship.source_type,
                     relationship.source, relationship.evidence
        """
        if connection is None:
            with self._connect() as opened_connection:
                rows = opened_connection.execute(sql, parameters).fetchall()
        else:
            rows = connection.execute(sql, parameters).fetchall()

        support_by_output: Dict[str, set[str]] = {}
        mane_select_outputs = set()
        for row in rows:
            support_by_output.setdefault(row["output_id"], set()).add(
                row["output_source"]
            )
            if (
                row["output_source"] == "mane_summary"
                and row["output_evidence"].strip().casefold() == "mane select"
            ):
                mane_select_outputs.add(row["output_id"])

        selection_reason = "all direct source assertions"
        if policy == "corroborated":
            rows = [
                row
                for row in rows
                if len(support_by_output.get(row["output_id"], set())) >= 2
            ]
            selection_reason = "supported by at least two independent sources"
        elif policy == "mane-select":
            rows = [
                row
                for row in rows
                if row["output_source"] == "mane_summary"
                and row["output_evidence"].strip().casefold() == "mane select"
            ]
            selection_reason = "direct MANE Select assertion"
        elif policy == "authoritative":
            priorities = AUTHORITATIVE_SOURCES.get(target_type, ())
            selected_source = next(
                (
                    source
                    for source in priorities
                    if any(row["output_source"] == source for row in rows)
                ),
                None,
            )
            if selected_source:
                rows = [
                    row for row in rows if row["output_source"] == selected_source
                ]
                selection_reason = (
                    f"highest-priority available authority for {target_type}: "
                    f"{selected_source}"
                )
            else:
                selection_reason = (
                    f"no configured authority for {target_type}; "
                    "returned all direct source assertions"
                )
        elif policy == "direct-only":
            selection_reason = "direct source assertions only"

        target_count = len({row["output_id"] for row in rows})
        status = "unique" if target_count == 1 else "ambiguous"
        results = []
        for row in rows:
            supporting_sources = sorted(support_by_output.get(row["output_id"], set()))
            if row["output_id"] in mane_select_outputs:
                agreement_status = "canonical_mane"
            elif len(supporting_sources) >= 2:
                agreement_status = "corroborated"
            else:
                agreement_status = "single_source"
            results.append(
                MappingResult(
                input_id=raw_value,
                normalized_input=normalized,
                requested_input_type=requested_type,
                matched_input_type=row["matched_input_type"],
                output_id=row["output_id"],
                output_type=target_type,
                mapping_status=status,
                mapping_count=target_count,
                mapping_policy=policy,
                selection_reason=selection_reason,
                source_support_count=len(supporting_sources),
                supporting_sources=";".join(supporting_sources),
                agreement_status=agreement_status,
                input_source=row["input_source"],
                input_evidence=row["input_evidence"],
                output_source=row["output_source"],
                output_evidence=row["output_evidence"],
            )
            )
        return results

    def convert_single(
        self,
        gene_id: object,
        input_type: str,
        output_types: Sequence[str],
        *,
        symbol_policy: Optional[str] = None,
        mapping_policy: Optional[str] = None,
        include_retired: bool = False,
        connection: Optional[sqlite3.Connection] = None,
    ) -> Dict[str, List[str]]:
        """Return all distinct output values for each requested namespace."""
        converted: Dict[str, List[str]] = {}
        for output_type in output_types:
            canonical = canonical_id_type(output_type)
            converted[canonical] = sorted(
                {
                    result.output_id
                    for result in self.query(
                        gene_id,
                        input_type,
                        canonical,
                        symbol_policy=symbol_policy,
                        mapping_policy=mapping_policy,
                        include_retired=include_retired,
                        connection=connection,
                    )
                    if result.output_id is not None
                }
            )
        return converted

    def convert_batch(
        self,
        gene_ids: Iterable[object],
        input_type: str,
        output_types: Sequence[str],
        *,
        include_unmapped: bool = True,
        output_format: str = "wide",
        separator: str = ";",
        symbol_policy: Optional[str] = None,
        mapping_policy: Optional[str] = None,
        include_retired: bool = False,
        progress: bool = False,
    ) -> pd.DataFrame:
        """Convert identifiers to provenance-rich long or multi-value wide rows."""
        if output_format not in {"long", "wide"}:
            raise ValueError("output_format must be 'long' or 'wide'")
        canonical_input = canonical_id_type(input_type)
        canonical_outputs = [canonical_id_type(value) for value in output_types]
        rows: List[Dict[str, object]] = []

        total = len(gene_ids) if hasattr(gene_ids, "__len__") else None
        inputs = track(
            gene_ids,
            enabled=progress,
            description="Converting identifiers",
            total=total,
            unit="ids",
        )
        with self.open_connection() as connection:
            for input_row, input_value in enumerate(inputs):
                raw_value = "" if pd.isna(input_value) else str(input_value).strip()
                if output_format == "long":
                    mapped = False
                    for output_type in canonical_outputs:
                        results = self.query(
                            raw_value,
                            canonical_input,
                            output_type,
                            symbol_policy=symbol_policy,
                            mapping_policy=mapping_policy,
                            include_retired=include_retired,
                            connection=connection,
                        )
                        if results:
                            mapped = True
                            for result in results:
                                result_row = result.as_dict()
                                result_row["input_row"] = input_row
                                rows.append(result_row)
                        elif include_unmapped:
                            result_row = MappingResult(
                                input_id=raw_value,
                                normalized_input=normalize_identifier(
                                    canonical_input, raw_value
                                ),
                                requested_input_type=canonical_input,
                                matched_input_type=canonical_input,
                                output_id=None,
                                output_type=output_type,
                                mapping_status="unmapped",
                                mapping_count=0,
                                mapping_policy=canonical_mapping_policy(
                                    mapping_policy or self.mapping_policy
                                ),
                                selection_reason="no direct mapping matched the policy",
                                agreement_status="unmapped",
                            ).as_dict()
                            result_row["input_row"] = input_row
                            rows.append(result_row)
                    if not mapped and not include_unmapped:
                        continue
                else:
                    row: Dict[str, object] = {
                        "input_row": input_row,
                        f"input_{canonical_input}": raw_value,
                    }
                    any_mapping = False
                    for output_type in canonical_outputs:
                        results = self.query(
                            raw_value,
                            canonical_input,
                            output_type,
                            symbol_policy=symbol_policy,
                            mapping_policy=mapping_policy,
                            include_retired=include_retired,
                            connection=connection,
                        )
                        values = sorted(
                            {result.output_id for result in results if result.output_id}
                        )
                        any_mapping = any_mapping or bool(values)
                        row[output_type] = separator.join(values)
                        row[f"{output_type}_mapping_count"] = len(values)
                        row[f"{output_type}_mapping_status"] = (
                            "unmapped"
                            if not values
                            else "unique"
                            if len(values) == 1
                            else "ambiguous"
                        )
                    if include_unmapped or any_mapping:
                        rows.append(row)
        return pd.DataFrame(rows)

    def get_statistics(self) -> Dict[str, object]:
        with self._connect() as connection:
            return {
                "database": str(self.db_path),
                "mapping_rows": connection.execute(
                    "SELECT COUNT(*) FROM mappings"
                ).fetchone()[0],
                "relationship_rows": connection.execute(
                    "SELECT COUNT(*) FROM relationships"
                ).fetchone()[0],
                "groups": connection.execute(
                    "SELECT COUNT(DISTINCT group_id) FROM mappings"
                ).fetchone()[0],
                "id_types": dict(
                    connection.execute(
                        "SELECT id_type, COUNT(DISTINCT id_value) "
                        "FROM mappings GROUP BY id_type ORDER BY id_type"
                    )
                ),
            }
