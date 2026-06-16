"""Build a deterministic SQLite identifier index from a GeneMapKit snapshot."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from .downloader import ID_TYPES, sha256sum
from genemapkit.utils.progress import track


LOGGER = logging.getLogger(__name__)
BUILD_SCHEMA_VERSION = "2.1.0"
HUMAN_TAXID = "9606"

ANCHOR_TYPES = {
    "hgnc_id",
    "entrez_id",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_protein_id",
    "refseq_rna_id",
    "refseq_protein_id",
    "uniprot_id",
    "ccds_id",
}

LABEL_TYPES = set(ID_TYPES) - ANCHOR_TYPES
Identifier = Tuple[str, str]
Record = Tuple[List[Identifier], str]
Parser = Callable[[Path], Iterator[Record]]

GENE_TYPES = {
    "symbol",
    "alias_symbol",
    "previous_symbol",
    "hgnc_id",
    "entrez_id",
    "ensembl_gene_id",
    "ccds_id",
    "omim_id",
    "ucsc_id",
    "vega_id",
    "ena_id",
}
TRANSCRIPT_TYPES = {"ensembl_transcript_id", "refseq_rna_id"}
PROTEIN_TYPES = {"ensembl_protein_id", "refseq_protein_id", "uniprot_id", "pdb_id"}
ANNOTATION_TYPES = {"go_id"}
PRECISE_TRANSCRIPT_PROTEIN_SOURCES = {
    "ensembl_gtf",
    "mane_summary",
    "ensembl_entrez",
    "ensembl_refseq",
    "ensembl_uniprot",
    "ncbi_gene2ensembl",
    "uniprot_human",
}

ENSEMBL_VERSION = re.compile(r"^(ENS[A-Z]*[GTP]\d+)\.\d+$")
REFSEQ_VERSION = re.compile(r"^([A-Z]{2}_\d+)\.\d+$")
CCDS_VERSION = re.compile(r"^(CCDS\d+)\.\d+$", re.IGNORECASE)
GTF_ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')


def normalize_identifier(id_type: str, value: object) -> str:
    """Normalize an identifier for lookup while preserving original values."""
    normalized = str(value).strip()
    if not normalized or normalized == "-":
        return ""
    if id_type.startswith("ensembl_"):
        match = ENSEMBL_VERSION.match(normalized)
        return match.group(1) if match else normalized
    if id_type.startswith("refseq_"):
        match = REFSEQ_VERSION.match(normalized)
        return match.group(1) if match else normalized
    if id_type == "ccds_id":
        match = CCDS_VERSION.match(normalized)
        return (match.group(1) if match else normalized).upper()
    if id_type in {"hgnc_id", "uniprot_id", "go_id", "pdb_id"}:
        return normalized.upper()
    return normalized


def classify_refseq(value: str) -> Optional[str]:
    accession = normalize_identifier("refseq_rna_id", value)
    prefix = accession.split("_", 1)[0]
    if prefix in {"NM", "XM", "NR", "XR"}:
        return "refseq_rna_id"
    if prefix in {"NP", "XP", "YP", "AP", "WP"}:
        return "refseq_protein_id"
    return None


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def split_values(value: object, separators: str = r"[|;]") -> List[str]:
    if value is None:
        return []
    return [
        part.strip()
        for part in re.split(separators, str(value))
        if part.strip() and part.strip() != "-"
    ]


def clean_record(values: Iterable[Identifier]) -> List[Identifier]:
    return sorted(
        {
            (id_type, str(value).strip())
            for id_type, value in values
            if id_type in ID_TYPES and str(value).strip() not in {"", "-"}
        }
    )


def relationship_type(left_type: str, right_type: str) -> str:
    levels = {
        "gene" if value in GENE_TYPES
        else "transcript" if value in TRANSCRIPT_TYPES
        else "protein" if value in PROTEIN_TYPES
        else "annotation" if value in ANNOTATION_TYPES
        else "other"
        for value in (left_type, right_type)
    }
    if levels == {"gene"}:
        return "gene_equivalent"
    if levels == {"gene", "transcript"}:
        return "gene_has_transcript"
    if levels == {"gene", "protein"}:
        return "gene_has_protein"
    if levels == {"gene", "annotation"}:
        return "gene_has_annotation"
    if levels == {"transcript"}:
        return "transcript_equivalent"
    if levels == {"protein"}:
        return "protein_equivalent"
    if levels == {"transcript", "protein"}:
        return "transcript_encodes_protein"
    return "associated"


def direct_relationships(
    source: str, record: Iterable[Identifier]
) -> List[Tuple[str, str, str, str, str]]:
    """Return direct source-asserted relationships without gene-group transitivity."""
    identifiers = sorted(
        {
            (id_type, normalize_identifier(id_type, value))
            for id_type, value in record
            if normalize_identifier(id_type, value)
        }
    )
    relationships = set()
    for index, left in enumerate(identifiers):
        for right in identifiers[index + 1 :]:
            left_type, left_value = left
            right_type, right_value = right
            if left_type == right_type and not (
                source == "ncbi_gene_history" and left_type == "entrez_id"
            ):
                continue
            if source == "uniprot_human" and "uniprot_id" not in {
                left_type,
                right_type,
            }:
                continue
            left_level = (
                "gene" if left_type in GENE_TYPES
                else "transcript" if left_type in TRANSCRIPT_TYPES
                else "protein" if left_type in PROTEIN_TYPES
                else "annotation" if left_type in ANNOTATION_TYPES
                else "other"
            )
            right_level = (
                "gene" if right_type in GENE_TYPES
                else "transcript" if right_type in TRANSCRIPT_TYPES
                else "protein" if right_type in PROTEIN_TYPES
                else "annotation" if right_type in ANNOTATION_TYPES
                else "other"
            )
            levels = {left_level, right_level}
            allowed = (
                left_level == right_level
                or "gene" in levels
                or (
                    levels == {"transcript", "protein"}
                    and source in PRECISE_TRANSCRIPT_PROTEIN_SOURCES
                )
            )
            if not allowed:
                continue
            rel_type = relationship_type(left_type, right_type)
            relationships.add(
                (left_type, left_value, right_type, right_value, rel_type)
            )
            relationships.add(
                (right_type, right_value, left_type, left_value, rel_type)
            )
    return sorted(relationships)


def parse_hgnc(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"hgnc_id", "symbol"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"HGNC is missing required columns: {sorted(required)}")
        for row in reader:
            values: List[Identifier] = [
                ("hgnc_id", row.get("hgnc_id", "")),
                ("symbol", row.get("symbol", "")),
                ("entrez_id", row.get("entrez_id", "")),
                ("ensembl_gene_id", row.get("ensembl_gene_id", "")),
                ("ucsc_id", row.get("ucsc_id", "")),
                ("vega_id", row.get("vega_id", "")),
            ]
            for column, id_type in (
                ("alias_symbol", "alias_symbol"),
                ("prev_symbol", "previous_symbol"),
                ("refseq_accession", "refseq_rna_id"),
                ("uniprot_ids", "uniprot_id"),
                ("ccds_id", "ccds_id"),
                ("omim_id", "omim_id"),
                ("ena", "ena_id"),
            ):
                values.extend((id_type, value) for value in split_values(row.get(column)))
            yield clean_record(values), "authority"


def parse_ncbi_gene_info(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            values: List[Identifier] = [
                ("entrez_id", row.get("GeneID", "")),
                ("symbol", row.get("Symbol", "")),
            ]
            values.extend(
                ("alias_symbol", value)
                for value in split_values(row.get("Synonyms"), r"\|")
            )
            for xref in split_values(row.get("dbXrefs"), r"\|"):
                if xref.startswith("Ensembl:"):
                    values.append(("ensembl_gene_id", xref.removeprefix("Ensembl:")))
                elif xref.startswith("HGNC:"):
                    values.append(("hgnc_id", xref.removeprefix("HGNC:")))
                elif xref.startswith("MIM:"):
                    values.append(("omim_id", xref.removeprefix("MIM:")))
            yield clean_record(values), "authority"


def parse_ensembl_gtf(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in {"gene", "transcript", "CDS"}:
                continue
            attributes = dict(GTF_ATTRIBUTE.findall(fields[8]))
            values = [
                ("ensembl_gene_id", attributes.get("gene_id", "")),
                ("ensembl_transcript_id", attributes.get("transcript_id", "")),
                ("ensembl_protein_id", attributes.get("protein_id", "")),
                ("symbol", attributes.get("gene_name", "")),
            ]
            yield clean_record(values), "annotation"


def parse_uniprot_selected(path: Path) -> Iterator[Record]:
    columns = {
        "uniprot_id": 0,
        "entrez_id": 2,
        "refseq": 3,
        "pdb_id": 5,
        "go_id": 6,
        "omim_id": 13,
        "ensembl_gene_id": 18,
        "ensembl_transcript_id": 19,
        "ensembl_protein_id": 20,
    }
    with open_text(path) as handle:
        for line in handle:
            row = line.rstrip("\n").split("\t")
            if len(row) <= columns["ensembl_protein_id"]:
                continue
            values: List[Identifier] = [("uniprot_id", row[columns["uniprot_id"]])]
            for id_type in (
                "entrez_id",
                "pdb_id",
                "go_id",
                "omim_id",
                "ensembl_gene_id",
                "ensembl_transcript_id",
                "ensembl_protein_id",
            ):
                values.extend(
                    (id_type, value)
                    for value in split_values(row[columns[id_type]], r";\s*")
                )
            for accession in split_values(row[columns["refseq"]], r";\s*"):
                id_type = classify_refseq(accession)
                if id_type:
                    values.append((id_type, accession))
            yield clean_record(values), "cross_reference"


def parse_mane(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            values = [
                ("entrez_id", row.get("NCBI_GeneID", "").removeprefix("GeneID:")),
                ("ensembl_gene_id", row.get("Ensembl_Gene", "")),
                ("hgnc_id", row.get("HGNC_ID", "")),
                ("symbol", row.get("symbol", "")),
                ("refseq_rna_id", row.get("RefSeq_nuc", "")),
                ("refseq_protein_id", row.get("RefSeq_prot", "")),
                ("ensembl_transcript_id", row.get("Ensembl_nuc", "")),
                ("ensembl_protein_id", row.get("Ensembl_prot", "")),
            ]
            yield clean_record(values), row.get("MANE_status", "MANE")


def parse_ensembl_xref(path: Path, target: str) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            xref = row.get("xref", "")
            if target == "entrez_id" and (
                row.get("db_name") != "EntrezGene" or not xref.isdigit()
            ):
                continue
            values: List[Identifier] = [
                ("ensembl_gene_id", row.get("gene_stable_id", "")),
                ("ensembl_transcript_id", row.get("transcript_stable_id", "")),
                ("ensembl_protein_id", row.get("protein_stable_id", "")),
            ]
            if target == "refseq":
                id_type = classify_refseq(xref)
                if id_type:
                    values.append((id_type, xref))
            else:
                values.append((target, xref))
            evidence = row.get("info_type", "cross_reference")
            yield clean_record(values), evidence


def collect_rejected_records(
    source: str, path: Path, max_records: Optional[int], max_examples: int = 10
) -> Dict[str, object]:
    """Count known source rows rejected because their namespace is incompatible."""
    if source != "ensembl_entrez":
        return {"count": 0, "reasons": {}, "examples": []}
    count = 0
    reasons: Dict[str, int] = defaultdict(int)
    examples = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for index, row in enumerate(reader, start=1):
            if max_records and index > max_records:
                break
            xref = row.get("xref", "")
            if row.get("db_name") == "EntrezGene" and xref.isdigit():
                continue
            reason = "not_numeric_entrezgene"
            count += 1
            reasons[reason] += 1
            if len(examples) < max_examples:
                examples.append(
                    {
                        "xref": xref,
                        "db_name": row.get("db_name", ""),
                        "info_type": row.get("info_type", ""),
                    }
                )
    return {"count": count, "reasons": dict(reasons), "examples": examples}


def parse_ccds(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            values = [
                ("ccds_id", row.get("ccds_id", "")),
                ("entrez_id", row.get("gene_id", "")),
                ("symbol", row.get("gene", "")),
            ]
            yield clean_record(values), row.get("ccds_status", "CCDS")


def parse_gene_history(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("tax_id") != HUMAN_TAXID:
                continue
            values = [
                ("entrez_id", row.get("GeneID", "")),
                ("entrez_id", row.get("Discontinued_GeneID", "")),
                ("previous_symbol", row.get("Discontinued_Symbol", "")),
            ]
            yield clean_record(values), "retired"


def parse_gene2ensembl(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("tax_id") != HUMAN_TAXID:
                continue
            values = [
                ("entrez_id", row.get("GeneID", "")),
                ("ensembl_gene_id", row.get("Ensembl_gene_identifier", "")),
                ("refseq_rna_id", row.get("RNA_nucleotide_accession.version", "")),
                ("ensembl_transcript_id", row.get("Ensembl_rna_identifier", "")),
                ("refseq_protein_id", row.get("protein_accession.version", "")),
                ("ensembl_protein_id", row.get("Ensembl_protein_identifier", "")),
            ]
            yield clean_record(values), "validation"


def parse_gene2refseq(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("tax_id") != HUMAN_TAXID:
                continue
            values = [
                ("entrez_id", row.get("GeneID", "")),
                ("refseq_rna_id", row.get("RNA_nucleotide_accession.version", "")),
                ("refseq_protein_id", row.get("protein_accession.version", "")),
            ]
            yield clean_record(values), "validation"


def parse_gene2go(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames:
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            if row.get("tax_id") != HUMAN_TAXID:
                continue
            yield clean_record(
                [("entrez_id", row.get("GeneID", "")), ("go_id", row.get("GO_ID", ""))]
            ), row.get("Evidence", "GO")


def parse_sifts(path: Path) -> Iterator[Record]:
    with open_text(path) as handle:
        header = None
        for line in handle:
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            yield clean_record(
                [
                    ("pdb_id", row.get("PDB", "")),
                    ("uniprot_id", row.get("SP_PRIMARY", "")),
                ]
            ), "SIFTS"


PARSERS: Dict[str, Parser] = {
    "hgnc": parse_hgnc,
    "ncbi_gene_info": parse_ncbi_gene_info,
    "ensembl_gtf": parse_ensembl_gtf,
    "uniprot_human": parse_uniprot_selected,
    "mane_summary": parse_mane,
    "ensembl_entrez": lambda path: parse_ensembl_xref(path, "entrez_id"),
    "ensembl_refseq": lambda path: parse_ensembl_xref(path, "refseq"),
    "ensembl_uniprot": lambda path: parse_ensembl_xref(path, "uniprot_id"),
    "ccds": parse_ccds,
    "ncbi_gene_history": parse_gene_history,
    "ncbi_gene2ensembl": parse_gene2ensembl,
    "ncbi_gene2refseq": parse_gene2refseq,
    "ncbi_gene2go": parse_gene2go,
    "sifts_pdb": parse_sifts,
}


class UnionFind:
    """Deterministic union-find over normalized stable identifiers."""

    def __init__(self):
        self.parent: Dict[Identifier, Identifier] = {}

    def find(self, item: Identifier) -> Identifier:
        self.parent.setdefault(item, item)
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != root:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: Identifier, right: Identifier) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            smaller, larger = sorted((left_root, right_root))
            self.parent[larger] = smaller


def load_manifest(raw_dir: Path) -> Dict:
    path = raw_dir / "manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Snapshot manifest not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def select_sources(
    raw_dir: Path,
    manifest: Dict,
    sources: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
) -> Dict[str, Path]:
    databases = manifest.get("databases", {})
    selected = sorted(sources) if sources else sorted(databases)
    excluded = set(exclude or [])
    unknown = sorted(set(selected) - set(PARSERS))
    if unknown:
        raise ValueError(f"No builder parser is registered for: {', '.join(unknown)}")
    result = {}
    for source in selected:
        if source in excluded:
            continue
        filename = databases.get(source, {}).get("filename")
        if filename and (raw_dir / filename).exists():
            result[source] = raw_dir / filename
    return result


def verify_snapshot(
    raw_dir: Path, manifest: Dict, selected: Dict[str, Path], strict: bool = True
) -> None:
    databases = manifest.get("databases", {})
    failures = []
    for source, path in selected.items():
        expected = databases[source].get("sha256")
        if not expected:
            failures.append(f"{source}: manifest has no SHA-256")
        elif sha256sum(path).lower() != expected.lower():
            failures.append(f"{source}: SHA-256 mismatch")
        expected_bytes = databases[source].get("bytes")
        if expected_bytes and path.stat().st_size != int(expected_bytes):
            failures.append(f"{source}: file-size mismatch")
    if failures:
        message = "Snapshot verification failed:\n- " + "\n- ".join(failures)
        if strict:
            raise ValueError(message)
        LOGGER.warning(message)


def iter_limited(parser: Parser, path: Path, max_records: Optional[int]):
    for index, record in enumerate(parser(path), start=1):
        if max_records and index > max_records:
            break
        yield record


def prepare_database(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode = DELETE;
        PRAGMA synchronous = NORMAL;
        PRAGMA temp_store = MEMORY;
        CREATE TABLE mappings (
            group_id INTEGER NOT NULL,
            id_type TEXT NOT NULL,
            id_value TEXT NOT NULL,
            original_value TEXT NOT NULL,
            source TEXT NOT NULL,
            evidence TEXT NOT NULL,
            UNIQUE(group_id, id_type, id_value, original_value, source, evidence)
        );
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE source_summary (
            source TEXT PRIMARY KEY,
            records INTEGER NOT NULL,
            records_without_anchor INTEGER NOT NULL,
            mapping_rows INTEGER NOT NULL,
            relationship_rows INTEGER NOT NULL
        );
        CREATE TABLE relationships (
            source_type TEXT NOT NULL,
            source_value TEXT NOT NULL,
            target_type TEXT NOT NULL,
            target_value TEXT NOT NULL,
            relationship_type TEXT NOT NULL,
            source TEXT NOT NULL,
            evidence TEXT NOT NULL,
            UNIQUE(source_type, source_value, target_type, target_value,
                   relationship_type, source, evidence)
        );
        """
    )


def build_database(
    raw_dir: Path,
    output: Path,
    *,
    sources: Optional[Sequence[str]] = None,
    exclude: Optional[Sequence[str]] = None,
    verify: bool = True,
    strict: bool = True,
    force: bool = False,
    max_records: Optional[int] = None,
    report_path: Optional[Path] = None,
    progress: bool = False,
) -> Dict:
    """Build the normalized any-to-any SQLite identifier index."""
    raw_dir, output = Path(raw_dir), Path(output)
    if max_records is not None and max_records < 1:
        raise ValueError("--max-records must be a positive integer")
    if output.exists() and not force:
        raise FileExistsError(f"Output already exists: {output}. Use --force to replace it.")
    manifest = load_manifest(raw_dir)
    selected = select_sources(raw_dir, manifest, sources=sources, exclude=exclude)
    requested = set(sources or manifest.get("databases", {})) - set(exclude or [])
    missing = sorted(requested - set(selected))
    if missing and strict:
        raise FileNotFoundError(f"Selected snapshot sources are missing: {', '.join(missing)}")
    if not selected:
        raise ValueError("No parsable source files were selected")
    if verify:
        verify_snapshot(raw_dir, manifest, selected, strict=strict)

    rejected_records = {
        source: rejection
        for source, path in selected.items()
        if (rejection := collect_rejected_records(source, path, max_records))["count"]
    }
    union_find = UnionFind()
    first_pass_counts: Dict[str, int] = {}
    no_anchor_counts: Dict[str, int] = defaultdict(int)

    for source, path in track(
        selected.items(),
        enabled=progress,
        description="Pass 1: grouping",
        total=len(selected),
        unit="sources",
    ):
        count = 0
        LOGGER.info("Grouping identifiers from %s", source)
        for record, _ in iter_limited(PARSERS[source], path, max_records):
            normalized = sorted(
                {
                    (id_type, normalize_identifier(id_type, value))
                    for id_type, value in record
                    if normalize_identifier(id_type, value)
                }
            )
            anchors = [item for item in normalized if item[0] in ANCHOR_TYPES]
            if not anchors:
                no_anchor_counts[source] += 1
            else:
                for anchor in anchors:
                    union_find.find(anchor)
                for anchor in anchors[1:]:
                    union_find.union(anchors[0], anchor)
            count += 1
        first_pass_counts[source] = count
        if count == 0 and strict:
            raise ValueError(f"Source produced zero records: {source}")

    roots = sorted({union_find.find(anchor) for anchor in union_find.parent})
    group_ids = {root: index for index, root in enumerate(roots, start=1)}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".part")
    if temporary.exists():
        temporary.unlink()
    connection = sqlite3.connect(temporary)
    prepare_database(connection)
    source_rows: Dict[str, int] = defaultdict(int)
    relationship_rows: Dict[str, int] = defaultdict(int)

    try:
        for source, path in track(
            selected.items(),
            enabled=progress,
            description="Pass 2: writing",
            total=len(selected),
            unit="sources",
        ):
            LOGGER.info("Writing mappings from %s", source)
            buffer = []
            relationship_buffer = []
            for record, evidence in iter_limited(PARSERS[source], path, max_records):
                anchors = [
                    (id_type, normalize_identifier(id_type, value))
                    for id_type, value in record
                    if id_type in ANCHOR_TYPES and normalize_identifier(id_type, value)
                ]
                if not anchors:
                    continue
                group_id = group_ids[union_find.find(anchors[0])]
                for id_type, original in record:
                    normalized = normalize_identifier(id_type, original)
                    if normalized:
                        buffer.append(
                            (group_id, id_type, normalized, original, source, evidence)
                        )
                relationship_buffer.extend(
                    (*relationship, source, evidence)
                    for relationship in direct_relationships(source, record)
                )
                if len(buffer) >= 50_000:
                    before = connection.total_changes
                    connection.executemany(
                        "INSERT OR IGNORE INTO mappings VALUES (?,?,?,?,?,?)", buffer
                    )
                    source_rows[source] += connection.total_changes - before
                    buffer.clear()
                if len(relationship_buffer) >= 50_000:
                    before = connection.total_changes
                    connection.executemany(
                        "INSERT OR IGNORE INTO relationships VALUES (?,?,?,?,?,?,?)",
                        relationship_buffer,
                    )
                    relationship_rows[source] += connection.total_changes - before
                    relationship_buffer.clear()
            if buffer:
                before = connection.total_changes
                connection.executemany(
                    "INSERT OR IGNORE INTO mappings VALUES (?,?,?,?,?,?)", buffer
                )
                source_rows[source] += connection.total_changes - before
            if relationship_buffer:
                before = connection.total_changes
                connection.executemany(
                    "INSERT OR IGNORE INTO relationships VALUES (?,?,?,?,?,?,?)",
                    relationship_buffer,
                )
                relationship_rows[source] += connection.total_changes - before
            connection.execute(
                "INSERT INTO source_summary VALUES (?,?,?,?,?)",
                (
                    source,
                    first_pass_counts[source],
                    no_anchor_counts[source],
                    source_rows[source],
                    relationship_rows[source],
                ),
            )
            connection.commit()

        metadata = {
            "build_schema_version": BUILD_SCHEMA_VERSION,
            "snapshot_schema_version": manifest.get("schema_version", ""),
            "snapshot_created_utc": manifest.get("created_utc", ""),
            "organism": manifest.get("organism", ""),
            "taxid": manifest.get("taxid", ""),
            "ensembl_release": manifest.get("ensembl_release", ""),
            "selected_sources": sorted(selected),
            "max_records_per_source": max_records or "",
            "complete_build": max_records is None,
            "relationship_model": "direct_source_assertions",
        }
        connection.executemany(
            "INSERT INTO metadata VALUES (?,?)",
            [(key, json.dumps(value) if isinstance(value, (list, bool)) else str(value))
             for key, value in metadata.items()],
        )
        LOGGER.info("Creating indexes and running ANALYZE/VACUUM (this may take several minutes)...")
        connection.executescript(
            """
            CREATE INDEX idx_lookup ON mappings(id_type, id_value, group_id);
            CREATE INDEX idx_group_type ON mappings(group_id, id_type, id_value);
            CREATE INDEX idx_source ON mappings(source);
            CREATE INDEX idx_relationship_lookup
              ON relationships(source_type, source_value, target_type, target_value);
            CREATE INDEX idx_relationship_target
              ON relationships(target_type, target_value, source_type, source_value);
            CREATE INDEX idx_relationship_source ON relationships(source);
            ANALYZE;
            VACUUM;
            """
        )
        total_rows = connection.execute("SELECT COUNT(*) FROM mappings").fetchone()[0]
        total_relationships = connection.execute(
            "SELECT COUNT(*) FROM relationships"
        ).fetchone()[0]
        type_counts = dict(
            connection.execute(
                "SELECT id_type, COUNT(DISTINCT id_value) FROM mappings GROUP BY id_type"
            ).fetchall()
        )
        connection.commit()
    finally:
        connection.close()

    if output.exists():
        output.unlink()
    temporary.replace(output)
    report = {
        "database": str(output),
        "database_sha256": sha256sum(output),
        "build_schema_version": BUILD_SCHEMA_VERSION,
        "groups": len(group_ids),
        "mapping_rows": total_rows,
        "relationship_rows": total_relationships,
        "id_type_counts": type_counts,
        "source_records": first_pass_counts,
        "source_mapping_rows": dict(source_rows),
        "source_relationship_rows": dict(relationship_rows),
        "records_without_anchor": dict(no_anchor_counts),
        "rejected_records": rejected_records,
        "missing_sources": missing,
        "verified_snapshot": verify,
        "complete_build": max_records is None,
    }
    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def convert_identifier(
    database: Path, value: str, input_type: str, output_type: str
) -> List[str]:
    """Query all target identifiers connected to an input identifier."""
    normalized = normalize_identifier(input_type, value)
    with sqlite3.connect(database) as connection:
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        if "relationships" in tables:
            return [
                row[0]
                for row in connection.execute(
                    """
                    SELECT DISTINCT target_value
                    FROM relationships
                    WHERE source_type = ? AND source_value = ? AND target_type = ?
                    ORDER BY target_value
                    """,
                    (input_type, normalized, output_type),
                )
            ]
        return [
            row[0]
            for row in connection.execute(
                """
                SELECT DISTINCT target.id_value
                FROM mappings AS source
                JOIN mappings AS target ON target.group_id = source.group_id
                WHERE source.id_type = ? AND source.id_value = ?
                  AND target.id_type = ?
                ORDER BY target.id_value
                """,
                (input_type, normalized, output_type),
            )
        ]
