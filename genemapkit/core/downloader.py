"""Reproducible raw database snapshot downloader for GeneMapKit.

This module is Stage 0 of the GeneMapKit data pipeline. It downloads immutable
upstream files into ``data/raw`` and records provenance in ``manifest.json``.
It does not parse, decompress, or combine source files.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import json
import logging
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import requests
from tqdm import tqdm


SNAPSHOT_SCHEMA_VERSION = "1.0.0"
DEFAULT_ENSEMBL_RELEASE = 116
CHUNK_SIZE = 1 << 16
TIERS = ("core", "crossref", "extended")
VALIDATION_SOURCES = ("ncbi_gene2ensembl", "ncbi_gene2refseq")

# Canonical namespaces used by source selection. Legacy names are normalized
# by ID_TYPE_ALIASES so existing callers continue to work.
ID_TYPES = (
    "symbol",
    "alias_symbol",
    "previous_symbol",
    "hgnc_id",
    "entrez_id",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "ensembl_protein_id",
    "refseq_rna_id",
    "refseq_protein_id",
    "uniprot_id",
    "ccds_id",
    "omim_id",
    "ucsc_id",
    "vega_id",
    "ena_id",
    "pdb_id",
    "go_id",
)

ID_TYPE_ALIASES = {
    "gene_symbol": "symbol",
    "aliases": "alias_symbol",
    "prev_symbol": "previous_symbol",
    "refseq_id": "refseq_rna_id",
}

# The default tiers intentionally exclude large all-species and niche sources.
# ``--all`` includes the extended tier; conversion-specific selection downloads
# only sources that add coverage for the requested namespaces.
SOURCE_REGISTRY: Dict[str, Dict] = {
    "hgnc": {
        "tier": "core",
        "priority": 1,
        "url": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt",
        "filename": "hgnc_complete_set.txt",
        "description": "HGNC complete human nomenclature and cross-reference set",
        "provides": [
            "symbol", "alias_symbol", "previous_symbol", "hgnc_id", "entrez_id",
            "ensembl_gene_id", "refseq_rna_id", "uniprot_id", "ccds_id",
            "omim_id", "ucsc_id", "vega_id", "ena_id",
        ],
    },
    "ncbi_gene_info": {
        "tier": "core",
        "priority": 1,
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        "filename": "Homo_sapiens.gene_info.gz",
        "description": "NCBI human Gene Info",
        "provides": [
            "symbol", "alias_symbol", "entrez_id", "ensembl_gene_id",
            "hgnc_id", "omim_id",
        ],
    },
    "ensembl_gtf": {
        "tier": "core",
        "priority": 1,
        "url": "https://ftp.ensembl.org/pub/release-{release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{release}.gtf.gz",
        "filename": "Homo_sapiens.GRCh38.{release}.gtf.gz",
        "description": "Ensembl GRCh38 genes and transcripts",
        "provides": ["symbol", "ensembl_gene_id", "ensembl_transcript_id"],
    },
    "uniprot_human": {
        "tier": "core",
        "priority": 1,
        "url": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz",
        "filename": "HUMAN_9606_idmapping_selected.tab.gz",
        "description": "UniProt human protein identifier mappings",
        "provides": [
            "symbol", "entrez_id", "ensembl_gene_id", "ensembl_transcript_id",
            "ensembl_protein_id", "refseq_protein_id", "uniprot_id", "pdb_id",
        ],
    },
    "mane_summary": {
        "tier": "core",
        "priority": 1,
        "url": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.5.summary.txt.gz",
        "filename": "MANE.GRCh38.v1.5.summary.txt.gz",
        "description": "MANE v1.5 high-confidence Ensembl/RefSeq transcript pairs",
        "provides": [
            "symbol", "hgnc_id", "entrez_id", "ensembl_gene_id",
            "ensembl_transcript_id", "ensembl_protein_id", "refseq_rna_id",
            "refseq_protein_id",
        ],
    },
    "ensembl_entrez": {
        "tier": "crossref",
        "priority": 2,
        "url": "https://ftp.ensembl.org/pub/release-{release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{release}.entrez.tsv.gz",
        "filename": "Homo_sapiens.GRCh38.{release}.entrez.tsv.gz",
        "description": "Ensembl to Entrez cross-reference table",
        "provides": [
            "symbol", "entrez_id", "ensembl_gene_id", "ensembl_transcript_id",
            "ensembl_protein_id",
        ],
    },
    "ensembl_refseq": {
        "tier": "crossref",
        "priority": 2,
        "url": "https://ftp.ensembl.org/pub/release-{release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{release}.refseq.tsv.gz",
        "filename": "Homo_sapiens.GRCh38.{release}.refseq.tsv.gz",
        "description": "Ensembl to RefSeq cross-reference table",
        "provides": [
            "ensembl_gene_id", "ensembl_transcript_id", "ensembl_protein_id",
            "refseq_rna_id", "refseq_protein_id",
        ],
    },
    "ensembl_uniprot": {
        "tier": "crossref",
        "priority": 2,
        "url": "https://ftp.ensembl.org/pub/release-{release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{release}.uniprot.tsv.gz",
        "filename": "Homo_sapiens.GRCh38.{release}.uniprot.tsv.gz",
        "description": "Ensembl to UniProt cross-reference table",
        "provides": [
            "ensembl_gene_id", "ensembl_transcript_id", "ensembl_protein_id",
            "uniprot_id",
        ],
    },
    "ccds": {
        "tier": "crossref",
        "priority": 2,
        "url": "https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt",
        "filename": "CCDS.current.txt",
        "description": "Current human Consensus CDS mappings",
        "provides": ["symbol", "entrez_id", "ccds_id"],
    },
    "ncbi_gene_history": {
        "tier": "crossref",
        "priority": 2,
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz",
        "filename": "gene_history.gz",
        "description": "NCBI discontinued and replaced Gene IDs; filter to taxid 9606",
        "provides": ["symbol", "entrez_id"],
        "taxon_filter": True,
    },
    "ncbi_gene2ensembl": {
        "tier": "extended",
        "priority": 3,
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz",
        "filename": "gene2ensembl.gz",
        "description": "NCBI Gene to Ensembl mappings; all species, filter to taxid 9606",
        "provides": [
            "entrez_id", "ensembl_gene_id", "ensembl_transcript_id",
            "ensembl_protein_id",
        ],
        "taxon_filter": True,
    },
    "ncbi_gene2refseq": {
        "tier": "extended",
        "priority": 3,
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz",
        "filename": "gene2refseq.gz",
        "description": "NCBI Gene to RefSeq mappings; all species, filter to taxid 9606",
        "provides": ["entrez_id", "refseq_rna_id", "refseq_protein_id"],
        "taxon_filter": True,
        "large": True,
    },
    "ncbi_gene2go": {
        "tier": "extended",
        "priority": 4,
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",
        "filename": "gene2go.gz",
        "description": "NCBI Gene to GO mappings; all species, filter to taxid 9606",
        "provides": ["entrez_id", "go_id"],
        "taxon_filter": True,
        "large": True,
    },
    "sifts_pdb": {
        "tier": "extended",
        "priority": 4,
        "url": "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz",
        "filename": "pdb_chain_uniprot.tsv.gz",
        "description": "SIFTS PDB chain to UniProt mappings",
        "provides": ["pdb_id", "uniprot_id"],
    },
}


def _normalize_id_type(id_type: str) -> str:
    return ID_TYPE_ALIASES.get(id_type, id_type)


def _validate_id_types(id_types: Iterable[str]) -> None:
    invalid = sorted(set(id_types) - set(ID_TYPES))
    if invalid:
        raise ValueError(
            f"Unsupported identifier namespace(s): {', '.join(invalid)}. "
            f"Supported namespaces: {', '.join(ID_TYPES)}"
        )


def sha256sum(path: Path) -> str:
    """Return the SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


class DatabaseDownloader:
    """Download raw source snapshots and record their provenance."""

    SOURCES = SOURCE_REGISTRY

    def __init__(
        self,
        data_dir: str = "data/raw",
        logger: Optional[logging.Logger] = None,
        ensembl_release: int = DEFAULT_ENSEMBL_RELEASE,
        timeout: int = 60,
        retries: int = 4,
    ):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logger or logging.getLogger(__name__)
        self.ensembl_release = ensembl_release
        self.timeout = timeout
        self.retries = retries

    def resolve_source(self, source: str) -> Dict:
        """Resolve release placeholders for a registry entry."""
        if source not in self.SOURCES:
            raise ValueError(
                f"Unknown source '{source}'. Available sources: {', '.join(self.SOURCES)}"
            )
        info = dict(self.SOURCES[source])
        info["url"] = info["url"].format(release=self.ensembl_release)
        info["filename"] = info["filename"].format(release=self.ensembl_release)
        return info

    def sources_for_tiers(self, tiers: Sequence[str]) -> List[str]:
        """Return sources in deterministic tier/priority order."""
        invalid = set(tiers) - set(TIERS)
        if invalid:
            raise ValueError(f"Unknown tiers: {sorted(invalid)}")
        return sorted(
            (name for name, info in self.SOURCES.items() if info["tier"] in tiers),
            key=lambda name: (
                TIERS.index(self.SOURCES[name]["tier"]),
                self.SOURCES[name].get("priority", 99),
                name,
            ),
        )

    def get_recommended_sources_for_conversion(
        self, input_type: str, output_types: List[str]
    ) -> List[str]:
        """Return default-tier sources that directly support a conversion."""
        input_id = _normalize_id_type(input_type)
        outputs = {_normalize_id_type(value) for value in output_types}
        _validate_id_types({input_id, *outputs})
        available = self.sources_for_tiers(("core", "crossref"))
        return [
            source
            for source in available
            if input_id in self.SOURCES[source]["provides"]
            and outputs.intersection(self.SOURCES[source]["provides"])
        ]

    def get_minimal_sources_for_conversion(
        self, input_type: str, output_types: List[str]
    ) -> List[str]:
        """Greedily select a small source set covering requested namespaces.

        Every selected source must contain the input namespace and at least one
        requested output namespace. This avoids downloading sources that cannot
        directly support the requested conversion.
        """
        input_id = _normalize_id_type(input_type)
        remaining = {_normalize_id_type(value) for value in output_types}
        _validate_id_types({input_id, *remaining})
        selected: List[str] = []
        candidates = self.sources_for_tiers(TIERS)

        while remaining:
            ranked = []
            for source in candidates:
                provides = set(self.SOURCES[source]["provides"])
                covered = remaining.intersection(provides)
                if input_id in provides and covered:
                    ranked.append(
                        (
                            len(covered),
                            -self.SOURCES[source].get("priority", 99),
                            source,
                            covered,
                        )
                    )
            if not ranked:
                break
            _, _, source, covered = max(ranked)
            selected.append(source)
            candidates.remove(source)
            remaining -= covered

        if remaining:
            self.logger.warning(
                "No direct source covers conversion from %s to: %s",
                input_id,
                ", ".join(sorted(remaining)),
            )
        return selected

    def coverage_table(self, sources: Optional[Iterable[str]] = None) -> Dict[str, List[str]]:
        """Return identifier namespace to source coverage."""
        selected = list(sources) if sources is not None else list(self.SOURCES)
        coverage: Dict[str, List[str]] = defaultdict(list)
        for source in selected:
            for id_type in self.SOURCES[source]["provides"]:
                coverage[id_type].append(source)
        return {id_type: coverage.get(id_type, []) for id_type in ID_TYPES}

    def print_coverage(self, sources: Optional[Iterable[str]] = None) -> None:
        coverage = self.coverage_table(sources)
        width = max(map(len, ID_TYPES))
        self.logger.info("Identifier coverage:")
        for id_type, source_names in coverage.items():
            self.logger.info(
                "  %-*s : %s",
                width,
                id_type,
                ", ".join(source_names) if source_names else "(not covered)",
            )

    def _validate_download(self, path: Path) -> None:
        """Perform inexpensive integrity checks without modifying raw data."""
        if not path.exists() or path.stat().st_size == 0:
            raise ValueError(f"Downloaded file is empty: {path}")
        if path.suffix == ".gz":
            with gzip.open(path, "rb") as handle:
                handle.read(1)

    def _download_file(self, url: str, destination: Path) -> None:
        """Download atomically with retry and backoff."""
        partial = destination.with_name(destination.name + ".part")
        for attempt in range(1, self.retries + 1):
            try:
                with requests.get(url, stream=True, timeout=self.timeout) as response:
                    response.raise_for_status()
                    total = int(response.headers.get("content-length", 0))
                    with partial.open("wb") as handle, tqdm(
                        total=total,
                        unit="iB",
                        unit_scale=True,
                        unit_divisor=1024,
                        desc=destination.name,
                        leave=False,
                    ) as progress:
                        for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                            if chunk:
                                handle.write(chunk)
                                progress.update(len(chunk))
                partial.replace(destination)
                self._validate_download(destination)
                return
            except Exception:
                if partial.exists():
                    partial.unlink()
                if attempt == self.retries:
                    raise
                self.logger.warning(
                    "Download attempt %d/%d failed for %s; retrying",
                    attempt,
                    self.retries,
                    destination.name,
                )
                time.sleep(3 * attempt)

    def _manifest_path(self) -> Path:
        return self.data_dir / "manifest.json"

    def _load_manifest(self) -> Dict:
        path = self._manifest_path()
        if path.exists():
            try:
                manifest = json.loads(path.read_text(encoding="utf-8"))
                existing_release = manifest.get("ensembl_release")
                if existing_release not in (None, self.ensembl_release):
                    raise ValueError(
                        f"Snapshot manifest pins Ensembl release {existing_release}, "
                        f"but release {self.ensembl_release} was requested. Use a "
                        "different snapshot directory."
                    )
                return manifest
            except (OSError, json.JSONDecodeError):
                self.logger.warning("Existing manifest is unreadable; rebuilding it")
        return {
            "schema_version": SNAPSHOT_SCHEMA_VERSION,
            "organism": "Homo sapiens",
            "taxid": 9606,
            "ensembl_release": self.ensembl_release,
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "databases": {},
        }

    def _write_manifest(self, manifest: Mapping) -> None:
        path = self._manifest_path()
        temporary = path.with_suffix(".json.part")
        temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)

    def download_source(self, source: str, force_download: bool = False) -> bool:
        """Download one source and update the snapshot manifest."""
        info = self.resolve_source(source)
        destination = self.data_dir / info["filename"]
        manifest = self._load_manifest()
        entry = {
            "url": info["url"],
            "filename": info["filename"],
            "tier": info["tier"],
            "description": info["description"],
            "provides": info["provides"],
            "taxon_filter_needed": info.get("taxon_filter", False),
            "large": info.get("large", False),
        }

        try:
            if destination.exists() and not force_download:
                self._validate_download(destination)
                self.logger.info("[skip] %s already exists", destination.name)
            else:
                self.logger.info("[get] %s", info["url"])
                self._download_file(info["url"], destination)
            entry.update(
                {
                    "status": "available",
                    "bytes": destination.stat().st_size,
                    "sha256": sha256sum(destination),
                    "verified_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                }
            )
            success = True
        except Exception as error:
            entry.update({"status": "failed", "error": str(error)})
            self.logger.error("Failed to download %s: %s", source, error)
            success = False

        manifest["ensembl_release"] = self.ensembl_release
        manifest["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
        manifest.setdefault("databases", {})[source] = entry
        self._write_manifest(manifest)
        return success

    def download_specific_sources(
        self, source_list: List[str], force_download: bool = False
    ) -> Dict[str, bool]:
        """Download an explicit source list."""
        ordered = list(dict.fromkeys(source_list))
        self.logger.info("Downloading %d source(s): %s", len(ordered), ", ".join(ordered))
        results: Dict[str, bool] = {}
        for index, source in enumerate(ordered, start=1):
            if source not in self.SOURCES:
                self.logger.error("Unknown source: %s", source)
                results[source] = False
                continue
            self.logger.info("[%d/%d] %s", index, len(ordered), source)
            results[source] = self.download_source(source, force_download)
        successful = sum(1 for value in results.values() if value)
        self.logger.info(
            "Download complete: %d/%d sources succeeded", successful, len(results)
        )
        return results

    def download_all(
        self, force_download: bool = False, include_optional: bool = True
    ) -> Dict[str, bool]:
        """Download all registry sources, or only core when optional is false."""
        tiers = TIERS if include_optional else ("core",)
        return self.download_specific_sources(
            self.sources_for_tiers(tiers), force_download=force_download
        )

    def download_default(self, force_download: bool = False) -> Dict[str, bool]:
        """Download the publication-oriented core and cross-reference tiers."""
        return self.download_specific_sources(
            self.sources_for_tiers(("core", "crossref")),
            force_download=force_download,
        )

    def validation_sources(self) -> List[str]:
        """Return independent sources used for benchmark validation."""
        return list(VALIDATION_SOURCES)

    def download_validation(
        self, force_download: bool = False, include_default: bool = True
    ) -> Dict[str, bool]:
        """Download independent benchmark-validation sources.

        By default this also ensures that the normal conversion snapshot is
        available. Validation files are not required for routine conversion.
        """
        sources = self.validation_sources()
        if include_default:
            sources = self.sources_for_tiers(("core", "crossref")) + sources
        return self.download_specific_sources(sources, force_download=force_download)

    def download_for_conversion(
        self,
        input_type: str,
        output_types: List[str],
        force_download: bool = False,
        minimal: bool = False,
    ) -> Dict[str, bool]:
        """Download only sources relevant to a requested conversion."""
        if minimal:
            sources = self.get_minimal_sources_for_conversion(input_type, output_types)
        else:
            sources = self.get_recommended_sources_for_conversion(input_type, output_types)
        self.logger.info("Selected sources: %s", ", ".join(sources) or "(none)")
        return self.download_specific_sources(sources, force_download)

    def get_file_info(self, source: str) -> Optional[Dict]:
        if source not in self.SOURCES:
            return None
        info = self.resolve_source(source)
        path = self.data_dir / info["filename"]
        if not path.exists():
            return None
        return {
            "filename": path.name,
            "size_mb": path.stat().st_size / (1024 * 1024),
            "description": info["description"],
            "sha256": sha256sum(path),
        }

    def list_downloaded_files(self) -> Dict[str, Dict]:
        return {
            source: info
            for source in self.SOURCES
            if (info := self.get_file_info(source)) is not None
        }
