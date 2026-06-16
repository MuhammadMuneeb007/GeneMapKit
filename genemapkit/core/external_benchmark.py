"""Cached multi-service external concordance benchmark."""

from __future__ import annotations

import csv
import json
import logging
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set
from urllib.parse import quote

import requests

from .builder import normalize_identifier
from .converter import GeneConverter, canonical_id_type, canonical_mapping_policy
from .external_validation import query_mygene_symbols
from .statistics import interval_fields
from genemapkit.utils.progress import track

LOGGER = logging.getLogger("genemapkit")


ProviderResults = Dict[str, Dict[str, Set[str]]]
ProviderFetcher = Callable[[Sequence[str], Sequence[str]], ProviderResults]

PROVIDER_OUTPUTS = {
    "mygene": {
        "ensembl_gene_id",
        "entrez_id",
        "refseq_rna_id",
        "refseq_protein_id",
        "uniprot_id",
    },
    "ensembl": {
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "ensembl_protein_id",
    },
    "ncbi": {"entrez_id"},
    "uniprot": {
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "ensembl_protein_id",
        "refseq_rna_id",
        "refseq_protein_id",
        "uniprot_id",
    },
    "gprofiler": {
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "ensembl_protein_id",
        "entrez_id",
        "refseq_rna_id",
        "refseq_protein_id",
        "uniprot_id",
    },
    "bridgedb": {
        "ensembl_gene_id",
        "ensembl_transcript_id",
        "ensembl_protein_id",
        "entrez_id",
        "refseq_rna_id",
        "refseq_protein_id",
        "uniprot_id",
    },
}
DEFAULT_OUTPUTS = (
    "ensembl_gene_id",
    "entrez_id",
    "refseq_rna_id",
    "refseq_protein_id",
    "uniprot_id",
)
DEFAULT_PROVIDERS = tuple(PROVIDER_OUTPUTS)
PROVIDER_ENDPOINTS = {
    "mygene": "https://mygene.info/v3/query",
    "ensembl": "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
    "ncbi": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "uniprot": "https://rest.uniprot.org/uniprotkb/search",
    "gprofiler": "https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
    "bridgedb": "https://webservice.bridgedb.org/Human/xrefsBatch",
}

GPROFILER_TARGETS = {
    "ensembl_gene_id": ("ENSG",),
    "ensembl_transcript_id": ("ENST",),
    "ensembl_protein_id": ("ENSP",),
    "entrez_id": ("ENTREZGENE_ACC",),
    "refseq_rna_id": ("REFSEQ_MRNA", "REFSEQ_NCRNA"),
    "refseq_protein_id": ("REFSEQ_PEPTIDE",),
    "uniprot_id": ("UNIPROTSWISSPROT", "UNIPROTSPTREMBL"),
}


def _empty_results(
    symbols: Sequence[str], outputs: Sequence[str]
) -> ProviderResults:
    return {symbol: {output: set() for output in outputs} for symbol in symbols}


def _add_normalized(
    target: Set[str], output_type: str, value: object
) -> None:
    normalized = normalize_identifier(output_type, value)
    if normalized:
        target.add(normalized)


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
            "%s returned HTTP %s; retrying in %.1f seconds",
            url,
            response.status_code,
            wait,
        )
        time.sleep(wait)
    return response


def query_ensembl_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 1000,
    timeout: float = 120.0,
    delay: float = 0.1,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> ProviderResults:
    """Query Ensembl REST in batches, expanding transcripts only when requested."""
    if batch_size < 1:
        raise ValueError("Ensembl batch_size must be a positive integer")
    outputs = [canonical_id_type(value) for value in output_types]
    results = _empty_results(symbols, outputs)
    client = session or requests.Session()
    expand = bool({"ensembl_transcript_id", "ensembl_protein_id"}.intersection(outputs))
    starts = range(0, len(symbols), batch_size)
    for start in track(
        starts,
        enabled=progress,
        description="Ensembl API",
        total=len(starts),
        unit="batches",
    ):
        batch = list(symbols[start : start + batch_size])
        response = _request_with_retries(
            client,
            "POST",
            PROVIDER_ENDPOINTS["ensembl"],
            max_retries=max_retries,
            params={"expand": 1 if expand else 0},
            json={"symbols": batch},
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": "GeneMapKit external concordance benchmark",
            },
            timeout=timeout,
        )
        response.raise_for_status()
        records = response.json()
        for symbol in batch:
            record = records.get(symbol)
            if not isinstance(record, Mapping):
                continue
            if str(record.get("display_name", "")).casefold() != symbol.casefold():
                continue
            if "ensembl_gene_id" in outputs:
                _add_normalized(
                    results[symbol]["ensembl_gene_id"],
                    "ensembl_gene_id",
                    record.get("id"),
                )
            for transcript in record.get("Transcript", []) or []:
                if "ensembl_transcript_id" in outputs:
                    _add_normalized(
                        results[symbol]["ensembl_transcript_id"],
                        "ensembl_transcript_id",
                        transcript.get("id"),
                    )
                translation = transcript.get("Translation") or {}
                if "ensembl_protein_id" in outputs:
                    _add_normalized(
                        results[symbol]["ensembl_protein_id"],
                        "ensembl_protein_id",
                        translation.get("id"),
                    )
        if delay:
            time.sleep(delay)
    return results


def query_ncbi_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 200,
    timeout: float = 60.0,
    delay: float = 0.35,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> ProviderResults:
    """Query NCBI E-utilities in attributed ESearch/ESummary batches."""
    if batch_size < 1:
        raise ValueError("NCBI batch_size must be a positive integer")
    outputs = [canonical_id_type(value) for value in output_types]
    results = _empty_results(symbols, outputs)
    client = session or requests.Session()
    starts = range(0, len(symbols), batch_size)
    for start in track(starts, enabled=progress, description="NCBI API", total=len(starts), unit="batches"):
        batch = list(symbols[start : start + batch_size])
        requested = {symbol.casefold(): symbol for symbol in batch}
        term = "(" + " OR ".join(f'"{symbol}"[Gene Name]' for symbol in batch) + ") AND 9606[Taxonomy ID]"
        response = _request_with_retries(
            client,
            "POST",
            PROVIDER_ENDPOINTS["ncbi"],
            max_retries=max_retries,
            data={
                "db": "gene",
                "term": term,
                "retmode": "json",
                "retmax": batch_size * 10,
                "tool": "GeneMapKit",
            },
            headers={"User-Agent": "GeneMapKit external concordance benchmark"},
            timeout=timeout,
        )
        response.raise_for_status()
        ids = response.json().get("esearchresult", {}).get("idlist", [])
        for id_start in range(0, len(ids), 500):
            id_batch = ids[id_start : id_start + 500]
            summary = _request_with_retries(
                client,
                "POST",
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                max_retries=max_retries,
                data={"db": "gene", "id": ",".join(id_batch), "retmode": "json", "tool": "GeneMapKit"},
                headers={"User-Agent": "GeneMapKit external concordance benchmark"},
                timeout=timeout,
            )
            summary.raise_for_status()
            records = summary.json().get("result", {})
            for value in id_batch:
                record = records.get(value, {})
                symbol = requested.get(str(record.get("name", "")).casefold())
                if symbol and "entrez_id" in outputs:
                    _add_normalized(results[symbol]["entrez_id"], "entrez_id", value)
        if delay:
            time.sleep(delay)
    return results


def query_uniprot_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 50,
    timeout: float = 120.0,
    delay: float = 0.1,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> ProviderResults:
    """Query UniProtKB using batched exact-gene OR searches."""
    if batch_size < 1:
        raise ValueError("UniProt batch_size must be a positive integer")
    outputs = [canonical_id_type(value) for value in output_types]
    results = _empty_results(symbols, outputs)
    client = session or requests.Session()
    starts = range(0, len(symbols), batch_size)
    for start in track(starts, enabled=progress, description="UniProt API", total=len(starts), unit="batches"):
        batch = list(symbols[start : start + batch_size])
        requested = {symbol.casefold(): symbol for symbol in batch}
        query = "(" + " OR ".join(f"gene_exact:{symbol}" for symbol in batch) + ") AND organism_id:9606"
        next_url = PROVIDER_ENDPOINTS["uniprot"]
        params = {
                "query": query,
                "format": "json",
                "size": 500,
                "fields": "accession,gene_names,xref_ensembl,xref_refseq",
            }
        while next_url:
            response = _request_with_retries(
                client, "GET", next_url, max_retries=max_retries, params=params,
                headers={"User-Agent": "GeneMapKit external concordance benchmark"}, timeout=timeout,
            )
            response.raise_for_status()
            for record in response.json().get("results", []):
                gene_names = {
                    str(gene.get("geneName", {}).get("value", "")).casefold()
                    for gene in record.get("genes", [])
                }
                matched_symbols = [requested[name] for name in gene_names if name in requested]
                for symbol in matched_symbols:
                    _add_uniprot_record(results[symbol], outputs, record)
            next_url = response.links.get("next", {}).get("url")
            params = None
        if delay:
            time.sleep(delay)
    return results


def _add_uniprot_record(
    symbol_results: Dict[str, Set[str]], outputs: Sequence[str], record: Mapping[str, object]
) -> None:
    if "uniprot_id" in outputs:
        _add_normalized(symbol_results["uniprot_id"], "uniprot_id", record.get("primaryAccession"))
    for xref in record.get("uniProtKBCrossReferences", []):
        values = [xref.get("id")]
        values.extend(
            item.get("value")
            for item in xref.get("properties", [])
            if isinstance(item, Mapping)
        )
        if xref.get("database") == "Ensembl":
            for value in values:
                if str(value).startswith("ENSG") and "ensembl_gene_id" in outputs:
                    _add_normalized(symbol_results["ensembl_gene_id"], "ensembl_gene_id", value)
                elif str(value).startswith("ENST") and "ensembl_transcript_id" in outputs:
                    _add_normalized(symbol_results["ensembl_transcript_id"], "ensembl_transcript_id", value)
                elif str(value).startswith("ENSP") and "ensembl_protein_id" in outputs:
                    _add_normalized(symbol_results["ensembl_protein_id"], "ensembl_protein_id", value)
        elif xref.get("database") == "RefSeq":
            for value in values:
                prefix = str(value).split("_", 1)[0]
                if prefix in {"NM", "NR", "XM", "XR"} and "refseq_rna_id" in outputs:
                    _add_normalized(symbol_results["refseq_rna_id"], "refseq_rna_id", value)
                elif prefix in {"NP", "XP", "YP", "WP"} and "refseq_protein_id" in outputs:
                    _add_normalized(symbol_results["refseq_protein_id"], "refseq_protein_id", value)


def query_gprofiler_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 1000,
    timeout: float = 120.0,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> ProviderResults:
    """Query g:Profiler g:Convert for exact human symbols."""
    outputs = [canonical_id_type(value) for value in output_types]
    results = _empty_results(symbols, outputs)
    client = session or requests.Session()
    if batch_size < 1:
        raise ValueError("g:Profiler batch_size must be a positive integer")
    requests_to_make = [
        (output, target, start)
        for output in outputs
        for target in GPROFILER_TARGETS[output]
        for start in range(0, len(symbols), batch_size)
    ]
    for output, target, start in track(
        requests_to_make,
        enabled=progress,
        description="g:Profiler API",
        total=len(requests_to_make),
        unit="queries",
    ):
        batch = list(symbols[start : start + batch_size])
        response = _request_with_retries(
            client, "POST",
            PROVIDER_ENDPOINTS["gprofiler"],
            max_retries=max_retries,
            json={"organism": "hsapiens", "query": batch, "target": target},
            headers={"User-Agent": "GeneMapKit external concordance benchmark"},
            timeout=timeout,
        )
        response.raise_for_status()
        for record in response.json().get("result", []):
            symbol = str(record.get("incoming", "")).strip()
            value = str(record.get("converted", "")).strip()
            if symbol not in results or value in {"", "None"}:
                continue
            if output == "uniprot_id":
                value = value.split(".", 1)[0]
            _add_normalized(results[symbol][output], output, value)
    return results


def query_bridgedb_symbols(
    symbols: Sequence[str],
    output_types: Sequence[str],
    *,
    batch_size: int = 50,
    timeout: float = 120.0,
    delay: float = 0.1,
    max_retries: int = 5,
    session: Optional[requests.Session] = None,
    progress: bool = False,
) -> ProviderResults:
    """Query BridgeDb Human xrefs using the batch endpoint."""
    if batch_size < 1:
        raise ValueError("BridgeDb batch_size must be a positive integer")
    outputs = [canonical_id_type(value) for value in output_types]
    results = _empty_results(symbols, outputs)
    client = session or requests.Session()
    starts = range(0, len(symbols), batch_size)
    for start in track(starts, enabled=progress, description="BridgeDb API", total=len(starts), unit="batches"):
        batch = list(symbols[start : start + batch_size])
        response = _request_with_retries(
            client, "POST",
            PROVIDER_ENDPOINTS["bridgedb"],
            max_retries=max_retries,
            data="".join(f"{symbol}\tHGNC\n" for symbol in batch),
            headers={
                "Accept": "text/plain",
                "Content-Type": "text/plain",
                "User-Agent": "GeneMapKit external concordance benchmark",
            },
            timeout=timeout,
        )
        response.raise_for_status()
        for line in response.text.splitlines():
            fields = line.split("\t", 2)
            if len(fields) != 3 or fields[0] not in results:
                continue
            symbol = fields[0]
            for item in fields[2].split(","):
                code, separator, value = item.partition(":")
                if not separator:
                    continue
                output = None
                if value.startswith("ENSG"):
                    output = "ensembl_gene_id"
                elif value.startswith("ENST"):
                    output = "ensembl_transcript_id"
                elif value.startswith("ENSP"):
                    output = "ensembl_protein_id"
                elif code == "Wg" and value.isdigit():
                    output = "entrez_id"
                elif code == "Q":
                    prefix = value.split("_", 1)[0]
                    output = (
                        "refseq_rna_id"
                        if prefix in {"NM", "NR", "XM", "XR"}
                        else "refseq_protein_id"
                        if prefix in {"NP", "XP", "YP", "WP"}
                        else None
                    )
                elif code in {"S", "Sp"}:
                    output = "uniprot_id"
                if output in outputs:
                    _add_normalized(results[symbol][output], output, value)
        if delay:
            time.sleep(delay)
    return results


PROVIDER_FETCHERS: Dict[str, ProviderFetcher] = {
    "mygene": query_mygene_symbols,
    "ensembl": query_ensembl_symbols,
    "ncbi": query_ncbi_symbols,
    "uniprot": query_uniprot_symbols,
    "gprofiler": query_gprofiler_symbols,
    "bridgedb": query_bridgedb_symbols,
}


def _load_cache(path: Path) -> ProviderResults:
    if not path.is_file():
        return {}
    raw = json.loads(path.read_text(encoding="utf-8"))
    return {
        symbol: {output: set(values) for output, values in outputs.items()}
        for symbol, outputs in raw.items()
    }


def _write_cache(path: Path, results: ProviderResults) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    serializable = {
        symbol: {output: sorted(values) for output, values in outputs.items()}
        for symbol, outputs in results.items()
    }
    path.write_text(json.dumps(serializable, indent=2) + "\n", encoding="utf-8")


def _fetch_with_cache(
    provider: str,
    symbols: Sequence[str],
    outputs: Sequence[str],
    cache_dir: Path,
    fetcher: ProviderFetcher,
    *,
    refresh_cache: bool,
    progress: bool = False,
) -> tuple[ProviderResults, Dict[str, object]]:
    cache_path = cache_dir / f"{provider}.json"
    cached = {} if refresh_cache else _load_cache(cache_path)
    missing = [symbol for symbol in symbols if symbol not in cached]
    started = time.perf_counter()
    if missing:
        if fetcher in PROVIDER_FETCHERS.values():
            fetched = fetcher(missing, outputs, progress=progress)
        else:
            fetched = fetcher(missing, outputs)
        for symbol in missing:
            cached[symbol] = fetched.get(
                symbol, {output: set() for output in outputs}
            )
        _write_cache(cache_path, cached)
    results = {
        symbol: {
            output: set(cached.get(symbol, {}).get(output, set()))
            for output in outputs
        }
        for symbol in symbols
    }
    return results, {
        "provider": provider,
        "requested_symbols": len(symbols),
        "queried_symbols": len(missing),
        "cached_symbols": len(symbols) - len(missing),
        "elapsed_seconds": round(time.perf_counter() - started, 6),
        "cache_used": bool(symbols) and not missing,
    }


def _set_status(gmk: Set[str], external: Set[str]) -> str:
    if not gmk and not external:
        return "both_unmapped"
    if gmk == external:
        return "exact_external_union"
    if gmk and external and gmk.issubset(external):
        return "supported_with_external_extras"
    if gmk and external and gmk.intersection(external):
        return "partial_support"
    if gmk and not external:
        return "genemapkit_only"
    if external and not gmk:
        return "external_only"
    return "conflict"


def benchmark_external_services(
    database: Path,
    symbols: Iterable[str],
    output_types: Sequence[str] = DEFAULT_OUTPUTS,
    *,
    providers: Sequence[str] = DEFAULT_PROVIDERS,
    mapping_policy: str = "all-supported",
    output_dir: Path = Path("results/external_benchmark"),
    refresh_cache: bool = False,
    fetchers: Optional[Mapping[str, ProviderFetcher]] = None,
    progress: bool = False,
) -> Dict[str, object]:
    """Compare final GeneMapKit mappings with multiple live external services."""
    unique_symbols = list(
        dict.fromkeys(str(value).strip() for value in symbols if str(value).strip())
    )
    if not unique_symbols:
        raise ValueError("No non-empty symbols were provided")
    outputs = list(dict.fromkeys(canonical_id_type(value) for value in output_types))
    selected_providers = list(dict.fromkeys(str(value).strip().lower() for value in providers))
    unknown = sorted(set(selected_providers) - set(PROVIDER_OUTPUTS))
    if unknown:
        raise ValueError(f"Unknown external provider(s): {', '.join(unknown)}")
    policy = canonical_mapping_policy(mapping_policy)
    converter = GeneConverter(database, mapping_policy=policy)
    output_dir = Path(output_dir)
    fetcher_map = dict(PROVIDER_FETCHERS)
    fetcher_map.update(fetchers or {})

    provider_results: Dict[str, ProviderResults] = {}
    provider_errors: Dict[str, str] = {}
    provider_timings: List[Dict[str, object]] = []
    for provider in track(
        selected_providers,
        enabled=progress,
        description="Querying external providers",
        total=len(selected_providers),
        unit="providers",
    ):
        supported_outputs = [value for value in outputs if value in PROVIDER_OUTPUTS[provider]]
        if not supported_outputs:
            continue
        LOGGER.info("External benchmark: querying %s for %d genes", provider, len(unique_symbols))
        try:
            provider_results[provider], timing = _fetch_with_cache(
                provider,
                unique_symbols,
                supported_outputs,
                output_dir / "cache",
                fetcher_map[provider],
                refresh_cache=refresh_cache,
                progress=progress,
            )
            provider_timings.append(timing)
        except (requests.RequestException, ValueError, KeyError, json.JSONDecodeError) as error:
            provider_errors[provider] = str(error)
            LOGGER.warning("External provider %s failed: %s", provider, error)

    details: List[Dict[str, object]] = []
    mapping_support: List[Dict[str, object]] = []
    summary_counts: Dict[str, Dict[str, int]] = {
        output: defaultdict(int) for output in outputs
    }
    provider_summary: Dict[tuple[str, str], Dict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )

    with converter.open_connection() as connection:
        for symbol in track(
            unique_symbols,
            enabled=progress,
            description="Comparing external mappings",
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
                available = [
                    provider
                    for provider in selected_providers
                    if provider in provider_results and output in PROVIDER_OUTPUTS[provider]
                ]
                external_sets = {
                    provider: set(provider_results[provider][symbol].get(output, set()))
                    for provider in available
                }
                external_union = (
                    set().union(*external_sets.values()) if external_sets else set()
                )
                common = gmk.intersection(external_union)
                status = (
                    _set_status(gmk, external_union)
                    if available
                    else "no_external_service"
                )
                counts = summary_counts[output]
                counts[status] += 1
                if available:
                    counts["evaluated_symbols"] += 1
                    counts["genemapkit_mappings"] += len(gmk)
                    counts["external_union_mappings"] += len(external_union)
                    counts["common_mappings"] += len(common)
                for provider, values in external_sets.items():
                    provider_common = gmk.intersection(values)
                    provider_counts = provider_summary[(provider, output)]
                    provider_counts["tested_symbols"] += 1
                    provider_counts["genemapkit_mappings"] += len(gmk)
                    provider_counts["provider_mappings"] += len(values)
                    provider_counts["common_mappings"] += len(provider_common)
                for mapping in sorted(gmk):
                    sources = sorted(
                        provider
                        for provider, values in external_sets.items()
                        if mapping in values
                    )
                    mapping_support.append(
                        {
                            "symbol": symbol,
                            "output_type": output,
                            "genemapkit_mapping": mapping,
                            "external_support_count": len(sources),
                            "supporting_external_services": ";".join(sources),
                            "support_status": (
                                "consensus"
                                if len(sources) >= 2
                                else "single_external_support"
                                if sources
                                else "no_external_service"
                                if not available
                                else "unsupported_externally"
                            ),
                        }
                    )
                detail = {
                    "symbol": symbol,
                    "output_type": output,
                    "mapping_policy": policy,
                    "status": status,
                    "available_external_services": ";".join(available),
                    "genemapkit": ";".join(sorted(gmk)),
                    "external_union": ";".join(sorted(external_union)),
                    "common": ";".join(sorted(common)),
                    "genemapkit_only": ";".join(sorted(gmk - external_union)),
                    "external_only": ";".join(sorted(external_union - gmk)),
                }
                for provider in selected_providers:
                    detail[provider] = (
                        ";".join(sorted(external_sets[provider]))
                        if provider in external_sets
                        else "NA"
                    )
                details.append(detail)

    summaries = []
    for output in outputs:
        counts = summary_counts[output]
        gmk_total = counts["genemapkit_mappings"]
        external_total = counts["external_union_mappings"]
        common = counts["common_mappings"]
        row = {
                "output_type": output,
                "requested_symbols": len(unique_symbols),
                "tested_symbols": counts["evaluated_symbols"],
                "exact_external_union": counts["exact_external_union"],
                "supported_with_external_extras": counts["supported_with_external_extras"],
                "partial_support": counts["partial_support"],
                "conflict": counts["conflict"],
                "genemapkit_only": counts["genemapkit_only"],
                "external_only": counts["external_only"],
                "both_unmapped": counts["both_unmapped"],
                "no_external_service": counts["no_external_service"],
                "genemapkit_mappings": gmk_total,
                "external_union_mappings": external_total,
                "common_mappings": common,
                "concordance_precision_pct": (
                    round(100.0 * common / gmk_total, 4) if gmk_total else None
                ),
                "concordance_recall_pct": (
                    round(100.0 * common / external_total, 4) if external_total else None
                ),
            }
        row.update(interval_fields("concordance_precision", common, gmk_total))
        row.update(interval_fields("concordance_recall", common, external_total))
        summaries.append(row)

    provider_summaries = []
    for (provider, output), counts in sorted(provider_summary.items()):
        common = counts["common_mappings"]
        gmk_total = counts["genemapkit_mappings"]
        provider_total = counts["provider_mappings"]
        row = {
                "provider": provider,
                "output_type": output,
                "tested_symbols": counts["tested_symbols"],
                "genemapkit_mappings": gmk_total,
                "provider_mappings": provider_total,
                "common_mappings": common,
                "concordance_precision_pct": (
                    round(100.0 * common / gmk_total, 4) if gmk_total else None
                ),
                "concordance_recall_pct": (
                    round(100.0 * common / provider_total, 4) if provider_total else None
                ),
            }
        row.update(interval_fields("concordance_precision", common, gmk_total))
        row.update(interval_fields("concordance_recall", common, provider_total))
        provider_summaries.append(row)

    report: Dict[str, object] = {
        "comparison_type": "multi_service_external_concordance_not_ground_truth",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "database": str(database),
        "database_metadata": converter.metadata(),
        "mapping_policy": policy,
        "input_type": "symbol",
        "output_types": outputs,
        "providers_requested": selected_providers,
        "providers_completed": sorted(provider_results),
        "provider_errors": provider_errors,
        "provider_timings": provider_timings,
        "provider_endpoints": {
            provider: PROVIDER_ENDPOINTS[provider] for provider in selected_providers
        },
        "tested_symbols": len(unique_symbols),
        "summaries": summaries,
        "provider_summaries": provider_summaries,
        "details": details,
        "mapping_support": mapping_support,
        "interpretation": (
            "External services integrate overlapping upstream databases. These metrics "
            "measure concordance and external support, not independent accuracy."
        ),
    }
    write_external_benchmark_report(output_dir, report)
    LOGGER.info("External benchmark complete: %d providers completed", len(provider_results))
    return report


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_external_benchmark_report(
    output_dir: Path, report: Mapping[str, object]
) -> None:
    """Write complete machine-readable and publication-friendly tables."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "external_benchmark.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    _write_csv(output_dir / "external_benchmark_summary.csv", report["summaries"])
    _write_csv(output_dir / "external_provider_summary.csv", report["provider_summaries"])
    _write_csv(output_dir / "external_provider_timings.csv", report["provider_timings"])
    _write_csv(output_dir / "external_benchmark_details.csv", report["details"])
    _write_csv(output_dir / "external_mapping_support.csv", report["mapping_support"])
    lines = [
        "# Multi-service External Concordance",
        "",
        "> External services overlap with GeneMapKit upstream sources. Report these "
        "results as concordance and external support, not independent accuracy.",
        "",
        f"- Mapping policy: `{report['mapping_policy']}`",
        f"- Tested symbols: `{report['tested_symbols']}`",
        f"- Completed providers: `{', '.join(report['providers_completed'])}`",
        "- Intervals: Wilson 95% confidence intervals.",
        "",
        "| Output type | Exact union | Supported + extras | Partial | Conflict | GMK-only | External-only | Precision (95% CI) | Recall (95% CI) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in report["summaries"]:
        lines.append(
            f"| `{row['output_type']}` | {row['exact_external_union']} | "
            f"{row['supported_with_external_extras']} | {row['partial_support']} | "
            f"{row['conflict']} | {row['genemapkit_only']} | {row['external_only']} | "
            f"{row['concordance_precision_pct']}% "
            f"[{row['concordance_precision_ci_low_pct']}, {row['concordance_precision_ci_high_pct']}] | "
            f"{row['concordance_recall_pct']}% "
            f"[{row['concordance_recall_ci_low_pct']}, {row['concordance_recall_ci_high_pct']}] |"
        )
    lines.extend(
        [
            "",
            "## Per-provider concordance",
            "",
            "| Provider | Output type | GeneMapKit mappings | Provider mappings | Common | Precision (95% CI) | Recall (95% CI) |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )
    for row in report["provider_summaries"]:
        lines.append(
            f"| `{row['provider']}` | `{row['output_type']}` | "
            f"{row['genemapkit_mappings']} | {row['provider_mappings']} | "
            f"{row['common_mappings']} | {row['concordance_precision_pct']}% "
            f"[{row['concordance_precision_ci_low_pct']}, {row['concordance_precision_ci_high_pct']}] | "
            f"{row['concordance_recall_pct']}% "
            f"[{row['concordance_recall_ci_low_pct']}, {row['concordance_recall_ci_high_pct']}] |"
        )
    lines.extend(
        [
            "",
            "## Provider query timings",
            "",
            "| Provider | Requested symbols | Queried symbols | Cached symbols | Seconds |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in report["provider_timings"]:
        lines.append(
            f"| `{row['provider']}` | {row['requested_symbols']} | "
            f"{row['queried_symbols']} | {row['cached_symbols']} | "
            f"{row['elapsed_seconds']} |"
        )
    if report["provider_errors"]:
        lines.extend(["", "## Provider errors", ""])
        for provider, error in report["provider_errors"].items():
            lines.append(f"- `{provider}`: {error}")
    (output_dir / "external_benchmark_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
