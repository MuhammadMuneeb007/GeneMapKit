"""Tests for external MyGene.info concordance reporting."""

import json
import logging
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

from genemapkit.core.builder import build_database
from genemapkit.core.converter import GeneConverter
from genemapkit.core.downloader import sha256sum
from genemapkit.core.external_benchmark import (
    benchmark_external_services,
    query_ensembl_symbols,
)
from genemapkit.core.external_validation import compare_mygene


class ExternalValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def tearDown(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def create_database(self, root: Path) -> Path:
        raw = root / "raw"
        raw.mkdir()
        hgnc = raw / "hgnc.tsv"
        hgnc.write_text(
            "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tentrez_id\t"
            "ensembl_gene_id\trefseq_accession\tuniprot_ids\tccds_id\t"
            "omim_id\tucsc_id\tvega_id\tena\n"
            "HGNC:11998\tTP53\t\t\t7157\tENSG00000141510.18\t"
            "NM_000546.6\tP04637\t\t\t\t\t\n",
            encoding="utf-8",
        )
        (raw / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": "1.0.0",
                    "organism": "Homo sapiens",
                    "taxid": 9606,
                    "ensembl_release": 116,
                    "databases": {
                        "hgnc": {
                            "filename": hgnc.name,
                            "bytes": hgnc.stat().st_size,
                            "sha256": sha256sum(hgnc),
                        }
                    },
                }
            ),
            encoding="utf-8",
        )
        database = root / "build" / "genemapkit.db"
        build_database(raw, database)
        return database

    def test_external_comparison_reports_set_concordance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = self.create_database(root)

            def fetcher(symbols, outputs):
                self.assertEqual(["TP53"], symbols)
                return {
                    "TP53": {
                        "ensembl_gene_id": {"ENSG00000141510"},
                        "entrez_id": {"7157"},
                        "refseq_rna_id": {"NM_000546", "NM_001126112"},
                    }
                }

            report = compare_mygene(
                database,
                ["TP53"],
                ["ensembl_gene_id", "entrez_id", "refseq_rna_id"],
                output_dir=root / "results",
                fetcher=fetcher,
            )

            summaries = {row["output_type"]: row for row in report["summaries"]}
            self.assertEqual(1, summaries["ensembl_gene_id"]["exact_set_inputs"])
            self.assertEqual(1, summaries["entrez_id"]["exact_set_inputs"])
            self.assertEqual(1, summaries["refseq_rna_id"]["overlap_inputs"])
            self.assertEqual(50.0, summaries["refseq_rna_id"]["concordance_recall_pct"])
            self.assertLess(
                summaries["ensembl_gene_id"]["concordance_precision_ci_low_pct"], 100.0
            )
            self.assertTrue((root / "results" / "mygene_concordance_summary.md").exists())
            self.assertEqual(
                "external_concordance_not_ground_truth", report["comparison_type"]
            )

    def test_multi_service_benchmark_reports_consensus_and_reuses_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = self.create_database(root)
            calls = defaultdict(int)

            def fetcher(provider, values):
                def fetch(symbols, outputs):
                    calls[provider] += 1
                    return {
                        "TP53": {
                            output: set(values.get(output, set())) for output in outputs
                        }
                    }

                return fetch

            fetchers = {
                "mygene": fetcher(
                    "mygene",
                    {
                        "ensembl_gene_id": {"ENSG00000141510"},
                        "entrez_id": {"7157"},
                        "uniprot_id": {"P04637"},
                    },
                ),
                "ensembl": fetcher(
                    "ensembl", {"ensembl_gene_id": {"ENSG00000141510"}}
                ),
                "ncbi": fetcher("ncbi", {"entrez_id": {"7157"}}),
                "uniprot": fetcher(
                    "uniprot",
                    {
                        "ensembl_gene_id": {"ENSG00000141510"},
                        "uniprot_id": {"P04637"},
                    },
                ),
                "gprofiler": fetcher(
                    "gprofiler",
                    {
                        "ensembl_gene_id": {"ENSG00000141510"},
                        "entrez_id": {"7157"},
                    },
                ),
                "bridgedb": fetcher(
                    "bridgedb",
                    {
                        "ensembl_gene_id": {"ENSG00000141510"},
                        "entrez_id": {"7157"},
                        "uniprot_id": {"P04637"},
                    },
                ),
            }
            output_dir = root / "external"
            report = benchmark_external_services(
                database,
                ["TP53"],
                ["ensembl_gene_id", "entrez_id", "uniprot_id"],
                providers=list(fetchers),
                output_dir=output_dir,
                fetchers=fetchers,
            )
            support = {
                (row["output_type"], row["genemapkit_mapping"]): row
                for row in report["mapping_support"]
            }
            self.assertEqual(
                5,
                support[("ensembl_gene_id", "ENSG00000141510")][
                    "external_support_count"
                ],
            )
            self.assertEqual(
                "consensus",
                support[("entrez_id", "7157")]["support_status"],
            )
            self.assertTrue((output_dir / "external_benchmark_summary.md").exists())
            self.assertTrue((output_dir / "external_mapping_support.csv").exists())
            self.assertIn(
                "concordance_precision_ci_low_pct", report["provider_summaries"][0]
            )

            def fail_fetcher(symbols, outputs):
                raise AssertionError("cached provider data should be reused")

            benchmark_external_services(
                database,
                ["TP53"],
                ["ensembl_gene_id", "entrez_id", "uniprot_id"],
                providers=list(fetchers),
                output_dir=output_dir,
                fetchers={provider: fail_fetcher for provider in fetchers},
            )
            self.assertEqual(
                {
                    "mygene": 1,
                    "ensembl": 1,
                    "ncbi": 1,
                    "uniprot": 1,
                    "gprofiler": 1,
                    "bridgedb": 1,
                },
                dict(calls),
            )

    def test_external_benchmark_reuses_one_sqlite_connection_for_local_compare(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = self.create_database(root)

            def fetch(symbols, outputs):
                return {
                    symbol: {
                        output: (
                            {"ENSG00000141510"}
                            if symbol == "TP53" and output == "ensembl_gene_id"
                            else set()
                        )
                        for output in outputs
                    }
                    for symbol in symbols
                }

            open_calls = 0
            original = GeneConverter.open_connection

            def counted_open(connection_self):
                nonlocal open_calls
                open_calls += 1
                return original(connection_self)

            with patch.object(GeneConverter, "open_connection", counted_open):
                benchmark_external_services(
                    database,
                    ["TP53", "MISSING1", "MISSING2"],
                    ["ensembl_gene_id", "entrez_id"],
                    providers=["mygene"],
                    output_dir=root / "external-single-connection",
                    fetchers={"mygene": fetch},
                )

            self.assertEqual(1, open_calls)

    def test_ensembl_queries_symbols_in_batches_without_unneeded_expansion(self):
        class Response:
            status_code = 200
            headers = {}

            def __init__(self, records):
                self.records = records

            def raise_for_status(self):
                return None

            def json(self):
                return self.records

        class Session:
            def __init__(self):
                self.calls = []

            def post(self, url, *, params, json, headers, timeout):
                self.calls.append(
                    {
                        "url": url,
                        "params": params,
                        "symbols": list(json["symbols"]),
                        "timeout": timeout,
                    }
                )
                return Response(
                    {
                        symbol: {
                            "display_name": symbol,
                            "id": f"ENSG{index:011d}",
                        }
                        for index, symbol in enumerate(json["symbols"], 1)
                    }
                )

            def request(self, method, url, **kwargs):
                self.assert_post(method)
                return self.post(url, **kwargs)

            @staticmethod
            def assert_post(method):
                if method != "POST":
                    raise AssertionError(f"unexpected method: {method}")

        session = Session()
        symbols = ["GENE1", "GENE2", "GENE3"]
        results = query_ensembl_symbols(
            symbols,
            ["ensembl_gene_id"],
            batch_size=2,
            delay=0,
            session=session,
        )

        self.assertEqual(2, len(session.calls))
        self.assertEqual([["GENE1", "GENE2"], ["GENE3"]], [
            call["symbols"] for call in session.calls
        ])
        self.assertEqual({0}, {call["params"]["expand"] for call in session.calls})
        self.assertEqual({"ENSG00000000001"}, results["GENE1"]["ensembl_gene_id"])


if __name__ == "__main__":
    unittest.main()
