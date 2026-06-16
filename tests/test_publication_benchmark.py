"""Tests for the combined publication benchmark."""

import json
import logging
import tempfile
import unittest
from pathlib import Path

from genemapkit.core.builder import build_database
from genemapkit.core.downloader import sha256sum
from genemapkit.core.publication_benchmark import (
    aggregate_held_out_reports,
    benchmark_local_performance,
    load_gtf_gene_coordinates,
    query_ensembl_coordinates,
)


class PublicationBenchmarkTests(unittest.TestCase):
    def setUp(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def tearDown(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def create_database(self, root: Path) -> Path:
        raw = root / "raw"
        raw.mkdir()
        hgnc = raw / "hgnc.tsv"
        hgnc.write_text(
            "hgnc_id\tsymbol\tentrez_id\tensembl_gene_id\trefseq_accession\t"
            "uniprot_ids\n"
            "HGNC:11998\tTP53\t7157\tENSG00000141510.18\tNM_000546.6\tP04637\n",
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

    def test_local_performance_reports_scaling_memory_and_offline_operation(self):
        with tempfile.TemporaryDirectory() as directory:
            database = self.create_database(Path(directory))
            report = benchmark_local_performance(
                database,
                ["TP53"],
                ["ensembl_gene_id", "entrez_id"],
                batch_sizes=[1, 3],
                repeats=1,
            )
            self.assertTrue(report["offline_operation_passed"])
            self.assertEqual([1, 3], [row["batch_size"] for row in report["scaling"]])
            self.assertGreater(report["scaling"][1]["identifiers_per_second"], 0)
            self.assertGreater(report["scaling"][1]["median_peak_process_rss_mb"], 0)

    def test_coordinate_parser_and_held_out_aggregation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / "genes.gtf"
            gtf.write_text(
                '1\tensembl\tgene\t100\t200\t.\t+\t.\t'
                'gene_id "ENSG00000141510.18"; gene_name "TP53";\n',
                encoding="utf-8",
            )
            coordinates = load_gtf_gene_coordinates(gtf, ["TP53", "MISSING"])
            self.assertEqual("ENSG00000141510", coordinates["TP53"][0]["ensembl_gene_id"])
            self.assertEqual([], coordinates["MISSING"])

            held_out = root / "held_out"
            held_out.mkdir()
            (held_out / "validation_metrics.json").write_text(
                json.dumps(
                    {
                        "database": "held_out.db",
                        "benchmarks": [
                            {"truth_source": "mane_summary", "independent": True},
                            {"truth_source": "mane_summary", "independent": False},
                        ],
                    }
                ),
                encoding="utf-8",
            )
            report = aggregate_held_out_reports([held_out, root / "missing"])
            self.assertEqual(1, len(report["independent_benchmarks"]))
            self.assertEqual("missing", report["sources"][1]["status"])

    def test_ensembl_coordinate_lookup_uses_batched_post_requests(self):
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

            def request(self, method, url, **kwargs):
                self.calls.append(
                    {
                        "method": method,
                        "url": url,
                        "params": kwargs["params"],
                        "symbols": list(kwargs["json"]["symbols"]),
                    }
                )
                records = {}
                for index, symbol in enumerate(kwargs["json"]["symbols"], 1):
                    records[symbol] = {
                        "id": f"ENSG{index:011d}",
                        "assembly_name": "GRCh38",
                        "seq_region_name": "17",
                        "start": index * 100,
                        "end": index * 100 + 50,
                        "strand": 1,
                    }
                return Response(records)

        with tempfile.TemporaryDirectory() as directory:
            session = Session()
            coordinates, timing = query_ensembl_coordinates(
                ["GENE1", "GENE2", "GENE3"],
                cache_path=Path(directory) / "coordinates.json",
                batch_size=2,
                delay=0,
                session=session,
            )

        self.assertEqual(2, len(session.calls))
        self.assertEqual([["GENE1", "GENE2"], ["GENE3"]], [
            call["symbols"] for call in session.calls
        ])
        self.assertEqual({"POST"}, {call["method"] for call in session.calls})
        self.assertEqual({0}, {call["params"]["expand"] for call in session.calls})
        self.assertEqual("ENSG00000000001", coordinates["GENE1"][0]["ensembl_gene_id"])
        self.assertEqual(3, timing["queried_symbols"])
        self.assertEqual(0, timing["cached_symbols"])


if __name__ == "__main__":
    unittest.main()
