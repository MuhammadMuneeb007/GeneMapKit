"""Tests for all-to-all conversion and path consistency benchmarking."""

import json
import logging
import tempfile
import unittest
from pathlib import Path

from genemapkit.core.builder import build_database
from genemapkit.core.downloader import sha256sum
from genemapkit.core.path_benchmark import benchmark_conversion_paths


class PathBenchmarkTests(unittest.TestCase):
    def setUp(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def tearDown(self) -> None:
        logging.getLogger("genemapkit").handlers.clear()

    def test_reports_coverage_round_trip_and_pivot_consistency(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw"
            raw.mkdir()
            hgnc = raw / "hgnc.tsv"
            hgnc.write_text(
                "hgnc_id\tsymbol\tentrez_id\tensembl_gene_id\tuniprot_ids\n"
                "HGNC:11998\tTP53\t7157\tENSG00000141510\tP04637\n",
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
            output = root / "results"
            report = benchmark_conversion_paths(
                database,
                id_types=["symbol", "entrez_id", "ensembl_gene_id", "uniprot_id"],
                pivot_type="entrez_id",
                path_sample=None,
                output_dir=output,
            )
            coverage = {
                (row["source_type"], row["target_type"]): row
                for row in report["conversion_coverage"]
            }
            self.assertEqual(
                100.0,
                coverage[("symbol", "ensembl_gene_id")]["conversion_success_pct"],
            )
            roundtrips = {
                (row["source_type"], row["intermediate_type"]): row
                for row in report["round_trip_consistency"]
            }
            self.assertEqual(
                100.0,
                roundtrips[("symbol", "entrez_id")]["exact_round_trip_pct"],
            )
            self.assertIn(
                "round_trip_recovery_ci_low_pct",
                roundtrips[("symbol", "entrez_id")],
            )
            paths = {
                (row["source_type"], row["target_type"]): row
                for row in report["pivot_path_consistency"]
            }
            self.assertEqual(
                100.0,
                paths[("symbol", "ensembl_gene_id")]["exact_nonempty_set_pct"],
            )
            self.assertTrue((output / "path_benchmark_summary.md").exists())


if __name__ == "__main__":
    unittest.main()
