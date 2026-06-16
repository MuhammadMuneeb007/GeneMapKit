"""Tests for relationship-level validation benchmarks."""

import json
import logging
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

from click.testing import CliRunner

from genemapkit.cli import cli
from genemapkit.core.builder import build_database
from genemapkit.core.downloader import sha256sum
from genemapkit.core.validation import validate_database


class ValidationTests(unittest.TestCase):
    def setUp(self):
        logging.getLogger("genemapkit").handlers.clear()

    def tearDown(self):
        logging.getLogger("genemapkit").handlers.clear()

    def create_snapshot_and_database(self, root: Path):
        raw = root / "raw"
        raw.mkdir()
        mane = raw / "mane.txt"
        mane.write_text(
            "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\t"
            "RefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\n"
            "GeneID:1\tENSG00000000001.1\tHGNC:1\tGENE1\tgene one\t"
            "NM_000001.1\tNP_000001.1\tENST00000000001.1\t"
            "ENSP00000000001.1\tMANE Select\n"
            "GeneID:1\tENSG00000000001.1\tHGNC:1\tGENE1\tgene one\t"
            "NM_000002.1\tNP_000002.1\tENST00000000002.1\t"
            "ENSP00000000002.1\tMANE Plus Clinical\n",
            encoding="utf-8",
        )
        gene2ensembl = raw / "gene2ensembl"
        gene2ensembl.write_text(
            "#tax_id\tGeneID\tEnsembl_gene_identifier\tRNA_nucleotide_accession.version\t"
            "Ensembl_rna_identifier\tprotein_accession.version\t"
            "Ensembl_protein_identifier\n"
            "9606\t1\tENSG00000000001.1\tNM_000001.1\tENST00000000001.1\t"
            "NP_000001.1\tENSP00000000001.1\n"
            "9606\t1\tENSG00000000001.1\tNM_000002.1\tENST00000000002.1\t"
            "NP_000002.1\tENSP00000000002.1\n",
            encoding="utf-8",
        )
        databases = {}
        for source, path in (
            ("mane_summary", mane),
            ("ncbi_gene2ensembl", gene2ensembl),
        ):
            databases[source] = {
                "filename": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256sum(path),
            }
        (raw / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": "1.0.0",
                    "organism": "Homo sapiens",
                    "taxid": 9606,
                    "ensembl_release": 116,
                    "databases": databases,
                }
            ),
            encoding="utf-8",
        )
        database = root / "build" / "genemapkit.db"
        build_database(raw, database, sources=["mane_summary"])
        return raw, database

    def test_independent_benchmark_confirms_precise_transcript_protein_relationships(self):
        with tempfile.TemporaryDirectory() as directory:
            raw, database = self.create_snapshot_and_database(Path(directory))

            report = validate_database(
                database, raw, sources=["ncbi_gene2ensembl"], sample=None
            )
            row = next(
                item
                for item in report["benchmarks"]
                if item["input_type"] == "ensembl_transcript_id"
                and item["output_type"] == "ensembl_protein_id"
            )

            self.assertTrue(row["independent"])
            self.assertEqual(100.0, row["precision_pct"])
            self.assertEqual(100.0, row["recall_pct"])
            self.assertEqual(0, row["overmapped_inputs"])
            self.assertLess(row["precision_ci_low_pct"], 100.0)
            self.assertEqual(100.0, row["precision_ci_high_pct"])

    def test_validation_reuses_one_sqlite_connection_per_conversion(self):
        with tempfile.TemporaryDirectory() as directory:
            raw, database = self.create_snapshot_and_database(Path(directory))

            from genemapkit.core.converter import GeneConverter

            original_connect = GeneConverter._connect
            opened = 0

            def counted_connect(converter):
                nonlocal opened
                opened += 1
                return original_connect(converter)

            with patch.object(GeneConverter, "_connect", counted_connect):
                report = validate_database(
                    database,
                    raw,
                    sources=["ncbi_gene2ensembl"],
                    sample=None,
                )

            self.assertTrue(report["benchmarks"])
            # Initialization, metadata, diagnostics, then one connection per conversion.
            self.assertLess(opened, 12)

    def test_require_held_out_skips_leaky_truth_source(self):
        with tempfile.TemporaryDirectory() as directory:
            raw, database = self.create_snapshot_and_database(Path(directory))

            report = validate_database(
                database,
                raw,
                sources=["mane_summary"],
                require_held_out=True,
            )

            self.assertEqual([], report["benchmarks"])
            self.assertEqual("mane_summary", report["skipped"][0]["source"])

    def test_validation_writes_reports_and_package_cli_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, database = self.create_snapshot_and_database(root)
            output = root / "validation"

            result = CliRunner().invoke(
                cli,
                [
                    "validate",
                    "--db",
                    str(database),
                    "--raw-dir",
                    str(raw),
                    "--source",
                    "ncbi_gene2ensembl",
                    "--sample",
                    "0",
                    "--policy",
                    "authoritative",
                    "--out",
                    str(output),
                ],
            )

            self.assertEqual(0, result.exit_code, result.output)
            self.assertTrue((output / "validation_metrics.json").exists())
            self.assertTrue((output / "validation_benchmarks.csv").exists())
            self.assertTrue((output / "validation_summary.md").exists())
            self.assertFalse((output / "validation_disagreements.csv").exists())
            report = json.loads((output / "validation_metrics.json").read_text())
            self.assertEqual("authoritative", report["mapping_policy"])
            self.assertEqual(
                {"authoritative"},
                {row["mapping_policy"] for row in report["benchmarks"]},
            )


if __name__ == "__main__":
    unittest.main()
