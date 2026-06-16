"""Fixture-based tests for the normalized Stage 1 database builder."""

import json
import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from genemapkit.cli import cli
from genemapkit.core.builder import build_database, convert_identifier
from genemapkit.core.downloader import sha256sum


class DatabaseBuilderTests(unittest.TestCase):
    def create_snapshot(self, root: Path) -> Path:
        raw = root / "raw"
        raw.mkdir()
        hgnc = raw / "hgnc_complete_set.txt"
        hgnc.write_text(
            "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tentrez_id\t"
            "ensembl_gene_id\trefseq_accession\tuniprot_ids\tccds_id\t"
            "omim_id\tucsc_id\tvega_id\tena\n"
            "HGNC:11998\tTP53\tP53\tTRP53\t7157\tENSG00000141510.18\t"
            "NM_000546.6\tP04637\tCCDS11118.1\t191170\t\t\t\n"
            "HGNC:5\tA1BG\tP53BP\t\t1\tENSG00000121410.14\t"
            "NM_130786.4\tP04217\tCCDS12976.1\t138670\t\t\t\n",
            encoding="utf-8",
        )
        mane = raw / "MANE.GRCh38.v1.5.summary.txt"
        mane.write_text(
            "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\t"
            "RefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\n"
            "GeneID:7157\tENSG00000141510.18\tHGNC:11998\tTP53\ttumor protein p53\t"
            "NM_000546.6\tNP_000537.3\tENST00000269305.9\t"
            "ENSP00000269305.4\tMANE Select\n",
            encoding="utf-8",
        )
        databases = {}
        for source, path in (("hgnc", hgnc), ("mane_summary", mane)):
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
        return raw

    def test_builds_version_normalized_gene_transcript_protein_index(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.create_snapshot(root)
            database = root / "build" / "genemapkit.db"
            report = root / "build" / "report.json"

            result = build_database(raw, database, report_path=report)

            self.assertTrue(database.exists())
            self.assertTrue(report.exists())
            self.assertTrue(result["complete_build"])
            self.assertEqual(
                ["ENSP00000269305"],
                convert_identifier(database, "7157", "entrez_id", "ensembl_protein_id"),
            )
            self.assertEqual(
                ["TP53"],
                convert_identifier(
                    database, "NM_000546.6", "refseq_rna_id", "symbol"
                ),
            )

    def test_alias_attaches_without_substring_matching(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "genemapkit.db"
            build_database(self.create_snapshot(root), database)

            self.assertEqual(
                ["TP53"],
                convert_identifier(database, "P53", "alias_symbol", "symbol"),
            )
            self.assertEqual(
                [],
                convert_identifier(database, "P5", "alias_symbol", "symbol"),
            )

    def test_direct_relationships_do_not_cross_transcripts_within_one_gene(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.create_snapshot(root)
            mane = raw / "MANE.GRCh38.v1.5.summary.txt"
            with mane.open("a", encoding="utf-8") as handle:
                handle.write(
                    "GeneID:7157\tENSG00000141510.18\tHGNC:11998\tTP53\t"
                    "tumor protein p53\tNM_999999.1\tNP_999999.1\t"
                    "ENST00000999999.1\tENSP00000999999.1\tMANE Plus Clinical\n"
                )
            manifest_path = raw / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["databases"]["mane_summary"].update(
                {"bytes": mane.stat().st_size, "sha256": sha256sum(mane)}
            )
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            database = root / "genemapkit.db"

            result = build_database(raw, database)

            self.assertGreater(result["source_relationship_rows"]["mane_summary"], 0)
            self.assertEqual(
                ["ENSP00000269305"],
                convert_identifier(
                    database,
                    "ENST00000269305",
                    "ensembl_transcript_id",
                    "ensembl_protein_id",
                ),
            )
            self.assertEqual(
                ["ENSP00000999999"],
                convert_identifier(
                    database,
                    "ENST00000999999",
                    "ensembl_transcript_id",
                    "ensembl_protein_id",
                ),
            )

    def test_ensembl_entrez_rejects_transcript_names_and_reports_them(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw"
            raw.mkdir()
            xref = raw / "ensembl_entrez.tsv"
            xref.write_text(
                "gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref\t"
                "db_name\tinfo_type\tsource_identity\txref_identity\tlinkage_type\n"
                "ENSG00000000001\tENST00000000001\tENSP00000000001\t1\t"
                "EntrezGene\tDEPENDENT\t-\t-\t-\n"
                "ENSG00000000002\tENST00000000002\tENSP00000000002\tGENE2-202\t"
                "EntrezGene_trans_name\tMISC\t-\t-\t-\n",
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
                            "ensembl_entrez": {
                                "filename": xref.name,
                                "bytes": xref.stat().st_size,
                                "sha256": sha256sum(xref),
                            }
                        },
                    }
                ),
                encoding="utf-8",
            )
            database = root / "genemapkit.db"

            report = build_database(raw, database)

            self.assertEqual(
                ["1"],
                convert_identifier(
                    database, "ENSG00000000001", "ensembl_gene_id", "entrez_id"
                ),
            )
            self.assertEqual(
                [],
                convert_identifier(
                    database, "ENSG00000000002", "ensembl_gene_id", "entrez_id"
                ),
            )
            self.assertEqual(
                1, report["rejected_records"]["ensembl_entrez"]["count"]
            )
            self.assertEqual(
                "GENE2-202",
                report["rejected_records"]["ensembl_entrez"]["examples"][0]["xref"],
            )

    def test_checksum_mismatch_fails_build(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.create_snapshot(root)
            (raw / "hgnc_complete_set.txt").write_text("changed", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                build_database(raw, root / "genemapkit.db")

    def test_existing_output_requires_force(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.create_snapshot(root)
            database = root / "genemapkit.db"
            build_database(raw, database)

            with self.assertRaisesRegex(FileExistsError, "Use --force"):
                build_database(raw, database)

    def test_package_cli_build_exposes_builder(self):
        runner = CliRunner()
        result_data = {"groups": 2, "mapping_rows": 5}
        with patch("genemapkit.cli.build_database", return_value=result_data) as build:
            result = runner.invoke(
                cli,
                [
                    "build",
                    "--raw-dir",
                    "snapshot",
                    "--out",
                    "index.db",
                    "--report",
                    "report.json",
                ],
            )
        logging.getLogger("genemapkit").handlers.clear()

        self.assertEqual(0, result.exit_code, result.output)
        build.assert_called_once()

    def test_same_snapshot_produces_same_database_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.create_snapshot(root)
            first = build_database(raw, root / "first.db")
            second = build_database(raw, root / "second.db")

            self.assertEqual(first["database_sha256"], second["database_sha256"])


if __name__ == "__main__":
    unittest.main()
