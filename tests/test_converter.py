"""Tests for the SQLite-backed Stage 2 converter and conversion CLI."""

import json
import logging
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from genemapkit.cli import cli
from genemapkit.core.builder import build_database
from genemapkit.core.converter import GeneConverter
from genemapkit.core.downloader import sha256sum


class SQLiteConverterTests(unittest.TestCase):
    def tearDown(self):
        logging.getLogger("genemapkit").handlers.clear()

    def create_database(self, root: Path) -> Path:
        raw = root / "raw"
        raw.mkdir()
        hgnc = raw / "hgnc_complete_set.txt"
        hgnc.write_text(
            "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tentrez_id\t"
            "ensembl_gene_id\trefseq_accession\tuniprot_ids\tccds_id\t"
            "omim_id\tucsc_id\tvega_id\tena\n"
            "HGNC:11998\tTP53\tP53|SHARED\tTRP53\t7157\t"
            "ENSG00000141510.18\tNM_000546.6\tP04637\tCCDS11118.1\t"
            "191170\t\t\t\n"
            "HGNC:5\tA1BG\tSHARED\t\t1\tENSG00000121410.14\t"
            "NM_130786.4\tP04217\tCCDS12976.1\t138670\t\t\t\n",
            encoding="utf-8",
        )
        history = raw / "gene_history"
        history.write_text(
            "#tax_id\tGeneID\tDiscontinued_GeneID\tDiscontinued_Symbol\t"
            "Discontinue_Date\n"
            "9606\t7157\t999\tOLDTP53\t20200101\n",
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
                        },
                        "ncbi_gene_history": {
                            "filename": history.name,
                            "bytes": history.stat().st_size,
                            "sha256": sha256sum(history),
                        }
                    },
                }
            ),
            encoding="utf-8",
        )
        database = root / "build" / "genemapkit.db"
        build_database(raw, database)
        return database

    def create_policy_database(self, root: Path) -> Path:
        raw = root / "raw"
        raw.mkdir()
        mane = raw / "mane.tsv"
        mane.write_text(
            "#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tRefSeq_nuc\t"
            "RefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\n"
            "GeneID:7157\tENSG00000141510\tHGNC:11998\tTP53\tNM_000001.1\t"
            "NP_000001.1\tENST00000269305.9\tENSP00000269305.4\tMANE Select\n"
            "GeneID:7157\tENSG00000141510\tHGNC:11998\tTP53\tNM_000003.1\t"
            "NP_000003.1\tENST00000269305.9\tENSP00000269305.4\t"
            "MANE Plus Clinical\n",
            encoding="utf-8",
        )
        ensembl_refseq = raw / "ensembl_refseq.tsv"
        ensembl_refseq.write_text(
            "gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref\t"
            "db_name\tinfo_type\n"
            "ENSG00000141510\tENST00000269305\tENSP00000269305\tNM_000001.1\t"
            "RefSeq_mRNA\tDIRECT\n"
            "ENSG00000141510\tENST00000269305\tENSP00000269305\tNM_000002.2\t"
            "RefSeq_mRNA\tDIRECT\n",
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
                        "mane_summary": {
                            "filename": mane.name,
                            "bytes": mane.stat().st_size,
                            "sha256": sha256sum(mane),
                        },
                        "ensembl_refseq": {
                            "filename": ensembl_refseq.name,
                            "bytes": ensembl_refseq.stat().st_size,
                            "sha256": sha256sum(ensembl_refseq),
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        database = root / "build" / "genemapkit.db"
        build_database(raw, database)
        return database

    def test_queries_sqlite_and_preserves_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            results = converter.query(7157, "entrez_id", "ensembl_gene_id")

            self.assertEqual({"ENSG00000141510"}, {row.output_id for row in results})
            self.assertEqual({"unique"}, {row.mapping_status for row in results})
            self.assertEqual({"hgnc"}, {row.output_source for row in results})

    def test_symbol_policy_and_exact_alias_matching(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            self.assertEqual(
                [],
                converter.query(
                    "P53", "symbol", "ensembl_gene_id", symbol_policy="approved"
                ),
            )
            self.assertEqual(
                {"ENSG00000141510"},
                {
                    row.output_id
                    for row in converter.query(
                        "P53", "symbol", "ensembl_gene_id", symbol_policy="all"
                    )
                },
            )
            self.assertEqual(
                [],
                converter.query(
                    "P5", "symbol", "ensembl_gene_id", symbol_policy="all"
                ),
            )

    def test_shared_alias_is_reported_as_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            results = converter.query("SHARED", "symbol", "entrez_id")

            self.assertEqual({"1", "7157"}, {row.output_id for row in results})
            self.assertEqual({"ambiguous"}, {row.mapping_status for row in results})
            self.assertEqual({2}, {row.mapping_count for row in results})
            self.assertEqual({"alias_symbol"}, {row.matched_input_type for row in results})

    def test_retired_targets_are_opt_in(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            current = converter.query("7157", "entrez_id", "entrez_id")
            with_retired = converter.query(
                "7157", "entrez_id", "entrez_id", include_retired=True
            )

            self.assertEqual(set(), {row.output_id for row in current})
            self.assertEqual({"999"}, {row.output_id for row in with_retired})

    def test_mapping_policies_select_direct_mane_and_authoritative_results(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_policy_database(Path(directory)))

            all_supported = converter.query(
                "ENST00000269305",
                "ensembl_transcript_id",
                "refseq_rna_id",
                mapping_policy="all-supported",
            )
            direct_only = converter.query(
                "ENST00000269305",
                "ensembl_transcript_id",
                "refseq_rna_id",
                mapping_policy="direct-only",
            )
            corroborated = converter.query(
                "ENST00000269305",
                "ensembl_transcript_id",
                "refseq_rna_id",
                mapping_policy="corroborated",
            )
            mane = converter.query(
                "ENST00000269305",
                "ensembl_transcript_id",
                "refseq_rna_id",
                mapping_policy="mane-select",
            )
            authoritative = converter.query(
                "ENST00000269305",
                "ensembl_transcript_id",
                "refseq_rna_id",
                mapping_policy="authoritative",
            )

            expected_all = {"NM_000001", "NM_000002", "NM_000003"}
            self.assertEqual(expected_all, {row.output_id for row in all_supported})
            self.assertEqual(expected_all, {row.output_id for row in direct_only})
            self.assertEqual({"NM_000001"}, {row.output_id for row in corroborated})
            self.assertEqual({"NM_000001"}, {row.output_id for row in mane})
            self.assertEqual(
                {"NM_000001", "NM_000003"},
                {row.output_id for row in authoritative},
            )
            self.assertEqual({"mane_summary"}, {row.output_source for row in mane})
            self.assertEqual({"mane-select"}, {row.mapping_policy for row in mane})
            self.assertIn("MANE", mane[0].selection_reason)
            self.assertEqual({2}, {row.source_support_count for row in corroborated})
            self.assertEqual(
                {"canonical_mane"}, {row.agreement_status for row in corroborated}
            )

    def test_invalid_mapping_policy_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            with self.assertRaisesRegex(ValueError, "Unknown mapping policy"):
                converter.query("TP53", "symbol", "entrez_id", mapping_policy="first")

    def test_batch_wide_preserves_multiple_values_and_long_preserves_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = GeneConverter(self.create_database(Path(directory)))

            wide = converter.convert_batch(
                ["SHARED"], "symbol", ["entrez_id"], output_format="wide"
            )
            long = converter.convert_batch(
                ["TP53"], "symbol", ["ensembl_gene_id"], output_format="long"
            )

            self.assertEqual("1;7157", wide.loc[0, "entrez_id"])
            self.assertEqual("ambiguous", wide.loc[0, "entrez_id_mapping_status"])
            self.assertIn("output_source", long.columns)
            self.assertIn("mapping_policy", long.columns)
            self.assertIn("selection_reason", long.columns)
            self.assertEqual({"hgnc"}, set(long["output_source"]))

    def test_package_cli_converts_batch_from_sqlite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = self.create_database(root)
            source = root / "input.csv"
            output = root / "output.csv"
            source.write_text("gene\nTP53\nSHARED\n", encoding="utf-8")

            result = CliRunner().invoke(
                cli,
                [
                    "convert",
                    str(source),
                    "--column",
                    "gene",
                    "--input-type",
                    "symbol",
                    "--output-type",
                    "entrez_id",
                    "--format",
                    "wide",
                    "--db",
                    str(database),
                    "--output",
                    str(output),
                ],
            )

            self.assertEqual(0, result.exit_code, result.output)
            converted = pd.read_csv(output, dtype=str)
            self.assertEqual("7157", converted.loc[0, "entrez_id"])
            self.assertEqual("1;7157", converted.loc[1, "entrez_id"])

    def test_package_cli_convert_one_exposes_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            database = self.create_database(Path(directory))

            result = CliRunner().invoke(
                cli,
                [
                    "convert-one",
                    "TP53",
                    "--from",
                    "symbol",
                    "--to",
                    "ensembl_gene_id",
                    "--db",
                    str(database),
                ],
            )

            self.assertEqual(0, result.exit_code, result.output)
            self.assertIn("ENSG00000141510", result.output)
            self.assertIn("source=hgnc", result.output)

    def test_package_cli_convert_one_accepts_mapping_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            database = self.create_policy_database(Path(directory))

            result = CliRunner().invoke(
                cli,
                [
                    "convert-one",
                    "ENST00000269305",
                    "--from",
                    "ensembl_transcript_id",
                    "--to",
                    "refseq_rna_id",
                    "--policy",
                    "mane-select",
                    "--db",
                    str(database),
                ],
            )

            self.assertEqual(0, result.exit_code, result.output)
            self.assertIn("NM_000001", result.output)
            self.assertNotIn("NM_000002", result.output)
            self.assertIn("policy=mane-select", result.output)


if __name__ == "__main__":
    unittest.main()
