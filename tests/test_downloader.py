"""Focused tests for Stage 0 source selection and manifest generation."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from download_databases import build_parser
from genemapkit.cli import cli
from genemapkit.core.downloader import DatabaseDownloader


class DatabaseDownloaderTests(unittest.TestCase):
    def test_default_excludes_extended_sources(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        default_sources = downloader.sources_for_tiers(("core", "crossref"))
        extended_sources = downloader.sources_for_tiers(("extended",))

        self.assertIn("hgnc", default_sources)
        self.assertIn("ensembl_refseq", default_sources)
        self.assertNotIn("ncbi_gene2go", default_sources)
        self.assertTrue(set(default_sources).isdisjoint(extended_sources))

    def test_validation_adds_independent_mapping_sources(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        validation_sources = downloader.validation_sources()

        self.assertEqual(
            ["ncbi_gene2ensembl", "ncbi_gene2refseq"], validation_sources
        )
        self.assertNotIn("ncbi_gene2go", validation_sources)
        self.assertNotIn("sifts_pdb", validation_sources)

    def test_download_validation_includes_default_and_validation_sources(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        with patch.object(
            downloader, "download_specific_sources", return_value={}
        ) as download:
            downloader.download_validation(include_default=True)

        selected = download.call_args.args[0]
        self.assertIn("hgnc", selected)
        self.assertIn("ensembl_refseq", selected)
        self.assertIn("ncbi_gene2ensembl", selected)
        self.assertIn("ncbi_gene2refseq", selected)
        self.assertNotIn("ncbi_gene2go", selected)

    def test_standalone_cli_accepts_all_with_validation(self):
        args = build_parser().parse_args(["--all", "--validation"])
        self.assertTrue(args.all)
        self.assertTrue(args.validation)

    def test_package_cli_all_with_validation_downloads_all_sources(self):
        runner = CliRunner()
        expected = {source: True for source in DatabaseDownloader.SOURCES}
        with patch.object(
            DatabaseDownloader, "download_all", return_value=expected
        ) as download_all:
            result = runner.invoke(cli, ["download", "--all", "--validation"])

        self.assertEqual(0, result.exit_code, result.output)
        download_all.assert_called_once_with(
            force_download=False, include_optional=True
        )

    def test_minimal_selection_requires_direct_input_output_coverage(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        sources = downloader.get_minimal_sources_for_conversion(
            "ensembl_transcript_id", ["refseq_id"]
        )

        self.assertEqual(["mane_summary"], sources)

    def test_legacy_refseq_namespace_is_normalized(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        sources = downloader.get_recommended_sources_for_conversion(
            "refseq_id", ["ensembl_transcript_id"]
        )

        self.assertIn("mane_summary", sources)
        self.assertIn("ensembl_refseq", sources)
        self.assertNotIn("hgnc", sources)

    def test_download_writes_manifest_and_skips_existing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            raw_dir = Path(directory)
            source_path = raw_dir / "hgnc_complete_set.txt"
            source_path.write_text("hgnc_id\tsymbol\nHGNC:11998\tTP53\n", encoding="utf-8")
            downloader = DatabaseDownloader(data_dir=directory)

            with patch.object(downloader, "_download_file") as download:
                self.assertTrue(downloader.download_source("hgnc"))
                download.assert_not_called()

            manifest = json.loads((raw_dir / "manifest.json").read_text(encoding="utf-8"))
            entry = manifest["databases"]["hgnc"]
            self.assertEqual("available", entry["status"])
            self.assertEqual(source_path.stat().st_size, entry["bytes"])
            self.assertEqual(64, len(entry["sha256"]))

    def test_snapshot_cannot_mix_ensembl_releases(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = {
                "schema_version": "1.0.0",
                "ensembl_release": 111,
                "databases": {},
            }
            Path(directory, "manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            downloader = DatabaseDownloader(data_dir=directory, ensembl_release=116)

            with self.assertRaisesRegex(ValueError, "different snapshot directory"):
                downloader._load_manifest()

    def test_unknown_identifier_namespace_is_rejected(self):
        downloader = DatabaseDownloader(data_dir=tempfile.mkdtemp())
        with self.assertRaisesRegex(ValueError, "Unsupported identifier"):
            downloader.get_minimal_sources_for_conversion("banana_id", ["symbol"])


if __name__ == "__main__":
    unittest.main()
