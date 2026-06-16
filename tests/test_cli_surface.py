"""Regression tests for the unified GeneMapKit CLI command surface."""

import tempfile
import unittest
import logging
from pathlib import Path
from unittest.mock import patch

from click.testing import CliRunner

from genemapkit.cli import cli
from genemapkit.core.downloader import DatabaseDownloader


class CliSurfaceTests(unittest.TestCase):
    def tearDown(self):
        logging.getLogger("genemapkit").handlers.clear()

    def test_every_command_help_page_loads(self):
        runner = CliRunner()
        commands = [
            "download",
            "download-for",
            "build",
            "validate",
            "compare-mygene",
            "benchmark-external",
            "benchmark-publication",
            "benchmark-paths",
            "convert-one",
            "convert",
            "info",
            "status",
        ]
        for command in commands:
            with self.subTest(command=command):
                result = runner.invoke(cli, [command, "--help"])
                self.assertEqual(0, result.exit_code, result.output)

    def test_status_accepts_raw_directory_aliases(self):
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as directory:
            raw_dir = Path(directory)
            (raw_dir / "hgnc_complete_set.txt").write_text(
                "hgnc_id\tsymbol\n", encoding="utf-8"
            )
            for option in ("--raw-dir", "--out-dir"):
                with self.subTest(option=option):
                    result = runner.invoke(cli, ["status", option, str(raw_dir)])
                    self.assertEqual(0, result.exit_code, result.output)
                    self.assertIn("available  hgnc", result.output)
                    self.assertIn("ncbi_gene2ensembl [validation]", result.output)

    def test_info_output_is_portable_and_current(self):
        result = CliRunner().invoke(cli, ["info"])
        self.assertEqual(0, result.exit_code, result.output)
        self.assertIn("refseq_rna_id", result.output)
        self.assertIn("refseq_protein_id", result.output)
        self.assertNotIn("refseq_id\n", result.output)

    def test_download_accepts_raw_directory_aliases(self):
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as directory:
            for option in ("--out-dir", "--raw-dir", "--out"):
                with self.subTest(option=option):
                    with patch.object(
                        DatabaseDownloader, "download_default", return_value={}
                    ):
                        result = runner.invoke(cli, ["download", option, directory])
                    self.assertEqual(0, result.exit_code, result.output)


if __name__ == "__main__":
    unittest.main()
