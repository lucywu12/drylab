#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony` package."""


import unittest

from codon_harmony.util import codon_table_downloader


class TestCodon_tools(unittest.TestCase):
    """Tests for `codon_tools` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.args_to_parse = [
            "--input",
            "misc/INPUT_LIST.fasta",
            "--output",
            "out.fasta",
        ]

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_save_codon_table_to_disk(self):
        """Test dowloading codon tables for offline use."""
        taxid = "9606"
        outfile = "human_codon_table.json"
        codon_table_downloader.save_codon_table_to_disk(taxid, outfile)

        try:
            codon_table_downloader.save_codon_table_to_disk("Notta species", outfile)
        except ValueError as ve:
            assert ve.args[0].startswith('"Notta species" is not a valid host id.')
        else:
            assert False

        codon_table_downloader.save_codon_table_to_disk(taxid, outfile.split(".")[0])
