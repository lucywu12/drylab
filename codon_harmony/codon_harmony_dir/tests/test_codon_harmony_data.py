#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_tools.data` module."""


import unittest

from codon_harmony import data


class TestCodon_tools_data(unittest.TestCase):
    """Tests for `codon_tools.data` module."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_species_taxid_map(self):
        """Test mapping from species to NCBI taxonomy id."""
        taxid = data._tax_id_from_species("homo_sapiens")
        assert taxid == 9606

        taxid = data._tax_id_from_species("homo sapiens")
        assert taxid == 9606

        taxid = data._tax_id_from_species("homo+sapiens")
        assert taxid == 9606

        taxid = data._tax_id_from_species("Homo sapiens")
        assert taxid == 9606

        taxid = data._tax_id_from_species("Escherichia coli B")
        assert taxid == 37762

        taxid = data._tax_id_from_species("Escherichia coli B str. REL606")
        assert taxid == 413997

        taxid = data._tax_id_from_species("Escherichia coli K12")
        assert taxid == 83333

        taxid = data._tax_id_from_species("Saccharomyces cerevisiae")
        assert taxid == 4932

        try:
            taxid = data._tax_id_from_species("Notta species")
        except ValueError as ve:
            assert ve.args[0].startswith('"Notta species" is not a valid host id.')
        else:
            assert False

    def test_codon_usage_from_id(self):
        """Test grabbing codon usage table from species or NCBI taxonomy id."""
        from_species = data.codon_tables("Homo sapiens")
        from_taxid = data.codon_tables(9606)
        assert from_species == from_taxid
        assert from_taxid["UUU"] == 0.46

        try:
            data.codon_tables(123456789101112)
        except ValueError as ve:
            assert ve.args[0].startswith('"123456789101112" is not a valid host id.')
        else:
            assert False

    def test_codon_usage_from_local_table(self):
        """Test loading a JSON-formatted codon usage table."""
        table_path = "misc/human_codon_table.json"
        from_file = data.codon_tables("irrelevant", table_path=table_path)
        assert from_file["UUU"] == 0.46

    def test_restriction_enzymes(self):
        """Test creating RestrictionEnzymesBatch from a list of restriction enzymes"""
        list_of_res_enz = ["XhoI", "HpaI", "PstI", "EcoRV", "NcoI", "BamHI"]
        reb = data.RestrictionEnzymes(list_of_res_enz)
        for res_enz in list_of_res_enz:
            assert res_enz in reb
        assert "NdeI" not in reb
