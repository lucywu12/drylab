#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony.util.seq` module."""


import unittest
import codon_harmony.util


class TestCodon_harmony_util_seq(unittest.TestCase):
    """Tests for `codon_harmony.util.seq` module."""

    def setUp(self):
        """Set up test fixtures, if any."""
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        self.test_aa = Seq("TESTTESTTEST", IUPAC.protein)
        self.test_dna = Seq(
            "ACCGAGTCTACCACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna
        )

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_back_translate(self):
        """Test `codon_harmony.util.seq` -- unoptimized reverse translation"""
        assert self.test_aa.back_translate() == self.test_dna

        try:
            self.test_dna.back_translate()
        except ValueError as ve:
            assert ve.args[0] == "Nucleic acids cannot be back translated!"
        else:
            assert False
