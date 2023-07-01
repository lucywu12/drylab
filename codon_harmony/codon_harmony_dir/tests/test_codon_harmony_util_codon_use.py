#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony.util.codon_use` module."""


import unittest
from codon_harmony.util import codon_use


class TestCodon_harmony_util_codon_use(unittest.TestCase):
    """Tests for `codon_harmony.util.codon_use` module."""

    def setUp(self):
        """Set up test fixtures, if any."""
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        self.test_aa = Seq("TESTTESTTEST", IUPAC.protein)
        self.test_dna = Seq(
            "ACCGAGTCTACCACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna
        )
        self.known_counts = {"ACC": 6, "GAG": 3, "TCT": 3}

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_codon_count(self):
        """Test `codon_harmony.util.codon_use.codon_count`"""

        cdn_count = codon_use.count_codons(self.test_dna)
        for codon, count in cdn_count.items():
            if codon in self.known_counts:
                assert self.known_counts[codon] == count
            else:
                assert not count

    def test_calc_profile(self):
        """Test `codon_harmony.util.codon_use.calc_profile`"""
        profile = codon_use.calc_profile(codon_use.count_codons(self.test_dna))
        for codon, frequency in profile.items():
            if codon in self.known_counts:
                assert frequency == 1.0 and type(frequency) == float
            else:
                assert frequency == 0

    def test_calc_cra(self):
        """Test `codon_harmony.util.codon_use.calc_codon_relative_adaptiveness`"""
        cai = codon_use.calc_codon_relative_adaptiveness(
            codon_use.count_codons(self.test_dna)
        )
        for codon, relative_adaptiveness in cai.index.items():
            if codon in self.known_counts:
                assert (
                    relative_adaptiveness == 1.0
                    and type(relative_adaptiveness) == float
                )
            else:
                assert relative_adaptiveness == 0

    def test_load_host_table(self):
        """Test `codon_harmony.util.codon_use._load_host_table`"""
        profile = codon_use._load_host_table(9606)

        human_profile = {
            "TTT": 0.46,
            "TTC": 0.54,
            "TTA": 0.079,
            "TTG": 0.129,
            "CTT": 0.129,
            "CTC": 0.198,
            "CTA": 0.069,
            "CTG": 0.396,
            "ATT": 0.36,
            "ATC": 0.47,
            "ATA": 0.17,
            "ATG": 1.0,
            "GTT": 0.18,
            "GTC": 0.24,
            "GTA": 0.12,
            "GTG": 0.46,
            "TAT": 0.44,
            "TAC": 0.56,
            "TAA": 0.297,
            "TAG": 0.238,
            "CAT": 0.42,
            "CAC": 0.58,
            "CAA": 0.27,
            "CAG": 0.73,
            "AAT": 0.47,
            "AAC": 0.53,
            "AAA": 0.43,
            "AAG": 0.57,
            "GAT": 0.46,
            "GAC": 0.54,
            "GAA": 0.42,
            "GAG": 0.58,
            "TCT": 0.19,
            "TCC": 0.22,
            "TCA": 0.15,
            "TCG": 0.05,
            "CCT": 0.29,
            "CCC": 0.32,
            "CCA": 0.28,
            "CCG": 0.11,
            "ACT": 0.25,
            "ACC": 0.36,
            "ACA": 0.28,
            "ACG": 0.11,
            "GCT": 0.267,
            "GCC": 0.396,
            "GCA": 0.228,
            "GCG": 0.109,
            "TGT": 0.46,
            "TGC": 0.54,
            "TGA": 0.465,
            "TGG": 1.0,
            "CGT": 0.081,
            "CGC": 0.182,
            "CGA": 0.111,
            "CGG": 0.202,
            "AGT": 0.15,
            "AGC": 0.24,
            "AGA": 0.212,
            "AGG": 0.212,
            "GGT": 0.16,
            "GGC": 0.34,
            "GGA": 0.25,
            "GGG": 0.25,
        }

        for codon, frequency in profile.items():
            self.assertAlmostEqual(frequency, human_profile[codon], places=3)

    def test_process_host_table(self):
        """Test `codon_harmony.util.codon_use.process_host_table`"""
        raw_table = codon_use._load_host_table(9606)
        processed_table = codon_use.process_host_table(
            9606, threshold=0.0, table_path=None
        )
        self.assertAlmostEqual(processed_table, raw_table)

        thresholds = [0.1, 0.2]
        for thresh in thresholds:
            rare_codons = [codon for codon, freq in raw_table.items() if freq < thresh]
            processed_table = codon_use.process_host_table(
                9606, threshold=thresh, table_path=None
            )
            for codon in rare_codons:
                assert not processed_table[codon]

    def test_host_codon_usage(self):
        """Test `codon_harmony.util.codon_use.host_codon_usage`"""
        processed_table = codon_use.process_host_table(
            9606, threshold=0.1, table_path=None
        )
        codon_use_by_aa, host_profile, cra = codon_use.host_codon_usage(
            9606, threshold=0.1, table_path=None
        )
        for AA, codon_freqs in codon_use_by_aa.items():
            self.assertAlmostEqual(sum(codon_freqs[-1]), 1.0)

            for codon, freq in zip(*codon_freqs):
                assert host_profile[codon] == freq
        self.assertAlmostEqual(host_profile, processed_table)
        self.assertAlmostEqual(
            cra.index, codon_use.calc_codon_relative_adaptiveness(host_profile).index
        )

    def test_process_host_table_w_logging(self):
        """Test `codon_harmony.util.codon_use.process_host_table` with logging"""
        from codon_harmony.util import logging

        with self.assertLogs("codon_harmony.util.codon_use", level="DETAIL") as cm:
            processed_table = codon_use.process_host_table(
                9606, threshold=0.0, table_path=None
            )

        assert (
            "DETAIL:codon_harmony.util.codon_use:Pre-threshold host table:" in cm.output
        )
        assert "DETAIL:codon_harmony.util.codon_use:TGT: 0.46" in cm.output
