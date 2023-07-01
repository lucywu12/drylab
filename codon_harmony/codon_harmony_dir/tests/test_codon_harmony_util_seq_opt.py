#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony.util.seq_opt` module."""


import unittest
from codon_harmony.util import seq_opt

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class TestCodon_harmony_util_seq_opt(unittest.TestCase):
    """Tests for `codon_harmony.util.seq_opt` module."""

    def setUp(self):
        """Set up test fixtures, if any."""
        self.test_aa = Seq("TESTTESTTEST", IUPAC.protein)
        self.test_dna = Seq(
            "ACCGAGTCTACCACCGAGTCTACCACCGAGTCTACC", IUPAC.unambiguous_dna
        )

        from codon_harmony.util import codon_use

        self.codon_use_table, self.host_profile, _ = codon_use.host_codon_usage(
            9606, 0.0
        )
        self.codons_count = codon_use.count_codons(self.test_dna)

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_mutate_codon(self):
        """Test `codon_harmony.util.seq_opt.mutate_codon`"""
        # Phe -- two codons
        codon_in = "TTT"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out == "TTC"

        # Leu -- six codns
        codon_in = "TTA"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out in ["TTG", "CTT", "CTC", "CTA", "CTG"]

        # Met -- one codon
        codon_in = "ATG"
        codon_out = seq_opt.mutate_codon(codon_in, self.codon_use_table)
        assert codon_out == codon_in

    def test_resample_codons(self):
        """Test `codon_harmony.util.seq_opt.resample_codons`"""
        test_dna_local = seq_opt.resample_codons(self.test_dna, self.codon_use_table)
        assert test_dna_local.translate() == self.test_dna.translate()
        assert test_dna_local.translate() == self.test_aa

    def test_compare_profiles(self):
        """Test `codon_harmony.util.seq_opt.compare_profiles`"""
        known_counts = {"ACC": 6, "GAG": 3, "TCT": 3}
        table, diff = seq_opt.compare_profiles(
            self.codons_count, self.host_profile, 1.0
        )
        self.assertAlmostEqual(diff, 0.5769, places=4)

        for codon, count_dict in table.items():
            if codon in known_counts:
                assert count_dict["input_count"] == known_counts[codon]
            else:
                assert not count_dict["input_count"]
            assert (
                count_dict["ideal_usage_abs"] + count_dict["difference_abs"]
                == count_dict["input_count"]
            )

        relevant_codons = {
            "Thr": ["ACC", "ACT", "ACG", "ACA"],
            "Glu": ["GAA", "GAG"],
            "Ser": ["TCA", "TCC", "TCG", "TCT", "AGT", "AGC"],
        }
        for _, codons in relevant_codons.items():
            assert not sum(table[codon]["difference"] for codon in codons)

        # modify the count to trigger a different code path
        codon_count_local = self.codons_count.copy()
        codon_count_local["ACG"] += 1
        codon_count_local["ACC"] -= 2
        table, diff = seq_opt.compare_profiles(
            codon_count_local, self.host_profile, 0.9
        )
        self.assertAlmostEqual(diff, 0.4545, places=4)
        for _, codons in relevant_codons.items():
            assert not sum(table[codon]["difference"] for codon in codons)

    def test_harmonize_codon_use_with_host(self):
        """Test `codon_harmony.util.seq_opt.harmonize_codon_use_with_host`"""
        mutation_profile, _ = seq_opt.compare_profiles(
            self.codons_count, self.host_profile, 1.0
        )
        test_dna_local = seq_opt.harmonize_codon_use_with_host(
            self.test_dna, mutation_profile
        )
        test_dna_local = str(test_dna_local)
        computed_counts = {}
        for i in range(len(test_dna_local) // 3):
            codon = test_dna_local[3 * i : 3 * (i + 1)]
            if codon not in computed_counts:
                computed_counts[codon] = 0
            computed_counts[codon] += 1

        for codon, count in computed_counts.items():
            assert mutation_profile[codon]["ideal_usage_abs"] == count

    def test_resample_codons_and_enforce_host_profile(self):
        """Test `codon_harmony.util.seq_opt.resample_codons_and_enforce_host_profile`"""
        test_dna_local = seq_opt.resample_codons_and_enforce_host_profile(
            self.test_dna, self.codon_use_table, self.host_profile, 1.0
        )
        assert test_dna_local.translate() == self.test_aa

    def test_gc_scan(self):
        """Test `codon_harmony.util.seq_opt.gc_scan`"""
        from codon_harmony.data import GCParams

        def _check_gc(gc):
            from Bio.SeqUtils import GC

            test_dna_local = seq_opt.gc_scan(self.test_dna, self.codon_use_table, gc)
            try:
                codon_window = gc.window_size // 3
            except TypeError:
                # this is triggered for gc3 -- "x0.5" -> (0.5 * len(dna)) // 3
                codon_window = int(float(gc.window_size[1:]) * len(test_dna_local)) // 3

            for i in range(0, len(test_dna_local) // 3):
                window = slice(
                    i * 3,
                    (i + codon_window) * 3
                    if (i + codon_window) * 3 < len(test_dna_local)
                    else len(test_dna_local),
                )
                gc_percent = GC(test_dna_local[window]) / 100
                assert gc_percent >= gc.low and gc_percent <= gc.high

            return test_dna_local

        gc1 = GCParams(name="test1", window_size=6, low=0.1, high=0.8)
        dna_should_be_the_same = _check_gc(gc1)
        assert dna_should_be_the_same == self.test_dna

        gc2 = GCParams(name="test2", window_size=6, low=0.1, high=0.5)
        dna_should_be_different = _check_gc(gc2)
        assert dna_should_be_different != self.test_dna

        gc3 = GCParams(name="test3", window_size="x0.5", low=0.1, high=0.5)
        dna_should_be_different = _check_gc(gc3)
        assert dna_should_be_different != self.test_dna

    def test_remove_restriction_sites(self):
        """Test `codon_harmony.util.seq_opt.remove_restriction_sites`"""
        from Bio import Restriction

        def _look_for_site(site, re_name, should_match=False):
            dna = Seq(site + str(self.test_dna), IUPAC.unambiguous_dna)
            proc_dna = seq_opt.remove_restriction_sites(
                dna,
                self.codon_use_table,
                Restriction.RestrictionBatch([Restriction.AllEnzymes.get(re_name)]),
            )
            assert (dna == proc_dna) == should_match

        # self.test_dna has the site for neither NdeI nor XhoI
        # prepend the site, confirm that it is replaced.
        _look_for_site("CATATG", "NdeI")
        _look_for_site("CTCGAG", "XhoI")
        # prepend the NdeI site, search against XhoI -- dna should be unchanged
        _look_for_site("CATATG", "XhoI", should_match=True)

    def test_remove_start_sites(self):
        """Test `codon_harmony.util.seq_opt.remove_start_sites`"""
        from codon_harmony.data import RibosomeBindingSites

        start = "TTG"
        dna = Seq(
            "CATATGGGGGGCTCGAGA" + start + str(self.test_dna), IUPAC.unambiguous_dna
        )
        rbs = {"test": "GGGGG"}
        proc_dna = seq_opt.remove_start_sites(dna, self.codon_use_table, rbs)

        assert dna != proc_dna
        assert dna[5:10] == rbs["test"]
        assert proc_dna[5:10] != rbs["test"]

        proc_dna = seq_opt.remove_start_sites(self.test_dna, self.codon_use_table, rbs)
        assert proc_dna == self.test_dna

    def test_remove_repeating_sequences(self):
        """Test `codon_harmony.util.seq_opt.remove_repeating_sequences`"""
        # test_dna contains repeats -- should be flagged here
        proc_dna = seq_opt.remove_repeating_sequences(
            self.test_dna, self.codon_use_table, 9
        )
        assert proc_dna != self.test_dna

        # test for three consecutive identical codons
        test_dna = Seq("TTGTTGTTG", IUPAC.unambiguous_dna)
        proc_dna = seq_opt.remove_repeating_sequences(test_dna, self.codon_use_table, 9)
        assert proc_dna != test_dna

    def test_remove_local_homopolymers(self):
        """Test `codon_harmony.util.seq_opt.remove_local_homopolymers`"""
        proc_dna = seq_opt.remove_local_homopolymers(
            self.test_dna, self.codon_use_table
        )
        assert proc_dna == self.test_dna

        # add a short stretch with >4 repeated NTs -- should be flagged
        test_dna = Seq("GTTTTT" + str(self.test_dna), IUPAC.unambiguous_dna)
        proc_dna = seq_opt.remove_local_homopolymers(test_dna, self.codon_use_table)
        assert proc_dna != test_dna

    def test_remove_hairpins(self):
        """Test `codon_harmony.util.seq_opt.remove_hairpins`"""
        proc_dna = seq_opt.remove_hairpins(self.test_dna, self.codon_use_table)
        assert proc_dna == self.test_dna

        stem_seq = Seq("CATATGGGGGGCTCGAGA", IUPAC.unambiguous_dna)
        test_dna = Seq(
            str(stem_seq) + "TTGTTGTTGA" + str(stem_seq.reverse_complement()),
            IUPAC.unambiguous_dna,
        )
        proc_dna = seq_opt.remove_hairpins(test_dna, self.codon_use_table)
        assert proc_dna != self.test_dna

    def test_remove_splice_sites(self):
        """Test `codon_harmony.uti.seq_opt.remove_splice_sites`"""
        known_donors = [
            "CGGTGAAT",
            "AAGTACGG",
            "CTGTGAGT",
            "TAGTGTGT",
            "TGGTCAGT",
            "AGGTGAAA",
            "GGGTGAGC",
            "AGGTACCT",
            "TGGTTAGA",
        ]

        known_acceptors = [
            "CTGCTACAGC",
            "GACCCACAGG",
            "CTGATTCAGG",
            "CTCTTGCAGT",
            "TCTGGTTAGA",
        ]

        dna_sequence = Seq(
            "ATGACAGCCAAGATGGAGACCACATTCTATGACGACGCCCTGAATGCTTCCTTCCTTCCCAGTGA"
            + "GTCTGGACCATACGGCTATTCTAATCCCAAAATCCTGAAGCAGTCCATGACACTGAACCTCGC"
            + "CGACCCCGTGGGCTCCCTTAAACCGCATCTCAGAGCTAAGAATTCCGACCTGCTTACTTCTCC"
            + "TGATGTGGGACTCCTGAAGCTCGCCAGCCCCGAGCTTGAACGGCTGATTATTCAGAGCTCCAA"
            + "TGGACACATCACCACCACTCCAACTCCCACACAGTTCCTGTGTCCCAAGAACGTTACCGACGA"
            + "GCAGGAGGGCTTCGCCGAAGGCTTTGTGCGCGCACTGGCTGAGCTTCACAGTCAAAATACTCT"
            + "GCCAAGCGTCACTTCTGCCGCCCAGCCTGTGAATGGAGCCGGGATGGTGGCACCTGCTGTTGC"
            + "TTCTGTTGCCGGCGGCTCCGGGTCTGGAGGCTTCAGTGCTTCTCTCCACAGTGAGCCGCCTGT"
            + "GTACGCTAACCTGTCCAACTTTAATCCTGGCGCTCTTAGCAGTGGCGGAGGCGCTCCTAGCTA"
            + "CGGGGCCGCAGGCCTCGCTTTTCCTGCCCAGCCACAGCAACAGCAGCAACCACCACATCATCT"
            + "GCCTCAGCAGATGCCAGTCCAACACCCACGCCTCCAGGCTCTGAAGGAGGAGCCTCAGACAGT"
            + "CCCTGAAATGCCTGGAGAAACCCCTCCACTCTCCCCAATCGATATGGAGAGCCAGGAACGGAT"
            + "CAAGGCAGAACGCAAGCGGATGAGGAATAGGATCGCCGCATCCAAGTGCCGCAAACGCAAGCT"
            + "GGAGCGGATTGCCCGGCTGGAGGAGAAGGTTAAAACCCTGAAGGCTCAAAACAGCGAGCTGGC"
            + "CTCTACCGCTAATATGCTTCGCGAGCAGGTTGCACAGCTCAAGCAGAAAGTGATGAACCACGT"
            + "CAACTCTGGATGCCAGCTCATGCTGACACAGCAGCTGCAAACATTCGGCAGTGGGGCAACCAA"
            + "CTTTTCCCTGCTGAAACAGGCCGGGGATGTTGAGGAGAACCCTGGTCCA",
            IUPAC.unambiguous_dna,
        )
        dna_sequence = seq_opt.remove_splice_sites(dna_sequence, self.codon_use_table)
        dna_sequence = str(dna_sequence)

        for site in known_donors + known_acceptors:
            assert site not in dna_sequence
