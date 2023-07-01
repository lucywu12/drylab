#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `codon_harmony.util.logging` module."""


import unittest
from codon_harmony.util import logging


class TestCodon_harmony_util_logging(unittest.TestCase):
    """Tests for `codon_harmony.util.logging` module."""

    def setUp(self):
        """Set up test fixtures, if any."""

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_log_levels(self):
        """Test `codon_harmony.util.logging`"""
        from codon_harmony.util import logging

        with self.assertLogs("test", level="DETAIL") as cm:
            logging.getLogger("test").detail("test detail")
            logging.getLogger("test").output("test output")

        assert cm.output == ["DETAIL:test:test detail", "OUTPUT:test:test output"]
