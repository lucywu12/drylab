=============
Codon Harmony
=============


.. image:: https://img.shields.io/pypi/v/codon_harmony.svg
        :target: https://pypi.python.org/pypi/codon_harmony
        
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
        :target: https://opensource.org/licenses/MIT
        :alt: MIT License

.. image:: https://img.shields.io/travis/lyellbio/codon-harmony.svg
        :target: https://travis-ci.org/lyellbio/codon-harmony

.. image:: https://readthedocs.org/projects/codon-harmony/badge/?version=latest
        :target: https://codon-harmony.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation status

.. image:: https://codecov.io/gh/lyellbio/codon-harmony/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/lyellbio/codon-harmony
        :alt: Coverage report

.. image:: https://pyup.io/repos/github/lyellbio/codon-harmony/shield.svg
     :target: https://pyup.io/repos/github/lyellbio/codon-harmony/
     :alt: Updates

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
     :target: https://github.com/ambv/black
     :alt: Code style: black


Amino acid reverse translation and DNA optimization tool based on species-specific codon-use distributions.
Species-specifc data can be found on the `Codon Usage Database`_ using the `NCBI Taxonomy database`_ id (e.g. 413997) or the organism's Latin name (e.g. *Escherichia coli* B). Mapping species names to Taxonomy IDs can be done here_.

.. _`Codon Usage Database`: http://www.kazusa.or.jp/codon
.. _`NCBI Taxonomy database`: http://www.ncbi.nlm.nih.gov/taxonomy
.. _here: https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi

* Documentation: https://codon-harmony.readthedocs.io


Features
--------

1. Reverse translates input amino acid sequence to DNA.
2. Calculates the host's per-AA codon usage profile – codons used less than a specified threshold (defaults to 10%) are dropped.
3. Compares the reverse-translated DNA sequence to the host profile, determines which codons are overused/underused.
4. Stochastically mutates codons according to host profile.
5. Ranks sequences by codon adaptation index relative to host
6. Processes DNA to remove unwanted features:

   * high GC content within a sliding window and across the entire sequence
   * unwanted restriction sites
   * alternate start positions (GA-rich regions 18 bp upstream of ATG/GTG/TTG)
   * 3-consecutive identical codons and 9-mer repeat chunks
   * areas with more than 4 (variable) consecutive identical bps ("local homopolymers")
   * RNA hairpins, detected by looking for 10-mers with reverse complements (including wobble bases) in the sequence
   * RNA splice sites, detected by similarity to consensus donor and acceptor site sequences
   
The process is repeated from step 3 for a specified number of cycles (defaults to 1000) OR until the per-AA codon profile of current DNA and host profile matches (within tolerance).

Future work
-----------

- More advanced RNA-structure removal

  * CONTRAfold_ – overkill for now
  * nupack_ – overkill for now

.. _CONTRAfold: http://contra.stanford.edu/contrafold/
.. _nupack: http://nupack.org
