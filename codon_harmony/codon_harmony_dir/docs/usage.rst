=====
Usage
=====

.. argparse::
   :module: codon_harmony.codon_harmony
   :func: get_parser
   :prog: codon_harmony

Executing Codon Harmony as a script
-----------------------------------

    python codon_harmony/codon_harmony.py --input misc/INPUT_LIST.fasta --output out.fasta

To get started, create a conda environment from the ``environment.yml`` file::

    conda env create -f environment.yml

contents of ``misc/INPUT_LIST.fasta``:

.. code-block:: text

  >test_sequence1|can be optimized with `max_relax` set to 0.1
  HHHHHHHHHH
  >test_sequence2|cannot be optimized with `max_relax` set to 0.1
  ACDEFGHIKLMNPQRSTVWY
  >test_sequence3|can be optimized with `max_relax` set to 0.1, has extreme GC content
  FFFFFFFFFFFF

Using Codon Harmony in a project
--------------------------------

.. code-block:: python

  import codon_harmony
  codon_harmony.runner()

The ``runner`` function will handle parsing all command line arguments.
