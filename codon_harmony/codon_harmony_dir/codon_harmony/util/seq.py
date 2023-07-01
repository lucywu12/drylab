import Bio

from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, Seq
from Bio.SeqUtils import CodonUsage, seq3


def back_translate(self):
    """Return the DNA sequence from an amino acid sequence by creating a new Seq object.
    The first codon in the synonymous codons list is always chosen for each amino acid;
    codon optimization is required after back translation.

    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import IUPAC
    >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
    >>> my_protein
    Seq('MAIVMGR', IUPACProtein())
    >>> my_protein.back_translate()
    Seq('ATGGCCATTGTAATGGGCCGCTG', IUPACUnambiguousDNA())

    Trying to back-transcribe a DNA or RNA sequence raises an
    exception:

    >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUG", IUPAC.unambiguous_rna)
    >>> messenger_rna.back_translate()
    Traceback (most recent call last):
    ...
    ValueError: Nucleic acids cannot be back translated!
    """
    base = Bio.Alphabet._get_base_alphabet(self.alphabet)
    if not isinstance(base, Bio.Alphabet.ProteinAlphabet):
        raise ValueError("Nucleic acids cannot be back translated!")

    # always use the first codon in the synonymous codons list for each AA
    return Seq(
        "".join(
            [
                CodonUsage.SynonymousCodons[seq3(AA).upper()][0]
                for AA in str(self).upper()
            ]
        ),
        IUPAC.unambiguous_dna,
    )


Seq.back_translate = back_translate
