import sys
from Bio import Seq


class SequenceContainer:

    def __init__(self, sequence):
        setattr(self, "start_fitness", sys.maxsize)
        setattr(self, "hairpins_fitness", sys.maxsize)
        assert isinstance(sequence, Seq.Seq)
        setattr(self, "sequence", sequence)
        setattr(self, "gc_fitness", sys.maxsize)
        setattr(self, "homo_fitness", sys.maxsize)
        setattr(self, "host_fitness", sys.maxsize)
        setattr(self, "repeats_fitness", sys.maxsize)
        setattr(self, "restriction_fitness", sys.maxsize)
        setattr(self, "splice_fitness", sys.maxsize)
        setattr(self, "dominated", False)
        setattr(self, "density", sys.maxsize)