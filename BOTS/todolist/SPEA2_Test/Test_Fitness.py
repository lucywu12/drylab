from todolist.SPEA2.Fitness_Functions import eval_gc_content, eval_hairpins, eval_homopolymers, eval_host, eval_repeats, \
    eval_restriction_sites, eval_splice_sites, eval_start_sites
from Bio.Alphabet import IUPAC
from Bio import Seq
from todolist.SPEA2.Bio_Structures import GC_content, RestrictionEnzymes
from todolist.SPEA2.Sequence_Container import SequenceContainer


def test_gc(sequence, gc):
    # assert eval_gc_content(sequence, gc) is 5
    print("GC", eval_gc_content(sequence, gc))


def test_hair(sequence):
    # assert eval_hairpins(sequence) is 1
    print("hair", eval_hairpins(sequence))


def test_homo(sequence):
    # assert eval_homopolymers(sequence) is 5
    print("homo", eval_homopolymers(sequence))


def test_host(sequence, ancestor_sequence):
    # assert eval_host(sequence, ancestor_sequence) is 7
    print("host", eval_host(sequence, ancestor_sequence))


def test_repeats(sequence):
    # assert eval_repeats(sequence) is 10
    print("repeats", eval_repeats(sequence))


def test_rest(sequence, restriction_sites):
    # assert eval_restriction_sites(sequence, restriction_sites)
    print("rest", eval_restriction_sites(sequence, restriction_sites))


def test_splice(sequence):
    # assert eval_splice_sites(sequence) is 1
    print("splice", eval_splice_sites(sequence))


def test_start(sequence):
    # assert eval_start_sites(sequence) is 1
    print("start", eval_start_sites(sequence))


ancestor_sequence = Seq.Seq(
    "TACGATATATGGGGATATATCGTATACGATATATACGATATATACGATATA",
    IUPAC.unambiguous_dna)
sequence = SequenceContainer(ancestor_sequence)
gc = GC_content[2]
test_gc(sequence, gc)
test_hair(sequence)
test_homo(sequence)
test_repeats(sequence)
rest = RestrictionEnzymes("NdeI XhoI HpaI PstI EcoRV NcoI BamHI".split())
test_rest(sequence, rest)
test_splice(sequence)
test_start(sequence)
sequence2 = SequenceContainer(
    Seq.Seq("TACGATATATGGGGAAATATCGTATACGATATA", IUPAC.unambiguous_dna))
test_host(sequence2, ancestor_sequence)
