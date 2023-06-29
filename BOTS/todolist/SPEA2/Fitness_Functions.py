import re
from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import GC
from todolist.SPEA2.Sequence_Container import SequenceContainer
from todolist.SPEA2.Bio_Structures import RibosomeBindingSites, splice_acceptors, splice_donors
import math


def eval_host(individual, ancestor_sequence):
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    return sum(1 for a, b in zip(sequence, ancestor_sequence) if a != b)


def eval_restriction_sites(individual, restrict_sites):
    """
    TODO: Make it remove rest sites
    """
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    # check unwanted restriction sites
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=sequence)
    result = analysis.full()
    # score the sequence based on the number of restriction sites
    score = 0
    for enz, cuts in result.items():
        for cut in cuts:
            score += 1
    return score


def eval_start_sites(individual, ribosome_binding_sites=RibosomeBindingSites, table_name="Standard"):
    """
    TODO: Make it remove start sites
    """
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(sequence))
    ]

    # None found
    if not len(start_codon_positions):
        return 0

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = sequence.tomutable()

    score = 0
    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for ii in range(2):
                    codon_idx = slice((codon_pos + ii) * 3, (codon_pos + ii + 1) * 3)
                    score += 1

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start: rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1
    return score


def eval_repeats(individual, window_size=10):
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    mutable_seq = sequence.tomutable()

    score = 0
    # iterate by codon, but map back to sequence-based indices
    for i in range(len(mutable_seq) // 3 - codon_window + 1):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
        )

        # check if the segment is found in the full sequence
        non_overlapping_matches = re.findall(
            str(mutable_seq[window]), str(mutable_seq)
        )

        if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
            # print(non_overlapping_matches)
            score += 1
    return score


def eval_homopolymers(individual, homopolymer_threshold=4):
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    seq = str(sequence.tomutable())

    score = 0

    # look at each (n_codons * 3)-mer
    idx = 0
    while idx < len(seq) - 4:
        current_count = 1
        current_letter = seq[idx]
        idx += 1
        while idx < len(seq) and seq[idx] is current_letter:
            current_count += 1
            if current_count > homopolymer_threshold:
                score += 1
            idx += 1
    return score


def eval_splice_sites(individual):
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")

    def _pass_back_matches(list_of_sites, curr_dna):
        dna = str(curr_dna)
        sites = set(m for expr in list_of_sites for m in re.finditer(expr, dna))
        try:
            sites.remove(None)
        except KeyError:
            pass
        # remove redundancy
        sites = set((site.span(), site[0]) for site in sites)
        codon_bounds = [
            (s[0][0] // 3, -(-s[0][1] // 3)) for s in sorted(sites, key=lambda x: x[0])
        ]
        return codon_bounds

    def _get_splice_sites(curr_dna):
        donor_sites = _pass_back_matches(splice_donors, curr_dna)
        acceptor_sites = _pass_back_matches(splice_acceptors, curr_dna)
        return set(donor_sites + acceptor_sites)

    return len(_get_splice_sites(sequence.tomutable))


def eval_gc_content(individual, gc):
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    window_size = gc.window_size  # tuples are immutable
    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(sequence))

    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(len(mutable_seq)):
        window = slice(
            i,
            (i + window_size)
            if (i + window_size) < len(mutable_seq)
            else len(mutable_seq),
        )
        gc_percent = GC(mutable_seq[window]) / 100

        if gc_percent > gc.high:
            score += (gc_percent - gc.high) * 100
        if window.stop is len(mutable_seq):
            break
    return round(score)


def eval_hairpins(individual, stem_length=10):
    """
    :param individual:
    :param stem_length:
    :return:
    """
    assert isinstance(individual, SequenceContainer)
    sequence = getattr(individual, "sequence")
    mutable_seq = sequence.tomutable()
    score = 0
    for i in range(0, len(mutable_seq) - stem_length, 3):
        stem_seq = mutable_seq[i: i + stem_length].toseq()
        # include wobble base pairing for G-[CT]
        hairpin_pattern = "".join(
            [nt if nt != "C" else "[CT]" for nt in stem_seq.reverse_complement()]
        )
        for hairpin in re.finditer(hairpin_pattern, str(mutable_seq)):
            if not (math.fabs(hairpin.start() - i) < stem_length or stem_length + 2 < hairpin.end() - hairpin.start()):
                score += 1
    return score
