import logging
import random
import re

from itertools import product

from Bio.Alphabet import IUPAC

from Bio.Data import CodonTable
from Bio.Restriction import Analysis
from Bio.SeqUtils import CodonUsage, GC, seq3

from . import Seq, MutableSeq, codon_use

logger = logging.getLogger(__name__)


def mutate_codon(codon_in, codon_use_table):
    """Select a synonymous codon in accordance with the frequency of use
    in the host organism.

    Args:
        codon_in (Bio.Seq.Seq): A single codon.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.

    Returns:
        Bio.Seq.Seq: A new codon.
    """
    AA = seq3(CodonTable.standard_dna_table.forward_table[str(codon_in)]).upper()

    synonymous_codons, codon_use_freq = codon_use_table[AA]
    if len(synonymous_codons) == 1:
        return codon_in

    # pick new codon
    codon_out = codon_in
    while codon_in == codon_out:
        codon_out = random.choices(synonymous_codons, codon_use_freq).pop()

    logger.detail("Mutating {} codon from {} to {}".format(AA, codon_in, codon_out))

    return codon_out


def resample_codons(dna_sequence, codon_use_table):
    """Generate a new DNA sequence by swapping synonymous codons.
    Codons are selected in accordance with their frequency of occurrence in
    the host organism.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    resampled_dna = "".join(
        [
            random.choices(*codon_use_table[seq3(AA).upper()]).pop()
            for AA in dna_sequence.translate()
        ]
    )

    return Seq(resampled_dna, IUPAC.unambiguous_dna)


def compare_profiles(codons_count, host_profile, relax):
    """Compute the deviation from the expected codon usage based on a host
    codon usage profile.

    Note:
        The `relax` parameter uniformly increases the host codon usage that
        is used to estimate the number of times each codon should appear in
        the sequence. These values are rounded and then iteratively adjusted
        to be consistent with the length of the sequence of interest.
        Increasing this parameter further distorts the codon use distribution
        from the host.

    Args:
        codons_count (dict{str, int}): A dictionary with each codon as
            keys and the number of times it appears in a gene as values.
        host_profile (dict{str, foat}): A dictionary with each codon as keys
            and the frequency of its use in the host organism as values.
        relax (float): The maximum deviation from the host profile to tolerate.

    Returns:
        dict{str, dict{str, int}}, float: A dictionary with each codon as keys,
        and dictionaries of the difference between the observed and expected
        codon usage.

        The number of mutations per residue that are needed to make
        the sequence match the host codon usage.
    """
    logger.info("Comparing codon usage profiles")
    table = {}
    # loop AAs
    for AA, synonymous_codons in CodonUsage.SynonymousCodons.items():
        logger.detail("Codon use for " + AA)
        temp_table = {}

        # calculate total usage of codon in input
        tot_usage = int(sum([codons_count[codon] for codon in synonymous_codons]))

        # calculate ideal usage of codon in host
        tot_ideal = 0
        for codon in synonymous_codons:
            ideal_usage_abs = int(round(host_profile[codon] * tot_usage, 0))
            ideal_usage = int(round(host_profile[codon] * relax * tot_usage, 0))
            logger.detail("{}: {}".format(codon, ideal_usage))
            tot_ideal += ideal_usage
            temp_table[codon] = {
                "input_count": codons_count[codon],
                "input_perc": codons_count[codon],
                "ideal_usage_abs": ideal_usage_abs,
                "ideal_usage": ideal_usage,
                "host_perc": host_profile[codon],
            }

        # account for rounding issues and relaxation of the host profile
        # by adjusting the most- and least-used codons as necessary

        # if the calculated usage exceeds the actual usage, subtract one from
        # the "ideal_use" of the least-used host codon that appears in the sequence
        while tot_ideal > tot_usage:
            codon = min(
                [c for c, d in temp_table.items() if d["ideal_usage"] > 0],
                key=lambda codon: temp_table[codon]["host_perc"],
            )
            temp_table[codon]["ideal_usage"] -= 1
            tot_ideal -= 1

        # if the calculated usage is less than the actual usage, add one to
        # the "ideal_use" of the most-used host codon
        while tot_ideal < tot_usage:
            codon = max(temp_table, key=lambda codon: temp_table[codon]["host_perc"])
            temp_table[codon]["ideal_usage"] += 1
            tot_ideal += 1

        # populate return table
        for codon, usage in temp_table.items():
            table[codon] = {
                "input_count": usage["input_count"],
                "ideal_usage_abs": usage["ideal_usage_abs"],
                "difference": usage["input_count"] - usage["ideal_usage"],
                "difference_abs": usage["input_count"] - usage["ideal_usage_abs"],
            }

    # calculate difference value
    number_of_residues, diff_total = 0, 0
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        for codon in synonymous_codons:
            number_of_residues += table[codon]["ideal_usage_abs"]
            diff_total += abs(table[codon]["difference_abs"])

    diff_total /= 2  # convert to the number of mutations needed
    diff = diff_total / number_of_residues
    logger.info(
        "Current number of deviations from host codon usage profile: "
        + "{} over {} codons".format(int(diff_total), number_of_residues)
    )
    logger.info("Current deviation from host codon usage profile: {:.0%}".format(diff))

    return table, diff


def harmonize_codon_use_with_host(dna_sequence, mutation_profile):
    """Adjust the codon usage in the DNA sequence to be consistent with
    the host profile.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        mutation_profile (dict{str, dict{str, int}}): A dictionary
            with each codon as keys, and dictionaries of the difference
            between the observed and expected codon usage.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("Harmonizing codon useage with host codon usage profile")

    mutable_seq = dna_sequence.tomutable()
    for _, synonymous_codons in CodonUsage.SynonymousCodons.items():
        # get the index of relevant codons in the sequence
        mutation_table = {}
        for codon in synonymous_codons:
            mutation_table[codon] = {
                "difference": mutation_profile[codon]["difference"]
            }

            pos_list = []
            for i in range(0, len(mutable_seq), 3):
                codon_idx = slice(i, i + 3)
                if mutable_seq[codon_idx] == codon:
                    pos_list.append(codon_idx)

            mutation_table[codon]["pos"] = pos_list

        # check if this AA even needs to be adjsuted
        tot_diff = sum(
            abs(mutation_table[codon]["difference"]) for codon in synonymous_codons
        )

        while tot_diff > 0:
            # randomly select a pair of codons to adjust
            codon_to_remove = random.choice(
                [c for c, d in mutation_table.items() if d["difference"] > 0]
            )
            codon_to_add = random.choice(
                [c for c, d in mutation_table.items() if d["difference"] < 0]
            )

            # randomly select the position to update
            codon_idx = random.choice(mutation_table[codon_to_remove]["pos"])

            # remove from sequence
            mutation_table[codon_to_remove]["pos"].remove(codon_idx)
            mutation_table[codon_to_remove]["difference"] -= 1

            # add to sequence
            mutation_table[codon_to_add]["pos"].append(codon_idx)
            mutation_table[codon_to_add]["difference"] += 1

            mutable_seq[codon_idx] = codon_to_add

            # update difference
            tot_diff = sum(
                abs(mutation_table[codon]["difference"]) for codon in synonymous_codons
            )
        # mutation_table now has difference = 0 for all codons

    return mutable_seq.toseq()


def resample_codons_and_enforce_host_profile(
    dna_sequence, codon_use_table, host_profile, relax
):
    """Generate a new DNA sequence by swapping synonymous codons.
    Codons are selected in accordance with their frequency of occurrence in
    the host organism and adjust the codon usage in the DNA sequence to
    match the host profile.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymous codon is used.
        host_profile (dict{str, foat}): A dictionary with each codon as keys
            and the frequency of its use in the host organism as values.
        relax (float): The maximum deviation from the host profile to tolerate.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("Relax coefficient: {}".format(relax))
    dna_sequence = resample_codons(dna_sequence, codon_use_table)

    # measure the deviation from the host profile and adjust accordingly
    mutation_table, _ = compare_profiles(
        codon_use.count_codons(dna_sequence), host_profile, relax
    )

    return harmonize_codon_use_with_host(dna_sequence, mutation_table)


def gc_scan(dna_sequence, codon_use_table, gc):
    """Scan across a sequence and replace codons to acheive a desired GC
    content within the window.

    Note:
        The following fields of the `GCParams` type are used in this
        function:

        * **window_size** (`int`) – Size of sliding window (in nucelotides) to
          examine for GC content. Window sizes can also be expressed as
          factors of the length of `dna_sequence` by passing a string
          that begins with "x" (e.g. "x0.5").

        * **low** (`float`) – Minimum GC content in window.

        * **high** (`float`) – Maximum GC content in window.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        gc (GCParams): A `namedtuple` with fields for name, window_size,
            minimum and maximum GC content.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info(
        "GC content scan -- window_size: {} nucleotides, threshold: {} < x < {}".format(
            gc.window_size, gc.low, gc.high
        )
    )

    window_size = gc.window_size  # tuples are immutable
    # some windows may be expressed as function of the sequence length
    if isinstance(window_size, str) and window_size.startswith("x"):
        window_size = int(float(window_size[1:]) * len(dna_sequence))

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3 + 1
    mutable_seq = dna_sequence.tomutable()

    # iterate by codon, but map back to sequence-based indices
    for i in range(len(mutable_seq) // 3):
        window = slice(
            i * 3,
            (i + codon_window) * 3
            if (i + codon_window) * 3 < len(mutable_seq)
            else len(mutable_seq),
        )
        logger.debug("Current segment: {}".format(mutable_seq[window]))

        gc_percent = GC(mutable_seq[window]) / 100
        count = 0  # counter to prevent infinite loop
        # check gc_percent of current segment
        while (
            gc_percent < gc.low or gc_percent > gc.high
        ) and count < codon_window * 2:
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice((i * 3) + position, ((i + 1) * 3) + position)

            init_codon = mutable_seq[codon_idx]
            new_codon = mutate_codon(init_codon, codon_use_table)

            if (GC(new_codon) < GC(init_codon) and gc_percent > gc.high) or (
                GC(new_codon) > GC(init_codon) and gc_percent < gc.low
            ):
                mutable_seq[codon_idx] = new_codon
                logger.debug("Mutating position: {}".format(position))
                gc_percent = GC(mutable_seq[window]) / 100

            count += 1

    return mutable_seq.toseq()


def remove_restriction_sites(dna_sequence, codon_use_table, restrict_sites):
    """Identify and remove seuences recognized by a set of restriction
    enzymes.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        restrict_sites (Bio.Restriction.RestrictionBatch): RestrictionBatch
            instance configured with the input restriction enzymes.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """

    logger.info("Removing restriction sites")

    # check each unwanted restriction site
    analysis = Analysis(restrictionbatch=restrict_sites, sequence=dna_sequence)
    result = analysis.full()

    mutable_seq = dna_sequence.tomutable()
    for enz, cuts in result.items():
        for cut in cuts:
            logger.info(
                "Restriction enzyme ({}) cut site detected at position {}.".format(
                    str(enz), cuts
                )
            )
            # map sequence position to codon position
            # subtract 1 from `cut` to convert from sequence to string indexing
            codon_pos, offset = divmod((cut - 1) - (enz.size // 2), 3)

            # ensure the whole codon we mutate is being recognized by the restriction enzyme
            if offset:
                codon_pos += 1
            codon_idx = slice(codon_pos * 3, (codon_pos + 1) * 3)

            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

    return mutable_seq.toseq()


def remove_start_sites(
    dna_sequence, codon_use_table, ribosome_binding_sites, table_name="Standard"
):
    """Identify and remove alternate start sites using a supplied set of
    ribosome binding sites and a codon table name.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        ribosome_binding_sites (dict{str, str}): A dictionary with named
            ribosome binding sites as keys and the corresponding sequences
            as values.
        table_name (str, optional): Name of a registered NCBI table. See
            `Bio.Data.CodonTable.unambiguous_dna_by_name.keys()` for
            options. Defaults to "Standard".

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    codon_table = CodonTable.unambiguous_dna_by_name[table_name]
    logger.info(
        "Removing alternate start sites: {}".format(", ".join(codon_table.start_codons))
    )

    # find all start codon sites (xTG)
    start_codon_positions = [
        m.start()
        for start_codon in codon_table.start_codons
        for m in re.finditer(start_codon, str(dna_sequence))
    ]

    if not len(start_codon_positions):
        logger.info("No start codons found in sequence")
        return dna_sequence

    logger.info(
        "Found {} start codons. Checking for upstream RBSs...".format(
            len(start_codon_positions)
        )
    )

    # check each start site for RBS
    # 18 base pairs upstream of each xTG, ignore 3 bp closest to xTG
    _rbs_offset = 18
    rbs_positions = [
        pos - _rbs_offset for pos in start_codon_positions if pos >= _rbs_offset
    ]
    mutable_seq = dna_sequence.tomutable()

    for rbs_start in rbs_positions:
        # ignore 3 bp closest to xTG
        rbs_stop = rbs_start + _rbs_offset - 3
        rbs_query_seq = str(mutable_seq[rbs_start:rbs_stop])

        logger.detail(
            'checking for sequence "{}" in "{}"'.format(
                rbs_query_seq, mutable_seq[rbs_stop : rbs_stop + 6]
            )
        )

        # check each unwanted RBS in each potential fragment
        for rbs, site in ribosome_binding_sites.items():
            logger.detail('checking for start site "{}" in "{}"'.format(rbs, site))
            search = rbs_query_seq.find(site)

            count = 0  # counter to prevent infinite loop
            while search != -1 and count < 10:
                # mutate residues if site is found
                codon_pos = (search + rbs_start) // 3
                for i in range(2):
                    codon_idx = slice((codon_pos + i) * 3, (codon_pos + i + 1) * 3)
                    new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
                    mutable_seq[codon_idx] = new_codon

                # reset sequence and search again
                rbs_query_seq = str(mutable_seq[rbs_start : rbs_stop + 3])
                search = rbs_query_seq.find(site)
                count += 1

    return mutable_seq.toseq()


def remove_repeating_sequences(dna_sequence, codon_use_table, window_size):
    """Identify and remove repeating sequences of codons or groups of
    codons within a DNA sequence.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        window_size (int): Size the window (in nucleotides) to examine.
            Window sizes are adjusted down to the nearest multiple of 3 so
            windows only contain complete codons.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("Scanning for repeated stretches of {} nucleotides".format(window_size))

    def _mutate_and_keep_looping(mutable_seq, window, offset):
        num_mutations = random.randint(1, 2)
        logger.debug("Mutating {} codons".format(num_mutations))
        for _ in range(num_mutations):
            position = random.randrange(0, len(mutable_seq[window]), 3)
            codon_idx = slice(offset + position, (offset + 3) + position)
            new_codon = mutate_codon(mutable_seq[codon_idx], codon_use_table)
            mutable_seq[codon_idx] = new_codon

        return True

    # iterate across overlapping chunks of complete codons
    codon_window = window_size // 3
    mutable_seq = dna_sequence.tomutable()

    current_cycle = 0  # prevent infinite loops (caused by poly-TRP or poly-MET)
    keep_looping = True
    # `keep_looping` if any mutation is made,
    # i.e. continue until both checks pass without mutations
    while keep_looping and (current_cycle < (codon_window * 10)):

        keep_looping = False

        # iterate by codon, but map back to sequence-based indices
        for i in range(len(mutable_seq) // 3):
            window = slice(
                i * 3,
                (i + codon_window) * 3
                if (i + codon_window) * 3 < len(mutable_seq)
                else len(mutable_seq),
            )

            # make each mutable codon immutable so it can be hashed later
            codons = [
                str(mutable_seq[window][i : i + 3])
                for i in range(0, len(mutable_seq[window]), 3)
            ]

            # check if all codons in the window are identical
            if len(set(codons)) == 1:
                logger.detail("All codons in window are identical: {}".format(codons))
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

            # check if the segment is found in the full sequence
            non_overlapping_matches = re.findall(
                str(mutable_seq[window]), str(mutable_seq)
            )
            if len(non_overlapping_matches) > 1 and len(mutable_seq[window]) > 3:
                logger.debug("Current window is repeated in the sequence")
                keep_looping = _mutate_and_keep_looping(mutable_seq, window, (i * 3))

        current_cycle += 1

    return mutable_seq.toseq()


def remove_local_homopolymers(
    dna_sequence, codon_use_table, n_codons=2, homopolymer_threshold=4
):
    """Identify and remove consecutive stretches of the same nucleotides
    using a sliding window of a fixed number of codons.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        n_codons (int, optional): Size of window (in codons) to examine.
            Defaults to 2.
        homopolymer_threshold (int): number of consecutive nucleotide
            repeats allowed. Defaults to 4.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    logger.info("Detecting and removing local homopolymers")
    mutable_seq = dna_sequence.tomutable()

    # look at each (n_codons * 3)-mer
    keep_looping = True
    while keep_looping:
        for i in range(0, len(mutable_seq), 3):
            window = slice(
                i,
                i + (n_codons * 3)
                if i + (n_codons * 3) < len(mutable_seq)
                else len(mutable_seq),
            )

            seq = str(mutable_seq[window])
            nt_counts = {letter: seq.count(letter) for letter in set(seq)}
            letter = max(nt_counts, key=lambda letter: nt_counts[letter])

            if nt_counts[letter] <= homopolymer_threshold:
                keep_looping = False
                continue

            logger.detail("Homopolymer ({}) detected at position {}".format(seq, i))
            logger.detail("{}, count={}".format(letter, nt_counts[letter]))

            for j in range(n_codons):
                codon_idx = slice(i + (j * 3), i + ((j + 1) * 3))
                mutable_seq[codon_idx] = mutate_codon(
                    mutable_seq[codon_idx], codon_use_table
                )
            keep_looping = True

    return mutable_seq.toseq()


def remove_hairpins(dna_sequence, codon_use_table, stem_length=10):
    """Identify and remove stretches of the equence that can form hairpins.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.
        stem_length (int, optional): Length of hairpin stem to detect.
            Defaults to 10.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """
    mutable_seq = dna_sequence.tomutable()
    for i in range(0, len(mutable_seq), stem_length):
        stem_seq = mutable_seq[i : i + stem_length].toseq()
        # include wobble base pairing for G-[CT]
        hairpin_pattern = "".join(
            [nt if nt != "C" else "[CT]" for nt in stem_seq.reverse_complement()]
        )
        for hairpin in re.finditer(hairpin_pattern, str(mutable_seq)):
            # floor start, ceil end
            pos = random.randint(hairpin.start() // 3, -(-hairpin.end() // 3))

            # don't run off the end of the sequence
            if (pos + 1) * 3 > len(mutable_seq):
                # account for the ceil AND shift back to the last complete codon
                pos -= 2

            codon_idx = slice(pos * 3, (pos + 1) * 3)
            mutable_seq[codon_idx] = mutate_codon(
                mutable_seq[codon_idx], codon_use_table
            )

    return mutable_seq.toseq()


# consensus donor seq is "GGGTRAGT"
# below is all possible versions with the first GT fixed, and at
# least 2 other NTs from the consensus seq
_splice_donors = [
    re.compile(r"GGGT\wAGT", re.UNICODE),
    re.compile(r"\wGGT\w[AT]GT", re.UNICODE),
    re.compile(r"G\wGT\wAGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]AGT", re.UNICODE),
    re.compile(r"GGGT[AG]\w[AG]T", re.UNICODE),
    re.compile(r"GGGT[AG]\wG\w", re.UNICODE),
    re.compile(r"GGGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"GGGT[AG]AG\w", re.UNICODE),
    re.compile(r"GGGT[AG]A[AG]T", re.UNICODE),
    re.compile(r"GGGT[AG]\wGT", re.UNICODE),
    # re.compile(r"\wGGT[AG]A[AG]\w", re.UNICODE), # redundant with below
    re.compile(r"\wGGT[AG][AT][ATG]\w", re.UNICODE),
    re.compile(r"\wGGT[AG]\wGT", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]\w", re.UNICODE),
    re.compile(r"G\wGT[AG]AG\w", re.UNICODE),
    re.compile(r"G\wGT[AG]A[AG]T", re.UNICODE),
    re.compile(r"G\wGT[AG]\wGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]AG\w", re.UNICODE),
    re.compile(r"\w\wGT[AG]\wGT", re.UNICODE),
    re.compile(r"\w\wGT[AG]\wG\w", re.UNICODE),
]

# consensus branch seq is "YTRAC"
# ignore branch points (for now) because they are small
# and occur 20-50 NTs upstream of acceptor -- not specific enough

# consensus acceptor seq is "YYYYYNCAGG"
# below are all sequences ending in NCAGG, NNAGG and NCAGN
# where at least 3 of the 5 upstream NTs are pyrimidines (Y, [TC])
_splice_acceptors = [
    re.compile(r"[TC][TC][TC]\w\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC][TC]\w[TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w[TC][TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC][TC][TC]\w[ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC][TC]\w[TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w[TC]\w[TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"\w\w[TC][TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w[TC]\w[TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC]\w\w[TC][TC][ATCG]CAG\w", re.UNICODE),
    re.compile(r"[TC][TC]\w\w[TC][ATCG]CAG\w", re.UNICODE),
]


def remove_splice_sites(dna_sequence, codon_use_table):
    """Identify and remove RNA splice sites within a DNA sequence.

    Args:
        dna_sequence (Bio.Seq.Seq): A read-only representation of
            the DNA sequence.
        codon_use_table (dict{str, list[list, list]}): A dictionary with
            each amino acid three-letter code as keys, and a list of two
            lists as values. The first list is the synonymous codons that
            encode the amino acid, the second is the frequency with which
            each synonymouscodon is used.

    Returns:
        Bio.Seq.Seq: A read-only representation of the new DNA sequence.
    """

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
        donor_sites = _pass_back_matches(_splice_donors, curr_dna)
        acceptor_sites = _pass_back_matches(_splice_acceptors, curr_dna)
        return set(donor_sites + acceptor_sites)

    mutable_seq = dna_sequence.tomutable()

    keep_looping = True
    n_times_through_unchanged = 0
    prev_seq = ""
    while keep_looping:
        # look for donor and acceptor seqs
        splice_sites = _get_splice_sites(mutable_seq)

        logger.info("Removing RNA splice site donors and acceptors.")
        for site_bounds in splice_sites:
            site_removed = False

            codon_list = list(range(site_bounds[0], site_bounds[-1] + 1))
            random.shuffle(codon_list)

            indices = [slice(cdn * 3, (cdn + 1) * 3) for cdn in codon_list]
            init_codons = [mutable_seq[idx] for idx in indices if str(mutable_seq[idx])]

            AAs = [
                seq3(CodonTable.standard_dna_table.forward_table[str(init)]).upper()
                for init in init_codons
            ]
            synonymous_codons = [codon_use_table[AA][0] for AA in AAs]
            substitutions = list(product(*synonymous_codons))
            random.shuffle(substitutions)

            for substitution in substitutions:
                for idx, codon in zip(indices, substitution):
                    mutable_seq[idx] = codon
                curr_sites = _get_splice_sites(mutable_seq)
                if site_bounds in curr_sites or len(curr_sites) >= len(splice_sites):
                    # put the starting sequence back
                    for idx, init_codon in zip(indices, init_codons):
                        mutable_seq[idx] = init_codon
                else:
                    site_removed = True
                    logger.detail("Removed site ({})!".format(site_bounds))
                    logger.detail(
                        "Remaining sites:\n" + "\n".join([str(s) for s in curr_sites])
                    )
                    break

            if site_removed:
                break

        remaining_sites = _get_splice_sites(mutable_seq)
        n_times_through_unchanged += int(str(mutable_seq) == prev_seq)
        prev_seq = str(mutable_seq)

        if not len(remaining_sites) or n_times_through_unchanged == 5:
            keep_looping = False

    return mutable_seq.toseq()
