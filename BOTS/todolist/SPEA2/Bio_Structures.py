import re
import json
from collections import namedtuple

import python_codon_tables as pct
from Bio import Entrez, Restriction


# Species-specifc data can be found on the Codon Usage Database
# (http://www.kazusa.or.jp/codon) using the NCBI Taxonomy database
# (http://www.ncbi.nlm.nih.gov/taxonomy) id (e.g. 413997)
# or the organism's Latin name (e.g. Escherichia coli B).
# Mapping species names to Taxonomy IDs can be done here:
_tax_id_url = "https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi"

# consensus donor seq is "GGGTRAGT"
# below is all possible versions with the first GT fixed, and at
# least 2 other NTs from the consensus seq
splice_donors = [
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
splice_acceptors = [
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


RibosomeBindingSites = {
    "rbs_0": "GGGGG",
    "rbs_1": "GGGGA",
    "rbs_2": "GGGAG",
    "rbs_3": "GGGAA",
    "rbs_4": "GGAGG",
    "rbs_5": "GGAGA",
    "rbs_6": "GGAAG",
    "rbs_7": "GGAAA",
    "rbs_8": "GAGGG",
    "rbs_9": "GAGGA",
    "rbs_10": "GAGAG",
    "rbs_11": "GAGAA",
    "rbs_12": "GAAGG",
    "rbs_13": "GAAGA",
    "rbs_14": "GAAAG",
    "rbs_15": "GAAAA",
    "rbs_16": "AGGGG",
    "rbs_17": "AGGGA",
    "rbs_18": "AGGAG",
    "rbs_19": "AGGAA",
    "rbs_20": "AGAGG",
    "rbs_21": "AGAGA",
    "rbs_22": "AGAAG",
    "rbs_23": "AGAAA",
    "rbs_24": "AAGGG",
    "rbs_25": "AAGGA",
    "rbs_26": "AAGAG",
    "rbs_27": "AAGAA",
    "rbs_28": "AAAGG",
    "rbs_29": "AAAGA",
    "rbs_30": "AAAAG",
    "rbs_31": "AAAAA",
}


def codon_tables(taxid, table_path=None):
    """Download the codon use table for the given species and return it as
    a dictionary.
    Returns:
        int: The NCBI taxonomy ID for the supplied species.
    Args:
        taxid (int): NCBI taxonomy ID for the desrired species.
        table_path (str): Defaults to None. Path to a JSON-formatted file representing the
            codon usage to consider. If None, the table is fetched from the internet.
    Raises:
        ValueError: If the NCBI taxonomy ID is not associated with a codon
        usage table, raise a ``ValueError`` informing the user and directing
        them to the NCBI Taxonomy Browser.
    Returns:
        dict{str, float}: A dictionary with codons as keys and the frequency
        that the codon is used to encode its amino acid as values.
    """
    if table_path is None:
        try:
            taxid = int(taxid)
        except ValueError as exc:
            taxid = _tax_id_from_species(taxid, exc)
        codon_table_by_aa = pct.download_codons_table(taxid)
    else:
        # load table from disk -- JSON format
        with open(table_path, "r") as table:
            codon_table_by_aa = json.load(table)

    return_dict = {}
    for _, codon_dict in codon_table_by_aa.items():
        for codon, frequency in codon_dict.items():
            return_dict[codon] = frequency
    if not return_dict:
        raise ValueError(
            '"{}" is not a valid host id. '.format(taxid)
            + "Supported hosts (Latin and NCBI taxonomy IDs) can be found at "
            + _tax_id_url
        )

    return return_dict


def _tax_id_from_species(species, exc=None):
    """Map the name of a species from a string to the NCBI taxonomy ID and
    return it.
    Args:
        species (str): Name of the species to map.
    Raises:
        ValueError: If the NCBI taxonomy ID cannot be determined, raise a ``ValueError``
        informing the user and directing them to the NCBI Taxonomy Browser.
    Returns:
        int: The NCBI taxonomy ID for the supplied species.
    """
    search_species = species.replace(" ", "+").replace("_", "+").strip()
    handle = Entrez.esearch(term=search_species, db="taxonomy")
    record = Entrez.read(handle)
    taxid = int(record["IdList"].pop())
    return taxid


class GCParams(namedtuple("GCParams", "name window_size low high")):
    """High and low values for GC-content within a specified window size.
    Attributes:
        name (str): Name of the parameter set.
        window_size (int): Number of nucleotides over which the GC content
            will be calculated.
        low (float): The minimum fraction of GC in the window.
        high (float): The maximum fraction of GC in the window.
    """


GC_content = [
    GCParams("IDT", 20, 0.15, 0.90),
    GCParams("twist", 50, 0.15, 0.80),
    GCParams("IDT_long", 100, 0.28, 0.68),
    GCParams("twist_long", "x3", 0.3, 0.65),
]


def RestrictionEnzymes(restriction_enzymes):
    """Create a RestrictionBatch instance to search for sites for a supplied
    list of restriction enzymes.
    Args:
        restriction_enzymes (list[str], optional): List of restriction
            enzymes to consider. Defaults to ["NdeI", "XhoI", "HpaI", "PstI",
            "EcoRV", "NcoI", "BamHI"].
    Returns:
        Bio.Restriction.Restriction.RestrictionBatch: RestrictionBatch instance
        configured with the input restriction enzymes.
    """
    return Restriction.RestrictionBatch(
        [Restriction.AllEnzymes.get(enz) for enz in restriction_enzymes]
    )
