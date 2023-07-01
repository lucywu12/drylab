import json
import python_codon_tables as pct
from . import logging
from ..data import _tax_id_from_species


logger = logging.getLogger(__name__)


def save_codon_table_to_disk(taxid, outfile):
    """Download a codon usage table and save it to disk for future use.

    Args:
        taxid (str): Name of the species to map.
        outfile (str): Full path to the file to write to store the codon
            usage information.
    """
    try:
        taxid = int(taxid)
    except ValueError as exc:
        taxid = _tax_id_from_species(taxid, exc)
    logger.info("Downloading host table for NCBI taxonomy ID {}...".format(taxid))
    codon_table_by_aa = pct.download_codons_table(taxid)
    if not outfile.endswith(".json"):
        logger.warning(
            'Output file "{}" does not have a JSON file extension.'.format(outfile)
        )
    with open(outfile, "w") as of:
        json.dump(codon_table_by_aa, of)
