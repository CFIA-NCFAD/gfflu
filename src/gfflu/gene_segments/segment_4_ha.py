import logging
import sys

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gfflu.peptides import Peptides

logger = logging.getLogger(__name__)


def naive_HA_signal_peptide_search(ha_aa: str) -> int:
    """Naive search for HA signal peptide in HA protein sequence

    The search is done by comparing the first X amino acids of the HA protein
    sequence to the X amino acids of the HA signal peptide sequences where X is
    the length of the HA signal peptide sequence.

    Returns length of best matching signal peptide if found, otherwise 0.
    """
    min_mismatch = sys.maxsize
    best_seq = ""
    best_seqid = None
    for seqid, seq in Peptides.ha_signal_peptides:
        mismatch = 0
        # stop at the end of the shorter AA sequences
        # the shorter sequence should be the HA signal peptide sequence
        stop = min(len(seq), len(ha_aa))
        for i in range(stop):
            if seq[i] != ha_aa[i]:
                mismatch += 1
        if mismatch < min_mismatch:
            min_mismatch = mismatch
            best_seq = seq
            best_seqid = seqid
    logger.info(f"HA signal peptide match: mismatch: {min_mismatch}")
    logger.info(f"HA signal peptide match: seqid: {best_seqid}")
    logger.info(f"HA signal peptide match: query : {best_seq}")
    logger.info(f"HA signal peptide match: target: {ha_aa[:len(best_seq)]}")
    return len(best_seq)


def find_HA_signal_peptide(seqrec: SeqRecord, ha_gene: SeqFeature) -> SeqFeature:
    """Find HA signal peptide in FASTA file for the HA gene segment (segment 4)"""
    ha_gene_start = ha_gene.location.start
    ha_gene_end = ha_gene.location.end
    # Assuming the HA gene is on the forward strand, the HA gene is extracted
    # from the FASTA file nucleotide sequence, and the HA gene is translated
    # into amino acids based on the expected HA gene start and end positions
    ha_aa = str(seqrec.seq[ha_gene_start:ha_gene_end].translate())
    # A naive search is performed to find the best matching signal peptide
    # sequence in the HA signal peptide database
    sigp_aa_len = naive_HA_signal_peptide_search(ha_aa)
    sigp_nt_len = sigp_aa_len * 3
    sigp_start = ha_gene_start
    sigp_end = ha_gene_start + sigp_nt_len
    logger.info(f"HA signal peptide: {sigp_start}-{sigp_end}")
    return SeqFeature(
        FeatureLocation(sigp_start, sigp_end),
        type="signal_peptide_region_of_CDS",
        strand=1,
        qualifiers={
            "ID": ["signal_peptide-HA"],
            "Parent": ["cds-HA"],
        },
    )
