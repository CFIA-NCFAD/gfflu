import logging
from pathlib import Path
from typing import Optional

import polars as pl
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def find_PB1_F2(fasta: Path, blastx: Path) -> Optional[SeqFeature]:
    """Find PB1-F2 protein in FASTA file using BLASTX results

    Tabular BLASTX results are parsed for PB1-F2 entries. If found, the first entry is used to determine the start
    position of the PB1-F2 protein in the FASTA file. The protein sequence is translated and the first stop codon is
    used to determine the length of the protein. The end position is calculated from the start position and the length.
    The start and end positions are used to create a SeqFeature object.

    Parameters
    ----------
    fasta : Path
        Path to FASTA file
    blastx : Path
        Path to BLASTX results file

    Returns
    -------
    Optional[SeqFeature]
        SeqFeature object for PB1-F2 protein or None if no PB1-F2 protein was found in the BLASTX results
    """
    df = pl.read_csv(
        blastx,
        separator="\t",
        comment_char="#",
        has_header=False,
        new_columns=[
            "query",
            "subject",
            "pct_id",
            "aln_len",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )
    df_pb1_f2 = df.filter(df["subject"].str.starts_with("PB1-F2")).sort("bitscore", descending=True)
    if df_pb1_f2.shape[0] == 0:
        return None
    pb1_f2_start = df_pb1_f2["qstart"][0]
    pb1_f2_start = pb1_f2_start - 1
    pb1_f2_subjectid = df_pb1_f2["subject"][0]
    seqrec: SeqRecord = next(iter(SeqIO.parse(fasta, "fasta")))
    seq_from_pb1_f2_start = seqrec.seq[pb1_f2_start:]
    pb1_f2_aa = str(seq_from_pb1_f2_start.translate())
    # find the first stop codon and add 1 to get the length of the PB1-F2 protein including stop codon
    pb1_f2_aa_len = next((i for i, x in enumerate(pb1_f2_aa) if x == "*"), 0) + 1
    if pb1_f2_aa_len <= 1:
        logger.warning(
            f"PB1-F2 protein length is {pb1_f2_aa_len}. Found no stop codon. PB1-F2 BLASTX results: {df_pb1_f2}"
        )
        return None
    pb1_f2_nt_len = pb1_f2_aa_len * 3
    pb1_f2_end = pb1_f2_start + pb1_f2_nt_len
    logger.info(f"PB1-F2: {pb1_f2_start}-{pb1_f2_end}")
    feature = SeqFeature(
        FeatureLocation(pb1_f2_start, pb1_f2_end),
        type="gene",
        strand=1,
        qualifiers={
            "ID": ["gene-PB1-F2"],
            "Target": [pb1_f2_subjectid],
            "gene": ["PB1-F2"],
            "gene_biotype": ["protein_coding"],
        },
    )
    feature.sub_features = [
        SeqFeature(
            FeatureLocation(pb1_f2_start, pb1_f2_end),
            type="CDS",
            strand=1,
            qualifiers={
                "ID": ["cds-PB1-F2"],
                "Parent": ["gene-PB1-F2"],
                "Target": [pb1_f2_subjectid],
                "gene": ["PB1-F2"],
                "product": ["PB1-F2 protein"],
            },
        ),
    ]
    return feature
