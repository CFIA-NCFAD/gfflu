import logging
import os
from pathlib import Path
from typing import Iterable, List, Optional

from BCBio import GFF
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def check_gff(gff_iterator: Iterable[SeqRecord], molecule_type: str = "DNA"):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.

    source: https://github.com/chapmanb/bcbb/blob/81442c07173aaa59cbc1be33dd4a6e0eee87c5f6/gff/Scripts/gff/gff_to_genbank.py#L43
    """
    for rec in gff_iterator:
        if "molecule_type" not in rec.annotations:
            rec.annotations["molecule_type"] = molecule_type
        yield flatten_features(rec)


def flatten_features(rec: SeqRecord) -> SeqRecord:
    """Make sub_features in an input rec flat for output.

    GenBank does not handle nested features, so we want to make
    everything top level.

    source: https://github.com/chapmanb/bcbb/blob/81442c07173aaa59cbc1be33dd4a6e0eee87c5f6/gff/Scripts/gff/gff_to_genbank.py#L52
    """
    out = []
    for f in rec.features:
        cur = [f]
        while cur:
            nextf = []
            for curf in cur:
                out.append(curf)
                if not hasattr(curf, "sub_features"):
                    continue
                sub_features: Optional[List[SeqFeature]] = curf.sub_features
                if sub_features:
                    for subf in sub_features:
                        if subf.type == "CDS":
                            translation = rec.seq[subf.location.start : subf.location.end].translate()
                            if translation[-1] == "*":
                                translation = translation[:-1]
                            subf.qualifiers["translation"] = [translation]
                        elif subf.type == "signal_peptide_region_of_CDS":
                            subf.type = "sig_peptide"
                        elif subf.type == "mature_protein_region_of_CDS":
                            subf.type = "mat_peptide"
                    nextf.extend(sub_features)
            cur = nextf
    rec.features = out
    return rec


def read_gff(gff: Path) -> Optional[SeqRecord]:
    recs = list(GFF.parse(gff))
    if not recs:
        logger.warning(f"Empty GFF file: {gff}")
        return None
    if len(recs) > 1:
        logger.warning(f"More than one record ({len(recs)}) in GFF file '{gff}'")
    return recs[0]


def write_gff(gff_recs: List[SeqRecord], out_file: os.PathLike) -> None:
    """Write out GFF file"""
    with open(out_file, "w") as out_handle:
        GFF.write(gff_recs, out_handle)


def read_fasta(fasta: Path) -> SeqRecord:
    # Assuming there should be at least 1 sequence in the FASTA file
    recs = list(SeqIO.parse(fasta, "fasta"))
    if not recs:
        msg = f"Empty FASTA file: {fasta}"
        raise ValueError(msg)
    if len(recs) > 1:
        logger.warning(f"More than one record ({len(recs)}) in FASTA file '{fasta}'")
    return recs[0]
