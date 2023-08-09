import logging
import os
from pathlib import Path
from typing import Iterable, List, Optional, Union

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqIO.InsdcIO import _insdc_location_string
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
    for feature in rec.features:
        features = [feature]
        while features:
            next_features = []
            for current_feature in features:
                out.append(current_feature)
                if not hasattr(current_feature, "sub_features"):
                    continue
                new_features = []
                sub_features: Optional[List[SeqFeature]] = current_feature.sub_features
                if sub_features:
                    cds_subfeatures = [x for x in sub_features if x.type == "CDS"]
                    if len(cds_subfeatures) > 1:
                        subf = cds_subfeatures[0]
                        locations = [x.location for x in cds_subfeatures]
                        new_location = CompoundLocation(locations)
                        subf.location = new_location
                        subf.qualifiers["translation"] = [get_translation(rec.seq, new_location)]
                        new_features.append(subf)
                    else:
                        subf = cds_subfeatures[0]
                        subf.qualifiers["translation"] = [get_translation(rec.seq, subf.location)]
                        new_features.append(subf)
                    for subf in sub_features:
                        if subf.type == "signal_peptide_region_of_CDS":
                            subf.type = "sig_peptide"
                        elif subf.type == "mature_protein_region_of_CDS":
                            subf.type = "mat_peptide"
                        else:
                            continue
                        subf.qualifiers["translation"] = [get_translation(rec.seq, subf.location)]
                        try:
                            subf.qualifiers["gene"] = [feature.qualifiers["gene"][0]]
                        except KeyError:
                            pass
                        try:
                            subf.qualifiers["product"] = [feature.qualifiers["product"][0]]
                        except KeyError:
                            pass
                        new_features.append(subf)
                    next_features.extend(new_features)
            features = next_features
    rec.features = out
    return rec


def get_translation(seq: Seq, location: Union[FeatureLocation, CompoundLocation]) -> str:
    """Get translation of sequence from start to end"""
    translation = location.extract(seq).translate()
    if translation[-1] == "*":
        translation = translation[:-1]
    return translation


def write_aa_fasta(flat_rec: SeqRecord, prefix: str, out_file: os.PathLike) -> None:
    """Write out FASTA file with amino acid sequences"""
    with open(out_file, "w") as out_handle:
        for feature in flat_rec.features:
            if feature.type not in {"CDS", "mat_peptide", "sig_peptide"}:
                continue
            translation = get_translation(flat_rec.seq, feature.location)
            try:
                gene = feature.qualifiers["gene"][0]
            except KeyError:
                gene = ""
            try:
                product = feature.qualifiers["product"][0]
            except KeyError:
                product = ""
            location_str = _insdc_location_string(location=feature.location, rec_length=len(flat_rec.seq))
            out_handle.write(f">{prefix}|{feature.type}|{gene}|{location_str} {product}\n{translation}\n")


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
