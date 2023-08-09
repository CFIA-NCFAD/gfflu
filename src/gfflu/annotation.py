import logging
from pathlib import Path
from typing import List, Mapping, Optional, Tuple

from _operator import itemgetter
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gfflu.blastx import run_blastx
from gfflu.gene_segments.segment_2_pb1 import find_PB1_F2
from gfflu.gene_segments.segment_4_ha import find_HA_signal_peptide
from gfflu.io import read_fasta, read_gff
from gfflu.miniprot import (
    get_feature_product,
    get_feature_qualifier_value,
    get_feature_type,
    get_features_by_gene,
    get_gene,
    run_miniprot,
)

logger = logging.getLogger(__name__)


def run_annotation(fasta: Path, outdir: Path, prefix: str) -> SeqRecord:
    seqrec = read_fasta(fasta)
    if len(seqrec) == 0:
        return seqrec
    miniprot_gff = run_miniprot(outdir, fasta, prefix)
    logger.info(f"Reading {miniprot_gff}")
    gffrec = read_gff(miniprot_gff)
    if gffrec is None or len(gffrec) == 0:
        return seqrec
    blastx_b6 = run_blastx(outdir, fasta, prefix)
    gffrec.seq = next(iter(SeqIO.parse(fasta, "fasta"))).seq
    # remove annotations so that SnpEff doesn't have any issues with the GFF file
    gffrec.annotations = {}
    gene_features = get_features_by_gene(gffrec)
    top_features = select_top_scoring_features_by_gene(gene_features)
    reannotate_rec(gffrec, top_features)
    # HA signal peptide is too small for Miniprot to detect
    if any(gene == "HA" for (gene, ftype, product) in gene_features):
        # get CDS feature for HA
        cds_feature = next(x for x in top_features for y in x.sub_features if y.type == "CDS")
        cds_feature.sub_features.append(find_HA_signal_peptide(seqrec, cds_feature))
    # merge sub_features by gene
    top_features = merge_subfeatures_by_gene(top_features)
    # PB1-F2 may be too divergent for Miniprot to detect (< 50% AA identity for included sequences)
    # TODO: supplement data/iav-annotation.faa with PB1-F2 from other strains (PK: 2023-06-05)
    if any(gene == "PB1" for (gene, ftype, product) in gene_features):
        pb1_f2_feature: Optional[SeqFeature] = find_PB1_F2(fasta, blastx_b6)
        if pb1_f2_feature:
            top_features.append(pb1_f2_feature)
    # sort sub_features by start position
    for f in top_features:
        f.sub_features.sort(key=lambda x: x.location.start)
    # sort features by start position
    top_features.sort(key=lambda x: x.location.start)
    logger.debug(f"top_features: {top_features}")
    gffrec.features = top_features
    return gffrec


def merge_subfeatures_by_gene(features: List[SeqFeature]) -> List[SeqFeature]:
    """Merge subfeatures by gene

    For example, if there are two CDS features for PB1, merge them into one CDS feature.
    """
    gene_features = {}
    f: SeqFeature
    # add CDS features to gene_features first
    logging.info(f"Processing {len(features)} features")
    for f in features:
        if f.sub_features[0].type == "CDS":
            logging.debug(f"Processing CDS {f.sub_features}")
            gene = get_gene(get_feature_qualifier_value(f, "Target"))
            if gene not in gene_features:
                gene_features[gene] = f
            else:
                logging.warning(f"Duplicate gene feature for {gene}: {f}")
    logging.debug(f"gene_features: {gene_features}")
    # add other features to gene_features
    for f in features:
        if f.sub_features[0].type in {"CDS", "gene"}:
            continue
        logging.debug(f"Processing non-CDS {f.sub_features}")
        gene = get_gene(get_feature_qualifier_value(f, "Target"))
        logging.debug(f"Adding {f} to {gene}")
        if gene not in gene_features:
            gene_features[gene] = f
        else:
            logging.debug(f"Adding {f.sub_features} to {gene_features[gene].sub_features}")
            gene_features[gene].sub_features += f.sub_features
    return list(gene_features.values())


def add_qualifiers(feature: SeqFeature, gene: str, ftype: str, product: str) -> None:
    """Add qualifiers to feature

    Need to add qualifiers to features and subfeatures.

    For gene type features:
      - ID: should be "gene-" prefix plus gene name (e.g. gene-NS2)
      - Name: should be gene name (e.g. NS2)
      - gene: should be gene name (e.g. NS2)
      - gene_biotype: should be "protein_coding"

    For CDS type features:
      - ID: should be "cds" prefix plus gene name (e.g. cds-NS2)
      - Parent: should be "gene" prefix plus gene name (e.g. gene-NS2)
      - gene: should be gene name (e.g. NS2)
      - product: should match what is typically used for the gene, e.g. for gene NS2, the product should be
        "nonstructural protein 2"

    """
    if feature.type == "gene":
        feature.qualifiers["ID"] = [f"gene-{gene}"]
        feature.qualifiers["Name"] = [gene]
        feature.qualifiers["gene"] = [gene]
        feature.qualifiers["gene_biotype"] = ["protein_coding"]
    elif ftype == "CDS" and feature.type == "CDS":
        feature.qualifiers["ID"] = [f"cds-{gene}"]
        feature.qualifiers["Parent"] = [f"gene-{gene}"]
        feature.qualifiers["gene"] = [gene]
        feature.qualifiers["product"] = [product]
    elif ftype == "signal_peptide_region_of_CDS":
        feature.type = "signal_peptide_region_of_CDS"
        feature.qualifiers["ID"] = [f"signal_peptide-{gene}"]
        feature.qualifiers["Parent"] = [f"cds-{gene}"]
    elif ftype == "mature_protein_region_of_CDS":
        feature.type = "mature_protein_region_of_CDS"
        feature.qualifiers["ID"] = [f"mature_protein-{gene}"]
        feature.qualifiers["Parent"] = [f"cds-{gene}"]
        feature.qualifiers["product"] = [product]


def reannotate_rec(rec: SeqRecord, top_features: List[SeqFeature]) -> None:
    """Reannotate SeqRecord with top features

    This function is used to reannotate a SeqRecord with the top features,
    i.e. the highest scoring features from a Miniprot GFF file.

    Miniprot outputs stop_codon subfeatures, which are removed from the final
    GFF file and the CDS length is adjusted accordingly.

    """
    rec.features = top_features
    for f in rec.features:
        target = get_feature_qualifier_value(f, "Target")
        gene = get_gene(target)
        ftype = get_feature_type(target)
        product = get_feature_product(target)
        product = product.replace("_", " ")
        logger.info(f"Processing {gene}|{ftype}|{product}|{f.location.start}|{f.location.end}")
        add_qualifiers(f, gene, ftype, product)

        new_sub_features = []
        for subfeat in f.sub_features:
            # skip stop codon subfeatures
            if subfeat.type == "stop_codon":
                continue

            add_qualifiers(subfeat, gene, ftype, product)
            new_sub_features.append(subfeat)
            logger.debug(f"  {subfeat}")

        f.sub_features = new_sub_features
        cds_subfeatures = [x for x in f.sub_features if x.type == "CDS"]
        if not cds_subfeatures:
            continue
        cds_subfeatures.sort(key=lambda x: x.location.start)
        gene_length_minus_stop = f.location.end - 3
        for cds_subfeature in cds_subfeatures:
            if cds_subfeature.location.end == gene_length_minus_stop:
                cds_subfeature.location = FeatureLocation(
                    cds_subfeature.location.start,
                    f.location.end,
                )


def select_top_scoring_features_by_gene(
    gene_features: Mapping[Tuple[str, str, str], List[Tuple[int, SeqFeature]]]
) -> List[SeqFeature]:
    top_features: List[SeqFeature] = []
    for features in gene_features.values():
        features.sort(key=itemgetter(0), reverse=True)
        top_feature: SeqFeature = features[0][1]
        if top_feature.sub_features[-1].type != "stop_codon":
            logging.info(
                f"Top feature at {top_feature.location} does not have a stop codon subfeature: "
                f"{top_feature.sub_features}"
            )
            top_score = features[0][0]
            for score, f in features:
                if f == top_feature:
                    continue
                if score == top_score:
                    logging.info(f"Multiple features with top score {top_score} at {f.location}")
                    if f.sub_features[-1].type == "stop_codon":
                        logging.info(
                            f"Selecting feature at {f.location} as top feature since it has a stop codon subfeature"
                        )
                        top_feature = f
                        break
        top_features.append(top_feature)
    if logger.isEnabledFor(logging.DEBUG):
        for f in top_features:
            logging.debug(f"Top feature: {f}")
            for subf in f.sub_features:
                logging.debug(f"  Top subfeature: {subf}")
    return top_features
