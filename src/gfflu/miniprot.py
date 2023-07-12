import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List, Mapping, Tuple

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from gfflu.peptides import iav_faa

logger = logging.getLogger(__name__)


def get_gene(target_value: str) -> str:
    """Get gene name from Miniprot target value

    >>> get_gene("NS2|CDS|nonstructural protein 2|seg8prot2A")
    'NS2'
    """
    return target_value.split("|")[0]


def get_feature_type(target_value: str) -> str:
    """Get feature type from Miniprot target value

    >>> get_feature_type("NS2|CDS|nonstructural protein 2|seg8prot2A")
    'CDS'
    """
    return target_value.split("|")[1]


def get_features_by_gene(gff_rec: SeqRecord) -> Mapping[Tuple[str, str, str], List[Tuple[int, SeqFeature]]]:
    """Read GFF file and return dict of (gene,type,product) -> features

    Within the Miniprot Target value will be the original sequence ID which will
    contain the gene name, feature type and product name. This is used to group
    unique features so that they can be sorted by score. The top scoring feature
    is then selected for each unique group. These will be used to generate the
    final GFF file.
    """
    gene_feature = defaultdict(list)
    for f in gff_rec.features:
        target = get_feature_qualifier_value(f, "Target")
        gene = get_gene(target)
        ftype = get_feature_type(target)
        product = get_feature_product(target)
        score = int(get_feature_qualifier_value(f, "score"))
        if f.type == "mRNA":
            f.type = "gene"
            gene_feature[(gene, ftype, product)].append((score, f))
    return gene_feature


def get_feature_product(target_value: str) -> str:
    """Get feature product from Miniprot target value

    >>> get_feature_product("NS2|CDS|nonstructural protein 2|seg8prot2A")
    'nonstructural protein 2'
    """
    return target_value.split("|")[2]


def get_feature_qualifier_value(feature: SeqFeature, key: str) -> str:
    """Get qualifier value from feature

    >>> get_feature_qualifier_value(feature, "Target")
    'NS2|CDS|nonstructural protein 2|seg8prot2A'
    """
    return feature.qualifiers[key][0]


def run_miniprot(outdir: Path, fasta: Path, prefix: str) -> Path:
    """Run Miniprot on FASTA file with Influenza A virus peptide sequences"""
    miniprot_out = outdir / f"{prefix}.miniprot.gff"
    logger.info(f"Running Miniprot on {fasta}")
    command = [
        "miniprot",
        "--gff",
        str(fasta.resolve().absolute()),
        str(iav_faa.resolve().absolute()),
    ]
    logger.info(f"Running command: $ {' '.join(command)}")
    subprocess.run(command, shell=False, check=True, stdout=miniprot_out.open("w"))  # noqa: S603
    return miniprot_out
