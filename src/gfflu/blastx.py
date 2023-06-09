import importlib
import importlib.resources
import logging
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)


def run_blastx(outdir: Path, fasta: Path) -> Path:
    """Run BLASTX on FASTA file with Influenza A virus peptide sequences"""
    blastx_out = outdir / fasta.with_suffix(".blastx.tsv").name
    logger.info(f"Running BLASTX on {fasta}")
    iav_faa = importlib.resources.path("gfflu.data", "iav-annotation.faa")
    command = [
        "blastx",
        "-word_size",
        "3",
        "-evalue",
        "0.1",
        "-db",
        str(iav_faa),
        "-query",
        str(fasta.resolve().absolute()),
        "-outfmt",
        "6",
        "-out",
        str(blastx_out.resolve().absolute()),
    ]
    logger.info(f"Running command: $ {' '.join(command)}")
    subprocess.run(command, shell=False, check=True)  # noqa: S603
    return blastx_out
