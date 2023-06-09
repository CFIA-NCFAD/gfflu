#!/usr/bin/env python
import logging
import sys
from pathlib import Path
from typing import Optional

import typer
from Bio import SeqIO

from gfflu.__about__ import __version__
from gfflu.annotation import run_annotation
from gfflu.io import check_gff, write_gff

app = typer.Typer()

logger = logging.getLogger(__name__)


def init_logging(level: int):
    from rich.logging import RichHandler
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=level,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )


def version_callback(value: bool):
    if value:
        typer.echo(f"gfflu version {__version__}")
        raise typer.Exit()


@app.command(
    epilog=f"gfflu version {__version__}; Python {sys.version_info.major}.{sys.version_info.minor}."
    f"{sys.version_info.micro}"
)
def main(
    fasta: Path = typer.Argument(..., exists=True, dir_okay=False),
    outdir: Path = typer.Option(Path("gfflu-outdir"), "--outdir", "-o", help="Output directory"),
    force: bool = typer.Option(False, "--force", "-f", is_flag=True),
    verbose: bool = typer.Option(False, "--verbose", "-v", is_flag=True),
    version: Optional[bool] = typer.Option(  # noqa: ARG001
        None,
        "--version",
        "-V",
        callback=version_callback,
        help=f"Print 'gfflu version {__version__}' and exit",
    ),
):
    """Annotate Influenza A virus sequences using Miniprot and BLASTX


    The Miniprot GFF for a particular reference sequence gene segment will have multiple
    annotations for the same gene. This script will select the top scoring annotation for
    each gene and write out a new GFF file that can be used with SnpEff.

    """
    init_logging(logging.DEBUG if verbose else logging.INFO)
    logger.info(f"Running gfflu version {__version__}")
    logger.info(f"Output directory: {outdir}")
    logger.info(f"FASTA file: {fasta}")
    logger.info(f"Overwrite output?: {force}")
    logger.info(f"Verbose?: {verbose}")
    if not outdir.exists():
        outdir.mkdir(parents=True)
    elif force:
        logger.warning(f"Output directory '{outdir}' exists, overwriting")
    else:
        logger.error(f"Output directory '{outdir}' exists, exiting")
        raise typer.Exit(1)

    rec = run_annotation(fasta, outdir)

    output_gff = outdir / f"{fasta.with_suffix('.gff').name}"

    logger.info(f"Writing {output_gff}")
    write_gff([rec], output_gff)
    gbk_path = outdir / f"{fasta.with_suffix('.gbk').name}"
    logger.info(f"Writing {gbk_path}")
    SeqIO.write(check_gff([rec], "DNA"), gbk_path, "genbank")


if __name__ == "__main__":
    app()
