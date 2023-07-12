from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from gfflu.annotation import run_annotation
from gfflu.io import check_gff

datadir = Path(__file__).parent / "data"


def test_gff_to_genbank():
    input_fasta = datadir / "Segment_4_HA.MH201222.fasta"
    with TemporaryDirectory(prefix="gfflu-") as tmpdir:
        gffrec = run_annotation(fasta=input_fasta, outdir=Path(tmpdir), prefix=input_fasta.stem)
        assert gffrec is not None
        assert len(gffrec.features) == 1
        feature_ha = gffrec.features[0]
        assert feature_ha.type == "gene"
        assert feature_ha.qualifiers["gene"] == ["HA"]
        assert len(feature_ha.sub_features) == 4
        from io import StringIO

        sio = StringIO()
        gbkrec = check_gff([gffrec], "DNA")
        SeqIO.write(gbkrec, sio, "genbank")
        sio.seek(0)
        genbank = sio.read()
        assert genbank.startswith("LOCUS")
        assert "sig_peptide" in genbank
        assert "mat_peptide" in genbank
