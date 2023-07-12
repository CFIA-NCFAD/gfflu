from pathlib import Path  # noqa: I001
from tempfile import TemporaryDirectory

import gfflu.io
from gfflu.miniprot import run_miniprot, get_feature_qualifier_value, get_gene

datadir = Path(__file__).parent / "data"
input_fasta = datadir / "Segment_1_PB2.MH201221.fasta"


def test_miniprot():
    with TemporaryDirectory(prefix="gfflu_") as tmpdir:
        outdir = Path(tmpdir)
        miniprot_out = run_miniprot(outdir, input_fasta, input_fasta.stem)
        assert miniprot_out.exists()
        assert miniprot_out.stat().st_size > 0
        assert miniprot_out.name.endswith(".miniprot.gff")
        assert miniprot_out.name == "Segment_1_PB2.MH201221.miniprot.gff"
        gffrec = gfflu.io.read_gff(miniprot_out)
        assert len(gffrec.features) > 1
        assert gffrec.features[0].type == "mRNA"
        target_value = get_feature_qualifier_value(gffrec.features[0], "Target")
        assert get_gene(target_value) == "PB2"
